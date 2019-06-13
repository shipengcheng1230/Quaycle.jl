export gen_alloc

## Allocation for inplace operations
abstract type AbstractAllocation{dim} end
abstract type TractionRateAllocation{dim} <: AbstractAllocation{dim} end
abstract type StressRateAllocation{dim} <: AbstractAllocation{dim} end
abstract type CompositeAllocation{dim} <: AbstractAllocation{dim} end

struct TractionRateAllocMatrix{T<:AbstractVecOrMat{<:Real}} <: TractionRateAllocation{1}
    dims::Dims{1}
    dτ_dt::T # traction rate
    dμ_dv::T # derivative of friction over velocity
    dμ_dθ::T # derivative of friction over state
    relv::T # relative velocity
end

struct TractionRateAllocFFTConv{T<:AbstractArray{<:Real}, U<:AbstractArray{<:Complex}, P<:FFTW.Plan} <: TractionRateAllocation{2}
    dims::Dims{2}
    dτ_dt::T # traction rate of interest
    dμ_dv::T # derivative of friction over velocity
    dμ_dθ::T # derivative of friction over state
    relv::T # relative velocity including zero-padding
    relvnp::T # relative velocity excluding zero-padding area
    dτ_dt_dft::U # stress rate in discrete fourier domain
    relv_dft::U # relative velocity in discrete fourier domain
    dτ_dt_buffer::T # stress rate including zero-padding zone for fft
    pf::P # real-value-FFT forward operator
end

struct StressRateAllocMatrix{dim, T, V, I<:Integer} <: StressRateAllocation{dim} # unstructured mesh
    nume::I # number of unstructured mesh elements
    numϵ::I # number of strain components considered
    numσ::I # number of independent stress components
    relϵ::T # relative strain rate
    dσ′_dt::T # deviatoric stress rate
    dς′_dt::V # norm of deviatoric stress rate

    function StressRateAllocMatrix(nume::I, numϵ::I, numσ::I, relϵ::T, dσ′_dt::T, dς′_dt::V) where {I, T, V}
        if numσ == 6
            dim = 3
        elseif numσ == 3
            dim = 2 # plane stress, remain further implementation
        elseif numσ == 2
            dim = 1 # antiplane stress, remain further implementation
        else
            error("Number of independent stress components should either be 6 (3-D) or 3 (2-D plane) or 2 (2-D antiplane).")
        end
        new{dim, T, V, I}(nume, numϵ, numσ, relϵ, dσ′_dt, dς′_dt)
    end
end

struct ViscoelasticCompositeAlloc{dim, A, B} <: CompositeAllocation{dim}
    e::A # traction rate allocation
    v::B # stress rate allocation

    function ViscoelasticCompositeAlloc(e::TractionRateAllocation{N}, v::StressRateAllocMatrix) where N
        new{N, typeof(e), typeof(v)}(e, v)
    end
end

"Generate 1-D computation allocation for computing traction rate."
gen_alloc(nξ::Integer; T=Float64) = TractionRateAllocMatrix((nξ,), [Vector{T}(undef, nξ) for _ in 1: 4]...)

"Generate 2-D computation allocation for computing traction rate."
function gen_alloc(nx::I, nξ::I; T=Float64) where I <: Integer
    FFTW.set_num_threads(parameters["FFT"]["NUM_THREADS"])
    x1 = Matrix{T}(undef, 2 * nx - 1, nξ)
    p1 = plan_rfft(x1, 1, flags=parameters["FFT"]["FLAG"])
    dτ_dt =

    return TractionRateAllocFFTConv(
        (nx, nξ),
        [Matrix{T}(undef, nx, nξ) for _ in 1: 3]...,
        zeros(T, 2nx-1, nξ), zeros(T, nx, nξ), # for relative velocity, including zero
        [Matrix{Complex{T}}(undef, nx, nξ) for _ in 1: 2]...,
        Matrix{T}(undef, 2nx-1, nξ),
        p1)
end

"Generate 3-D computation allocation for computing stress rate."
function gen_alloc(nume::I, numϵ::I, numσ::I; T=Float64) where I<:Integer
    relϵ = Matrix{T}(undef, nume, numϵ)
    relϵ = [Vector{T}(undef, nume) for _ in 1: numϵ]
    σ′ = [Vector{T}(undef, nume) for _ in 1: numσ]
    ς′ = Vector{T}(undef, nume)
    return StressRateAllocMatrix(nume, numϵ, numσ, relϵ, σ′, ς′)
end

gen_alloc(mesh::LineOkadaMesh) = gen_alloc(mesh.nξ; T=typeof(mesh.Δξ))
gen_alloc(mesh::RectOkadaMesh) = gen_alloc(mesh.nx, mesh.nξ; T=typeof(mesh.Δx))
gen_alloc(me::SBarbotMeshEntity{3}, numϵ::Integer) = gen_alloc(length(me.tag), numϵ, 6; T=eltype(me.x1))
gen_alloc(mf::OkadaMesh, me::SBarbotMeshEntity{3}, numϵ::Integer) = ViscoelasticCompositeAlloc(gen_alloc(mf), gen_alloc(me, numϵ))

## traction & stress rate operators
@inline function relative_velocity!(alloc::TractionRateAllocMatrix, vpl::T, v::AbstractVector) where T
    @inbounds @fastmath @threads for i = 1: alloc.dims[1]
        alloc.relv[i] = v[i] - vpl
    end
end

@inline function relative_velocity!(alloc::TractionRateAllocFFTConv, vpl::T, v::AbstractMatrix) where T
    @inbounds @fastmath @threads for j = 1: alloc.dims[2]
        @simd for i = 1: alloc.dims[1]
            alloc.relv[i,j] = v[i,j] - vpl # there are zero paddings in `alloc.relv`
            alloc.relvnp[i,j] = alloc.relv[i,j] # copy-paste, useful for `LinearAlgebra.BLAS`
        end
    end
end

@inline function relative_strain!(alloc::StressRateAllocation, ϵ₀::AbstractVector, ϵ::AbstractVecOrMat)
    @inbounds @fastmath for j = 1: alloc.numϵ
        @threads for i = 1: alloc.nume
            alloc.relϵ[j][i] = ϵ[j][i] - ϵ₀[j]
        end
    end
end

@inline function deviatoric_stress!(alloc::StressRateAllocation{3})
    @inbounds @fastmath @threads for i = 1: alloc.nume
        σkk = (alloc.dσ′_dt[1][i] + alloc.dσ′_dt[4][i] + alloc.dσ′_dt[6][i]) / 3
        alloc.dσ′_dt[1][i] -= σkk
        alloc.dσ′_dt[4][i] -= σkk
        alloc.dσ′_dt[6][i] -= σkk
    end
end

@inline function stress_norm!(alloc::StressRateAllocation)
    fill!(alloc.dς′_dt, zero(eltype(alloc.dς′_dt)))
    @inbounds @fastmath @threads for i = 1: alloc.nume
        @simd for c = 1: alloc.numσ
            alloc.dς′_dt[i] += alloc.dσ′_dt[c][i]^2 # for higher precision use `hypot` or `norm`
        end
        alloc.dς′_dt[i] = sqrt(alloc.dς′_dt[i])
    end
end

"Traction rate within 1-D elastic plane."
@inline function dτ_dt!(gf::AbstractArray{T, 2}, alloc::TractionRateAllocMatrix) where T<:Number
    mul!(alloc.dτ_dt, gf, alloc.relv)
end

"Traction rate within 2-D elastic plane. Using FFT to convolving translational symmetric tensor."
@inline function dτ_dt!(gf::AbstractArray{T, 3}, alloc::TractionRateAllocFFTConv) where {T<:Complex, U<:Number}
    mul!(alloc.relv_dft, alloc.pf, alloc.relv)
    fill!(alloc.dτ_dt_dft, zero(T))
    @inbounds @fastmath @threads for j = 1: alloc.dims[2]
        for l = 1: alloc.dims[2]
            @simd for i = 1: alloc.dims[1]
                @inbounds alloc.dτ_dt_dft[i,j] += gf[i,j,l] * alloc.relv_dft[i,l]
            end
        end
    end
    ldiv!(alloc.dτ_dt_buffer, alloc.pf, alloc.dτ_dt_dft)
    @inbounds @fastmath @threads for j = 1: alloc.dims[2]
        @simd for i = 1: alloc.dims[1]
            @inbounds alloc.dτ_dt[i,j] = alloc.dτ_dt_buffer[i,j]
        end
    end
end

"Traction rate from (ℕ+1)-D inelastic entity to (ℕ)-D elastic entity."
@inline function dτ_dt!(gf::NTuple{N, <:AbstractMatrix{T}}, alloc::ViscoelasticCompositeAlloc) where {T, N}
    @inbounds for i = 1: N
        BLAS.gemv!('N', one(T), gf[i], alloc.v.relϵ[i], one(T), vec(alloc.e.dτ_dt))
    end
end

"Stress rate from (ℕ)-D elastic entity to (ℕ+1)-D inelastic entity."
@inline function dσ_dt!(gf::NTuple{N, <:AbstractMatrix{T}}, alloc::ViscoelasticCompositeAlloc) where {T, N}
    @inbounds for i = 1: N
        BLAS.gemv!('N', one(T), gf[i], vec(alloc.e.relvnp), zero(T), alloc.v.dσ′_dt[i])
    end
end

"Stress rate within inelastic entity."
@inline function dσ_dt!(gf::NTuple{N1, NTuple{N2, <:AbstractMatrix{T}}}, alloc::StressRateAllocMatrix) where {T, N1, N2}
    @inbounds for j = 1: N1, i = 1: N2
        BLAS.gemv!('N', one(T), gf[j][i], alloc.relϵ[j], one(T), alloc.dσ′_dt[i])
    end
end
