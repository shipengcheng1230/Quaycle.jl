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

struct StressRateAllocMatrix{dim, T<:AbstractMatrix, V<:AbstractVector, I<:Integer} <: StressRateAllocation{dim} # unstructured mesh
    nume::I # number of unstructured mesh elements
    numϵ::I # number of strain components considered
    numσ::I # number of independent stress components
    reldϵ::T # relative strain rate
    σ′::T # deviatoric stress, <matrix>
    ς′::V # norm of deviatoric stress, <vector>

    function StressRateAllocMatrix(nume::I, numϵ::I, numσ::I, reldϵ::T, σ′::T, ς′::V) where {I, T, V, M}
        if numσ == 6
            dim = 3
        elseif numσ == 3
            dim = 2 # plane stress, remain further implementation
        elseif numσ == 2
            dim = 1 # antiplane stress, remain further implementation
        else
            error("Number of independent stress components should either be 6 (3-D) or 3 (2-D plane) or 2 (2-D antiplane).")
        end
        new{dim, T, V, I}(nume, numϵ, numσ, reldϵ, σ′, ς′)
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
    reldϵ = Matrix{T}(undef, nume, numϵ)
    σ′ = Matrix{T}(undef, nume, numσ)
    ς′ = Vector{T}(undef, nume)
    return StressRateAllocMatrix(nume, numϵ, numσ, reldϵ, σ′, ς′)
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

@inline function relative_strain_rate!(alloc::StressRateAllocation, dϵ₀::AbstractVector, dϵ::AbstractVecOrMat)
    @inbounds @fastmath for j in 1: alloc.numϵ
        @threads for i = 1: alloc.nume
            alloc.reldϵ[i,j] = dϵ[i,j] - dϵ₀[j]
        end
    end
end

@inline function deviatoric_stress!(σ::AbstractVecOrMat, alloc::StressRateAllocation{3})
    BLAS.blascopy!(6alloc.nume, σ, 1, alloc.σ′, 1)
    @inbounds @fastmath @threads for i = 1: alloc.nume
        σkk = (σ[i,1] + σ[i,4] + σ[i,6]) / 3
        alloc.σ′[i,1] -= σkk
        alloc.σ′[i,4] -= σkk
        alloc.σ′[i,6] -= σkk
    end
    # could use https://github.com/Jutho/Strided.jl
    alloc.ς′ .= sqrt.(vec(sum(abs2, alloc.σ′; dims=2))) # for higher precision use `hypot` or `norm`
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
@inline function dτ_dt!(gf::AbstractMatrix{T}, alloc::ViscoelasticCompositeAlloc) where T
    BLAS.gemv!('N', one(T), gf, vec(alloc.v.reldϵ), one(T), vec(alloc.e.dτ_dt))
end

"Stress rate from (ℕ)-D elastic entity to (ℕ+1)-D inelastic entity."
@inline function dσ_dt!(dσ::AbstractVecOrMat, gf::AbstractMatrix{T}, alloc::TractionRateAllocation) where T
    BLAS.gemv!('N', one(T), gf, vec(alloc.relvnp), zero(T), vec(dσ))
end

"Stress rate within inelastic entity."
@inline function dσ_dt!(dσ::AbstractVecOrMat, gf::AbstractMatrix{T}, alloc::StressRateAllocMatrix) where T
    BLAS.gemv!('N', one(T), gf, vec(alloc.reldϵ), one(T), vec(dσ))
end
