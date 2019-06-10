export gen_alloc

## Allocation for inplace operations
abstract type AbstractAllocation{dim} end
abstract type TractionRateAllocation{dim} <: AbstractAllocation{dim} end
abstract type StressRateAllocation{dim} <: AbstractAllocation{dim} end
abstract type CompositeAllocation{dim} <: AbstractAllocation{dim} end

struct TractionRateAllocMatrix{T<:AbstractVecOrMat{<:Real}} <: TractionRateAllocation{1}
    dims::Dims{1}
    dτ_dt::T # traction rate
    dμ_dv::T
    dμ_dθ::T
    relv::T # relative velocity
end

struct TractionRateAllocFFTConv{T<:AbstractArray{<:Real}, U<:AbstractArray{<:Complex}, P<:FFTW.Plan} <: TractionRateAllocation{2}
    dims::Dims{2}
    dτ_dt::T # traction rate of interest
    dμ_dv::T
    dμ_dθ::T
    relv::T # relative velocity
    dτ_dt_dft::U # stress rate in discrete fourier domain
    relv_dft::U # relative velocity in discrete fourier domain
    dτ_dt_buffer::T # stress rate including zero-padding zone for fft
    pf::P # fft operator
end

struct StressRateAllocMatrix{T, I<:Integer} <: StressRateAllocation{-1} # unstructured mesh
    nume::I # number of unstructured elements
    numϵ::I # number of strain components
    relϵ::T # relative strain rate
    σ′::T # deviatoric stress
    σ′_norm::T # norm of deviatoric stress
end

struct ViscoelasticCompositeAlloc{dim, A, B} <: CompositeAllocation{dim}
    e::A # traction rate
    v::B # stress rate

    function ViscoelasticCompositeAlloc(e::TractionRateAllocation{N}, v::StressRateAllocMatrix) where N
        new{N, typeof(e), typeof(v)}(e, v)
    end
end

"Generate 1-D computation allocation for computing traction rate."
gen_alloc(nξ::Integer; T=Float64) = TractionRateAllocMatrix((nξ,), [Vector{T}(undef, nξ) for _ in 1: 4]...)

"Generate 2-D computation allocation for computing traction rate."
function gen_alloc(nx::I, nξ::I; T=Float64) where {I <: Integer}
    FFTW.set_num_threads(parameters["FFT"]["NUM_THREADS"])
    x1 = Matrix{T}(undef, 2 * nx - 1, nξ)
    p1 = plan_rfft(x1, 1, flags=parameters["FFT"]["FLAG"])

    return TractionRateAllocFFTConv(
        (nx, nξ),
        [Matrix{T}(undef, nx, nξ) for _ in 1: 3]...,
        zeros(T, 2nx-1, nξ), # for relative velocity
        [Matrix{Complex{T}}(undef, nx, nξ) for _ in 1: 2]...,
        Matrix{T}(undef, 2nx-1, nξ),
        p1)
end

"Generate 3-D computation allocation for computing stress rate."
function gen_alloc(nume::Integer, numϵ::Integer, T=Float64)
    relϵ = Matrix{T}(undef, nume, numϵ)
    σ′ = Matrix{T}(undef, nume, 6) # 6 independent components in 3-dimension space
    σ′_norm = Vector{T}(undef, nume)
    return StressRateAllocMatrix(nume, numϵ, relϵ, σ′, σ′_norm)
end

gen_alloc(mesh::LineOkadaMesh) = gen_alloc(mesh.nξ; T=typeof(mesh.Δξ))
gen_alloc(mesh::RectOkadaMesh) = gen_alloc(mesh.nx, mesh.nξ; T=typeof(mesh.Δx))
gen_alloc(me::SBarbotMeshEntity{3}, numϵ::Integer) = gen_alloc(length(me.tag), numϵ, eltype(me.x1))
gen_alloc(mf::OkadaMesh, me::SBarbotMeshEntity{3}, numϵ::Integer) = ViscoelasticCompositeAlloc(gen_alloc(mf), gen_alloc(me, numϵ))

## traction & stress rate operators
@inline function relative_velocity(alloc::TractionRateAllocMatrix, vpl::T, v::AbstractVector) where T
    @fastmath @threads for i = 1: alloc.dims[1]
        @inbounds alloc.relv[i] = v[i] - vpl
    end
end

@inline function relative_velocity(alloc::TractionRateAllocFFTConv, vpl::T, v::AbstractMatrix) where T
    @fastmath @threads for j = 1: alloc.dims[2]
        @simd for i = 1: alloc.dims[1]
            @inbounds alloc.relv[i,j] = v[i,j] - vpl # there are zero paddings in `alloc.relv`
        end
    end
end

@inline function relative_strain(alloc::StressRateAllocation, ϵ₀::AbstractVector, ϵ::AbstractVecOrMat)
    @inbounds @fastmath for j = 1: alloc.numϵ
        @threads for i = 1: alloc.nume
            alloc.relϵ[i,j] = ϵ[i,j] - ϵ₀[j]
        end
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
    @fastmath @threads for j = 1: alloc.dims[2]
        for l = 1: alloc.dims[2]
            @simd for i = 1: alloc.dims[1]
                @inbounds alloc.dτ_dt_dft[i,j] += gf[i,j,l] * alloc.relv_dft[i,l]
            end
        end
    end
    ldiv!(alloc.dτ_dt_buffer, alloc.pf, alloc.dτ_dt_dft)
    @fastmath @threads for j = 1: alloc.dims[2]
        @simd for i = 1: alloc.dims[1]
            @inbounds alloc.dτ_dt[i,j] = alloc.dτ_dt_buffer[i,j]
        end
    end
end

"Traction rate from 3-D inelastic volume to 2-D plane."
@inline function dτ_dt!(gf::AbstractArray, alloc::ViscoelasticCompositeAlloc)
    BLAS.gemm!('N', 'N', LinearAlgebra.I, gf, alloc.v.relϵ, alloc.e.dτ_dt)
end

"Stress rate from 2-D plane to 3-D inelastic volume."
@inline function dσ_dt!(gf::NTuple{6, <:AbstractMatrix}, alloc::ViscoelasticCompositeAlloc, ϵ₀::AbstractVector, ϵ::AbstractArray)
    for i = 1: 6

    end
end

"Stress rate within 3-D inelastic volume."
@inline function dσ_dt!(gf::NTuple{6, <:AbstractMatrix}, alloc::StressRateAllocMatrix, ϵ₀::AbstractVector, ϵ::AbstractArray)

end
