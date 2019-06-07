export gen_alloc

## Allocation for inplace operations

abstract type AbstractAllocation{dim} end
abstract type OkadaRSFAllocation{dim} <: AbstractAllocation{dim} end
abstract type OkadaSBarbotRSFAllocation{dim} <: AbstractAllocation{dim} end

struct OkadaRSFAllocMatrix{T<:AbstractVecOrMat{<:Real}} <: OkadaRSFAllocation{1}
    dims::Dims{1}
    dτ_dt::T # stress rate
    dμ_dv::T
    dμ_dθ::T
    relv::T # relative velocity
end

struct OkadaRSFAllocFFTConv{T<:AbstractArray{<:Real}, U<:AbstractArray{<:Complex}, P<:FFTW.Plan} <: OkadaRSFAllocation{2}
    dims::Dims{2}
    dτ_dt::T # stress rate of interest
    dμ_dv::T
    dμ_dθ::T
    relv::T # relative velocity
    dτ_dt_dft::U # stress rate in discrete fourier domain
    relv_dft::U # relative velocity in discrete fourier domain
    dτ_dt_buffer::T # stress rate including zero-padding zone for fft
    pf::P # fft operator
end

struct OkadaSBarbotAllocFFTConv{T} <: OkadaSBarbotRSFAllocation{3}
    dims::Dims
    okada::OkadaRSFAllocFFTConv # okada fft conv allocation
    relϵ::T # relative strain rate
    σ′::T # deviatoric stress
    σ′_norm::T # norm of deviatoric stress
end

gen_alloc(mesh::LineOkadaMesh) = gen_alloc(mesh.nξ; T=typeof(mesh.Δξ))
gen_alloc(mesh::RectOkadaMesh) = gen_alloc(mesh.nx, mesh.nξ; T=typeof(mesh.Δx))

"Generate 1-D computation allocation for computing stress rate."
gen_alloc(nξ::Integer; T=Float64) = OkadaRSFAllocMatrix((nξ,), [Vector{T}(undef, nξ) for _ in 1: 4]...)

"Generate 2-D computation allocation for computing stress rate."
function gen_alloc(nx::I, nξ::I; T=Float64) where {I <: Integer}
    FFTW.set_num_threads(parameters["FFT"]["NUM_THREADS"])
    x1 = Matrix{T}(undef, 2 * nx - 1, nξ)
    p1 = plan_rfft(x1, 1, flags=parameters["FFT"]["FLAG"])

    alloc = OkadaRSFAllocFFTConv(
        (nx, nξ),
        [Matrix{T}(undef, nx, nξ) for _ in 1: 3]...,
        zeros(T, 2nx-1, nξ), # for relative velocity
        [Matrix{Complex{T}}(undef, nx, nξ) for _ in 1: 2]...,
        Matrix{T}(undef, 2nx-1, nξ),
        p1)
end

## Stress rate operators
"Stress rate in 1-D elastic plane."
@inline function dτ_dt!(gf::AbstractArray{T, 2}, alloc::OkadaRSFAllocMatrix, vpl::T, v::AbstractVector) where T<:Number
    @fastmath @threads for i = 1: alloc.dims[1]
        @inbounds alloc.relv[i] = v[i] - vpl
    end
    mul!(alloc.dτ_dt, gf, alloc.relv)
end

"Stress rate in 2-D elastic plane. Using FFT to convolving translational symmetric tensor."
@inline function dτ_dt!(gf::AbstractArray{T, 3}, alloc::OkadaRSFAllocFFTConv, vpl::U, v::AbstractMatrix) where {T<:Complex, U<:Number}
    @fastmath @threads for j = 1: alloc.dims[2]
        @simd for i = 1: alloc.dims[1]
            @inbounds alloc.relv[i,j] = v[i,j] - vpl
        end
    end
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
