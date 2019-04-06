## operators of greens function

export gf_operator, gen_alloc

abstract type Allocation{dim} end
abstract type OkadaGFAllocation{dim} <: Allocation{dim} end


struct OkadaGFAllocMatrix{T<:AbstractVecOrMat{<:Real}} <: OkadaGFAllocation{1}
    dims::Dims{1}
    dτ_dt::T # stress rate
    relv::T # relative velocity
end


struct OkadaGFAllocFFTConv{T<:AbstractArray{<:Real}, U<:AbstractArray{<:Complex}, P<:Plan} <: OkadaGFAllocation{2}
    dims::Dims{2}
    dτ_dt::T # stress rate of interest
    relv::T # relative velocity
    dτ_dt_dft::U # stress rate in discrete fourier domain
    relv_dft::U # relative velocity in discrete fourier domain
    dτ_dt_buffer::T # stress rate including zero-padding zone for fft
    pf::P # fft operator
end


gen_alloc(gtype::Val{:okada}, mesh::SimpleLineGrid) = gen_alloc(gtype, mesh.nξ; T=typeof(mesh.Δξ))

gen_alloc(gtype::Val{:okada}, mesh::SimpleRectGrid) = gen_alloc(gtype, mesh.nx, mesh.nξ; T=typeof(mesh.Δx))

gen_alloc(gtype::Val{:okada}, nξ::Integer; T=Float64) = OkadaGFAllocMatrix((nξ,), [Vector{T}(undef, nξ) for _ in 1: 2]...)


function gen_alloc(::Val{:okada}, nx::I, nξ::I; T=Float64) where {I <: Integer}
    FFTW.set_num_threads(parameters["FFT"]["NUM_THREADS"])
    x1 = Matrix{T}(undef, 2 * nx - 1, nξ)
    p1 = plan_rfft(x1, 1, flags=parameters["FFT"]["FLAG"])

    alloc = OkadaGFAllocFFTConv(
        (nx, nξ),
        Matrix{T}(undef, nx, nξ),
        zeros(T, 2nx-1, nξ), # for relative velocity
        [Matrix{Complex{T}}(undef, nx, nξ) for _ in 1: 2]...,
        Matrix{T}(undef, 2nx-1, nξ),
        p1)
end


function gf_operator(::Val{:okada}, gf::AbstractArray, alloc::OkadaGFAllocation, vpl::T) where T
    f! = (alloc, v) -> dτ_dt!(gf, alloc, vpl, v)
    return f!
end


@inline function dτ_dt!(gf::AbstractArray{T, 2}, alloc::OkadaGFAllocMatrix, vpl::T, v::AbstractVector) where T<:Number
    @fastmath @inbounds @simd for i = 1: alloc.dims[1]
        alloc.relv[i] = vpl - v[i]
    end
    mul!(alloc.dτ_dt, gf, alloc.relv)
end


@inline function dτ_dt!(gf::AbstractArray{T, 3}, alloc::OkadaGFAllocFFTConv, vpl::U, v::AbstractMatrix) where {T<:Complex, U<:Number}
    @fastmath @inbounds for j = 1: alloc.dims[2]
        @simd for i = 1: alloc.dims[1]
            alloc.relv[i,j] = vpl - v[i,j]
        end
    end
    mul!(alloc.relv_dft, alloc.pf, alloc.relv)
    fill!(alloc.dτ_dt_dft, zero(T))

    @fastmath @inbounds for l = 1: alloc.dims[2]
        for j = 1: alloc.dims[2]
            @simd for i = 1: alloc.dims[1]
                alloc.dτ_dt_dft[i,j] += gf[i,j,l] * alloc.relv_dft[i,l]
            end
        end
    end

    ldiv!(alloc.dτ_dt_buffer, alloc.pf, alloc.dτ_dt_dft)

    @fastmath @inbounds for j = 1: alloc.dims[2]
        @simd for i = 1: alloc.dims[1]
            alloc.dτ_dt[i,j] = alloc.dτ_dt_buffer[i,j]
        end
    end
end
