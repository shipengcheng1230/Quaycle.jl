## static green's function
export greens_tensor
include("gfkernels/dc3d.jl")

## okada
@gen_shared_chunk_call okada_disp_gf_tensor

function okada_disp_gf_tensor(mesh::LineTopCenterMesh, λ::T, μ::T, ft::PlaneFault; kwargs...) where T
    st = SharedArray{T}(mesh.nξ, mesh.nξ)
    okada_disp_gf_tensor!(st, mesh, λ, μ, ft; kwargs...)
    return sdata(st)
end

function okada_disp_gf_tensor(mesh::RectTopCenterMesh, λ::T, μ::T, ft::PlaneFault; fourier_domain=true, kwargs...) where T
    st = SharedArray{T}(mesh.nx, mesh.nξ, mesh.nξ)
    okada_disp_gf_tensor!(st, mesh, λ, μ, ft; kwargs...)

    function __convert_to_fourier_domain__()
        FFTW.set_num_threads(parameters["FFT"]["NUM_THREADS"])
        x1 = zeros(T, 2 * mesh.nx - 1)
        p1 = plan_rfft(x1, flags=parameters["FFT"]["FLAG"])
        st_dft = Array{Complex{T}}(undef, mesh.nx, mesh.nξ, mesh.nξ)
        @inbounds for l = 1: mesh.nξ, j = 1: mesh.nξ
            # The most tricky part to ensure correct FFT
            # Ref -> (https://github.com/JuliaMatrices/ToeplitzMatrices.jl/blob/cbe29c344be8363f33eb17090121f8cff600b72e/src/ToeplitzMatrices.jl#L627)
            st_dft[:,j,l] .= p1 * [st[:,j,l]; reverse(st[2:end,j,l])]
        end
        return st_dft
    end

    fourier_domain ? __convert_to_fourier_domain__() : sdata(st)
end

function okada_disp_gf_tensor_chunk!(st::SharedArray{T, 2}, subs::AbstractArray, mesh::LineTopCenterMesh, λ::T, μ::T, ft::PlaneFault; ax_ratio::Real=12.5) where T
    ud = unit_dislocation(ft)
    ax = mesh.nξ * mesh.Δξ * ax_ratio .* [-one(T), one(T)]
    α = (λ + μ) / (λ + 2μ)
    @inbounds @simd for sub in subs
        i, j = sub[1], sub[2]
        u = dc3d_okada(mesh.x, mesh.y[i], mesh.z[i], α, mesh.dep, mesh.dip, ax, mesh.aξ[j], ud)
        st[i,j] = shear_traction(ft, u, λ, μ, mesh.dip)
    end
end

function okada_disp_gf_tensor_chunk!(st::SharedArray{T, 3}, subs::AbstractArray, mesh::RectTopCenterMesh, λ::T, μ::T, ft::PlaneFault; nrept::Integer=2, buffer_ratio::Integer=0) where T
    ud = unit_dislocation(ft)
    lrept = (buffer_ratio + one(T)) * (mesh.Δx * mesh.nx)
    u = Vector{T}(undef, 12)
    α = (λ + μ) / (λ + 2μ)
    for sub in subs
        i, j, l = sub[1], sub[2], sub[3]
        okada_gf_periodic_bc!(u, mesh.x[i], mesh.y[j], mesh.z[j], α, mesh.dep, mesh.dip, mesh.ax[1], mesh.aξ[l], ud, nrept, lrept)
        @inbounds st[i,j,l] = shear_traction(ft, u, λ, μ, mesh.dip)
    end
end

function okada_gf_periodic_bc!(u::AbstractVector{T}, x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T},
    nrept::Integer, lrept::T) where {T<:Real}
    fill!(u, zero(T))
    @fastmath @simd for i = -nrept: nrept
        u .+= dc3d_okada(x, y, z, α, depth, dip, al .+ i * lrept, aw, disl)
    end
end

@inline unit_dislocation(::DIPPING) = [0.0, 1.0, 0.0]
@inline unit_dislocation(::STRIKING) = [1.0, 0.0, 0.0]

@inline function shear_traction(::DIPPING, u::AbstractVector, λ::T, μ::T, dip::T) where T
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sind(2dip) + τyz * cosd(2dip))
end

@inline function shear_traction(::STRIKING, u::AbstractVector, λ::T, μ::T, dip::T) where T
    # Special case for `dip = 90.` against above, however, at x-y plane instead.
    μ * (u[5] + u[7])
end

##
abstract type OkadaGFAllocation{dim} <: AbstractAllocation{dim} end

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

gen_alloc(gtype::Val{:okada}, mesh::LineTopCenterMesh) = gen_alloc(gtype, mesh.nξ; T=typeof(mesh.Δξ))
gen_alloc(gtype::Val{:okada}, mesh::RectTopCenterMesh) = gen_alloc(gtype, mesh.nx, mesh.nξ; T=typeof(mesh.Δx))
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

@inline function dτ_dt!(gf::AbstractArray{T, 2}, alloc::OkadaGFAllocMatrix, vpl::T, v::AbstractVector) where T<:Number
    @fastmath @threads for i = 1: alloc.dims[1]
        @inbounds alloc.relv[i] = vpl - v[i]
    end
    mul!(alloc.dτ_dt, gf, alloc.relv)
end

@inline function dτ_dt!(gf::AbstractArray{T, 3}, alloc::OkadaGFAllocFFTConv, vpl::U, v::AbstractMatrix) where {T<:Complex, U<:Number}
    @fastmath @threads for j = 1: alloc.dims[2]
        @simd for i = 1: alloc.dims[1]
            @inbounds alloc.relv[i,j] = vpl - v[i,j]
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
