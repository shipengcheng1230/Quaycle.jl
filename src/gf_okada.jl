## static green's function
export okada_disp_gf_tensor, gen_alloc

## okada displacement green's function

@gen_shared_chunk_call okada_disp_gf_tensor false

"Okada green's function in 1-D elastic fault in [`LineOkadaMesh`](@ref)."
function okada_disp_gf_tensor(mesh::LineOkadaMesh, λ::T, μ::T, ft::PlaneFault; kwargs...) where T
    st = SharedArray{T}(mesh.nξ, mesh.nξ)
    okada_disp_gf_tensor!(st, mesh, λ, μ, ft; kwargs...)
    return sdata(st)
end

"Okada green's function in 2-D elastic fault in [`RectOkadaMesh`](@ref). Translational symmetry is considered."
function okada_disp_gf_tensor(mesh::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault; fourier_domain=true, kwargs...) where T
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

function okada_disp_gf_tensor_chunk!(st::SharedArray{T, 2}, subs::AbstractArray, mesh::LineOkadaMesh, λ::T, μ::T, ft::PlaneFault; ax_ratio::Real=12.5) where T
    ud = unit_dislocation(ft)
    ax = mesh.nξ * mesh.Δξ * ax_ratio .* [-one(T), one(T)]
    α = (λ + μ) / (λ + 2μ)
    pm = Progress(length(subs); dt=0.5, desc="Computing green's function ...", color=:light_cyan)
    @inbounds @simd for sub in subs
        i, j = sub[1], sub[2]
        u = dc3d(mesh.x, mesh.y[i], mesh.z[i], α, mesh.dep, mesh.dip, ax, mesh.aξ[j], ud)
        st[i,j] = shear_traction(ft, u, λ, μ, mesh.dip)
        next!(pm)
    end
end

function okada_disp_gf_tensor_chunk!(st::SharedArray{T, 3}, subs::AbstractArray, mesh::RectOkadaMesh, λ::T, μ::T, ft::PlaneFault; nrept::Integer=2, buffer_ratio::Integer=0) where T
    ud = unit_dislocation(ft)
    lrept = (buffer_ratio + one(T)) * (mesh.Δx * mesh.nx)
    u = Vector{T}(undef, 12)
    α = (λ + μ) / (λ + 2μ)
    pm = Progress(length(subs); dt=0.5, desc="Computing green's function ...", color=:light_cyan)
    for sub in subs
        i, j, l = sub[1], sub[2], sub[3]
        okada_gf_periodic_bc!(u, mesh.x[i], mesh.y[j], mesh.z[j], α, mesh.dep, mesh.dip, mesh.ax[1], mesh.aξ[l], ud, nrept, lrept)
        @inbounds st[i,j,l] = shear_traction(ft, u, λ, μ, mesh.dip)
        next!(pm)
    end
end

"""
Periodic summation of green's function.
- `lrept` represents jumping periodic block for simulating locked region at along-strike edge.
- `nrept` represents number of blocks to be summed.
"""
function okada_gf_periodic_bc!(u::AbstractVector{T}, x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T},
    nrept::Integer, lrept::T) where {T<:Real}
    fill!(u, zero(T))
    @fastmath @simd for i = -nrept: nrept
        u .+= dc3d(x, y, z, α, depth, dip, al .+ i * lrept, aw, disl)
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

"Allocation for computing stress rate in 1-D elastic fault."
struct OkadaGFAllocMatrix{T<:AbstractVecOrMat{<:Real}} <: OkadaGFAllocation{1}
    dims::Dims{1}
    dτ_dt::T # stress rate
    relv::T # relative velocity
end

"Allocation for computing stress rate in 2-D elastic fault."
struct OkadaGFAllocFFTConv{T<:AbstractArray{<:Real}, U<:AbstractArray{<:Complex}, P<:Plan} <: OkadaGFAllocation{2}
    dims::Dims{2}
    dτ_dt::T # stress rate of interest
    relv::T # relative velocity
    dτ_dt_dft::U # stress rate in discrete fourier domain
    relv_dft::U # relative velocity in discrete fourier domain
    dτ_dt_buffer::T # stress rate including zero-padding zone for fft
    pf::P # fft operator
end

gen_alloc(mesh::LineOkadaMesh) = gen_alloc(mesh.nξ; T=typeof(mesh.Δξ))
gen_alloc(mesh::RectOkadaMesh) = gen_alloc(mesh.nx, mesh.nξ; T=typeof(mesh.Δx))

"Generate 1-D computation allocation for computing stress rate."
gen_alloc(nξ::Integer; T=Float64) = OkadaGFAllocMatrix((nξ,), [Vector{T}(undef, nξ) for _ in 1: 2]...)

"Generate 2-D computation allocation for computing stress rate."
function gen_alloc(nx::I, nξ::I; T=Float64) where {I <: Integer}
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

"Stress rate in 1-D elastic plane."
@inline function dτ_dt!(gf::AbstractArray{T, 2}, alloc::OkadaGFAllocMatrix, vpl::T, v::AbstractVector) where T<:Number
    @fastmath @threads for i = 1: alloc.dims[1]
        @inbounds alloc.relv[i] = vpl - v[i]
    end
    mul!(alloc.dτ_dt, gf, alloc.relv)
end

"Stress rate in 2-D elastic plane. Using FFT to convolving translational symmetric tensor."
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

## okada strain green's function

"Compute 6 strain components from [`dc3d`](@ref) output."
function strain_components!(ϵ::T, u::T) where T<:AbstractVector
    ϵ[1] = u[4]
    ϵ[2] = (u[5] + u[7]) / 2
    ϵ[3] = (u[6] + u[10]) / 2
    ϵ[4] = u[8]
    ϵ[5] = (u[9] + u[11]) / 2
    ϵ[6] = u[12]
end

@gen_shared_chunk_call okada_strain_gf_tensor true

"Compute strain green's function from [`RectOkadaMesh`](@ref) to [`SBarbotTet4MeshEntity`](@ref) or [`SBarbotHex8MeshEntity`](@ref)"
function okada_strain_gf_tensor(mf::RectOkadaMesh, ma::SBarbotMeshEntity{3}, λ::T, μ::T, ft::PlaneFault, comp::AbstractVector; kwargs...) where T
    st = ntuple(_ -> SharedArray{T}(mf.nx * mf.nξ, length(ma.tag)), Val(length(comp)))
    okada_strain_gf_tensor!(st, mf, ma, λ, μ, ft, comp; kwargs...)
    return [sdata(x) for x in st]
end

function okada_strain_gf_tensor_chunk!(
    st::NTuple{N, <:SharedArray}, subs::AbstractArray, mf::RectOkadaMesh, ma::SBarbotMeshEntity{3},
    λ::T, μ::T, ft::PlaneFault, comp::AbstractVector; nrept::Integer=2, buffer_ratio::Integer=0
    ) where {T, N}
    ud = unit_dislocation(ft)
    lrept = (buffer_ratio + one(T)) * (mf.Δx * mf.nx)
    α = (λ + μ) / (λ + 2μ)
    u = Vector{T}(undef, 12)
    ϵ = Vector{T}(undef, 6)
    i2s = CartesianIndices((mf.nx, mf.nξ))

    @inbounds @fastmath @simd for sub in subs
        i, j = sub[1], sub[2] # index of fault, index of volume
        q = i2s[i] # return (ix, iξ)
        okada_gf_periodic_bc!(u, ma.x2[j], ma.x1[j], -ma.x3[j], α, mf.dep, mf.dip, mf.ax[q[1]], mf.aξ[q[2]], ud, nrept, lrept)
        strain_components!(ϵ, u)
        for (ic, c) in enumerate(comp)
            st[ic][i,j] = ϵ[c]
        end
    end
end
