## static green's function
export greens_function

@traitdef IsOnYZPlane{X}
@traitimpl IsOnYZPlane{DIPPING}

const KERNELDIR = joinpath(@__DIR__, "gfkernels")

for f in filter!(x -> endswith(x, ".jl"), readdir(KERNELDIR))
    include(abspath(joinpath(KERNELDIR, f)))
end

## okada
greens_function(::Val{:okada}, mesh::SimpleMesh, λ::Real, μ::Real, dip::Real, ft::FT; kwargs...) where {FT<:PlaneFault} = okada_greens_function(mesh, promote(λ, μ, dip)..., ft; kwargs...)

@inline @traitfn function shear_traction(::FT, u::AbstractVector, λ::T, μ::T, dip::T) where {T<:Real, FT<:PlaneFault; IsOnYZPlane{FT}}
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180))
end

@inline @traitfn function shear_traction(::FT, u::AbstractVector, λ::T, μ::T, dip::T) where {T<:Real, FT<:PlaneFault; !IsOnYZPlane{FT}}
    # Special case for `dip = 90.` against above, however, at x-y plane instead.
    μ * (u[5] + u[7])
end

function okada_greens_function(mesh::SimpleMesh{1}, λ::T, μ::T, dip::T, ft::PlaneFault; kwargs...) where T
    st = SharedArray{T}(mesh.nξ, mesh.nξ)
    okada_greens_function!(st, mesh, λ, μ, dip, ft; kwargs...)
    return sdata(st)
end

function okada_greens_function(mesh::SimpleMesh{2}, λ::T, μ::T, dip::T, ft::PlaneFault; fourier_domain=true, kwargs...) where T
    st = SharedArray{T}(mesh.nx, mesh.nξ, mesh.nξ)
    okada_greens_function!(st, mesh, λ, μ, dip, ft; kwargs...)

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

function okada_greens_function!(st::SharedArray, mesh::SimpleMesh, λ::T, μ::T, dip::T, ft::PlaneFault; kwargs...) where T
    @sync begin
        for p in procs(st)
            @async remotecall_wait(okada_gf_shared_chunk!, WorkerPool(workers()), st, mesh, λ, μ, dip, ft; kwargs...)
        end
    end
end

function okada_gf_shared_chunk!(st::SharedArray, mesh::SimpleMesh, λ::T, μ::T, dip::T, ft::PlaneFault; kwargs...) where T
    i2s = CartesianIndices(st)
    inds = localindices(st)
    subs = i2s[inds]
    okada_gf_chunk!(st, mesh, λ, μ, dip, ft, subs; kwargs...)
end

function okada_gf_chunk!(st::SharedArray{T, 2}, mesh::SimpleMesh{1}, λ::T, μ::T, dip::T, ft::PlaneFault, subs::AbstractArray; ax_ratio::Real=12.5) where T
    ud = unit_dislocation(ft)
    ax = mesh.nξ * mesh.Δξ * ax_ratio .* [-one(T), one(T)]
    α = (λ + μ) / (λ + 2μ)
    y, z = mesh.ξ .* cosd(dip), mesh.ξ .* sind(dip)
    @inbounds @simd for sub in subs
        i, j = sub[1], sub[2]
        u = dc3d_okada(zero(T), y[i], z[i], α, zero(T), dip, ax, mesh.aξ[j], ud)
        st[i,j] = shear_traction(ft, u, λ, μ, dip)
    end
end

function okada_gf_chunk!(st::SharedArray{T, 3}, mesh::SimpleMesh{2}, λ::T, μ::T, dip::T, ft::PlaneFault, subs::AbstractArray; nrept::Integer=2, buffer_ratio::Integer=1) where T
    ud = unit_dislocation(ft)
    lrept = (buffer_ratio + one(T)) * (mesh.Δx * mesh.nx)
    u = Vector{T}(undef, 12)
    α = (λ + μ) / (λ + 2μ)
    y, z = mesh.ξ .* cosd(dip), mesh.ξ .* sind(dip)
    for sub in subs
        i, j, l = sub[1], sub[2], sub[3]
        # For simple rectangular mesh, `depth` here are fixed at 0.
        okada_gf_periodic_bc!(u, mesh.x[i], y[j], z[j], α, zero(T), dip, mesh.ax[1], mesh.aξ[l], ud, nrept, lrept)
        st[i,j,l] = shear_traction(ft, u, λ, μ, dip)
    end
end

function okada_gf_periodic_bc!(u::AbstractVector{T}, x::T, y::T, z::T, α::T, depth::T, dip::T, al::AbstractVector{T}, aw::AbstractVector{T}, disl::AbstractVector{T},
    nrept::Integer, lrept::T) where {T<:Real}
    fill!(u, zero(T))
    @fastmath @simd for i = -nrept: nrept
        u .+= dc3d_okada(x, y, z, α, depth, dip, al .+ i * lrept, aw, disl)
    end
end

@inline @traitfn unit_dislocation(::FT) where {FT<:PlaneFault; IsOnYZPlane{FT}} = [0.0, 1.0, 0.0]
@inline @traitfn unit_dislocation(::FT) where {FT<:PlaneFault; !IsOnYZPlane{FT}} = [1.0, 0.0, 0.0]
