## static green's function
export greens_tensor

const KERNELDIR = joinpath(@__DIR__, "gfkernels")

for f in filter!(x -> endswith(x, ".jl"), readdir(KERNELDIR))
    include(abspath(joinpath(KERNELDIR, f)))
end

## okada
greens_tensor(::Val{:okada}, fault::BasicFaultSpace, λ::Real, μ::Real; kwargs...) where {FT<:PlaneFault} = okada_greens_tensor(fault, promote(λ, μ)...; kwargs...)

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

function okada_greens_tensor(fault::BasicFaultSpace{1}, λ::T, μ::T; kwargs...) where T
    st = SharedArray{T}(fault.mesh.nξ, fault.mesh.nξ)
    okada_greens_tensor!(st, fault, λ, μ; kwargs...)
    return sdata(st)
end

function okada_greens_tensor(fault::BasicFaultSpace{2}, λ::T, μ::T; fourier_domain=true, kwargs...) where T
    st = SharedArray{T}(fault.mesh.nx, fault.mesh.nξ, fault.mesh.nξ)
    okada_greens_tensor!(st, fault, λ, μ; kwargs...)

    function __convert_to_fourier_domain__()
        FFTW.set_num_threads(parameters["FFT"]["NUM_THREADS"])
        x1 = zeros(T, 2 * fault.mesh.nx - 1)
        p1 = plan_rfft(x1, flags=parameters["FFT"]["FLAG"])
        st_dft = Array{Complex{T}}(undef, fault.mesh.nx, fault.mesh.nξ, fault.mesh.nξ)
        @inbounds for l = 1: fault.mesh.nξ, j = 1: fault.mesh.nξ
            # The most tricky part to ensure correct FFT
            # Ref -> (https://github.com/JuliaMatrices/ToeplitzMatrices.jl/blob/cbe29c344be8363f33eb17090121f8cff600b72e/src/ToeplitzMatrices.jl#L627)
            st_dft[:,j,l] .= p1 * [st[:,j,l]; reverse(st[2:end,j,l])]
        end
        return st_dft
    end

    fourier_domain ? __convert_to_fourier_domain__() : sdata(st)
end

function okada_greens_tensor!(st::SharedArray, fault::BasicFaultSpace, λ::T, μ::T; kwargs...) where T
    @sync begin
        for p in procs(st)
            @async remotecall_wait(okada_gf_shared_chunk!, WorkerPool(workers()), st, fault, λ, μ; kwargs...)
        end
    end
end

function okada_gf_shared_chunk!(st::SharedArray, fault::BasicFaultSpace, λ::T, μ::T; kwargs...) where T
    i2s = CartesianIndices(st)
    inds = localindices(st)
    subs = i2s[inds]
    okada_gf_chunk!(st, fault, λ, μ, subs; kwargs...)
end

function okada_gf_chunk!(st::SharedArray{T, 2}, fault::BasicFaultSpace{1}, λ::T, μ::T, subs::AbstractArray; ax_ratio::Real=12.5) where T
    ud = unit_dislocation(fault.ft)
    ax = fault.mesh.nξ * fault.mesh.Δξ * ax_ratio .* [-one(T), one(T)]
    α = (λ + μ) / (λ + 2μ)
    @inbounds @simd for sub in subs
        i, j = sub[1], sub[2]
        u = dc3d_okada(fault.mesh.x, fault.mesh.y[i], fault.mesh.z[i], α, fault.mesh.dep, fault.mesh.dip, ax, fault.mesh.aξ[j], ud)
        st[i,j] = shear_traction(fault.ft, u, λ, μ, fault.mesh.dip)
    end
end

function okada_gf_chunk!(st::SharedArray{T, 3}, fault::BasicFaultSpace{2}, λ::T, μ::T, subs::AbstractArray; nrept::Integer=2, buffer_ratio::Integer=1) where T
    ud = unit_dislocation(fault.ft)
    lrept = (buffer_ratio + one(T)) * (fault.mesh.Δx * fault.mesh.nx)
    u = Vector{T}(undef, 12)
    α = (λ + μ) / (λ + 2μ)
    for sub in subs
        i, j, l = sub[1], sub[2], sub[3]
        # For simple rectangular mesh, `depth` here are fixed at 0.
        okada_gf_periodic_bc!(u, fault.mesh.x[i], fault.mesh.y[j], fault.mesh.z[j], α, fault.mesh.dep, fault.mesh.dip, fault.mesh.ax[1], fault.mesh.aξ[l], ud, nrept, lrept)
        @inbounds st[i,j,l] = shear_traction(fault.ft, u, λ, μ, fault.mesh.dip)
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
