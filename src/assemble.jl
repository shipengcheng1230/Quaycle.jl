## assemble the system derivative function

export assemble

function assemble(
    stype::Val{:okada}, fs::CentralSymmetryFS, sp::HomoFaultProperties, fp::RSFrictionalProperties,
    u0::AbstractArray, tspan::NTuple{2};
    kwargs...)
    gf = greens_function(stype, fs.mesh, sp.λ, sp.μ, fs.dip, fs.faulttype; kwargs...)
    alloc = gen_alloc(stype, fs.mesh)
    gfop! = gf_operator(stype, gf, alloc, sp.vpl)
    f! = (du, u, p, t) -> derivative_kernel!(du, u, p[1], p[2], alloc, gfop!)
    return ODEProblem(f!, u0, tspan, (fp, sp))
end

@inline function dv_dθ_dt!(::CForm,
    dv::T, dθ::T, dδ::T, v::T, θ::T, frp::RSFrictionalProperties, fap::HomoFaultProperties, alloc::OkadaGFAllocation
    ) where {T<:AbstractVecOrMat}
    @fastmath @inbounds @threads for i = 1: prod(alloc.dims)
        dδ[i] = v[i]
        dμ_dθ = frp.σ[i] * frp.b[i] / θ[i]
        dμ_dv = frp.σ[i] * frp.a[i] / v[i]
        dθ[i] = dθ_dt(frp.sel, v[i], θ[i], frp.L[i])
        dv[i] = dv_dt(alloc.dτ_dt[i], dμ_dv, dμ_dθ, dθ[i], fap.η)
    end
end

@inline function dv_dθ_dt!(::RForm,
    dv::T, dθ::T, dδ::T, v::T, θ::T, frp::RSFrictionalProperties, fap::HomoFaultProperties, alloc::OkadaGFAllocation
    ) where {T<:AbstractVecOrMat}
    @fastmath @inbounds @threads for i = 1: prod(alloc.dims)
        dδ[i] = v[i]
        ψ1 = exp((fap.f0 + frp.b[i] * log(fap.v0 * θ[i] / frp.L[i])) / frp.a[i]) / 2fap.v0
        ψ2 = frp.σ[i] * ψ1 / hypot(1, v[i] * ψ1)
        dμ_dv = frp.a[i] * ψ2
        dμ_dθ = frp.b[i] / θ[i] * v[i] * ψ2
        dθ[i] = dθ_dt(frp.sel, v[i], θ[i], frp.L[i])
        dv[i] = dv_dt(alloc.dτ_dt[i], dμ_dv, dμ_dθ, dθ[i], fap.η)
    end
end

function derivative_kernel!(
    du::AbstractArray{T}, u::AbstractArray{T}, frp::RSFrictionalProperties, fap::HomoFaultProperties,
    alloc::OkadaGFAllocation{dim}, gfop!::Function,
    ) where {T<:Real, dim}
    _ndim = dim + 1
    v = selectdim(u, _ndim, 1)
    θ = selectdim(u, _ndim, 2)
    dv = selectdim(du, _ndim, 1)
    dθ = selectdim(du, _ndim, 2)
    dδ = selectdim(du, _ndim, 3)
    clamp!(θ, zero(T), Inf)
    gfop!(alloc, v)
    dv_dθ_dt!(frp.flf, dv, dθ, dδ, v, θ, frp, fap, alloc)
end
