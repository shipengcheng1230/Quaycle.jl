## assemble the system derivative function

export assemble

## single degree of freedom

@inline dτ_dt(K::T, v::T, vpl::T) where {T<:Number} = K * (vpl - v)

@inline dv_dt(dτdt::T, dμdv::T, dμdθ::T, dθdt::T, η::T) where {T<:Number} = (dτdt - dμdθ * dθdt) / (dμdv + η)

@inline function dv_dθ_dt(::CForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, ::Vararg{T},
    ) where {T<:Number, SE<:StateEvolutionLaw}
    dμdθ = σ * b / θ
    dμdv = σ * a / v
    dθdt = dθ_dt(se, v, θ, L)
    dτdt = dτ_dt(k, v, vpl)
    return dv_dt(dτdt, dμdv, dμdθ, dθdt, η), dθdt
end

@inline function dv_dθ_dt(::RForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T,
    ) where {T<:Number, SE<:StateEvolutionLaw}
    ψ1 = exp((f0 + b * log(v0 * clamp(θ, zero(T), Inf) / L)) / a) / 2v0
    ψ2 = σ * ψ1 / sqrt(1 + (v * ψ1)^2)
    dμdv = a * ψ2
    dμdθ = b / θ * v * ψ2
    dθdt = dθ_dt(se, v, θ, L)
    dτdt = dτ_dt(k, v, vpl)
    return dv_dt(dτdt, dμdv, dμdθ, dθdt, η), dθdt
end

function assemble(p::SingleDofRSFProperties, u0::AbstractArray, tspan::NTuple;
    flf::FrictionLawForm=CForm(), se::StateEvolutionLaw=DieterichStateLaw()) where T<:Real
    (typeof(flf) == RForm && p.η ≈ 0) && @warn "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt."
    op! = (du, u, p, t) -> du .= dv_dθ_dt(flf, se, u[1], u[2], p.a, p.b, p.L, p.k, p.σ, p.η, p.vpl, p.f0, p.v0)
    return ODEProblem(op!, u0, tspan, p)
end

## elastic okada system

@inline function dv_dθ_dt!(::CForm, se::StateEvolutionLaw,
    dv::T, dθ::T, v::T, θ::T, p::ElasticRSFProperties, alloc::OkadaGFAllocation
    ) where {T<:AbstractVecOrMat}
    @fastmath @inbounds @threads for i = 1: prod(alloc.dims)
        dμ_dθ = p.σ[i] * p.b[i] / θ[i]
        dμ_dv = p.σ[i] * p.a[i] / v[i]
        dθ[i] = dθ_dt(se, v[i], θ[i], p.L[i])
        dv[i] = dv_dt(alloc.dτ_dt[i], dμ_dv, dμ_dθ, dθ[i], p.η)
    end
end

@inline function dv_dθ_dt!(::RForm, se::StateEvolutionLaw,
    dv::T, dθ::T, v::T, θ::T, p::ElasticRSFProperties, alloc::OkadaGFAllocation
    ) where {T<:AbstractVecOrMat}
    @fastmath @inbounds @threads for i = 1: prod(alloc.dims)
        ψ1 = exp((p.f0 + p.b[i] * log(p.v0 * θ[i] / p.L[i])) / p.a[i]) / 2p.v0
        ψ2 = p.σ[i] * ψ1 / hypot(1, v[i] * ψ1)
        dμ_dv = p.a[i] * ψ2
        dμ_dθ = p.b[i] / θ[i] * v[i] * ψ2
        dθ[i] = dθ_dt(se, v[i], θ[i], p.L[i])
        dv[i] = dv_dt(alloc.dτ_dt[i], dμ_dv, dμ_dθ, dθ[i], p.η)
    end
end

@generated function (∂u∂t)(du::AbstractArray{T}, u::AbstractArray{T}, p::ElasticRSFProperties, alloc::OkadaGFAllocation{N}, gf::AbstractArray, flf::FrictionLawForm, se::StateEvolutionLaw
    ) where {T, N}
    quote
        v = selectdim(u, $(N+1), 1)
        θ = selectdim(u, $(N+1), 2)
        dv = selectdim(du, $(N+1), 1)
        dθ = selectdim(du, $(N+1), 2)
        clamp!(θ, zero(T), Inf)
        dτ_dt!(gf, alloc, p.vpl, v)
        dv_dθ_dt!(flf, se, dv, dθ, v, θ, p, alloc)
    end
end

function assemble(
    fs::BasicFaultSpace, p::ElasticRSFProperties, u0::AbstractArray, tspan::NTuple{2};
    flf::FrictionLawForm=RForm(), se::StateEvolutionLaw=DieterichStateLaw(), savefilename="gf.h5", savename="gf", kwargs...
    )
    gf = okada_disp_gf_tensor(fs.mesh, p.λ, p.μ, fs.ft; kwargs...)
    h5write(savefilename, savename, gf)
    return assemble(gf, fs, p, u0, tspan; flf=flf, se=se)
end

function assemble(
    gf::AbstractArray, fs::BasicFaultSpace, p::ElasticRSFProperties, u0::AbstractArray, tspan::NTuple{2};
    flf::FrictionLawForm=RForm(), se::StateEvolutionLaw=DieterichStateLaw(),
    )
    alloc = gen_alloc(fs.mesh)
    f! = (du, u, p, t) -> ∂u∂t(du, u, p, alloc, gf, flf, se)
    return ODEProblem(f!, u0, tspan, p)
end
