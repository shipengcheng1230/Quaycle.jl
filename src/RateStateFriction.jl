module RateStateFriction

using DifferentialEquations
using Parameters

export StateEvolutionLaw,DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export FrictionLawForm, CForm, RForm
export dθ_dt, dμ_dt, dv_dθ_dt
export friction
export AbstractMaterialProfile, ZeroDimMaterialProfile, EarthquakeCycleProblem

abstract type StateEvolutionLaw end
abstract type DieterichStateLaw <: StateEvolutionLaw end
abstract type RuinaStateLaw <: StateEvolutionLaw end
abstract type PrzStateLaw <: StateEvolutionLaw end

abstract type FrictionLawForm end
abstract type CForm <: FrictionLawForm end
abstract type RForm <: FrictionLawForm end

function dθ_dt(::Type{DieterichStateLaw}, v::T, θ::T, L::T) where {T<:Number}
    1 - v * θ / L
end

function dθ_dt(::Type{RuinaStateLaw}, v::T, θ::T, L::T) where {T<:Number}
    x = v * θ / L
    -x * log(x)
end

function dθ_dt(::Type{PrzStateLaw}, v::T, θ::T, L::T) where {T<:Number}
    1 - (v * θ / 2L) ^ 2
end

function friction(::Type{CForm}, v::T, θ::T, L::T, a::T, b::T, f0=0.6, v0=1e-6) where {T<:Number}
    f0 + a * log(v / v0) + b * log(v0 * θ / L)
end

function friction(::Type{RForm}, v::T, θ::T, L::T, a::T, b::T, f0=0.6, v0=1e-6) where {T<:Number}
    a * asinh(v / 2v0 * exp((f0 + b * log(v0 * θ / L)) / a))
end

function dμ_dt(K::T, v::T, vpl::T) where {T<:Number}
    K * (vpl - v)
end

function dv_dt(dμdt::T, dμdv::T, dμdθ::T, dθdt::T, η::T) where {T<:Number}
    (dμdt - dμdθ * dθdt) / (dμdv + η)
end

function dv_dθ_dt(::Type{CForm}, ev::EV, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, ::Vararg{T},
    ) where {T<:Number, EV<:Type{<:StateEvolutionLaw}}
    dμdθ = σ * b / θ
    dμdv = σ * a / v
    dθdt = dθ_dt(ev, v, θ, L)
    dμdt = dμ_dt(k, v, vpl)
    return dv_dt(dμdt, dμdv, dμdθ, dθdt, η), dθdt
end

function dv_dθ_dt(::Type{RForm}, ev::EV, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T,
    ) where {T<:Number, EV<:Type{<:StateEvolutionLaw}}
    ψ1 = exp((f0 + b * log(v0 * θ / L)) / a) / 2v0
    ψ2 = σ * ψ1 / sqrt(1 + (v * ψ1)^2)
    dμ_dv = a * ψ2
    dμ_dθ = b / θ * v * ψ2
    dθdt = dθ_dt(ev, v, θ, L)
    dμdt = dμ_dt(k, v, vpl)
    return dv_dt(dμdt, dμdv, dμdθ, dθdt, η), dθdt
end

abstract type AbstractMaterialProfile{DIMS} end

@with_kw mutable struct ZeroDimMaterialProfile{EV, F, T} <: AbstractMaterialProfile{0} where {
    EV<:Type{<:StateEvolutionLaw}, F<:Type{<:FrictionLawForm}, T<:Number}
    ev::EV
    fform::F
    a::T
    b::T
    L::T
    k::T
    σ::T
    η::T
    vpl::T
    f0::T
    v0::T
    function ZeroDimMaterialProfile(ev::Type{EV}, fform::Type{F}, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T) where {
        EV<:StateEvolutionLaw, F<:FrictionLawForm, T<:Number}
        new{Type{EV}, Type{F}, T}(ev, fform, a, b, L, k, σ, η, vpl, f0, v0)
    end
end

dv_dθ_dt(p::ZeroDimMaterialProfile, v, θ) = dv_dθ_dt(p.fform, p.ev, v, θ, p.a, p.b, p.L, p.k, p.σ, p.η, p.vpl, p.f0, p.v0)
friction(p::ZeroDimMaterialProfile, v, θ) = friction(p.fform, v, θ, p.L, p.a, p.b, p.f0, p.v0)
friction(p::ZeroDimMaterialProfile, vθ::AbstractArray{T,1}) where {T<:Number} = friction(p, vθ[1], vθ[2])

function derivations!(du, u, p::ZeroDimMaterialProfile, t)
    du[1], du[2] = dv_dθ_dt(p, u[1], u[2])
end

function EarthquakeCycleProblem(p::ZeroDimMaterialProfile, u0, tspan)
    ODEProblem(derivations!, u0, tspan, p)
end

end # module
