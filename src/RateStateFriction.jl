module RateStateFriction

using DifferentialEquations

export StateEvolutionLaw,DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export FrictionLawForm, CForm, RForm
export dθ_dt, dμ_dt, dv_dθ_dt
export friction
export AbstractMaterialProfile, OneDofMP

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

function dv_dθ_dt(::Type{CForm}, v::T, θ::T, a::T, b::T, L::T, σ::T, η::T, dμdt::T, ev::EV,
    ) where {T<:Number, EV<:Type{<:StateEvolutionLaw}}
    dμdθ = σ * b / θ
    dμdv = σ * a / v
    dθdt = dθ_dt(ev, v, θ, L)
    return dv_dt(dμdt, dμdv, dμdθ, dθdt, η), dθdt
end

function dv_dθ_dt(::Type{RForm}, v::T, θ::T, a::T, b::T, L::T, σ::T, η::T, dμdt::T, ev::EV, f0=0.6, v0=1e-6,
    ) where {T<:Number, EV<:Type{<:StateEvolutionLaw}}
    ψ1 = exp((f0 + b * log(v0 * θ / L)) / a) / 2v0
    ψ2 = σ * ψ1 / sqrt(1 + (v * ψ1)^2)
    dμ_dv = a * ψ2
    dμ_dθ = b / θ * v * ψ2
    dθdt = dθ_dt(ev, v, θ, L)
    return dv_dt(dμdt, dμdv, dμdθ, dθdt, η), dθdt
end

abstract type AbstractMaterialProfile{DIMS} end

mutable struct OneDofMP{EV, T} <: AbstractMaterialProfile{1} where {EV<:Type{<:StateEvolutionLaw}, T<:Number}
    ev::EV
    a::T
    b::T
    L::T
    k::T
    vpl::T
    f0::T
    v0::T
    function OneDofMP(ev::Type{EV}, a::T, b::T, L::T, k::T, vpl::T, f0::T, v0::T) where {EV<:Type{<:StateEvolutionLaw}, T<:Number}
        new{EV, T}(ev, a, b, L, k, vpl, f0, v0)
    end
end

function ODE_vθ!(du, u, p, t)

end

function simulate(rsf::OneDofMP; kwargs...)

end

end # module
