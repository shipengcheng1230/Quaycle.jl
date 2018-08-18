module RateStateFriction

using DifferentialEquations: ODEProblem
using Parameters

export StateEvolutionLaw, DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export FrictionLawForm, CForm, RForm
export dθ_dt, dμ_dt, dv_dθ_dt
export friction
export AbstractMaterialProfile, ZeroDimMaterialProfile, EarthquakeCycleProblem

abstract type StateEvolutionLaw end
struct DieterichStateLaw <: StateEvolutionLaw end
struct RuinaStateLaw <: StateEvolutionLaw end
struct PrzStateLaw <: StateEvolutionLaw end


dθ_dt(::DieterichStateLaw, v::T, θ::T, L::T) where {T<:Number} = 1 - v * θ / L

function dθ_dt(::RuinaStateLaw, v::T, θ::T, L::T) where {T<:Number}
    x = v * θ / L
    -x * log(clamp(x, zero(T), Inf))
end

dθ_dt(::PrzStateLaw, v::T, θ::T, L::T) where {T<:Number} = 1 - (v * θ / 2L) ^ 2

abstract type FrictionLawForm end
struct CForm <: FrictionLawForm end
struct RForm <: FrictionLawForm end

function friction(::CForm, v::T, θ::T, L::T, a::T, b::T, f0::T, v0::T) where {T<:Number}
    f0 + a * log(v / v0) + b * log(v0 * θ / L)
end

function friction(::RForm, v::T, θ::T, L::T, a::T, b::T, f0::T, v0::T) where {T<:Number}
    a * asinh(v / 2v0 * exp((f0 + b * log(v0 * θ / L)) / a))
end

dμ_dt(K::T, v::T, vpl::T) where {T<:Number} = K * (vpl - v)

dv_dt(dμdt::T, dμdv::T, dμdθ::T, dθdt::T, η::T) where {T<:Number} = (dμdt - dμdθ * dθdt) / (dμdv + η)

function dv_dθ_dt(::CForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, ::Vararg{T},
    ) where {T<:Number, SE<:StateEvolutionLaw}
    dμdθ = σ * b / θ
    dμdv = σ * a / v
    dθdt = dθ_dt(se, v, θ, L)
    dμdt = dμ_dt(k, v, vpl)
    return dv_dt(dμdt, dμdv, dμdθ, dθdt, η), dθdt
end

function dv_dθ_dt(::RForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T,
    ) where {T<:Number, SE<:StateEvolutionLaw}
    ψ1 = exp((f0 + b * log(v0 * θ / L)) / a) / 2v0
    ψ2 = σ * ψ1 / sqrt(1 + (v * ψ1)^2)
    dμ_dv = a * ψ2
    dμ_dθ = b / θ * v * ψ2
    dθdt = dθ_dt(se, v, θ, L)
    dμdt = dμ_dt(k, v, vpl)
    return dv_dt(dμdt, dμdv, dμdθ, dθdt, η), dθdt
end

abstract type AbstractPhysicalProperties{DIMS} end

@with_kw mutable struct PhysicalProperties{DIM} <: AbstractPhysicalProperties{DIM}
    a # contrib from velocity
    b # contrib from state
    L # critical distance
    k # stiffness tensor
    σ # effective normal stress
    η # radiation damping
    vpl # plate rate
    f0 # ref. frictional coeff
    v0 # ref. velocity

    function PhysicalProperties(a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T) where {
        T<:Union{<:Number, AbstractArray{<:Number}}}
        dims = maximum([ndims(x) for x in (a, b, L, σ, η, vpl, f0, v0)])
        new{dims}(a, b, L, k, σ, η, vpl, f0, v0)
    end
end

function dv_dθ_dt(p::PhysicalProperties{0}, v, θ, se, fform)
    dv_dθ_dt(fform, se, v, θ, p.a, p.b, p.L, p.k, p.σ, p.η, p.vpl, p.f0, p.v0)
end

friction(p::PhysicalProperties{0}, v, θ; fform=CForm()) = friction(fform, v, θ, p.L, p.a, p.b, p.f0, p.v0)
friction(p::PhysicalProperties{0}, vθ::AbstractArray{T,1}; kwargs...) where {T<:Number} = friction(p, vθ[1], vθ[2]; kwargs...)

function derivations!(du, u, p::PhysicalProperties{0}, t, se::StateEvolutionLaw, fform::FrictionLawForm)
    du[1], du[2] = dv_dθ_dt(p, u[1], u[2], se, fform)
end

function EarthquakeCycleProblem(p::PhysicalProperties{0}, u0, tspan; se=DieterichStateLaw(), fform=CForm())
    f! = (du, u, p, t) -> derivations!(du, u, p, t, se, fform)
    ODEProblem(f!, u0, tspan, p)
end

end # module
