using DifferentialEquations: ODEProblem
using Parameters

export StateEvolutionLaw, DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export CForm, RForm
export friction
export MaterialProperties, EarthquakeCycleProblem

abstract type StateEvolutionLaw end

"`\\frac{\\mathrm{d}θ}{\\mathrm{d}t} = 1 - \\frac{v θ}{L}``"
struct DieterichStateLaw <: StateEvolutionLaw end

"``\\frac{\\mathrm{d}θ}{\\mathrm{d}t} = -\\frac{v θ}{L} * \\log{\\frac{v θ}{L}}``"
struct RuinaStateLaw <: StateEvolutionLaw end

"``\\frac{\\mathrm{d}θ}{\\mathrm{d}t} = 1 - (\\frac{v θ}{2L})^2``"
struct PrzStateLaw <: StateEvolutionLaw end

dθ_dt(::DieterichStateLaw, v::T, θ::T, L::T) where {T<:Number} = 1 - v * θ / L

function dθ_dt(::RuinaStateLaw, v::T, θ::T, L::T) where {T<:Number}
    x = v * θ / L
    -x * log(clamp(x, zero(T), Inf))
end

dθ_dt(::PrzStateLaw, v::T, θ::T, L::T) where {T<:Number} = 1 - (v * θ / 2L) ^ 2

abstract type FrictionLawForm end
struct CForm <: FrictionLawForm end # conventinal form
struct RForm <: FrictionLawForm end # regularized form

"""
    friction(::FrictionLawForm, v::T, θ::T, L::T, a::T, b::T, f0::T, v0::T) where {T<:Number}

Calculate friction given by the form of fomula as well as other necessary parameters.
"""
function friction(::CForm, v::T, θ::T, L::T, a::T, b::T, f0::T, v0::T) where {T<:Number}
    f0 + a * log(v / v0) + b * log(v0 * θ / L)
end

function friction(::RForm, v::T, θ::T, L::T, a::T, b::T, f0::T, v0::T) where {T<:Number}
    a * asinh(v / 2v0 * exp((f0 + b * log(v0 * θ / L)) / a))
end

dτ_dt(K::T, v::T, vpl::T) where {T<:Number} = K * (vpl - v)

dv_dt(dτdt::T, dμdv::T, dμdθ::T, dθdt::T, η::T) where {T<:Number} = (dτdt - dμdθ * dθdt) / (dμdv + η)

function dv_dθ_dt(::CForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, ::Vararg{T},
    ) where {T<:Number, SE<:StateEvolutionLaw}
    dμdθ = σ * b / θ
    dμdv = σ * a / v
    dθdt = dθ_dt(se, v, θ, L)
    dτdt = dτ_dt(k, v, vpl)
    return dv_dt(dτdt, dμdv, dμdθ, dθdt, η), dθdt
end

function dv_dθ_dt(::RForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T,
    ) where {T<:Number, SE<:StateEvolutionLaw}
    ψ1 = exp((f0 + b * log(v0 * clamp(θ, zero(T), Inf) / L)) / a) / 2v0
    ψ2 = σ * ψ1 / sqrt(1 + (v * ψ1)^2)
    dμdv = a * ψ2
    dμdθ = b / θ * v * ψ2
    dθdt = dθ_dt(se, v, θ, L)
    dτdt = dτ_dt(k, v, vpl)
    return dv_dt(dτdt, dμdv, dμdθ, dθdt, η), dθdt
end

abstract type AbstractMaterialProperties{DIMS} end

@with_kw mutable struct MaterialProperties{dim, T<:Number} <: AbstractMaterialProperties{0}
    a::T # contrib from velocity
    b::T # contrib from state
    L::T # critical distance
    k::T # stiffness tensor
    σ::T # effective normal stress
    η::T # radiation damping
    vpl::T # plate rate
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    function MaterialProperties(a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T) where {T<:Number}
        new{0, T}(a, b, L, k, σ, η, vpl, f0, v0)
    end
end

function dv_dθ_dt(p::MaterialProperties{0}, v, θ, se, fform)
    dv_dθ_dt(fform, se, v, θ, p.a, p.b, p.L, p.k, p.σ, p.η, p.vpl, p.f0, p.v0)
end

friction(p::MaterialProperties{0}, v, θ; fform=CForm()) = friction(fform, v, θ, p.L, p.a, p.b, p.f0, p.v0)

friction(p::MaterialProperties{0}, vθ::AbstractVecOrMat{T}; kwargs...) where {T<:Number} = friction(p, vθ[1], vθ[2]; kwargs...)

friction(p::MaterialProperties{0}, vθ::AbstractArray{T}) where {T<:AbstractVecOrMat} = friction.(Ref(p), vθ)

function derivations!(du, u, p::MaterialProperties{0}, t, se::StateEvolutionLaw, fform::FrictionLawForm)
    du[1], du[2] = dv_dθ_dt(p, u[1], u[2], se, fform)
end

function EarthquakeCycleProblem(p::MaterialProperties{0}, u0, tspan; se=DieterichStateLaw(), fform=CForm())
    (fform == RForm() && p.η ≈ 0.0) && @warn "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt."
    f! = (du, u, p, t) -> derivations!(du, u, p, t, se, fform)
    ODEProblem(f!, u0, tspan, p)
end
