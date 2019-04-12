## rate-and-state friction

export CForm, RForm
export DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export SingleDegreeSystem, assemble

abstract type StateEvolutionLaw end

"``\\frac{\\mathrm{d}θ}{\\mathrm{d}t} = 1 - \\frac{V θ}{L}``"
struct DieterichStateLaw <: StateEvolutionLaw end

"``\\frac{\\mathrm{d}θ}{\\mathrm{d}t} = -\\frac{V θ}{L} * \\log{\\frac{V θ}{L}}``"
struct RuinaStateLaw <: StateEvolutionLaw end

"``\\frac{\\mathrm{d}θ}{\\mathrm{d}t} = 1 - (\\frac{V θ}{2L})^2``"
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
- Conventional Form:
```math
f(V, θ) = f_0 + a \\ln{\\frac{V}{V_0}} + b \\ln{\\left(\\frac{V_0 θ}{L}\\right)}
```
- Regularized Form:
```math
f(V, θ) = a \\sinh ^{-1}{\\left[\\frac{V}{2V_0} \\exp{\\frac{f_0 + b \\ln{\\left(V_0 θ/L\\right)}}{a}}\\right]}
```
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

@with_kw mutable struct SingleDegreeSystem{T<:Number}
    a::T # contrib from velocity
    b::T # contrib from state
    L::T # critical distance
    k::T # stiffness tensor
    σ::T # effective normal stress
    η::T # radiation damping
    vpl::T # plate rate
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity
end

function dv_dθ_dt(p::SingleDegreeSystem, v, θ, se, fform)
    dv_dθ_dt(fform, se, v, θ, p.a, p.b, p.L, p.k, p.σ, p.η, p.vpl, p.f0, p.v0)
end

friction(p::SingleDegreeSystem, v, θ; fform=CForm()) = friction(fform, v, θ, p.L, p.a, p.b, p.f0, p.v0)

friction(p::SingleDegreeSystem, vθ::AbstractVecOrMat{T}; kwargs...) where {T<:Number} = friction(p, vθ[1], vθ[2]; kwargs...)

friction(p::SingleDegreeSystem, vθ::AbstractArray{T}; kwargs...) where {T<:AbstractVecOrMat} = friction.(Ref(p), vθ; kwargs...)

function derivations!(du, u, p::SingleDegreeSystem, t, se::StateEvolutionLaw, fform::FrictionLawForm)
    du[1], du[2] = dv_dθ_dt(p, u[1], u[2], se, fform)
end

function assemble(p::SingleDegreeSystem, u0, tspan; se=DieterichStateLaw(), fform=CForm())
    (fform == RForm() && p.η ≈ 0.0) && @warn "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt."
    f! = (du, u, p, t) -> derivations!(du, u, p, t, se, fform)
    ODEProblem(f!, u0, tspan, p)
end
