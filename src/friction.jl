## rate-and-state friction
export CForm, RForm
export StateEvolutionLaw, DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export friction

"""
Currently support:
- [`DieterichStateLaw`](@ref)
- [`RuinaStateLaw`](@ref)
- [`PrzStateLaw`](@ref)
"""
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
"Conventional form, see [`friction`](@ref)"
struct CForm <: FrictionLawForm end # conventinal form
"Regularized form, see [`friction`](@ref)"
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
function friction(::CForm, v::T, θ::T, a::T, b::T, L::T, f0::T, v0::T) where T
    f0 + a * log(v / v0) + b * log(v0 * θ / L)
end

function friction(::RForm, v::T, θ::T, a::T, b::T, L::T, f0::T, v0::T) where T
    a * asinh(v / 2v0 * exp((f0 + b * log(v0 * θ / L)) / a))
end
