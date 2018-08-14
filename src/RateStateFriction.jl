module RateStateFriction

export StateEvolutionLaw
export DieterichState, RuinaState, PrzState, NagataState
export CompositeState
export dθ_dt, dμ_dt, dμ_dv, dμ_dθ
export friction

abstract type StateEvolutionLaw end
abstract type DieterichState <: StateEvolutionLaw end
abstract type RuinaState <: StateEvolutionLaw end
abstract type PrzState <: StateEvolutionLaw end
abstract type NagataState <: StateEvolutionLaw end

const CompositeState = AbstractArray{Type{StateEvolutionLaw},N} where {N}

function dθ_dt(::Type{DieterichState}, v::T, θ::T, L::T, _...) where {T<:Number}
    1. - v * θ / L
end

function dθ_dt(::Type{RuinaState}, v::T, θ::T, L::T, _...) where {T<:Number}
    x = v * θ / L
    -x * log(x)
end

function dθ_dt(::Type{PrzState}, v::T, θ::T, L::T, _...) where {T<:Number}
    1. - (v * θ / 2L) ^ 2
end

function dθ_dt(::Type{NagataState}, v::T, θ::T, L::T, c::T, b::T, dμdt::T) where {T<:Number}
    1. - (v * θ / L) - (c / b * θ * dμdt)
end

function dθ_dt(states::CompositeState, v::T, θ::T, L::T, c::T, b::T, dμdt::T) where {T<:Number}
    r = zeros(states)
    for i, s in enumerate(states)
        r[i] = dθ_dt(s, v, θ, L, c, b, dμdt)
    end
    r    
end

function friction(::Val{:conventional}, v, θ, L, a, b, f0=0.6, v0=1e-6)
    f0 + a * log(v / v0) + b * log(v0 * θ / L)
end

function friction(::Val{:regularized}, v, θ, L, a, b, f0=0.6, v0=1e-6)
    a * asinh(v / 2v0 * exp((f0 + b * log(v0 * θ / L)) / a))
end

function dμ_dt(K::T, v::T, vpl::T) where {T<:Number}
    K * (vpl - v)
end

function dv_dt(dμdt::T, dμdv::T, dμdθ::T, dθdt::T, η::T) where {T<:Number}
    (dμdt - dμdθ * dθdt) / (dμdv + η)
end

function dμ_dv(v::T, a::T) where {T<:Number}
    a / v
end

function dμ_dθ(θ::T, b::T) where {T<:Number}
    b / θ
end


end
