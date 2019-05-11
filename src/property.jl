## properties interface

export SingleDofRSFProperties, ElasticRSFProperties

import Base.fieldnames
import Base.==

abstract type AbstractProperties end

# https://github.com/jw3126/Setfield.jl
@with_kw mutable struct SingleDofRSFProperties{T<:Real} <: AbstractProperties
    a::T # contrib from velocity
    b::T # contrib from state
    L::T # critical distance
    k::T # stiffness tensor
    σ::T # effective normal stress
    η::T # radiation damping
    vpl::T # plate rate
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    @assert a > 0
    @assert b > 0
    @assert L > 0
    @assert k > 0
    @assert σ > 0
    @assert η ≥ 0
    @assert vpl > 0
    @assert f0 > 0
    @assert v0 > 0
end

@with_kw struct ElasticRSFProperties{T<:Real, U<:AbstractVecOrMat} <: AbstractProperties
    a::U # contrib from velocity
    b::U # contrib from state
    L::U # critical distance
    σ::U # effective normal stress
    λ::T # Lamé first constants
    μ::T # Lamé second constants
    η::T # radiation damping
    vpl::T # plate rate
    f0::T = 0.6 # ref. frictional coeff
    v0::T = 1e-6 # ref. velocity

    @assert size(a) == size(b)
    @assert size(b) == size(L)
    @assert size(L) == size(σ)
    @assert f0 > 0
    @assert v0 > 0
    @assert λ > 0
    @assert μ > 0
    @assert η > 0
    @assert vpl > 0
end

fieldnames(::Val{:SingleDofRSFProperties}) = ("a", "b", "L", "k", "σ", "η", "vpl", "f0", "v0")
fieldnames(::Val{:ElasticRSFProperties}) = ("a", "b", "L", "σ", "λ", "μ", "η", "vpl", "f0", "v0")

const __prop_names__ = (:SingleDofRSFProperties, :ElasticRSFProperties)

for nn in __prop_names__
    eval(
        quote
            fieldnames(p::$(nn)) = fieldnames(Val($(QuoteNode(nn))))
            description(p::$(nn)) = String($(QuoteNode(nn)))
        end
    )
end

function Base.:(==)(p1::P, p2::P) where P<:AbstractProperties
    reduce(&, [getfield(p1, Symbol(name)) == getfield(p2, Symbol(name)) for name in fieldnames(p1)])
end
