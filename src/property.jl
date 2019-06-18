## Property interface

export SingleDofRSFProperty, ElasticRSFProperty, DislocationCreepProperty, DiffusionCreepProperty,
    CompositePlasticDeformationProperty, ViscoelasticMaxwellProperty

import Base.fieldnames
import Base.==

abstract type AbstractProperty end
abstract type PlasticDeformationProperty <: AbstractProperty end

# https://github.com/jw3126/Setfield.jl
"""
System property for single degree of freedom under rate-state friction.

- `a`: contrib from velocity
- `b`: contrib from state
- `L`: critical distance
- `k`: spring stiffness
- `σ`: effective normal stress
- `η`: radiation damping
- `vpl`: plate rate
- `f0` = 0.6: ref. frictional coeff
- `v0` = 1e-6: ref. velocity
"""
@with_kw mutable struct SingleDofRSFProperty{T<:Real} <: AbstractProperty
    a::T # contrib from velocity
    b::T # contrib from state
    L::T # critical distance
    k::T # spring stiffness
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

"""
System property for multiple fault patches under rate-state friction.

- `a`: contrib from velocity
- `b`: contrib from state
- `L`: critical distance
- `σ`: effective normal stress
- `λ`: Lamé first constants
- `μ`: Lamé second constants
- `η`: radiation damping
- `vpl`: plate rate
- `f0` = 0.6: ref. frictional coeff
- `v0` = 1e-6: ref. velocity
"""
@with_kw struct ElasticRSFProperty{T<:Real, U<:AbstractVecOrMat} <: AbstractProperty
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

"""
Compose all three type of plastic deformation, see
    [(Kohlstedt & Hansen, 2015)](https://www.sciencedirect.com/science/article/pii/B9780444538024000427).
    Each field is the overall equivalent factor not dependent on stress.
"""
@with_kw struct CompositePlasticDeformationProperty{U, I, V} <: PlasticDeformationProperty
    disl::U
    n::I
    diff::U
    peie::U
    ϵref::U
    ϵind::V

    @assert size(disl) == size(n)
    @assert size(n) == size(diff)
    @assert size(diff) == size(peie)
    @assert size(ϵref) == size(ϵind)
end

"""
System properties for plastic deformation of dislocation creep.
"""
@with_kw struct DislocationCreepProperty{V<:AbstractVector} <: PlasticDeformationProperty
    A::V # prefactor
    n::V # power law stress exponent
    fH₂0::V # water content
    r::V # water fugacity exponent
    α::V # melting constant
    ϕ::V # melting fraction
    Q::V # activation energy
    P::V # pressure
    Ω::V # activation volume
    T::V # temperature
end

"""
System properties for plastic deformation of diffusion creep.
"""
@with_kw struct DiffusionCreepProperty{V<:AbstractVector} <: PlasticDeformationProperty
    A::V # prefactor
    d::V # grain size
    m::V # grain size exponent
    fH₂0::V # water content
    r::V # water fugacity exponent
    α::V # melting constant
    ϕ::V # melting fraction
    Q::V # activation energy
    P::V # pressure
    Ω::V # activation volume
    T::V # temperature
end

"""
Composite property for viscoelastic rheology of maxwell representation.
"""
struct ViscoelasticMaxwellProperty{T1, T2} <: AbstractProperty
    pe::T1
    pv::T2

    function ViscoelasticMaxwellProperty(pe::ElasticRSFProperty, pv::CompositePlasticDeformationProperty)
        new{typeof(pe), typeof(pv)}(pe, pv)
    end
end

composite_factor(pv::DislocationCreepProperty) = @. pv.A * pv.fH₂0^(pv.r) * exp(pv.α * pv.ϕ) * exp(-(pv.Q + pv.P * pv.Ω) / 𝙍 / pv.T)
composite_factor(pv::DiffusionCreepProperty) = @. pv.A * pv.d^(-pv.m) * pv.fH₂0^(pv.r) * exp(pv.α * pv.ϕ) * exp(-(pv.Q + pv.P * pv.Ω) / 𝙍 / pv.T)

function ViscoelasticMaxwellProperty(pe::ElasticRSFProperty{T}, ϵref, ϵind, pvs...) where T
    @assert length(pvs) ≤ 3 "Received more than 3 types of plastic deformation mechanisms."
    disl, diff, peie, n = [zeros(T, size(pvs[1].A)) for _ in 1: 4]
    for pv in pvs
        if isa(pv, DislocationCreepProperty)
            disl .= composite_factor(pv)
            n .= pv.n
        end
        if isa(pv, DiffusionCreepProperty)
            diff .= composite_factor(pv)
        end
    end
    ViscoelasticMaxwellProperty(pe, CompositePlasticDeformationProperty(disl, n, diff, peie, ϵref, ϵind))
end

const prop_field_names = Dict(
    :SingleDofRSFProperty => ("a", "b", "L", "k", "σ", "η", "vpl", "f0", "v0"),
    :ElasticRSFProperty => ("a", "b", "L", "σ", "λ", "μ", "η", "vpl", "f0", "v0"),
    :DislocationCreepProperty => ("A", "n", "fH₂0", "r", "α", "ϕ", "Q", "P", "Ω", "T"),
    :DiffusionCreepProperty => ("A", "d", "m", "fH₂0", "r", "α", "ϕ", "Q", "P", "Ω", "T"),
    :ViscoelasticMaxwellProperty => ("pe", "pv"),
    :CompositePlasticDeformationProperty => ("disl", "n", "diff", "peie", "ϵref", "ϵind"),
    )

for (nn, fn) in prop_field_names
    @eval begin
        fieldnames(p::$(nn)) = $(fn)
        description(p::$(nn)) = String($(QuoteNode(nn)))
    end
end

function Base.:(==)(p1::P, p2::P) where P<:AbstractProperty
    reduce(&, [getfield(p1, Symbol(name)) == getfield(p2, Symbol(name)) for name in fieldnames(p1)])
end

## shortcut function
friction(flf::FrictionLawForm, v::T, θ::T, p::SingleDofRSFProperty) where T = friction(flf, v, θ, p.a, p.b, p.L, p.f0, p.v0)
friction(flf::FrictionLawForm, u::AbstractVecOrMat{T}, p::SingleDofRSFProperty) where T<:Real = friction(flf, u[1], u[2], p)
friction(flf::FrictionLawForm, u::AbstractArray{T}, p::SingleDofRSFProperty) where T<:AbstractVecOrMat = friction.(Ref(flf), u, Ref(p))
