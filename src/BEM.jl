module BEM

include(joinpath(@__DIR__, "Fault.jl"))
include(joinpath(@__DIR__, "RateStateFriction.jl"))
include(joinpath(@__DIR__, "dc3d.jl"))

using SimpleTraits
using Parameters
using ToeplitzMatrices
using TensorOperations

using .RateStateFriction
using .Fault

export AbstractBEMGrids, BEMGrid

"""
    Symmetric properties in stiffness tensor w.r.t. different fault types.
"""
@traitdef IsTranslationalSymmetric{ftype}
@traitdef IsReflectiveSymmetric{ftype}
@traitdef IsReflectiveASymmetric{ftype}

@traitimpl IsTranslationalSymmetric{NormalFault}
@traitimpl IsTranslationalSymmetric{ThrustFault}
@traitimpl IsTranslationalSymmetric{StrikeSlipFault}

@traitimpl IsReflectiveSymmetric{NormalFault}
@traitimpl IsReflectiveSymmetric{ThrustFault}

@traitimpl IsReflectiveASymmetric{StrikeSlipFault}

abstract type AbstractBEMGrids{dim, isuniform} end

struct BEMGrid{N, U<:AbstractVector} <: AbstractBEMGrids{N, true}
    x::U
    ξ::U
    Δx
    Δξ
    nx
    nξ

    function BEMGrid(x::U, ξ::U, Δx, Δξ, nx, nξ) where {U<:AbstractVector}
        N = isempty(x) ? 1 : 2
        N == 2 ? new{N, U}(x, ξ, Δx, Δξ, nx, nξ) : new{N, U}(x, ξ, nothing, Δξ, nothing, nξ)
    end
end

BEMGrid(ξ, Δξ, nξ) = BEMGrid([], ξ, nothing, Δξ, nothing, nξ)

function discretize(fa::PlaneFaultDomain{ftype, 1}, Δξ) where {ftype}
    ξ, nξ = _divide_segment(Val(:halfspace), fa[:ξ], Δξ)
    BEMGrid(ξ, Δξ, nξ)
end

function discretize(fa::PlaneFaultDomain{ftype, 2}, Δx, Δξ) where {ftype}
    ξ, nξ = _divide_segment(Val(:halfspace), fa[:ξ], Δξ)
    x, nx = _divide_segment(Val(:symmatzero), fa[:x], Δx)
    BEMGrid(x, ξ, Δx, Δξ, nx, nξ)
end

function _divide_segment(::Val{:halfspace}, x::T, Δx::T) where {T<:Number}
    xi = collect(range(-x+Δx, stop=zero(T), step=Δx)) .- Δx/2
    return xi, length(xi)
end

function _divide_segment(::Val{:symmatzero}, x::T, Δx::T) where {T<:Number}
    xi = collect(range(-x/2 + Δx/2, stop=x/2 - Δx/2, step=Δx))
    return xi, length(xi)
end

function stiffness_tensor()

end

end # module
