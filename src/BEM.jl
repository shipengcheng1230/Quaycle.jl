module BEM

include(joinpath(@__DIR__, "Fault.jl"))
include(joinpath(@__DIR__, "RateStateFriction.jl"))
include(joinpath(@__DIR__, "dc3d.jl"))

using SimpleTraits
using Parameters
using ToeplitzMatrices
using TensorOperations
using Threads

using .RateStateFriction
using .Fault

export AbstractBEMGrids, BEMGrid

"Symmetric properties in stiffness tensor w.r.t. different fault types."
@traitdef IsTranslationalSymmetrical{ftype}
@traitdef IsReflectiveSymmetrical{ftype}
@traitdef IsReflectiveASymmetrical{ftype}

@traitimpl IsTranslationalSymmetrical{NormalFault}
@traitimpl IsTranslationalSymmetrical{ThrustFault}
@traitimpl IsTranslationalSymmetrical{StrikeSlipFault}

@traitimpl IsReflectiveSymmetrical{NormalFault}
@traitimpl IsReflectiveSymmetrical{ThrustFault}

@traitimpl IsReflectiveASymmetrical{StrikeSlipFault}

abstract type AbstractBEMGrids{dim, isuniform} end

"Along-strike will be places on x-axis while along-downdip on yz-plane, w.r.t Okada's dc3d coordinates."
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

"""
    shear_stress(::Type{<:PlaneFault}, u12, λ, μ, dip)

Calculate the shear traction on the fault plane w.r.t. fault types.

# Arguments
- `u12::AbstractArray{<:Number, 1}`: the output from dc3d_okada
- `λ::Number`: Lamé's first parameter
- `μ::Number`: shear modulus
- `dip::Number`: plane dip angle

# Reference
- [plane stress](https://nnakata.oucreate.com/page/Teaching_files/GEOPHYS130/GEOPHYS130_notes_all.pdf)

"""
@inline @views @inbounds @traitfn function shear_stress(::Type{FT}, u12, λ, μ, dip) where {FT<:PlaneFault; IsReflectiveSymmetrical{FT}}
    σzz = (λ + 2μ) * u[12] + λ * u[4] + λ * u[8]
    σyy = (λ + 2μ) * u[8] + λ * u[4] + λ * u[12]
    τyz = μ * (u[11] + u[9])
    -((σzz - σyy)/2 * sinpi(2dip/180) + τyz * cospi(2dip/180))
end

@inline @views @inbounds @traitfn function shear_stress(::Type{FT}, u12, λ, μ, dip) where {FT<:PlaneFault; IsReflectiveASymmetrical{FT}}
    μ * (u[5] + u[7])
end

function stiffness_element()

end


end # module
