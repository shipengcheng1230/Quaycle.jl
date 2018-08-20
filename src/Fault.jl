module Fault

using Base: getindex

export AbstractFault, PlaneFault, CurvedFault
export NormalFault, ThrustFault, StrikeSlipFault
export AbstractFaultDomain, PlaneFaultDomain
export fault

abstract type AbstractFault end
abstract type PlaneFault <: AbstractFault end
abstract type CurvedFault <: AbstractFault end

abstract type NormalFault <: PlaneFault end
abstract type ThrustFault <: PlaneFault end
abstract type StrikeSlipFault <: PlaneFault end

abstract type AbstractFaultDomain{ftype, dim} end

struct PlaneFaultDomain{ftype, dim, T} <: AbstractFaultDomain{ftype, dim}
    dip::T
    span::NTuple{dim, T}

    function PlaneFaultDomain(ftype::Type{FT}, dip::T, span::U) where
        {FT<:PlaneFault, T<:Number, U<:NTuple{N, T}} where {N}
        zero(T) ≤ dip ≤ 90 * one(T) || error("Fautl dip angle ∈ [0, π/2].")
        any(@. span ≤ zero(T)) && error("Fault domain must be larger than zero.")
        (N == 1 || N == 2) || error("Fault domain should be given by along-downdip [and along-strike] length[s].")
        new{FT, N, T}(dip, span)
    end
end

fault(ftype::Type{<:PlaneFault}, dip, span) = PlaneFaultDomain(ftype, dip, span)
fault(ftype::Type{<:PlaneFault}, dip, span::Number) = PlaneFaultDomain(ftype, dip, tuple(span))
fault(ftype::Type{<:PlaneFault}, dip, span::AbstractVector{<:Number}) = PlaneFaultDomain(ftype, dip, tuple(span...))

getindex(fd::PlaneFaultDomain{ftype, 1}, name::Symbol) where {ftype} = name == :ξ ? fd.span[1] : error(
    "One dim fault domain only contains index of `ξ`.")
getindex(fd::PlaneFaultDomain{ftype, 2}, name::Symbol) where {ftype} = name == :ξ ? fd.span[2] : name == :x ? fd.span[1] : error(
    "Two dim fault domain only ontains index of `x` and `ξ`.")

end # module
