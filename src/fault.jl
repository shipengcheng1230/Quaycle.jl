import Base.getindex

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

        @assert(zero(T) ≤ dip ≤ 90 * one(T), "Fault dip angle: $dip received, must ∈ [0, 90].")
        @assert(all(@. span ≥ zero(T)), "Fault domain: $span received, must > 0.")
        @assert((N == 1 || N == 2), "Fault domain dim: $N received, must be `1` (along-downdip) or `2` (plus along-strike).")
        (FT == StrikeSlipFault && dip ≉ 90 * one(T)) && error("Dip angle: $dip received, for strike-slip faults must be `90`.")
        new{FT, N, T}(dip, span)
    end
end

"""
    fault(ftype::Type{<:PlaneFault}, dip, span)
Generate a fault given the fault type, dip angle and its spatial span.
## Arguments
- `ftype::Type{<:PlaneFault}`: type of plane fault
- `dip`: dip angle in degree
- `span`: spatial span of fault size
"""
fault(ftype::Type{<:PlaneFault}, dip::T, span::NTuple{N, T}) where {N, T} = PlaneFaultDomain(ftype, dip, span)
fault(ftype::Type{<:PlaneFault}, dip, span::Number) = fault(ftype, dip, tuple(span))
fault(ftype::Type{<:PlaneFault}, dip, span::AbstractVector{<:Number}) = fault(ftype, dip, tuple(span...))
fault(ftype::Type{<:PlaneFault}, dip::T, span::NTuple{N, T}) where {N, T<:Integer} = fault(ftype, Float64(dip), span)

function fault(ftype::Type{<:PlaneFault}, dip::T1, span::NTuple{N, T2}) where {N, T1, T2}
    pt = promote_type(T1, T2)
    fault(ftype, pt(dip), convert(NTuple{N, pt}, span))
end

fault(ftype::Type{StrikeSlipFault}, span) = fault(ftype, one(eltype(span)) * 90, span)
fault(ftype::Type{StrikeSlipFault}, span::Number) = fault(ftype, tuple(span))
fault(ftype::Type{StrikeSlipFault}, span::AbstractVector{<:Number}) = fault(ftype, tuple(span...))

getindex(fd::PlaneFaultDomain{ftype, 1}, name::Symbol) where {ftype} = name == :ξ ? fd.span[1] : error(
    "One dim fault domain only contains index of `ξ`, received `$name`.")
getindex(fd::PlaneFaultDomain{ftype, 2}, name::Symbol) where {ftype} = name == :ξ ? fd.span[2] : name == :x ? fd.span[1] : error(
    "Two dim fault domain only contains index of `x` and `ξ`., received `$name`.")
