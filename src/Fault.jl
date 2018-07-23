module Fault

abstract type FaultType end
abstract type PlaneFaultType <: FaultType end
abstract type CurvedFaultType <: FaultType end

abstract type NormalFaultType <: PlaneFaultType end
abstract type ThrustFaultType <: PlaneFaultType end
abstract type StrikeSlipFaultType <: PlaneFaultType end
abstract type MixedFaultType <: PlaneFaultType end

abstract type DiscretizedMethod end
abstract type UniformlyDiscretized <: DiscretizedMethod end
abstract type MeshedDiscretized <: DiscretizedMethod end

abstract type FaultDomain end

struct GenericFaultDomain{FT} <: FaultDomain where {FT <: FaultType}
    faulttype::FT
end

struct PlaneFaultDomain{FT, DIM, T, A} <: FaultDomain where {
    FT <: PlaneFaultType, DIM <: Int, T <: Number, A <: NTuple{DIM, T}}
    faulttype::FT
    dim::DIM
    dip::T
    span::A
end

function create_fault(ft::Type{FT}) where {FT <: FaultType}
    GenericFaultDomain(ft)
end

function create_fault(ft::Type{FT}, dim, dip, span) where {FT <: PlaneFaultType}
    PlaneFaultDomain(ft, dim, len, wd, dip)
end

function discretize_fault(ft::Type{FT}) where {FT <: FaultType}
    info("Cannot descretize generic fault type.", prefix="Fault: ")
end

function discretize_fault(ft::PlaneFault, nx::Vararg)

end

end
