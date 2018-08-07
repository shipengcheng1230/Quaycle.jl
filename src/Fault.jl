module Fault

abstract type AbstractFaultType end
abstract type PlaneAbstractFaultType <: AbstractFaultType end
abstract type CurvedAbstractFaultType <: AbstractFaultType end

abstract type NormalAbstractFaultType <: PlaneAbstractFaultType end
abstract type ThrustAbstractFaultType <: PlaneAbstractFaultType end
abstract type StrikeSlipAbstractFaultType <: PlaneAbstractFaultType end
abstract type MixedAbstractFaultType <: PlaneAbstractFaultType end

abstract type DiscretizedMethod end
abstract type UniformlyDiscretized <: DiscretizedMethod end
abstract type MeshedDiscretized <: DiscretizedMethod end

abstract type AbstractFaultDomain end

struct GenericFaultDomain{FT} <: AbstractFaultDomain where {FT <: AbstractFaultType}
    AbstractFaultType::FT
end

struct PlaneFaultDomain{FT, DIM, T, A} <: GenericFaultDomain where {
    FT <: PlaneAbstractFaultType, T <: Number, A <: NTuple{DIM, T}}
    AbstractFaultType::FT
    dip::T
    span::A
    dim::Int
end

function create_fault(ft::Type{FT}) where {FT <: AbstractFaultType}
    GenericFaultDomain(ft)
end

function create_fault(ft::Type{FT}, dim, dip, span) where {FT <: PlaneAbstractFaultType}
    PlaneFaultDomain(ft, dim, len, wd, dip)
end

function discretize_fault(ft::Type{FT}) where {FT <: AbstractFaultType}
    info("Cannot descretize generic fault type.", prefix="Fault: ")
end

function discretize_fault(ft::PlaneFault, nx::Vararg)

end

end
