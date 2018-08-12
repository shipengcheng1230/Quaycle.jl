module Fault

abstract type AbstractFaultType end
abstract type PlaneFaultType <: AbstractFaultType end
abstract type CurvedFaultType <: AbstractFaultType end

abstract type NormalFaultType <: PlaneFaultType end
abstract type ThrustFaultType <: PlaneFaultType end
abstract type StrikeSlipFaultType <: PlaneFaultType end
abstract type MixedAFaultType <: PlaneFaultType end

abstract type AbstractDiscretizedMethod end
abstract type FDM <: AbstractDiscretizedMethod end
abstract type FEM <: AbstractDiscretizedMethod end
abstract type BEM <: AbstractDiscretizedMethod end

abstract type AbstractFaultDomain end

struct GenericFaultDomain{FT} <: AbstractFaultDomain where {FT <: AbstractFaultType}
    AbstractFaultType::FT
end

struct PlaneFaultDomain{FT, DIM, T, A} <: GenericFaultDomain where {
    FT <: PlaneAbstractFaultType, T <: Number, A <: AbstractArray{T}}
    faulttype::FT
    dip::T
    span::A
end

struct CurvedFaultDomain{FT} <: GenericFaultDomain where {FT <: CurvedFaultType}
    faulttype::FT
end


end
