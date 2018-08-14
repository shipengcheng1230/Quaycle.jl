module Fault

export AbstractFaultType, PlaneFaultType, CurvedFaultType
export NormalFaultType, ThrustFaultType, StrikeSlipFaultType
export NormalFault, ThrustFault, StrikeSlipFault
export AbstractFaultDomain, PlaneFaultDomain
export AbstractGridType, FiniteDifferenceGrid, UniformFiniteDifferenceGrid
export fault, discretize

abstract type AbstractFaultType end
abstract type PlaneFaultType <: AbstractFaultType end
abstract type CurvedFaultType <: AbstractFaultType end

abstract type NormalFaultType <: PlaneFaultType end
abstract type ThrustFaultType <: PlaneFaultType end
abstract type StrikeSlipFaultType <: PlaneFaultType end

struct NormalFault <: NormalFaultType end
struct ThrustFault <: ThrustFaultType end
struct StrikeSlipFault <: StrikeSlipFaultType end

abstract type AbstractFaultDomain end

struct PlaneFaultDomain{FT,T,N} <: AbstractFaultDomain where {FT<:PlaneFaultType,T<:Number,N}
    fault_type::FT
    domain::NTuple{N,T} # https://discourse.julialang.org/t/why-is-it-not-matching-ntuple-n-t-when-specifying-n/13382
    dip::T

    function PlaneFaultDomain(ft::FT, dom::NTuple{N,T}, dip::T) where {FT,T<:Number,N}
        zero(T) ≤ dip ≤ one(T) * 90 || error("Fautl dip angle ∈ [0, π/2].")
        any(@. dom < zero(T)) && error("Fault domain must be larger than zero.")
        (N == 1 || N == 2) || error("Fault domain should be given by along-downdip [and along-strike] length[s].")
        new{FT,T,N}(ft, dom, dip)
    end
end

function fault(ft::FT, dom::NTuple{N,T}, dip::T) where {FT<:PlaneFaultType,N,T<:Number}
    PlaneFaultDomain(ft, dom, dip)
end

function fault(ft::FT, dom::T, dip::T) where {FT<:PlaneFaultType,T<:Number}
    PlaneFaultDomain(ft, tuple(dom), dip)
end

abstract type AbstractGridType{ndim} end
abstract type FiniteDifferenceGrid{ndim} <: AbstractGridType{ndim} end

struct OnFaultUniformFDGrid{E,T,Val(2)} <: FiniteDifferenceGrid{Val(2)} where {E<:Number,T<:AbstractArray{E}}
    x::T # along-strike
    ξ::T # along-downdip
    Δx::E
    Δξ::E
end


function discretize(fd::PlaneFaultDomain, nelm::NTuple{N,T}) where {N,T<:Integer}
    dims = length(fd.domain)
    dims == N || error("Grid numbers mismatch with fault domain.")
    any(@. nelm ≤ zero(T)) && error("Grid numbers must be larger than 0.")

    lξ, nξ = fd.domain[1], nelm[1]
    Δξ = lξ / nξ
    ξ = collect(range(zero(lξ), stop=-lξ+Δξ, length=nξ)) .- Δξ/2
    if N == 1
        OnFaultUniformFDGrid(missing, ξ, missing, Δξ)
    elseif N == 2
        lx, nx = fd.domain[2], nelm[2]
        Δx = lx / nx
        x = collect(range(-lx/2+Δx/2, stop=lx/2-Δx/2, length=nx))
        OnFaultUniformFDGrid(x, ξ, Δx, Δξ)
    end
end

function discretize(dom::NTuple{1}, nelm::NTuple{1,<:Integer})

end

function discretize(dom::NTuple{2}, nelm::NTuple{2,<:Integer})

end

end
