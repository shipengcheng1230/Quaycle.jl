export fault
export DIPPING, STRIKING

abstract type AbstractFault end
abstract type PlaneFault <: AbstractFault end

"Dipping, indicate dislocation occurs at y-z plane in [`dc3d`](@ref) use."
struct DIPPING <: PlaneFault end
"Striking, indicate dislocation occurs at x-direction in [`dc3d`](@ref) use."
struct STRIKING <: PlaneFault end

abstract type FaultSpace end

"Basic fault space encapsulating mesh and fault type"
struct BasicFaultSpace{dim, FT, M} <: FaultSpace
    mesh::M
    ft::FT

    function BasicFaultSpace(mesh::M, ft::FT, dim::Integer) where {M<:TopCenterPlaneMesh, FT<:PlaneFault}
        new{dim, FT, M}(mesh, ft)
    end
end

BasicFaultSpace(mesh::LineTopCenterMesh, ft::PlaneFault) = BasicFaultSpace(mesh, ft, 1)
BasicFaultSpace(mesh::RectTopCenterMesh, ft::PlaneFault) = BasicFaultSpace(mesh, ft, 2)

fault(::Val{:topcenter}, ft::FT, mesh::M) where {FT<:PlaneFault, M<:TopCenterPlaneMesh} = BasicFaultSpace(mesh, ft)

"Generate fault space encapsulating basic topcenter mesh [`LineTopCenterMesh`](@ref) and fault type `ft`."
function fault(ftype::Val{:topcenter}, ft::FT, ξ::T, Δξ::T, dip::T) where {T<:Real, FT<:PlaneFault}
    mesh = gen_mesh(ftype, ξ, Δξ, dip)
    return BasicFaultSpace(mesh, ft)
end

"Generate fault space encapsulating basic topcenter mesh [`RectTopCenterMesh`](@ref) and fault type `ft`."
function fault(ftype::Val{:topcenter}, ft::FT, x::T, ξ::T, Δx::T, Δξ::T, dip::T) where {T<:Real, FT<:PlaneFault}
    mesh = gen_mesh(ftype, x, ξ, Δx, Δξ, dip)
    return BasicFaultSpace(mesh, ft)
end

fault(fty::Val{:topcenter}, ft::STRIKING, ξ::T, Δξ::T) where T = fault(fty, ft, ξ, Δξ, 90*one(T))
fault(fty::Val{:topcenter}, ft::STRIKING, x::T, ξ::T, Δx::T, Δξ::T) where T = fault(fty, ft, x, ξ, Δx, Δξ, 90*one(T))
