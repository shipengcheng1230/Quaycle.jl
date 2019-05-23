export fault
export DIPPING, STRIKING

abstract type AbstractFault end
abstract type PlaneFault <: AbstractFault end

"Dipping, indicate dislocation occurs at y-z plane in [`dc3d`](@ref) use."
struct DIPPING <: PlaneFault end
"Striking, indicate dislocation occurs at x-direction in [`dc3d`](@ref) use."
struct STRIKING <: PlaneFault end

abstract type AbstractFaultSpace end

"Okada fault space encapsulating mesh and fault type"
struct OkadaFaultSpace{dim, FT, M} <: AbstractFaultSpace
    mesh::M
    ft::FT

    function OkadaFaultSpace(mesh::M, ft::FT, dim::Integer) where {M<:OkadaMesh, FT<:PlaneFault}
        new{dim, FT, M}(mesh, ft)
    end
end

OkadaFaultSpace(mesh::OkadaMesh{1}, ft::PlaneFault) = OkadaFaultSpace(mesh, ft, 1)
OkadaFaultSpace(mesh::OkadaMesh{2}, ft::PlaneFault) = OkadaFaultSpace(mesh, ft, 2)

"Generate fault space encapsulating [`LineOkadaMesh`](@ref) and fault type `ft`."
function fault(ftype::Val{:LineOkada}, ft::FT, ξ::T, Δξ::T, dip::T) where {T<:Real, FT<:PlaneFault}
    mesh = gen_mesh(ftype, ξ, Δξ, dip)
    return OkadaFaultSpace(mesh, ft)
end

"Generate fault space encapsulating [`RectOkadaMesh`](@ref) and fault type `ft`."
function fault(ftype::Val{:RectOkada}, ft::FT, x::T, ξ::T, Δx::T, Δξ::T, dip::T) where {T<:Real, FT<:PlaneFault}
    mesh = gen_mesh(ftype, x, ξ, Δx, Δξ, dip)
    return OkadaFaultSpace(mesh, ft)
end

fault(fty::Val{:LineOkada}, ft::STRIKING, ξ::T, Δξ::T) where T = fault(fty, ft, ξ, Δξ, 90*one(T))
fault(fty::Val{:RectOkada}, ft::STRIKING, x::T, ξ::T, Δx::T, Δξ::T) where T = fault(fty, ft, x, ξ, Δx, Δξ, 90*one(T))
