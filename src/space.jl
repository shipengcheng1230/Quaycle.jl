export compose
export DIPPING, STRIKING

abstract type AbstractFaultType end
abstract type PlaneFault <: AbstractFaultType end
abstract type FlatPlaneFault <: PlaneFault end

"Dipping, indicate dislocation occurs at y-z plane in [`dc3d`](@ref) use."
struct DIPPING <: FlatPlaneFault end
"Striking, indicate dislocation occurs at x-direction in [`dc3d`](@ref) use."
struct STRIKING <: FlatPlaneFault end

abstract type AbstractFaultSpace end

"Okada fault space encapsulating mesh and fault type"
struct PlaneFaultSpace{N, M} <: AbstractFaultSpace
    mesh::M

    function PlaneFaultSpace(mesh::AbstractMesh{N}) where N
        new{N, typeof(mesh)}(mesh)
    end
end

"Generate fault space encapsulating [`LineOkadaMesh`](@ref)"
function fault(ftype::Val{:LineOkada}, ξ::T, Δξ::T, dip::T) where T<:Real
    mesh = gen_mesh(ftype, ξ, Δξ, dip)
    return PlaneFaultSpace(mesh)
end

"Generate fault space encapsulating [`RectOkadaMesh`](@ref)."
function fault(ftype::Val{:RectOkada}, x::T, ξ::T, Δx::T, Δξ::T, dip::T) where T<:Real
    mesh = gen_mesh(ftype, x, ξ, Δx, Δξ, dip)
    return PlaneFaultSpace(mesh)
end

fault(mesh::OkadaMesh) = PlaneFaultSpace(mesh)

abstract type AbstractLithosphereAsthenosphereSpace end

struct LithAsthSpace{dim, M1, M2} <: AbstractLithosphereAsthenosphereSpace
    me::M1
    mv::M2

    function LithAsthSpace(me::M1, mv::M2, dim::Integer) where {M1<:AbstractMesh, M2<:AbstractMesh}
        new{dim, M1, M2}(me, mv)
    end
end

"""
    compose(me::AbstractMesh, me::SBarbotMeshEntity{N})

Create a composite space consisting of both fault mesh and asthenosphere mesh.

## Arguments
- `me::AbstractMesh`: fault mesh
- `me::SBarbotMeshEntity{N}`: asthenosphere mesh
"""
compose(me::AbstractMesh, mv::SBarbotMeshEntity{N}) where N = LithAsthSpace(me, mv, N)
compose(fa::PlaneFaultSpace, me::SBarbotMeshEntity{N}) where N = compose(fa.mesh, me)
