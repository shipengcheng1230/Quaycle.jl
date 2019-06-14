## mesh utilities

export gen_mesh

abstract type AbstractMesh{dim} end
abstract type StructuredMesh{dim} <: AbstractMesh{dim} end
abstract type RectLinearMesh{dim} <: StructuredMesh{dim} end
abstract type OkadaMesh{dim} <: RectLinearMesh{dim} end

"""
Generate a uniform line mesh in accordance with [`dc3d`](@ref) usage, i.e the line sits
at y-z plane, started from (0, 0, 0) and extended into negative half space.
"""
@with_kw struct LineOkadaMesh{T<:AbstractVector, U<:Real, I<:Integer, S<:AbstractVector} <: OkadaMesh{1}
    ξ::T # centroid along downdip
    Δξ::U
    nξ::I
    aξ::S
    x::U
    y::T
    z::T
    dep::U # fault origin depth
    dip::U # fault dipping angle

    @assert length(ξ) == nξ
    @assert length(aξ) == nξ
    @assert length(y) == nξ
    @assert length(z) == nξ
end

"""
Generate a uniform rectangular mesh in accordance with [`dc3d`](@ref) usage, i.e the rectangular sits
parallel to x-axis, top edge starts from z = 0 and centered at x = 0.
The geometry extends into negative half space and rotate around the pivot of (y=0, z=0).
"""
@with_kw struct RectOkadaMesh{T<:AbstractArray, U<:Real, I<:Integer, S<:AbstractArray} <: OkadaMesh{2}
    x::T # centroid along strike
    Δx::U
    nx::I
    ax::S
    ξ::T # centroid along downdip
    Δξ::U
    nξ::I
    aξ::S
    y::T
    z::T
    dep::U # fault origin depth
    dip::U # fault dipping angle

    @assert length(x) == nx
    @assert length(ax) == nx
    @assert length(ξ) == nξ
    @assert length(aξ) == nξ
    @assert length(y) == nξ
    @assert length(z) == nξ
end

"""
    gen_mesh(::Val{:LineOkada}, ξ::T, Δξ::T, dip::T)

Generate [`LineOkadaMesh`](@ref)

## Arguments
- `ξ`: downdip length
- `Δξ`: downdip interval
- `dip`: dipping angle
"""
function gen_mesh(::Val{:LineOkada}, ξ::T, Δξ::T, dip::T) where T
    ξ, nξ, aξ, y, z = mesh_downdip(ξ, Δξ, dip)
    return LineOkadaMesh(ξ, Δξ, nξ, aξ, zero(T), y, z, zero(T), dip)
end

"""
    gen_mesh(::Val{:RectOkada}, x::T, ξ::T, Δx::T, Δξ::T, dip::T)

Generate [`RectOkadaMesh`](@ref)

## Arguments
- `x`: along strike length
- `ξ`: downdip length
- `Δx`: along strike interval
- `Δξ`: downdip interval
- `dip`: dipping angle
"""
function gen_mesh(::Val{:RectOkada}, x::T, ξ::T, Δx::T, Δξ::T, dip::T) where T
    ξ, nξ, aξ, y, z = mesh_downdip(ξ, Δξ, dip)
    x, nx, ax = mesh_strike(x, Δx)
    return RectOkadaMesh(x, Δx, nx, ax, ξ, Δξ, nξ, aξ, y, z, zero(T), dip)
end

function mesh_downdip(ξ::T, Δξ::T, dip::T) where T
    ξi = range(zero(T), stop=-ξ+Δξ, step=-Δξ) .- Δξ/2 |> collect
    aξ = [[w - Δξ/2, w + Δξ/2] for w in ξi]
    y, z = ξi .* cosd(dip), ξi .* sind(dip)
    return ξi, length(ξi), aξ, y, z
end

function mesh_strike(x::T, Δx::T) where T
    xi = range(-x/2 + Δx/2, stop=x/2 - Δx/2, step=Δx) |> collect
    ax = [[w - Δx/2, w + Δx/2] for w in xi]
    return xi, length(xi), ax
end

## unstructured mesh entities
abstract type UnstructuredMesh{dim} <: AbstractMesh{dim} end
abstract type SBarbotMeshEntity{dim} <: UnstructuredMesh{dim} end

"Mesh entities of Tet4 for using strain-stress green's function."
@with_kw struct SBarbotTet4MeshEntity{P<:AbstractArray, Q<:AbstractVector, T<:AbstractVector} <: SBarbotMeshEntity{3}
    x1::P # centroid +y
    x2::P # centroid +x
    x3::P # centroid -z
    A::Q
    B::Q
    C::Q
    D::Q
    tag::T

    @assert minimum(x3) > 0
    @assert map(x -> x[3], A) |> minimum > 0
    @assert map(x -> x[3], B) |> minimum > 0
    @assert map(x -> x[3], C) |> minimum > 0
    @assert map(x -> x[3], D) |> minimum > 0
    @assert size(tag) == size(x1)
    @assert size(tag) == size(x2)
    @assert size(tag) == size(x3)
    @assert size(tag) == size(A)
    @assert size(tag) == size(B)
    @assert size(tag) == size(C)
    @assert size(tag) == size(D)
end

"Mesh entities of Hex8 for using strain-stress green's function."
@with_kw struct SBarbotHex8MeshEntity{P<:AbstractVector, A<:AbstractVector, U<:Number} <: SBarbotMeshEntity{3}
    x1::P # centroid +y
    x2::P # centroid +x
    x3::P # centroid -z
    q1::P # anchor +y
    q2::P # anchor +x
    q3::P # anchor -z
    L::P # length ←y→ when θ = 0°
    T::P # thickness ←x→ when θ = 0°
    W::P # width ←z→, always
    θ::U # strike, clockwise away from +y
    tag::A # element tag

    @assert minimum(x3) > 0
    @assert minimum(q3) ≥ 0
    @assert minimum(L) > 0
    @assert minimum(W) > 0
    @assert minimum(T) > 0
    @assert size(tag) == size(x1)
    @assert size(tag) == size(x2)
    @assert size(tag) == size(x3)
    @assert size(tag) == size(q1)
    @assert size(tag) == size(q2)
    @assert size(tag) == size(q3)
    @assert size(tag) == size(L)
    @assert size(tag) == size(T)
    @assert size(tag) == size(W)
end
