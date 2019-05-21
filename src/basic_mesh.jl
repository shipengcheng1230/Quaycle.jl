## mesh utilities

export gen_mesh

abstract type AbstractMesh{dim} end
abstract type StructuredMesh{dim} <: AbstractMesh{dim} end
abstract type BasicTransfiniteStructuredMesh{dim} <: StructuredMesh{dim} end
abstract type OkadaMesh{dim} <: BasicTransfiniteStructuredMesh{dim} end

"""
Generate a uniform line mesh in accordance with [`dc3d`](@ref) usage, i.e the line sits
at y-z plane, started from (0, 0, 0) and extended into negative half space.
"""
struct LineOkadaMesh{T<:AbstractVector, U<:Real, I<:Integer, S<:AbstractVector} <: OkadaMesh{1}
    ξ::T # along downdip
    Δξ::U
    nξ::I
    aξ::S
    x::U
    y::T
    z::T
    dep::U # fault origin depth
    dip::U # fault dipping angle
end

"""
Generate a uniform rectangular mesh in accordance with [`dc3d`](@ref) usage, i.e the rectangular sits
parallel to x-axis, top edge starts from z = 0 and centered at x = 0.
The geometry extends into negative half space and rotate around the pivot of (y=0, z=0).
"""
struct RectOkadaMesh{T<:AbstractArray, U<:Real, I<:Integer, S<:AbstractArray} <: OkadaMesh{2}
    x::T # along strike
    Δx::U
    nx::I
    ax::S
    ξ::T
    Δξ::U
    nξ::I
    aξ::S
    y::T
    z::T
    dep::U # fault origin depth
    dip::U # fault dipping angle
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

"Mesh entities of Hex8 for using strain-stress green's function."
struct SBarbotTet4MeshEntity <: UnstructuredMesh{3}

end

"Mesh entities of Tet4 for using strain-stress green's function."
struct SBarbotHex8MeshEntity{T, A, U} <: UnstructuredMesh{3}
    x1::T
    x2::T
    x3::T
    q1::T
    q2::T
    q3::T
    L::T
    T::T
    W::T
    tag::A
    θ::U
end
