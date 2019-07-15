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
    ξ, nξ, aξ, y, z = _equidist_mesh_downdip(ξ, Δξ, dip)
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
    ξ, nξ, aξ, y, z = _equidist_mesh_downdip(ξ, Δξ, dip)
    x, nx, ax = _equidist_mesh_strike(x, Δx)
    return RectOkadaMesh(x, Δx, nx, ax, ξ, Δξ, nξ, aξ, y, z, zero(T), dip)
end

function _equidist_mesh_downdip(ξ::T, Δξ::T, dip::T) where T
    ξi = range(zero(T), stop=-ξ+Δξ, step=-Δξ) .- Δξ/2 |> collect
    aξ = [[w - Δξ/2, w + Δξ/2] for w in ξi]
    y, z = ξi .* cosd(dip), ξi .* sind(dip)
    return ξi, length(ξi), aξ, y, z
end

function _equidist_mesh_strike(x::T, Δx::T) where T
    xi = range(-x/2 + Δx/2, stop=x/2 - Δx/2, step=Δx) |> collect
    ax = [[w - Δx/2, w + Δx/2] for w in xi]
    return xi, length(xi), ax
end

## unstructured mesh entities
abstract type UnstructuredMesh{dim} <: AbstractMesh{dim} end
abstract type SBarbotMeshEntity{dim} <: UnstructuredMesh{dim} end
abstract type TriangularMesh <: UnstructuredMesh{2} end

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
    @assert size(tag) ==
        size(x1) == size(x2) == size(x3) ==
        size(A) == size(B) == size(C) == size(D)
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
    @assert size(tag) ==
        size(x1) == size(x2) == size(x3) ==
        size(q1) == size(q2) == size(q3) ==
        size(L) == size(T) == size(W)
end

"Mesh entities of Tri3 for using dislocaiton-stress green's function"
@with_kw struct TDTri3MeshEntity{T<:AbstractVector, V<:AbstractVector, VI<:AbstractVector} <: TriangularMesh
    x::T # centroid x
    y::T # centroid y
    z::T # centroid z
    A::V
    B::V
    C::V
    ss::V # strike
    ds::V # dip
    ts::V # tensile
    tag::VI

    @assert size(x) == size(y) == size(z)
    @assert size(A) == size(B) == size(C) == size(ss) == size(ds) == size(ts)
    @assert mapreduce(x -> x[3], (a, b) -> a < b ? a : b, A) ≤ 0
    @assert mapreduce(x -> x[3], (a, b) -> a < b ? a : b, B) ≤ 0
    @assert mapreduce(x -> x[3], (a, b) -> a < b ? a : b, C) ≤ 0
end

## geometry
"""
This function calculates the strike and downdip direction given a triangular
    element in Local east, north, up (ENU) coordinates.
    Notice when the element is coplanar at ``z = 0`` it cannot tell such
    information. Also element vertices cannot be collinear.
    It's worth mentioning that the normal vector, as in our coordinate choice,
    shall point upwards and `A`, `B` and `C` are in couter-clockwise direction
    when viewing from above. If those points are coplanar horizontally, the normal
    vector shall point to negative-y-section.


The geometry vector is obtained through solving plane equation:
```math
ax + by + cz + d = 0
```

You may use [Reduce.jl](https://github.com/chakravala/Reduce.jl) to obtain the solution
    without populating a temporary matrix, for instance if ``d ≠ 0``:
```julia
using Reduce

R"(mat((A[1], A[2], A[3]), (B[1], B[2], B[3]), (C[1], C[2], C[3])))^(-1) * mat((-1), (-1), (-1))" |> rcall
```
"""
function triangle_geometric_vector!(A::V, B::V, C::V, ss::V, ds::V, ts::V; atol=1e-12) where V <: AbstractVector{T} where T
    ts .= normal_vector(A, B, C)
    if ts[3] < zero(T) || (ts[3] ≈ zero(T) && ts[2] > zero(T)) # point to hanging wall
        ts .*= -one(T)
        for i in eachindex(A)
            A[i], B[i] = B[i], A[i]
        end
    end
    if A[3] ≈ B[3] ≈ C[3]
        @warn "Coplanar at `z = const` where we cannot tell strike and downdip direction."
        ss .= NaN, NaN, NaN
        ds .= NaN, NaN, NaN
        ts .= zero(T), zero(T), one(T)
    elseif isapprox(ts[1], zero(T); atol=atol) # parallel -x
        ss .= one(T), zero(T), zero(T)
        ds .= cross(ts, ss)
    elseif isapprox(ts[2], zero(T); atol=atol) # parallel -y
        ss .= zero(T), one(T), zero(T)
        ds .= cross(ts, ss)
    elseif isapprox(ts[3], zero(T); atol=atol) # parallel -z
        ds .= zero(T), zero(T), one(T)
        ss .= cross(ds, ts)
    else
        # r = (A[2]*C[3] - C[2]*A[3])*B[1] - (B[2]*C[3] - C[2]*B[3])*A[1] - (A[2]*B[3] - B[2]*A[3])*C[1]
        a = - ((A[3] - C[3])*B[2] - (B[3] - C[3])*A[2] - (A[3] - B[3])*C[2]) # / r
        b = (A[3] - C[3])*B[1] - (B[3] - C[3])*A[1] - (A[3] - B[3])*C[1] # / r
        # c = - ((A[2] - C[2])*B[1] - (B[2] - C[2])*A[1] - (A[2] - B[2])*C[1]) / r

        ss[1], ss[2], ss[3] = a, -b, zero(T)
        if a < zero(T) ss .*= -one(T) end # to positive x-axis
        normalize!(ss)
        ds .= cross(ts, ss)
        normalize!(ds)
    end
end

@inline normal_vector(p1::T, p2::T, p3::T) where T<:AbstractVector = normalize!(cross(p2 - p1, p3 - p1))
