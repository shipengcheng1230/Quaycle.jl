# most of the macro are expected to use at top-level scope

export gen_gmsh_mesh

## mesh generator

"Code snippet for adding a line from (x, y, z) -> (x+dx, y+dy, z+dz)."
function geo_line(x::T, y::T, z::T, dx::T, dy::T, dz::T, reg::Integer) where T<:Real
    @addPoint begin
        x, y, z, 0.0, reg
        x+dx, y+dy, z+dz, 0.0, reg+1
    end
    @addLine begin
        reg, reg+1, reg+2
    end
    return reg+3
end

"Code snippet for adding a rectangular parallel to x-axis with topleft at (x, y, z), length of `dx` and width of `hypot(dy, dz)`."
function geo_rect_x(x::T, y::T, z::T, dx::T, dy::T, dz::T, reg::Integer) where T<:Real
    @addPoint begin
        x, y, z, 0.0, reg
        x+dx, y, z, 0.0, reg+1
        x+dx, y+dy, z+dz, 0.0, reg+2
        x, y+dy, z+dz, 0.0, reg+3
    end
    @addLine begin
        # top left -> top right -> bottom right (reverse) -> bottom left (reverse) -> top left
        reg, reg+1, reg+4
        reg+1, reg+2, reg+5
        reg+3, reg+2, reg+6
        reg, reg+3, reg+7
    end
    gmsh.model.geo.addCurveLoop([reg+4, reg+5, -reg-6, -reg-7], reg+8)
    gmsh.model.geo.addPlaneSurface([reg+8], reg+9)
    return reg+10
end

"Code snippet for adding [`LineOkadaMesh`](@ref)."
function geo_okada_line(y::T, z::T, nξ::I, reg::I) where {T<:Real, I<:Integer}
    _reg = geo_line(0.0, 0.0, 0.0, 0.0, y, z, reg)
    @setTransfiniteCurve begin
        _reg-1, nξ+1, "Progression", 1.0
    end
    return _reg
end

"Code snippet for adding [`RectOkadaMesh`](@ref)."
function geo_okada_rect(x::T, y::T, z::T, nx::I, nξ::I, reg::I) where {T<:Real, I<:Integer}
    _reg = geo_rect_x(-x/2, 0.0, 0.0, x, y, z, reg)
    gmsh.model.geo.mesh.setTransfiniteSurface(_reg-1, "Left", [_reg-10, _reg-9, _reg-8, _reg-7])
    @setTransfiniteCurve begin
        _reg-6, nx+1, "Progression", 1.0
        _reg-5, nξ+1, "Progression", 1.0
        _reg-4, nx+1, "Progression", 1.0
        _reg-3, nξ+1, "Progression", 1.0
    end
    return _reg
end

"Code snippet for adding box extruded from surface (x-y ↑ z)."
function geo_box_extruded_from_surfaceXY(llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, ny::I, rfx::T, rfy::T, rfzn::AbstractVector, rfzh::AbstractVector, reg::I=1) where {T<:Real, I<:Integer}
    _reg = geo_rect_x(llx, lly, llz, dx, dy, 0.0, reg)
    gmsh.model.geo.mesh.setTransfiniteSurface(_reg-1, "Left", [_reg-10, _reg-9, _reg-8, _reg-7])
    @setTransfiniteCurve begin
        _reg-6, nx+1, "Bump", rfx
        _reg-5, ny+1, "Bump", rfy
        _reg-4, nx+1, "Bump", rfx
        _reg-3, ny+1, "Bump", rfy
    end
    gmsh.model.geo.mesh.setRecombine(2, _reg-1)
    gmsh.model.geo.extrude([(2, _reg-1)], 0.0, 0.0, -dz, rfzn, rfzh, true)
    return _reg
end

"Generate equivalent [`LineTopCenterMesh`](@ref) via [Gmsh](http://gmsh.info/) buildin engine."
function gen_gmsh_mesh(::Val{:LineOkada}, ξ::T, Δξ::T, dip::T; filename::AbstractString="temp.msh", reg::Integer=1) where T
    @gmsh_do begin
        dy, dz = -ξ * cosd(dip), -ξ * sind(dip)
        nξ = round(Int, ξ / Δξ) # same as counting length of `range` in `mesh_downdip`
        _reg = geo_okada_line(dy, dz, nξ, reg)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(1)
        gmsh.write(filename)
        return _reg
    end
end

"Generate equivalent [`RectTopCenterMesh`](@ref) via [Gmsh](http://gmsh.info/) buildin engine."
function gen_gmsh_mesh(::Val{:RectOkada}, x::T, ξ::T, Δx::T, Δξ::T, dip::T; filename::AbstractString="temp.msh", reg::Integer=1) where T
    @gmsh_do begin
        y, z = -ξ * cosd(dip), -ξ * sind(dip)
        nξ = round(Int, ξ / Δξ) # same as counting length of `range` in `mesh_downdip`
        nx = round(Int, x / Δx) # same as counting length of `range` in `mesh_strike`
        _reg = geo_okada_rect(x, y, z, nx, nξ, reg)
        gmsh.model.geo.mesh.setRecombine(2, _reg-1)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(filename)
        return _reg
    end
end

"""
    gen_gmsh_mesh(::Val{:BoxHexByExtrude},
        llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, ny::I,
        rfx::T, rfy::T, rfzn::AbstractVector, rfzh::AbstractVector;
        filename::AbstractString="temp.msh") where {T, I}

Gernate a box for [Asthenosphere](https://en.wikipedia.org/wiki/Asthenosphere) using 8-node hexahedron
    elements (via setting transfinite curve).

## Arguments
- `llx`, `lly`, `llz`: coordinates of low-left corner on the top surface
- `dx`, `dy`, `dz`: x-, y-, z-extension
- `nx`, `ny`: number of cells along x-, y-axis
- `rfx`, `rfy`: refinement coefficients along x-, y-axis using ["Bump"] algorithm, please refer `gmsh.model.geo.mesh.setTransfiniteCurve`
- `rfzn`: number of cells along z-axis, please refer `numElements` in `gmsh.model.geo.extrude`
- `rfzh`: accumulated height of cells along z-axis, please refer `heights` in `gmsh.model.geo.extrude`
"""
function gen_gmsh_mesh(::Val{:BoxHexExtrudeFromSurface},
    llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, ny::I,
    rfx::T, rfy::T, rfzn::AbstractVector, rfzh::AbstractVector;
    filename::AbstractString="temp.msh", reg::Integer=1) where {T, I}

    @gmsh_do begin
        _reg = geo_box_extruded_from_surfaceXY(llx, lly, llz, dx, dy, dz, nx, ny, rfx, rfy, rfzn, rfzh, reg)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.write(filename)
        return _reg
    end
end

"Adding a combination of [`RectOkadaMesh`](@ref) and `BoxHexExtrudeFromSurface` mesh"
function gen_gmsh_mesh(
    mtype1::Val{:RectOkada}, x::T, ξ::T, Δx::T, Δξ::T, dip::T,
    mtype2::Val{:BoxHexExtrudeFromSurface}, llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, ny::I, rfx::T, rfy::T, rfzn::AbstractVector, rfzh::AbstractVector;
    filename="temp.msh", reg::Integer=1)

    y, z = -ξ * cosd(dip), -ξ * sind(dip)
    nξ = round(Int, ξ / Δξ) # same as counting length of `range` in `mesh_downdip`
    nx = round(Int, x / Δx) # same as counting length of `range` in `mesh_strike`
    _reg = geo_okada_rect(x, y, z, nx, nξ, reg)
    _reg2 = geo_box_extruded_from_surfaceXY(llx, lly, llz, dx, dy, dz, nx, ny, rfx, rfy, rfzn, rfzh, _reg)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.write(filename)
    return _reg2
end

## elements tag mapping

"Compute `[i,j] => tag` from [`RectOkadaMesh`](@ref) to unstructured mesh file."
function indice2tag(mesh::RectOkadaMesh, file::AbstractString)
    @gmsh_open file begin
        pos = Iterators.product(1: mesh.nx, 1: mesh.nξ)
        map(p -> gmsh.model.mesh.getElementByCoordinates(mesh.x[p[1]], mesh.y[p[2]], mesh.z[p[2]], 2)[1] |> Int, pos)
    end
end

"Compute `[i] => tag` from [`LineOkadaMesh`](@ref) to unstructured mesh file."
function indice2tag(mesh::LineOkadaMesh, file::AbstractString)
    @gmsh_open file begin
        pos = 1: mesh.nξ
        map(p -> gmsh.model.mesh.getElementByCoordinates(0.0, mesh.y[p], mesh.z[p], 1)[1] |> Int, pos)
    end
end

## mesh IO
