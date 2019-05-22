# most of the macro are expected to use at top-level scope

export gen_gmsh_mesh, read_gmsh_mesh, indice2tag

## code snippet for mesh generator

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
    gmsh.model.geo.mesh.setRecombine(2, _reg-1)
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
    outDimTags = gmsh.model.geo.extrude([(2, _reg-1)], 0.0, 0.0, -dz, rfzn, rfzh, true)
    return outDimTags
end

## concrete mesh generator

"Generate equivalent [`LineTopCenterMesh`](@ref) via [Gmsh](http://gmsh.info/) buildin engine."
function gen_gmsh_mesh(::Val{:LineOkada}, ξ::T, Δξ::T, dip::T; filename::AbstractString="temp.msh", reg::Integer=1) where T
    @gmsh_do begin
        dy, dz = -ξ * cosd(dip), -ξ * sind(dip)
        nξ = round(Int, ξ / Δξ) # same as counting length of `range` in `mesh_downdip`
        _reg = geo_okada_line(dy, dz, nξ, reg)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(1)
        gmsh.write(filename)
    end
end

"Generate equivalent [`RectTopCenterMesh`](@ref) via [Gmsh](http://gmsh.info/) buildin engine."
function gen_gmsh_mesh(::Val{:RectOkada}, x::T, ξ::T, Δx::T, Δξ::T, dip::T; filename::AbstractString="temp.msh", reg::Integer=1) where T
    @gmsh_do begin
        y, z = -ξ * cosd(dip), -ξ * sind(dip)
        nξ = round(Int, ξ / Δξ) # same as counting length of `range` in `mesh_downdip`
        nx = round(Int, x / Δx) # same as counting length of `range` in `mesh_strike`
        _reg = geo_okada_rect(x, y, z, nx, nξ, reg)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(filename)
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
    end
end

"Adding a combination of [`RectOkadaMesh`](@ref) and `BoxHexExtrudeFromSurface` mesh for asthenosphere."
function gen_gmsh_mesh(mf::RectOkadaMesh, ::Val{:BoxHexExtrudeFromSurface},
    llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, ny::I, rfx::T, rfy::T, rfzn::AbstractVector, rfzh::AbstractVector;
    filename="temp.msh", reg::Integer=1,
    faulttag=(1, "fault"), asthenospheretag=(1, "asthenosphere")) where {T, I}

    x, ξ = mf.Δx * mf.nx, mf.Δξ * mf.nξ
    y, z = -ξ * cosd(mf.dip), -ξ * sind(mf.dip)
    @gmsh_do begin
        _reg = geo_okada_rect(x, y, z, mf.nx, mf.nξ, reg)
        gmsh.model.addPhysicalGroup(2, [_reg-1], faulttag[1])
        gmsh.model.setPhysicalName(2, faulttag...)
        # type of `_reg2` is tuple, each is `(dim, tag)`
        _reg2 = geo_box_extruded_from_surfaceXY(llx, lly, llz, dx, dy, dz, nx, ny, rfx, rfy, rfzn, rfzh, _reg)
        # it is assumed that the box is the first added volume entity
        volumetag = _reg2[findfirst(x -> x[1] == 3, _reg2)][2]
        gmsh.model.addPhysicalGroup(3, [volumetag], asthenospheretag[1])
        gmsh.model.setPhysicalName(3, asthenospheretag...)
        @addOption begin
            "Mesh.SaveAll", 1 # the mesh is incorrect without this
        end
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.write(filename)
    end
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

macro check_and_get_mesh_entity(ecode)
    esc(quote
        nodes = gmsh.model.mesh.getNodes()
        @assert nodes[1][1] == 1 && nodes[1][end] == length(nodes[1]) "Number of nodes tags are not continuous."
        volumetag = phytag ≥ 0 ? gmsh.model.getEntitiesForPhysicalGroup(3, phytag) : phytag
        es = gmsh.model.mesh.getElements(3, volumetag)
        @assert length(es[1]) == 1 "Got more than one element type."
        etag = es[2][1]
        numelements = length(etag)
        etype = es[1][1]
        @assert etype == $(ecode) "Got element type $(gmsh.model.mesh.getElementProperties(etype)[1]), should be $(gmsh.model.mesh.getElementProperties($(ecode))[1])."
        numnodes = gmsh.model.mesh.getElementProperties(etype)[4]
        centers = gmsh.model.mesh.getBarycenters(es[1][1], volumetag, 0, 1)
        x2, x1, x3 = centers[1: 3: end], centers[2: 3: end], -centers[3: 3: end]
    end)
end

function read_gmsh_mesh(::Val{:SBarbotHex8}, f::AbstractString; phytag::Integer=1::Integer, rotate::Number=90.0)
    @gmsh_open f begin
        @check_and_get_mesh_entity(5)
        q1, q2, q3, L, T, W = [Vector{Float64}(undef, numelements) for _ in 1: 6]
        @inbounds @fastmath @simd for i in 1: numelements
            # in Gmsh Hex8, vertex-1 and vertex-7 are the volume diagonal
            ntag1 = es[3][1][numnodes*i-numnodes+1]
            ntag7 = es[3][1][numnodes*i-numnodes+7]
            p1x, p1y, p1z = nodes[2][3*ntag1-2], nodes[2][3*ntag1-1], nodes[2][3*ntag1]
            p7x, p7y, p7z = nodes[2][3*ntag7-2], nodes[2][3*ntag7-1], nodes[2][3*ntag7]
            L[i] = abs(p7x - p1x)
            W[i] = abs(p7z - p1z)
            T[i] = abs(p7y - p1y)
            q1[i] = x1[i]
            q2[i] = x2[i] - L[i]/2
            q3[i] = x3[i] - W[i]/2
        end
        return SBarbotHex8MeshEntity(x1, x2, x3, q1, q2, q3, L, T, W, etag, rotate)
    end
end

function read_gmsh_mesh(::Val{:SBarbotTet4}, f::AbstractString; phytag=1::Integer)
    @gmsh_open f begin
        @check_and_get_mesh_entity(4)
        A, B, C, D = [[Vector{Float64}(undef, 3) for _ in 1: numelements] for _ in 1: 4]
        for i in 1: numelements
            ta, tb, tc, td = selectdim(es[3][1], 1, numnodes*i-numnodes+1: numnodes*i-numnodes+4)
            A[i] .= selectdim(nodes[2], 1, 3*ta-2: 3*ta)
            B[i] .= selectdim(nodes[2], 1, 3*tb-2: 3*tb)
            C[i] .= selectdim(nodes[2], 1, 3*tc-2: 3*tc)
            D[i] .= selectdim(nodes[2], 1, 3*td-2: 3*td)
        end
        return SBarbotTet4MeshEntity(x1, x2, x3, A, B, C, D, etag)
    end
end
