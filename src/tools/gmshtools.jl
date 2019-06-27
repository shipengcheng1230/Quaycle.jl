# most of the macro are expected to use at top-level scope
export gen_gmsh_mesh, read_gmsh_mesh, gmsh_vtk_output_cache

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
"""
    gen_gmsh_mesh(::Val{:LineOkada}, ξ::T, Δξ::T, dip::T;
        filename::AbstractString="temp.msh", reg::Integer=1) where T

Generate equivalent [`LineOkadaMesh`](@ref) via [Gmsh](http://gmsh.info/) buildin engine.

## Extra Arguments
- `filename::AbstractString="temp.msh"`: name of the generated mesh file. The file ext will be automatically handled by Gmsh.
- `reg::Integer=1`: the starting tag for entity of any dimension

The rest arguments stay the same as [`gen_mesh`](@ref).
"""
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

"""
    gen_gmsh_mesh(::Val{:RectOkada}, x::T, ξ::T, Δx::T, Δξ::T, dip::T;
        filename::AbstractString="temp.msh", reg::Integer=1) where T

Generate equivalent [`RectOkadaMesh`](@ref) via [Gmsh](http://gmsh.info/) buildin engine.

## Extra Arguments
- `filename::AbstractString="temp.msh"`: name of the generated mesh file. The file ext will be automatically handled by Gmsh.
- `reg::Integer=1`: the starting tag for entity of any dimension

The rest arguments stay the same as [`gen_mesh`](@ref).
"""
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
    gen_gmsh_mesh(mf::OkadaMesh; kwargs...)

Generate an equivalent unstructured mesh as `mf::OkadaMesh`

## Arguments
- `mf::OkadaMesh`: the structured mesh
- `kwargs...`: stay the same as other methods for [`gen_gmsh_mesh`](@ref)
"""
gen_gmsh_mesh(mf::LineOkadaMesh; kwargs...) = gen_gmsh_mesh(Val(:LineOkada), mf.nξ * mf.Δξ, mf.Δξ, mf.dip; kwargs...)
gen_gmsh_mesh(mf::RectOkadaMesh; kwargs...) = gen_gmsh_mesh(Val(:RectOkada), mf.nx * mf.Δx, mf.nξ * mf.Δξ, mf.Δx, mf.Δξ, mf.dip; kwargs...)

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
- `rfx`, `rfy`: refinement coefficients along x-, y-axis using **Bump** algorithm, please refer `gmsh.model.geo.mesh.setTransfiniteCurve`
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

"""
    gen_gmsh_mesh(mf::RectOkadaMesh,
        ::Val{:BoxHexByExtrude},
        llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, ny::I,
        rfx::T, rfy::T, rfzn::AbstractVector, rfzh::AbstractVector;
        filename::AbstractString="temp.msh") where {T, I}
Generate a mesh combinating [`RectOkadaMesh`](@ref) and `BoxHexExtrudeFromSurface` mesh for asthenosphere.
    The first argument is the corresponding [`RectOkadaMesh`](@ref), the rest ones stay the same.
"""
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

tag2linearindice(tags::AbstractVecOrMat) = Dict(tags[k] => k for k in 1: length(tags))

## mesh IO
macro check_and_get_mesh_entity(ecode)
    esc(quote
        nodes = gmsh.model.mesh.getNodes()
        @assert nodes[1][1] == 1 && nodes[1][end] == length(nodes[1]) "Number of nodes tags are not continuous."
        volumetag = phytag ≥ 0 ? gmsh.model.getEntitiesForPhysicalGroup(3, phytag) : phytag
        @assert length(volumetag) == 1 "Multiple entities associated with physical group $(phytag), please distinguish them!"
        es = gmsh.model.mesh.getElements(3, volumetag[1])
        @assert length(es[1]) == 1 "Got more than one element type."
        etag = es[2][1]
        numelements = length(etag)
        etype = es[1][1]
        @assert etype == $(ecode) "Got element type $(gmsh.model.mesh.getElementProperties(etype)[1]), should be $(gmsh.model.mesh.getElementProperties($(ecode))[1])."
        numnodes = gmsh.model.mesh.getElementProperties(etype)[4]
        centers = gmsh.model.mesh.getBarycenters(es[1][1], volumetag[1], 0, 1)
        x2, x1, x3 = centers[1: 3: end], centers[2: 3: end], -centers[3: 3: end]
    end)
end

"""
    read_gmsh_mesh(::Val{:SBarbotHex8}, f::AbstractString;
        phytag::Integer=-1, rotate::Number=0.0, reverse=false, check=false)

Read the mesh and construct mesh entity infomation for SBarbot Hex8 Green's function use.

## Arguments
- `f`: mesh file name
- `phytag`: physical tag for targeting volume entity. If smaller than `0`, retrieve all elements in all 3-dimensional entities. If in
    this case, your mesh must contain only one element type.
- `rotate`: the angle of strike direction, see [`sbarbot_disp_hex8!`](@ref). If your meshing box isn't parallel to x, y-axis, your must
    provide your strike angle manually. By default, the strike angle is zero
- `reverse`: if `true`, reverse the along-x, y-node tag during read. By default, 1→4 in x-axis, 1→2 in y-axis, 1→5 in z-axis
- `check`: if `true`, check that number of distinctive `q1` equals that of `x1`, same for `q2` and `x2` at orthogonal direction,
    which should hold for transfinite mesh.

## Notice
- This function can only be used for Hex8 element with each element lying parallel to z-axis.

- The `check` procedure is not complete for arbitrary strike angle (0 < θ < 90). The user should take a close
    look on the node ordering for one element to ensure the x-, y-extent are correctly resolved by change `reverse` accordingly.
"""
function read_gmsh_mesh(::Val{:SBarbotHex8}, f::AbstractString; phytag::Integer=-1, rotate::Number=0.0, reverse=false, check=false)
    @gmsh_open f begin
        @check_and_get_mesh_entity(5)
        q1, q2, q3, L, T, W = [Vector{Float64}(undef, numelements) for _ in 1: 6]
        @inbounds @fastmath @simd for i in 1: numelements
            ntag1 = es[3][1][numnodes*i-numnodes+1]
            ntag2 = es[3][1][numnodes*i-numnodes+2]
            ntag4 = es[3][1][numnodes*i-numnodes+4]
            ntag5 = es[3][1][numnodes*i-numnodes+5]
            reverse && begin ntag2, ntag4 = ntag4, ntag2 end
            p1x, p1y, p1z = nodes[2][3*ntag1-2], nodes[2][3*ntag1-1], nodes[2][3*ntag1]
            p2x, p2y = nodes[2][3*ntag2-2], nodes[2][3*ntag2-1]
            p4x, p4y = nodes[2][3*ntag4-2], nodes[2][3*ntag4-1]
            p5z = nodes[2][3*ntag5]
            T[i] = hypot(p1x - p4x, p1y - p4y)
            W[i] = abs(p1z - p5z)
            L[i] = hypot(p1x - p2x, p1y - p2y)
            q1[i] = x1[i] - L[i] / 2 * cosd(rotate)
            q2[i] = x2[i] - L[i] / 2 * sind(rotate)
            q3[i] = x3[i] - W[i] / 2
        end
        if check
            f = x -> round(x; digits=3)
            if rotate ≈ 0
                @assert unique(f, x1) |> length == unique(f, q1) |> length "Unmatched `q1` and `x1`, please reset keyword `reverse` to opposite."
            elseif rotate ≈ 90
                @assert unique(f, x2) |> length == unique(f, q2) |> length "Unmatched `q2` and `x2`, please reset keyword `reverse` to opposite."
            end
        end
        SBarbotHex8MeshEntity(x1, x2, x3, q1, q2, q3, L, T, W, rotate, etag)
    end
end

"""
    read_gmsh_mesh(::Val{:SBarbotTet4}, f::AbstractString; phytag::Integer=-1)

Read the mesh and construct mesh entity infomation for SBarbot Tet4 Green's function use.

## Arguments
- `f`: mesh file name
- `phytag`: physical tag for targeting volume entity. If smaller than `0`, retrieve all elements in all 3-dimensional entities. If in
    this case, your mesh must contain only one element type.
"""
function read_gmsh_mesh(::Val{:SBarbotTet4}, f::AbstractString; phytag::Integer=-1)
    @gmsh_open f begin
        @check_and_get_mesh_entity(4)
        A, B, C, D = [[Vector{Float64}(undef, 3) for _ in 1: numelements] for _ in 1: 4]
        @inbounds @fastmath @simd for i in 1: numelements
            ta, tb, tc, td = selectdim(es[3][1], 1, numnodes*i-numnodes+1: numnodes*i-numnodes+4)
            A[i][1] = nodes[2][3*ta-1]
            A[i][2] = nodes[2][3*ta-2]
            A[i][3] = -nodes[2][3*ta]
            B[i][1] = nodes[2][3*tb-1]
            B[i][2] = nodes[2][3*tb-2]
            B[i][3] = -nodes[2][3*tb]
            C[i][1] = nodes[2][3*tc-1]
            C[i][2] = nodes[2][3*tc-2]
            C[i][3] = -nodes[2][3*tc]
            D[i][1] = nodes[2][3*td-1]
            D[i][2] = nodes[2][3*td-2]
            D[i][3] = -nodes[2][3*td]
        end
        SBarbotTet4MeshEntity(x1, x2, x3, A, B, C, D, etag)
    end
end

## paraview
function gmsh_write_vtk_cache(file, phydim, phytag)
    @gmsh_open file begin
        nodes = gmsh.model.mesh.getNodes()
        pts = reshape(nodes[2], 3, :)
        if length(gmsh.model.getPhysicalGroups()) == 0 # no physical group assigned
            entag = -1
        else
            entag = gmsh.model.getEntitiesForPhysicalGroup(phydim, phytag)
            @assert length(entag) == 1 "Multiple entities associated with physical group $(phytag), please distinguish them!"
        end
        es = gmsh.model.mesh.getElements(phydim, entag[1])
        @assert length(es[1][1]) == 1 "More than one type of elements found!"
        etag = es[2][1]
        celltype = gmshcelltype2vtkcelltype[es[1][1]]
        nnode = gmsh.model.mesh.getElementProperties(es[1][1])[4]
        nume = length(es[2][1])
        conn = reshape(es[3][1], Int(nnode), :)
        cells = [MeshCell(celltype, view(conn, :, i)) for i in 1: nume]
        return (pts, cells, etag)
    end
end

"""
    gmsh_vtk_output_cache(file::AbstractString, mf::OkadaMesh{N}, phytag::Integer=-1, datatype=Float64) where N

Create cache of structured [`OkadaMesh`](@ref) for VTK output, which handles the data mapping from structured data
    to unstructured mesh. It's worth mention that currently [WriteVTK](https://github.com/jipolanco/WriteVTK.jl)
    cannot write inclined plane in 3D space. As a workround, it seeks transfering from Gmsh unstructured (transfinite)
    mesh.

## Arguments
- `file::AbstractString`: mesh file containing fault mesh (transfinite)
- `mf::OkadaMesh{N}`: equivalent structured mesh, must match unstructured mesh in the file above
- `phytag::Integer=-1`: physical group tag associated with fault mesh in the mesh file. If smaller than 0,
    retrieve all entities in physical group whose dimension is determined by `mf`. If in this case, only one
    such entity, assumed to be the fault, shall exist.
- `datatype=Float64`: data type for temporary array storing the mapped data
"""
function gmsh_vtk_output_cache(file::AbstractString, mf::OkadaMesh{N}, phytag::Integer=-1, datatype=Float64) where N
    tag = indice2tag(mf, file)
    tagmap = tag2linearindice(tag)
    pts, cells, etag = gmsh_write_vtk_cache(file, N, phytag)
    nume = length(cells)
    dat = Vector{datatype}(undef, nume)
    lidx = Base.OneTo(nume)
    VTKStructuredScalarConversionCache(tagmap, etag, dat, cells, pts, lidx)
end
"""
    gmsh_vtk_output_cache(file::AbstractString, phydim::I, phytag::I) where I<:Integer

Create cache of unstructured mesh for VTK output.

## Arguments
- `file::AbstractString`: mesh file
- `phydim`: physical group dimension, which you will querry
- `phytag`: physical group tag associated with `phydim`. If smaller than 0,
    retrieve all entities in physical group whose dimension is `phydim`. If in this case, only one
    such entity, binded with that physical group, shall exist. If you would like to write multi-block data,
    create VTK output caches for each physical group.
"""
function gmsh_vtk_output_cache(file::AbstractString, phydim::I, phytag::I) where I<:Integer
    pts, cells, _ = gmsh_write_vtk_cache(file, phydim, phytag)
    VTKUnStructuredCache(cells, pts)
end
