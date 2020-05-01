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

"Code snippet for adding a rectangular parallel to x-axis with bottom left at (x, y, z), length of `dx` and width of `hypot(dy, dz)`."
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
    gen_gmsh_mesh(::Val{:InPlaneX},
        llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, nv::I;
        rfxstr::S="Bump", rfx::T=1.0,
        rfvstr::S="Progression", rfv::T=1.0,
        filename::AbstractString="temp.msh", reg::Integer=1) where {T, I, S}

Generate an rectangle with one side parallel to YZ plane.

## Arguments
- `llx`, `lly`, `llz`: coordinates of low-left corner on the top surface
- `dx`, `dy`, `dz`: x-, y-, z-extension
- `nx`, `nv`: number of cells along x-, vertical-axis

## KWARGS
- `rfxstr`, `rfx`, `rfvstr`, `rfv`: the transfinite mesh parameters, see Gmsh docs.
"""
function gen_gmsh_mesh(::Val{:InPlaneX},
    llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, nv::I;
    rfxstr::S="Bump", rfx::T=1.0,
    rfvstr::S="Progression", rfv::T=1.0,
    filename::AbstractString="temp.msh", reg::Integer=1) where {T, I, S}

    @gmsh_do begin
        _reg = geo_rect_x(llx, lly, llz, dx, dy, dz, reg::Integer)
        gmsh.model.geo.mesh.setTransfiniteSurface(_reg-1, "Left", [_reg-10, _reg-9, _reg-8, _reg-7])
        @setTransfiniteCurve begin
            _reg-6, nx+1, rfxstr, rfx
            _reg-5, nv+1, rfvstr, rfv
            _reg-4, nx+1, rfxstr, rfx
            _reg-3, nv+1, rfvstr, rfv
        end
        gmsh.model.geo.mesh.setRecombine(2, _reg-1)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(filename)
    end
end

"""
    gen_gmsh_mesh(mf::RectOkadaMesh, ::Val{:InPlaneX},
        llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, nv::I;
        rfxstr::S="Bump", rfx::T=1.0,
        rfvstr::S="Progression", rfv::T=1.0,
        filename="temp.msh", reg::Integer=1,
        faulttag=(1, "fault"), asthenospheretag=(2, "asthenosphere")) where {T, I, S}

Generate an rectangle with one side parallel to YZ plane along with a Rect Okada mesh.
"""
function gen_gmsh_mesh(mf::RectOkadaMesh, ::Val{:InPlaneX},
    llx::T, lly::T, llz::T, dx::T, dy::T, dz::T, nx::I, nv::I;
    rfxstr::S="Bump", rfx::T=1.0,
    rfvstr::S="Progression", rfv::T=1.0,
    filename="temp.msh", reg::Integer=1,
    faulttag=(1, "fault"), asthenospheretag=(2, "asthenosphere")) where {T, I, S}

    x, ξ = mf.Δx * mf.nx, mf.Δξ * mf.nξ
    y, z = -ξ * cosd(mf.dip), -ξ * sind(mf.dip)

    @gmsh_do begin
        _reg = geo_okada_rect(x, y, z, mf.nx, mf.nξ, reg)
        gmsh.model.addPhysicalGroup(2, [_reg-1], faulttag[1])
        gmsh.model.setPhysicalName(2, faulttag...)
        gmsh.model.geo.synchronize()
        _reg2 = geo_rect_x(llx, lly, llz, dx, dy, dz, _reg)
        gmsh.model.geo.mesh.setTransfiniteSurface(_reg2-1, "Left", [_reg2-10, _reg2-9, _reg2-8, _reg2-7])
        @setTransfiniteCurve begin
            _reg2-6, nx+1, rfxstr, rfx
            _reg2-5, nv+1, rfvstr, rfv
            _reg2-4, nx+1, rfxstr, rfx
            _reg2-3, nv+1, rfvstr, rfv
        end
        gmsh.model.geo.mesh.setRecombine(2, _reg2-1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [_reg2-1], asthenospheretag[1])
        gmsh.model.setPhysicalName(2, asthenospheretag...)
        @addOption begin
            "Mesh.SaveAll", 1 # the mesh is incorrect without this
        end
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
        # it is assumed that the box is the first added 3D volume entity
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
        if VersionNumber(gmsh.GMSH_API_VERSION_MAJOR, gmsh.GMSH_API_VERSION_MINOR) < v"4.3"
            map(p -> gmsh.model.mesh.getElementByCoordinates(mesh.x[p[1]], mesh.y[p[2]], mesh.z[p[2]])[1] |> Int, pos)
        else
            map(p -> gmsh.model.mesh.getElementByCoordinates(mesh.x[p[1]], mesh.y[p[2]], mesh.z[p[2]], 2)[1] |> Int, pos)
        end
    end
end

"Compute `[i] => tag` from [`LineOkadaMesh`](@ref) to unstructured mesh file."
function indice2tag(mesh::LineOkadaMesh, file::AbstractString)
    @gmsh_open file begin
        if VersionNumber(gmsh.GMSH_API_VERSION_MAJOR, gmsh.GMSH_API_VERSION_MINOR) < v"4.3"
            map(p -> gmsh.model.mesh.getElementByCoordinates(0.0, mesh.y[p], mesh.z[p])[1] |> Int, eachindex(mesh.ξ))
        else
            map(p -> gmsh.model.mesh.getElementByCoordinates(0.0, mesh.y[p], mesh.z[p], 1)[1] |> Int, eachindex(mesh.ξ))
        end
    end
end

tag2linearindice(tags::AbstractVecOrMat) = Dict(tags[k] => k for k in 1: length(tags))

## mesh IO
_get_element_node_num(etype::Integer) = _get_element_properties(etype)[4]
_get_element_dim(etype::Integer) = _get_element_properties(etype)[2]
_check_element_type(given::Integer, target::Integer) =
    @assert given == target "Got element type $(_get_element_properties(given)[1]), should be $(_get_element_properties(target)[1])."

function _get_element_properties(etype::Integer)
    @gmsh_do begin
        gmsh.model.mesh.getElementProperties(etype)
    end
end

function _get_nodes(f)
    @gmsh_open f begin
        nodes = gmsh.model.mesh.getNodes()
        @assert nodes[1][1] == 1 && nodes[1][end] == length(nodes[1]) "Number of nodes tags are not continuous."
        return nodes
    end
end

function _get_centers(f, etype, entag::Number)
    @gmsh_open f begin
        gmsh.model.mesh.getBarycenters(etype, entag, 0, 1)
    end
end

_get_centers(f, etype, entags::AbstractVector) = mapreduce(x -> _get_centers(f, etype, x), vcat, entags)

function _get_all_elements_in_physical_group(f, pdim::Integer, ptag::Integer)
    @gmsh_open f begin
        entag = ptag > 0 ? gmsh.model.getEntitiesForPhysicalGroup(pdim, ptag) : ptag
        isa(entag, Integer) && return gmsh.model.mesh.getElements(pdim, entag)
        @assert length(entag) ≥ 1 "No entity found on physical group $(ptag) of dimension $(pdim)."
        es = gmsh.model.mesh.getElements(pdim, entag[1])
        etype = deepcopy(es[1])
        etags = deepcopy(es[2])
        econn = deepcopy(es[3])
        for i ∈ 2: length(entag)
            _es = gmsh.model.mesh.getElements(pdim, entag[i])
            append!(etype, _es[1])
            foreach(append!, etags, _es[2])
            foreach(append!, econn, _es[3])
        end
        return etype, etags, econn
    end
end

function _get_entity_tags_in_physical_group(f, pdim::Integer, ptag::Integer)
    ptag < 0 ? ptag :
        @gmsh_open f begin
            gmsh.model.getEntitiesForPhysicalGroup(pdim, ptag)
        end
end

macro _check_and_get_mesh_entity(f, phytag, ecode)
    esc(quote
        nodes = _get_nodes(f)
        edim = _get_element_dim($(ecode))
        es = _get_all_elements_in_physical_group(f, edim, phytag)
        @assert length(unique(es[1])) == 1 "More than one element type found."
        _check_element_type(es[1][1], $(ecode))
        numelements = length(es[2][1])
        numnodes = _get_element_node_num(es[1][1])
        entag = _get_entity_tags_in_physical_group(f, edim, phytag)
        centers = _get_centers(f, es[1][1], entag)
    end)
end

"""
    read_gmsh_mesh(::Val{:InPlaneX}, f::AbstractString; phytag::Integer=2, inxz::Bool=true)

Read the mesh and construct mesh entity infomation for SBarbot Quad4 in-plane (x-z, no y) Green's function use.
"""
function read_gmsh_mesh(::Val{:InPlaneX}, f::AbstractString; phytag::Integer=2, inxz::Bool=true, rotate::Real=0.0, reverse::Bool=false)
    !inxz && error("InPlaneX denote quad4 at X-Z plane without Y extent.")
    @_check_and_get_mesh_entity(f, phytag, 3)

    x2, x3 = centers[1: 3: end], -centers[3: 3: end] # assume y at 0
    q2, q3, T, W = [Vector{Float64}(undef, numelements) for _ in 1: 4]
    x1 = zeros(eltype(x2), size(x2))

    @inbounds @fastmath @simd for i in 1: numelements
        ntag1 = es[3][1][numnodes*i-numnodes+1]
        ntag2 = es[3][1][numnodes*i-numnodes+2]
        ntag3 = es[3][1][numnodes*i-numnodes+3]
        reverse && begin ntag1, ntag3 = ntag3, ntag1 end
        p1x, p1z = nodes[2][3*ntag1-2], nodes[2][3*ntag1]
        p2x, p2z = nodes[2][3*ntag2-2], nodes[2][3*ntag2]
        p3x, p3z = nodes[2][3*ntag3-2], nodes[2][3*ntag3]

        W[i] = hypot(p1x - p2x, p1z - p2z)
        T[i] = hypot(p2x - p3x, p2z - p3z)
        q2[i] = x2[i] - W[i] / 2 * cosd(rotate)
        q3[i] = x3[i] - W[i] / 2 * sind(rotate)
    end
    SBarbotQuad4InPlaneMeshEntity(x1, x2, x3, q2, q3, T, W, rotate, es[2][1])
end

"""
    read_gmsh_mesh(::Val{:SBarbotHex8}, f::AbstractString;
        phytag::Integer=-1, rotate::Number=0.0, reverse=false, check=false)

Read the mesh and construct mesh entity infomation for SBarbot Hex8 Green's function use.

## Arguments
- `f`: mesh file name
- `phytag`: physical tag for targeting volume entity. If smaller than `0`, retrieve all elements in all 3-dimensional entities. If in
    this case, your mesh must contain only one element type, which should be Hex8
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
    @_check_and_get_mesh_entity(f, phytag, 5)
    x2, x1, x3 = centers[1: 3: end], centers[2: 3: end], -centers[3: 3: end]
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
    SBarbotHex8MeshEntity(x1, x2, x3, q1, q2, q3, L, T, W, rotate, es[2][1])
end

"""
    read_gmsh_mesh(::Val{:SBarbotTet4}, f::AbstractString; phytag::Integer=-1)

Read the mesh and construct mesh entity infomation for SBarbot Tet4 Green's function use.

## Arguments
- `f`: mesh file name
- `phytag`: physical tag for targeting volume entity. If smaller than `0`, retrieve all elements in all 3-dimensional entities. If in
    this case, your mesh must contain only one element type, which should be Tet4.
"""
function read_gmsh_mesh(::Val{:SBarbotTet4}, f::AbstractString; phytag::Integer=-1)
    @_check_and_get_mesh_entity(f, phytag, 4)
    x2, x1, x3 = centers[1: 3: end], centers[2: 3: end], -centers[3: 3: end]
    A, B, C, D = [[Vector{Float64}(undef, 3) for _ in 1: numelements] for _ in 1: 4]
    @inbounds @fastmath @simd for i in 1: numelements
        ta = es[3][1][(i-1)*numnodes+1]
        tb = es[3][1][(i-1)*numnodes+2]
        tc = es[3][1][(i-1)*numnodes+3]
        td = es[3][1][(i-1)*numnodes+4]
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
    SBarbotTet4MeshEntity(x1, x2, x3, A, B, C, D, es[2][1])
end

"""
    read_gmsh_mesh(::Val{:TDTri3}, f::AbstractString; phytag::Integer=-1)

Read the mesh and construct mesh entity infomation for triangular Green's function use.

## Arguments
- `f`: mesh file name
- `phytag`: physical tag for targeting volume entity. If smaller than `0`, retrieve all elements in all 2-dimensional entities. If in
    this case, your mesh must contain only one element type, which should be Tri3.
- `atol`: absolute tolerance, by default `1e-12`, to determine whether the triangle is parallel to axis.
    If this reading procedure does not resolve the slip direction correctly, try to lower this value.
"""
function read_gmsh_mesh(::Val{:TDTri3}, f::AbstractString; phytag::Integer=-1, atol::Real=1e-12)
    @_check_and_get_mesh_entity(f, phytag, 2)
    x, y, z = centers[1: 3: end], centers[2: 3: end], centers[3: 3: end]
    A, B, C, ss, ds, ts = [[Vector{Float64}(undef, 3) for _ in 1: numelements] for _ in 1: 6]
    @inbounds @fastmath @simd for i ∈ 1: numelements
        tag1 = es[3][1][(i-1)*numnodes+1]
        tag2 = es[3][1][(i-1)*numnodes+2]
        tag3 = es[3][1][(i-1)*numnodes+3]
        A[i][1], A[i][2], A[i][3] = nodes[2][3tag1-2], nodes[2][3tag1-1], nodes[2][3tag1]
        B[i][1], B[i][2], B[i][3] = nodes[2][3tag2-2], nodes[2][3tag2-1], nodes[2][3tag2]
        C[i][1], C[i][2], C[i][3] = nodes[2][3tag3-2], nodes[2][3tag3-1], nodes[2][3tag3]
        triangle_geometric_vector!(A[i], B[i], C[i], ss[i], ds[i], ts[i]; atol=atol)
    end
    TDTri3MeshEntity(x, y, z, A, B, C, ss, ds, ts, es[2][1])
end

## paraview
function gmsh_write_vtk_cache(file, phydim, phytag)
    nodes = _get_nodes(file)
    es = _get_all_elements_in_physical_group(file, phydim, phytag)
    @assert length(unique(es[1])) == 1 "More than one element type found."
    nnode = _get_element_node_num(es[1][1])
    etag = es[2][1]
    celltype = gmshcelltype2vtkcelltype[es[1][1]]
    nume = length(es[2][1])
    conn = reshape(es[3][1], Int(nnode), :)
    cells = [MeshCell(celltype, view(conn, :, i)) for i in 1: nume]
    return (reshape(nodes[2], 3, :), cells, etag)
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
    pts, cells, etag = gmsh_write_vtk_cache(file, phydim, phytag)
    VTKUnStructuredCache(cells, pts, etag)
end
