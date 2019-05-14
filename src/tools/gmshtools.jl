export gen_gmsh_mesh

"Generate equivalent [`LineTopCenterMesh`](@ref) via [Gmsh](http://gmsh.info/)."
function gen_gmsh_mesh(::Val{:topcenter}, ξ::T, Δξ::T, dip::T; filename::AbstractString="temp.msh") where T
    @gmsh_do begin
        model = gmsh.model
        factory = model.geo
        y, z = -ξ * cosd(dip), -ξ * sind(dip)
        # same as counting length of `range` in `mesh_downdip`
        nξ = ceil(Int, ξ / Δξ)
        @addPoint begin
            0.0, 0.0, 0.0, 1
            0.0, y, z, 2
        end
        @addLine begin
            1, 2, 1
        end
        @setTransfiniteCurve begin
            1, nξ+1, "Progression", 1.0
        end
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(1)
        gmsh.write(filename)
    end
end

"Generate equivalent [`RectTopCenterMesh`](@ref) via [Gmsh](http://gmsh.info/)."
function gen_gmsh_mesh(::Val{:topcenter}, x::T, ξ::T, Δx::T, Δξ::T, dip::T; filename::AbstractString="temp.msh") where T
    @gmsh_do begin
        model = gmsh.model
        factory = model.geo
        y, z = -ξ * cosd(dip), -ξ * sind(dip)
        # same as counting length of `range` in `mesh_downdip`
        nξ = round(Int, ξ / Δξ)
        # same as counting length of `range` in `mesh_strike`
        nx = round(Int, x / Δx)
        @addPoint begin
            -x/2, 0.0, 0.0, 1
            x/2, 0.0, 0.0, 2
            x/2, y, z, 3
            -x/2, y, z, 4
        end
        @addLine begin
            # top left -> top right -> bottom right (reverse) -> bottom left (reverse) -> top left
            1, 2, 1
            2, 3, 2
            4, 3, 3
            1, 4, 4
        end
        factory.addCurveLoop([1, 2, -3, -4], 1)
        factory.addPlaneSurface([1], 1)
        factory.mesh.setTransfiniteSurface(1, "Left", [1, 2, 3, 4])
        @setTransfiniteCurve begin
            1, nx+1, "Progression", 1.0
            2, nξ+1, "Progression", 1.0
            3, nx+1, "Progression", 1.0
            4, nξ+1, "Progression", 1.0
        end
        factory.mesh.setRecombine(2, 1)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(filename)
    end
end
