using Test
using GmshTools
using LinearAlgebra

# This file currently isn't included in `runtests.jl` but tested locally on MacOS

@testset "Gmsh Okada Line" begin
    filename = tempname() * ".msh"
    gen_gmsh_mesh(Val(:LineOkada), 100.0, 20.0, 45.0; filename=filename)
    els1d = @gmsh_open filename begin
        gmsh.model.mesh.getElements(1)
    end
    @test els1d[1][1] == 1 # 2-node line
    @test length(els1d[2][1]) == round(Int, 100.0 / 20.0)
    rm(filename)
end

@testset "Gmsh Okada Rect" begin
    filename = tempname() * ".msh"
    gen_gmsh_mesh(Val(:RectOkada), 100.0, 50.0, 20.0, 10.0, 45.0; filename=filename)
    els2d = @gmsh_open filename begin
        gmsh.model.mesh.getElements(2)
    end
    @test els2d[1][1] == 3 # 4-node quadrangle
    @test length(els2d[2][1]) == round(Int, 100.0 / 20.0 * 50.0 / 10.0)
    rm(filename)
end

@testset "Gmsh Box by Extrude" begin # may fail due to empty mesh, rerun it
    rfzn = ones(Float64, 10)
    rfzh = accumulate((x, y) -> x * y, fill(2.0, size(rfzn))) |> cumsum
    normalize!(rfzh, Inf)

    filename = tempname() * ".msh"
    gen_gmsh_mesh(Val(:BoxHexExtrudeFromSurface), -100.0, -100.0, 0.0, 200.0, 200.0, 100.0, 10, 10, 1.1, 5.0, rfzn, rfzh; filename=filename)
    el3d = @gmsh_open filename begin
        gmsh.model.mesh.getElements(3)
    end
    @test el3d[1][1] == 5 # 8-node hexahedron
    @test length(el3d[2][1]) == 10 * 10 * 10
    rm(filename)
end
