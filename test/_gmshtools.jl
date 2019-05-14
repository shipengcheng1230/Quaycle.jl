using Test
using GmshTools

# This file currently isn't included in `runtests.jl`

gen_gmsh_mesh(Val(:topcenter), 100.0, 50.0, 13.0, 7.0, 45.0)

@testset "Gmsh Top Center Line" begin
    filename = tempname() * ".msh"
    gen_gmsh_mesh(Val(:topcenter), 100.0, 20.0, 45.0; filename=filename)
    els1d = @gmsh_open filename begin
        gmsh.model.mesh.getElements(1)
    end
    @test els1d[1][1] == 1 # 2-node line
    @test length(els1d[2][1]) == round(Int, 100.0 / 20.0)
    rm(filename)
end

@testset "Gmsh Top Center Line" begin
    filename = tempname() * ".msh"
    gen_gmsh_mesh(Val(:topcenter), 100.0, 50.0, 20.0, 10.0, 45.0; filename=filename)
    els2d = @gmsh_open filename begin
        gmsh.model.mesh.getElements(2)
    end
    @test els2d[1][1] == 3 # 4-node quadrangle
    @test length(els2d[2][1]) == round(Int, 100.0 / 20.0 * 50.0 / 10.0)
    rm(filename)
end
