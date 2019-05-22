using Test
using GmshTools
using LinearAlgebra

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

@testset "Gmsh Box by Extrude" begin
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

@testset "indice to tag" begin
    @testset "LineOkada" begin
        fname = tempname() * ".msh"
        gen_gmsh_mesh(Val(:LineOkada), 100.0, 10.0, 90.0; filename=fname)
        mesh = gen_mesh(Val(:LineOkada), 100.0, 10.0, 90.0)
        i2t = indice2tag(mesh, fname)
        @gmsh_open fname begin
            i = rand(1: mesh.nξ)
            @test gmsh.model.mesh.getElementByCoordinates(0.0, mesh.y[i], mesh.z[i], 1)[1] == i2t[i]
        end
    end
    @testset "RectOkada" begin
        fname = tempname() * ".msh"
        gen_gmsh_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 90.0; filename=fname)
        mesh = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 90.0)
        i2t = indice2tag(mesh, fname)
        @gmsh_open fname begin
            i, j = rand(1: mesh.nx), rand(1: mesh.nξ)
            @test gmsh.model.mesh.getElementByCoordinates(mesh.x[i], mesh.y[j], mesh.z[j], 2)[1] == i2t[i,j]
        end
    end
end

@testset "Combined Okada Rect and Box" begin
    fname = tempname() * ".msh"
    mf = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 90.0)
    rfzn = ones(Float64, 10)
    rfzh = accumulate((x, y) -> x * y, fill(1.5, size(rfzn))) |> cumsum
    normalize!(rfzh, Inf)
    gen_gmsh_mesh(mf, Val(:BoxHexExtrudeFromSurface), -50.0, -50.0, -60.0, 100.0, 100.0, 100.0, 10, 10, 4.0, 1.0, rfzn, rfzh; filename=fname)
    els = @gmsh_open fname begin
        gmsh.model.mesh.getElements(2)
    end
    @gmsh_open fname begin
        entag = gmsh.model.getEntitiesForPhysicalGroup(2, 1)
        @test gmsh.model.mesh.getElements(2, entag[1])[2][1] |> length == 5 * 10
    end
    @gmsh_open fname begin
        entag = gmsh.model.getEntitiesForPhysicalGroup(3, 1)
        @test gmsh.model.mesh.getElements(3, entag[1])[2][1] |> length == 10^3
    end
    rm(fname)
end

@testset "SBarbot Hex8 mesh entities" begin
    fname = tempname() * ".msh"
    rfzn = ones(Int, 7)
    rfzh = accumulate((x, y) -> x * y, fill(2.0, length(rfzn))) |> cumsum
    normalize!(rfzh, Inf)
    gen_gmsh_mesh(Val(:BoxHexExtrudeFromSurface), -50.0, -50.0, -60.0, 100.0, 100.0, 100.0, 20, 25, 4.0, 5.0, rfzn, rfzh; filename=fname)
    mc = read_gmsh_mesh(Val(:SBarbotHex8), fname; phytag=-1)
    fround = x -> round(x; digits=3)
    unique(fround, mc.L) |> length == 20 / 2
    unique(fround, mc.W) |> length == 7
    unique(fround, mc.T) |> length == 25 ÷ 2 + 1
    rm(fname)
end