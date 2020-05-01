using Test
using GmshTools
using Gmsh_SDK_jll
using LinearAlgebra
using Logging

using Quaycle: indice2tag, geo_okada_rect, geo_box_extruded_from_surfaceXY,
    _get_all_elements_in_physical_group, _get_entity_tags_in_physical_group,
    _get_centers, _get_element_node_num

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
            @test gmsh.model.mesh.getElementByCoordinates(0.0, mesh.y[i], mesh.z[i])[1] == i2t[i]
        end
        rm(fname)
    end
    @testset "RectOkada" begin
        fname = tempname() * ".msh"
        gen_gmsh_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 90.0; filename=fname)
        mesh = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 90.0)
        i2t = indice2tag(mesh, fname)
        @gmsh_open fname begin
            i, j = rand(1: mesh.nx), rand(1: mesh.nξ)
            @test gmsh.model.mesh.getElementByCoordinates(mesh.x[i], mesh.y[j], mesh.z[j])[1] == i2t[i,j]
        end
        rm(fname)
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
    mf = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 90.0)
    nx, ny, nz = 20, 25, 7
    fname = tempname() * ".msh"
    rfzn = ones(Int, nz)
    rfzh = accumulate((x, y) -> x * y, fill(2.0, length(rfzn))) |> cumsum
    normalize!(rfzh, Inf)
    gen_gmsh_mesh(mf, Val(:BoxHexExtrudeFromSurface), -50.0, -50.0, -60.0, 100.0, 100.0, 100.0, nx, ny, 4.0, 5.0, rfzn, rfzh;
        filename=fname, faulttag=(100, "fault"), asthenospheretag=(999, "asthenosphere"))
    mc = read_gmsh_mesh(Val(:SBarbotHex8), fname; phytag=999, reverse=false)
    fround = x -> round(x; digits=3)
    @test unique(fround, mc.L) |> length == ceil(ny / 2)
    @test unique(fround, mc.W) |> length == nz
    @test unique(fround, mc.T) |> length == ceil(nx / 2)
    @test unique(fround, mc.q1) |> length == ny
    @test unique(fround, mc.q2) |> length == nx
    @test unique(fround, mc.q3) |> length == nz

    # check unmatched read
    mc = read_gmsh_mesh(Val(:SBarbotHex8), fname; phytag=-1, reverse=true)
    @test length(unique(fround, mc.q1)) == length(unique(fround, mc.L)) * length(unique(fround, mc.x1))
    @test length(unique(fround, mc.q2)) == length(unique(fround, mc.x2))

    mc = read_gmsh_mesh(Val(:SBarbotHex8), fname; phytag=-1, reverse=false, rotate=90.0)
    @test length(unique(fround, mc.q2)) == length(unique(fround, mc.x2)) * length(unique(fround, mc.L))
    @test length(unique(fround, mc.q1)) == length(unique(fround, mc.x1))

    @test_throws AssertionError begin
        read_gmsh_mesh(Val(:SBarbotHex8), fname; phytag=-1, reverse=true, check=true)
    end
    @test_throws AssertionError begin
        read_gmsh_mesh(Val(:SBarbotHex8), fname; rotate=90.0, phytag=-1, reverse=false, check=true)
    end
    rm(fname)
end

@testset "read Tri3 mesh" begin
    filename = tempname() * ".msh"
    @testset "perpendicular to z" begin
        @gmsh_do begin
            @addPoint begin
                0.0, 0.0, 0.0, 0.0, 1
                1.0, 0.0, 0.0, 0.0, 2
                1.0, 1.0, 0.0, 0.0, 3
                0.0, 1.0, 0.0, 0.0, 4
            end
            @addLine begin
                1, 2, 1
                2, 3, 2
                3, 4, 3
                4, 1, 4
            end
            gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
            gmsh.model.geo.addPlaneSurface([1], 1)
            gmsh.model.addPhysicalGroup(2, [1], 99)
            gmsh.model.setPhysicalName(2, 99, "FAULT")
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)
            gmsh.write(filename)
        end

        mf = with_logger(NullLogger()) do
            read_gmsh_mesh(Val(:TDTri3), filename; phytag=99)
        end
        @test map(i ->
            cross(mf.B[i] - mf.A[i], mf.C[i] - mf.A[i]) × mf.ts[i] |> norm ≈ 0 &&
            all(isnan.(mf.ss[i])) && all(isnan.(mf.ds[i])) && mf.ts[i] ≈ [0, 0, 1], eachindex(mf.tag)) |> all
    end
    @testset "parallel to z" begin
        @gmsh_do begin
            @addPoint begin
                0.0, 0.0, 0.0, 0.0, 1
                1.0, 0.0, 0.0, 0.0, 2
                1.0, 0.0, -1.0, 0.0, 3
                0.0, 0.0, -1.0, 0.0, 4
            end
            @addLine begin
                1, 2, 1
                2, 3, 2
                3, 4, 3
                4, 1, 4
            end
            gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
            gmsh.model.geo.addPlaneSurface([1], 1)
            gmsh.model.addPhysicalGroup(2, [1], 99)
            gmsh.model.setPhysicalName(2, 99, "FAULT")
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)
            gmsh.write(filename)
        end

        mf = read_gmsh_mesh(Val(:TDTri3), filename; phytag=99)
        @test map(i ->
            cross(mf.B[i] - mf.A[i], mf.C[i] - mf.A[i]) × mf.ts[i] |> norm ≈ 0 &&
            mf.ss[i] ≈ [1, 0, 0] && mf.ds[i] ≈ [0, 0, 1] && mf.ts[i] ≈ [0, -1, 0], eachindex(mf.tag)) |> all
    end
    @testset "arbitrary plane" begin
        l, w = 1.5, 0.75
        dip, strike = 30.0, 10.0
        depth = 2.0

        p1 = [l/2 * sind(strike), l/2 * cosd(strike), -depth]
        p2 = [-p1[1], -p1[2], -2.0]
        p4 = [p1[1] + w * cosd(dip) * cosd(strike), p1[2] - w * cosd(dip) * sind(strike), p1[3] - w * sind(dip)]
        p3 = [p4[1] - l * sind(strike), p4[2] - l * cosd(strike), p4[3]]

        # true solution
        A = hcat(p1, p2, p3)'
        a, b, c = A \ [-1; -1; -1]
        ts = (p2 - p1) × (p3 - p1) |> normalize!
        ss = normalize([a, -b, 0.0])
        ds = ts × ss |> normalize

        @gmsh_do begin
            @addPoint begin
                p1..., 0.0, 1
                p2..., 0.0, 2
                p3..., 0.0, 3
                p4..., 0.0, 4
            end
            @addLine begin
                1, 2, 1
                2, 3, 2
                3, 4, 3
                4, 1, 4
            end
            gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
            gmsh.model.geo.addPlaneSurface([1], 1)
            gmsh.model.addPhysicalGroup(2, [1], 99)
            gmsh.model.setPhysicalName(2, 99, "FAULT")
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)
            gmsh.write(filename)
        end

        mf = read_gmsh_mesh(Val(:TDTri3), filename; phytag=99)
        @test map(i -> mf.ts[i] ≈ ts && mf.ss[i] ≈ ss && mf.ds[i] ≈ ds, 1: length(mf.tag)) |> all

        # reverse the curve loop
        @gmsh_do begin
            @addPoint begin
                p1..., 0.0, 1
                p2..., 0.0, 2
                p3..., 0.0, 3
                p4..., 0.0, 4
            end
            @addLine begin
                1, 2, 1
                2, 3, 2
                3, 4, 3
                4, 1, 4
            end
            gmsh.model.geo.addCurveLoop([4, 3, 2, 1], 1)
            gmsh.model.geo.addPlaneSurface([1], 1)
            gmsh.model.addPhysicalGroup(2, [1], 99)
            gmsh.model.setPhysicalName(2, 99, "FAULT")
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)
            gmsh.write(filename)
        end

        mf = read_gmsh_mesh(Val(:TDTri3), filename; phytag=99)
        @test map(i -> mf.ts[i] ≈ ts && mf.ss[i] ≈ ss && mf.ds[i] ≈ ds, eachindex(mf.tag)) |> all
    end
    rm(filename)
end

@testset "multiple entities within one physical group" begin
    filename = tempname() * ".msh"
    mf1 = gen_mesh(Val(:RectOkada), 80.0e3, 10.0e3, 5e3, 5e3, 90.0)
    mf2 = gen_mesh(Val(:RectOkada), 80.0e3, 10.0e3, 5e3, 5e3, 30.0)
    rfzn = ones(5)
    rfzh = accumulate((x, y) -> x * y, fill(1.0, length(rfzn))) |> cumsum |> x -> normalize(x, Inf)
    llx1, lly1, llz1 =  -50.0e3, -10.0e3, -10.0e3
    llx2, lly2, llz2 =  -50.0e3, -10.0e3, -50.0e3
    dx, dy, dz = 100.0e3, 20.0e3, 10e3
    rfx, rfy = 1.0, 1.0
    nx, ny = 25, 10
    faulttag = (100, "fault")
    asthenospheretag = (999, "asthenosphere")

    @gmsh_do begin
        x, ξ = mf1.nx * mf1.Δx,  mf1.nξ * mf1.Δξ
        _reg1 = geo_okada_rect(x, ξ * cosd(mf1.dip), ξ * sind(mf1.dip), mf1.nx, mf1.nξ, 1)
        x, ξ = mf2.nx * mf2.Δx,  mf2.nξ * mf2.Δξ
        _reg2 = geo_okada_rect(x, ξ * cosd(mf2.dip), ξ * sind(mf2.dip), mf2.nx, mf2.nξ, _reg1)
        gmsh.model.addPhysicalGroup(2, [_reg1-1, _reg2-1], faulttag[1])
        gmsh.model.setPhysicalName(2, faulttag...)
        _reg3 = geo_box_extruded_from_surfaceXY(llx1, lly1, llz1, dx, dy, dz, nx, ny, rfx, rfy, rfzn, rfzh, _reg2)
        volumetag1 = _reg3[findfirst(x -> x[1] == 3, _reg3)][2]
        _reg4 = maximum([x[2] for x in _reg3]) + 1
        _reg5 = geo_box_extruded_from_surfaceXY(llx2, lly2, llz2, dx, dy, dz, nx, ny, rfx, rfy, rfzn, rfzh, _reg4)
        volumetag2 = _reg5[findfirst(x -> x[1] == 3, _reg5)][2]
        gmsh.model.addPhysicalGroup(3, [volumetag1, volumetag2], asthenospheretag[1])
        gmsh.model.setPhysicalName(3, asthenospheretag...)
        @addOption begin
            "Mesh.SaveAll", 1 # the mesh is incorrect without this
        end
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.write(filename)
    end

    e1 = _get_all_elements_in_physical_group(filename, 2, faulttag[1])
    entag1=  _get_entity_tags_in_physical_group(filename, 2, faulttag[1])
    c1 = _get_centers(filename, e1[1][1], entag1)
    nnode1 = _get_element_node_num(e1[1][1])

    e2 = _get_all_elements_in_physical_group(filename, 3, asthenospheretag[1])
    entag2=  _get_entity_tags_in_physical_group(filename, 3, asthenospheretag[1])
    c2 = _get_centers(filename, e2[1][1], entag2)
    nnode2 = _get_element_node_num(e2[1][1])

    @test length(e1[2][1]) == mf1.nx * mf1.nξ + mf2.nx * mf2.nξ
    @test length(e1[3][1]) == length(e1[2][1]) * nnode1
    @test length(c1) == length(e1[2][1]) * 3

    @test length(e2[2][1]) == nx * ny * length(rfzn) * 2
    @test length(e2[3][1]) == length(e2[2][1]) * nnode2
    @test length(c2) == length(e2[2][1]) * 3

    rm(filename)
end

@testset "paraview cache" begin
    tmp = tempname() * ".msh"

    mf = gen_mesh(Val(:LineOkada), 80.0e3, 1e3, 45.0)
    gen_gmsh_mesh(mf; filename=tmp)
    vcache = gmsh_vtk_output_cache(tmp, mf, -1)
    @test vcache.dat |> length == mf.nξ
    @test vcache.tag |> length == mf.nξ

    mf = gen_mesh(Val(:RectOkada), 80.0e3, 10.0e3, 1e3, 1e3, 30.0)
    gen_gmsh_mesh(mf; filename=tmp)
    vcache = gmsh_vtk_output_cache(tmp, mf, -1)
    @test vcache.dat |> length == mf.nx * mf.nξ
    @test size(vcache.pts, 2) == (mf.nx + 1) * (mf.nξ + 1)

    vcache = gmsh_vtk_output_cache(tmp, 2, -1)
    @test vcache.cells |> length == mf.nx * mf.nξ
    @test size(vcache.pts, 2) == (mf.nx + 1) * (mf.nξ + 1)
    @test vcache.tag |> length == mf.nx * mf.nξ

    rm(tmp)
end

@testset "InPlaneX with Okada Rect" begin
    llx, lly, llz, dx, dy, dz, nx, nv = -40., 0.0, -50., 80., 0.0, 20., 8, 4
    rfxstr, rfvstr = "Bump", "Progression"
    rfx, rfv = 1.2, 1.2
    filename = "temp.msh"
    mf = gen_mesh(Val(:RectOkada), 100.0, 20.0, 5.0, 5.0, 90.0)

    gen_gmsh_mesh(mf, Val(:InPlaneX), llx, lly, llz, dx, dy, dz, nx, nv;
        rfxstr=rfxstr, rfx=rfx, rfvstr=rfvstr, rfv=rfv, filename=filename,
        faulttag=(15, "fault"), asthenospheretag=(23, "asthenosphere"))

    @gmsh_open filename begin
        phytag = gmsh.model.getPhysicalGroups(2)
        @test length(phytag) == 2
        @test phytag[1][1] == phytag[2][1] # 2D
        @test phytag[1][2] == 15
        @test phytag[2][2] == 23

        entag1 = gmsh.model.getEntitiesForPhysicalGroup(2, phytag[1][2])
        els1 = gmsh.model.mesh.getElements(2, entag1[1])
        @test els1[1][1] == 3
        @test length(els1[2][1]) == 20 * 4

        entag2 = gmsh.model.getEntitiesForPhysicalGroup(2, phytag[2][2])
        els2 = gmsh.model.mesh.getElements(2, entag2[1])
        @test els2[1][1] == 3
        @test length(els2[2][1]) == 8 * 4
    end
    rm(filename)
end

@testset "InPlaneX read gmsh" begin
    llx, lly, llz, dx, dy, dz, nx, nv = -40., 0.0, -50., 80., 0.0, 20., 5, 6
    rfxstr, rfvstr = "Bump", "Progression"
    rfx, rfv = 1.0, 1.0
    filename = tempname() * ".msh"
    mf = gen_mesh(Val(:RectOkada), 100.0, 20.0, 5.0, 5.0, 90.0)

    gen_gmsh_mesh(mf, Val(:InPlaneX), llx, lly, llz, dx, dy, dz, nx, nv;
        rfxstr=rfxstr, rfx=rfx, rfvstr=rfvstr, rfv=rfv, filename=filename,
        faulttag=(15, "fault"), asthenospheretag=(23, "asthenosphere"))

    ma = read_gmsh_mesh(Val(:InPlaneX), filename; phytag=23)
    @test unique(x -> round(x; digits=3), ma.W)[1] ≈ dx / nx
    @test unique(x -> round(x; digits=3), ma.T)[1] ≈ hypot(dy, dz) / nv
    rm(filename)
end
