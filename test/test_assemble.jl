using Test
using GmshTools

@testset "Okada Assemble" begin
    @testset "1D fault" begin
        mesh = gen_mesh(Val(:LineOkada), 10., 2.0, 45.0)
        gf = stress_greens_func(mesh, 1.0, 1.0, DIPPING())
        p = RateStateQuasiDynamicProperty([rand(mesh.nξ) for _ in 1: 4]..., rand(4)...)
        u0 = ArrayPartition([rand(mesh.nξ) for _ in 1: 2]...)
        prob = assemble(gf, p, u0, (0., 1.0); flf=CForm())
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end

    @testset "2D fault" begin
        mesh = gen_mesh(Val(:RectOkada), 10., 10., 2., 2., 90.)
        gf = stress_greens_func(mesh, 1.0, 1.0, STRIKING(); buffer_ratio=1.0)
        p = RateStateQuasiDynamicProperty([rand(mesh.nx, mesh.nξ) for _ in 1: 4]..., rand(4)...)
        u0 = ArrayPartition([rand(mesh.nx, mesh.nξ) for _ in 1: 2]...)
        prob = assemble(gf, p, u0, (0., 1.0))
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end

    @testset "triangular fault mesh" begin
        filename = tempname() * ".msh"
        @gmsh_do begin
            reg = JuEQ.geo_rect_x(-40e3, 0.0, -10e3, 80e3, 0.0, 10e3, 1)
            gmsh.model.addPhysicalGroup(2, [reg-1], 99)
            gmsh.model.setPhysicalName(2, 99, "FAULT")
            @addOption begin
                "Mesh.CharacteristicLengthMax", 10000.0
                "Mesh.CharacteristicLengthMin", 10000.0
            end
            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)
            gmsh.write(filename)
        end
        mesh = read_gmsh_mesh(Val(:TDTri3), filename; phytag=99)
        gf = stress_greens_func(mesh, 1.0, 1.0, DIPPING())
        p = RateStateQuasiDynamicProperty([rand(length(mesh.tag)) for _ in 1: 4]..., rand(4)...)
        u0 = ArrayPartition([rand(length(mesh.tag)) for _ in 1: 2]...)
        prob = assemble(gf, p, u0, (0.0, 1.0))
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
        rm(filename)
    end

    @testset "2D fault 3D asthenosphere" begin
        mf = gen_mesh(Val(:RectOkada), 10., 10., 2., 2., 90.)
        rfzn = ones(2)
        rfzh = [0.5, 1.0]
        tmp = tempname() * ".msh"
        gen_gmsh_mesh(mf,
            Val(:BoxHexExtrudeFromSurface), -5.0, -5.0, -10.0, 10.0, 10.0, 20.0, 5, 5, 1.0, 1.2, rfzn, rfzh;
            filename=tmp)
        ma = read_gmsh_mesh(Val(:SBarbotHex8), tmp; phytag=1)
        gg = compose_stress_greens_func(mf, ma, 1.0, 1.0, STRIKING(), (:xx, :xy, :xz))
        pe = RateStateQuasiDynamicProperty([rand(mf.nx, mf.nξ) for _ in 1: 4]..., rand(4)...)
        pdisl = DislocationCreepProperty([rand(length(ma.tag)) for _ in 1: 10]...)
        pdiff = DiffusionCreepProperty([rand(length(ma.tag)) for _ in 1: 11]...)
        pc = compose(pe, rand(3), pdisl, pdiff)
        v0 = rand(mf.nx, mf.nξ)
        θ0 = rand(mf.nx, mf.nξ)
        ϵ0 = rand(length(ma.tag), 3)
        σ0 = rand(length(ma.tag), 6)
        u0 = ArrayPartition(v0, θ0, ϵ0, σ0)
        prob = assemble(gg, pc, u0, (0., 1.0))
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
        rm(tmp)
    end
end
