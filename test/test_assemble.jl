using Test

@testset "Okada Assemble" begin
    @testset "1D fault" begin
        mesh = gen_mesh(Val(:LineOkada), 10., 2.0, 45.0)
        gf = okada_stress_gf_tensor(mesh, 1.0, 1.0, DIPPING())
        p = ElasticRSFProperty([rand(mesh.nξ) for _ in 1: 4]..., rand(6)...)
        u0 = ArrayPartition([rand(mesh.nξ) for _ in 1: 2]...)
        prob = assemble(mesh, gf, p, u0, (0., 1.0); flf=CForm())
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end

    @testset "2D fault" begin
        mesh = gen_mesh(Val(:RectOkada), 10., 10., 2., 2., 90.)
        gf = okada_stress_gf_tensor(mesh, 1.0, 1.0, STRIKING(); buffer_ratio=1.0)
        p = ElasticRSFProperty([rand(mesh.nx, mesh.nξ) for _ in 1: 4]..., rand(6)...)
        u0 = ArrayPartition([rand(mesh.nx, mesh.nξ) for _ in 1: 2]...)
        prob = assemble(mesh, gf, p, u0, (0., 1.0))
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end
end
