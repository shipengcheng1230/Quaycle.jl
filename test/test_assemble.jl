using Test

@testset "Okada Assemble" begin
    @testset "1D fault" begin
        fa = fault(Val(:topcenter), STRIKING(), 10., 2.0)
        p = ElasticRSFProperties([rand(fa.mesh.nξ) for _ in 1: 4]..., rand(6)...)
        u0 = rand(fa.mesh.nξ, 2)
        prob, = assemble(fa, p, u0, (0., 1.0))
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end

    @testset "2D fault" begin
        fa = fault(Val(:topcenter), STRIKING(), 10., 10., 2., 2.)
        p = ElasticRSFProperties([rand(fa.mesh.nx, fa.mesh.nξ) for _ in 1: 4]..., rand(6)...)
        u0 = rand(fa.mesh.nx, fa.mesh.nξ, 2)
        prob, = assemble(fa, p, u0, (0., 1.0); buffer_ratio=1)
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end
end
