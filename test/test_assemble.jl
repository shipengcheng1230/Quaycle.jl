using Test
using InteractiveUtils

@testset "Okada Assemble" begin

    @testset "1D fault" begin
        fa = fault(Val(:CSFS), STRIKING(), 10., 2.0)
        frprop = init_friction_prop(fa)
        faprop = init_fault_prop(rand(6)...)
        u0 = rand(fa.mesh.nξ, 3)
        prob = assemble(Val(:okada), fa, faprop, frprop, u0, (0., 1.0))
        du = similar(u0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end

    @testset "2D fault" begin
        fa = fault(Val(:CSFS), STRIKING(), (10., 10.), (2., 2.))
        frprop = init_friction_prop(fa)
        faprop = init_fault_prop(rand(6)...)
        u0 = rand(fa.mesh.nx, fa.mesh.nξ, 3)
        prob = assemble(Val(:okada), fa, faprop, frprop, u0, (0., 1.0), buffer_ratio=1)
        du = similar(u0)
        ret = @code_typed prob.f(du, u0, prob.p, 1.0)
        @inferred prob.f(du, u0, prob.p, 1.0)
    end

end
