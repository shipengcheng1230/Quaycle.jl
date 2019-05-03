using Test

@testset "trait function for unit dislocation" begin
    @test JuEQ.unit_dislocation(DIPPING()) == [0.0, 1.0, 0.0]
    @test JuEQ.unit_dislocation(STRIKING()) == [1.0, 0.0, 0.0]
end

mesh = gen_mesh(Val(:topcenter), 100.0, 100.0, 2.0, 2.0, 44.0)
greens_tensor(Val(:okada), mesh, 1.0, 1.0, DIPPING(); buffer_ratio=1)
