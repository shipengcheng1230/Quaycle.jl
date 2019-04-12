using Test

@testset "trait function for unit dislocation" begin
    @test unit_dislocation(DIPPING()) == [0.0, 1.0, 0.0]
    @test unit_dislocation(STRIKING()) == [1.0, 0.0, 0.0]
end
