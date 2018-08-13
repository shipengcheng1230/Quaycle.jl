using Test

@testset "Fault creation" begin
    @test_throws ErrorException("Fautl dip angle ∈ [0, π/2].") fault(NormalFault(), (1, 2), -1)
    @test_throws ErrorException("Fault domain must be larger than zero.") fault(ThrustFault(), (-1, 2), 2)
    @test_throws ErrorException("Fault domain should be given by along-downdip [and along-strike] length[s].") fault(StrikeSlipFault(), (1, 2, 3), 2)
end

@testset "Fault uniform FD grids" begin
    global fa = fault(ThrustFault(), (30., 60.), 10.)
    @testset "Generate error" begin
        @test_throws ErrorException("Grid numbers mismatch with fault domain.") discretize(fa, (64, 128, 20))
        @test_throws ErrorException("Grid numbers must be larger than 0.") discretize(fa, (-1, 128))
    end

    @testset "Interval corretion" begin
        global g = discretize(fa, (64, 128))
        @test abs(g.x[2] - g.x[1]) ≈ 60 / 128
        @test abs(g.ξ[2] - g.ξ[1]) ≈ 30 / 64
    end
end
