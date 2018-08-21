using Test

@testset "Fault creation" begin
    @test_throws ErrorException("Fautl dip angle ∈ [0, π/2].") fault(NormalFault, -1, (1, 2))
    @test_throws ErrorException("Fault domain must be larger than zero.") fault(ThrustFault, 2, (-1, 2))
    @test_throws ErrorException("Fault domain should be given by along-downdip [and along-strike] length[s].") fault(StrikeSlipFault, 2, (1, 2, 3))
end

@testset "Equivalent initialization" begin
    fd1 = fault(NormalFault, 1, (1,))
    fd2 = fault(NormalFault, 1, 1)
    fd3 = fault(NormalFault, 1, [1])
    @test fd1 == fd2
    @test fd2 == fd3
end

@testset "Get index" begin
    @testset "one dim" begin
        fd = fault(NormalFault, 1, 1)
        @test fd.span[1] == fd[:ξ]
    end
    @testset "two dim" begin
        fd = fault(NormalFault, 1, (1, 2))
        @test fd.span[1] == fd[:x]
        @test fd.span[2] == fd[:ξ]
    end
end
