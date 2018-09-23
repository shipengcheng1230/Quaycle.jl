@testset "Fault creation" begin
    @test_throws AssertionError("Fault dip angle: -1 received, must ∈ [0, 90].") fault(NormalFault, -1, (1, 2))
    @test_throws AssertionError("Fault domain: (-1, 2) received, must > 0.") fault(ThrustFault, 2, (-1, 2))
    @test_throws AssertionError("Fault domain dim: 3 received, must be `1` (along-downdip) or `2` (plus along-strike).") fault(StrikeSlipFault, 2, (1, 2, 3))
    @test_throws ErrorException("Dip angle: 10 received, for strike-slip faults must be `90`.") fault(StrikeSlipFault, 10, 1)
end

@testset "Equivalent initialization" begin
    @testset "Default construction" begin
        fd1 = fault(NormalFault, 1, (1,))
        fd2 = fault(NormalFault, 1, 1)
        fd3 = fault(NormalFault, 1, [1])
        @test fd1 == fd2
        @test fd2 == fd3
    end

    @testset "Vertical strike-slip fault" begin
        fd1 = fault(StrikeSlipFault, 1.)
        fd2 = fault(StrikeSlipFault, 90., 1.)
        @test fd1 == fd2
    end

    @testset "Heterogeneous input type" begin
        fd1 = fault(NormalFault, 1, (1.,))
        fd2 = fault(NormalFault, 1., (1.,))
        @test fd1 == fd2
    end
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
