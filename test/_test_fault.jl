using Test

@testset "Fault creation" begin
    @test_throws AssertionError("Fault dip angle: -1.0 received, must ∈ [0, 90].") fault(NormalFault, -1, (1, 2))
    @test_throws AssertionError("Fault domain: (-1.0, 2.0) received, must > 0.") fault(ThrustFault, 2, (-1, 2))
    @test_throws AssertionError("Fault domain dim: 3 received, must be `1` (along-downdip) or `2` (plus along-strike).") fault(StrikeSlipFault, 2, (1, 2, 3))
    @test_throws ErrorException("Dip angle: 10.0 received, for strike-slip faults must be `90`.") fault(StrikeSlipFault, 10, 1)
end

@testset "Equivalent initialization" begin
    @testset "Default construction" begin
        fa1 = fault(NormalFault, 1, (1,))
        fa2 = fault(NormalFault, 1, 1)
        fa3 = fault(NormalFault, 1, [1])
        @test fa1 == fa2
        @test fa2 == fa3
    end

    @testset "Vertical strike-slip fault" begin
        fa1 = fault(StrikeSlipFault, 1.)
        fa2 = fault(StrikeSlipFault, 90., 1.)
        @test fa1 == fa2
    end

    @testset "Heterogeneous input type" begin
        fa1 = fault(NormalFault, 1, (1.,))
        fa2 = fault(NormalFault, 1., (1.,))
        @test fa1 == fa2
    end
end

@testset "Get index" begin
    @testset "one dim" begin
        fa = fault(NormalFault, 1, 1)
        @test fa.span[1] == fa[:ξ]
    end
    @testset "two dim" begin
        fa = fault(NormalFault, 1, (1, 2))
        @test fa.span[1] == fa[:x]
        @test fa.span[2] == fa[:ξ]
    end
end
