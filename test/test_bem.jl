@testset "Fault uniform FD grids" begin
    fa = fault(ThrustFault, 10., (60., 30.))

    @testset "Interval corretion" begin
        gd = discretize(fa, fa[:x]/128, fa[:ξ]/64)
        @test abs(gd.x[2] - gd.x[1]) ≈ fa[:x] / 128
        @test abs(gd.ξ[2] - gd.ξ[1]) ≈ fa[:ξ] / 64
    end
end
