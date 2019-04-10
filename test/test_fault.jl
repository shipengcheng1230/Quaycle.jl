using Test

@testset "Fault Space" begin

    @testset "Line Fault" begin
        ft1 = fault(Val(:CSFS), STRIKING(), 90.0, 10.0, 1.0)
        ft2 = fault(Val(:CSFS), STRIKING(), 10, 1.0)
        @test ft1.mesh.aξ == ft2.mesh.aξ
    end

    @testset "Rect Fault" begin
        ft3 = fault(Val(:CSFS), STRIKING(), 90.0, (10.0, 5.0), (1.0, 0.5))
        ft4 = fault(Val(:CSFS), STRIKING(), 90.0, 10.0, 5.0, 1.0, 0.5)
        ft5 = fault(Val(:CSFS), STRIKING(), 90.0, 10, 5, 1, 0.5)
        ft6 = fault(Val(:CSFS), STRIKING(), 10.0, 5.0, 1.0, 0.5)
        @test ft3.mesh.aξ == ft4.mesh.aξ
        @test ft4.mesh.aξ == ft5.mesh.aξ
        @test ft5.mesh.aξ == ft6.mesh.aξ
    end

end
