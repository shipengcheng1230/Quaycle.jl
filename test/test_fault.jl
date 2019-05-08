using Test

@testset "Fault Space" begin

    @testset "Top Center Line Fault" begin
        ft = fault(Val(:topcenter), DIPPING(), 100.0, 2.0, 33.0)
        @test ft.ft == DIPPING()
        @test ft.mesh.nξ == 50
    end

    @testset "Top Center Rect Fault" begin
        ft = fault(Val(:topcenter), DIPPING(), 100.0, 50.0, 2.0, 2.0, 33.0)
        @test ft.ft == DIPPING()
        @test ft.mesh.nx == 50
        @test ft.mesh.nξ == 25
    end

    @testset "Top Center Striking Fault" begin
        ft = fault(Val(:topcenter), STRIKING(), 100.0, 2.0)
        @test ft.mesh.dip == 90.0

        ft = fault(Val(:topcenter), STRIKING(), 100.0, 50.0, 2.0, 2.0)
        @test ft.mesh.dip == 90.0
    end

end
