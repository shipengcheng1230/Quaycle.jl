using Test

@testset "Okada Fault Space" begin

    @testset "Okada Line Fault" begin
        ft = fault(Val(:LineOkada), DIPPING(), 100.0, 2.0, 33.0)
        @test ft.ft == DIPPING()
        @test ft.mesh.nξ == 50
    end

    @testset "Okada Rect Fault" begin
        ft = fault(Val(:RectOkada), DIPPING(), 100.0, 50.0, 2.0, 2.0, 33.0)
        @test ft.ft == DIPPING()
        @test ft.mesh.nx == 50
        @test ft.mesh.nξ == 25
    end

    @testset "Okada Striking Fault" begin
        ft = fault(Val(:LineOkada), STRIKING(), 100.0, 2.0)
        @test ft.mesh.dip == 90.0

        ft = fault(Val(:RectOkada), STRIKING(), 100.0, 50.0, 2.0, 2.0)
        @test ft.mesh.dip == 90.0
    end

end
