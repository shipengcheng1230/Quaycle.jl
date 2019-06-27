using Test
using JuEQ: fault, SBarbotTet4MeshEntity, SBarbotHex8MeshEntity

@testset "Okada Fault Space" begin

    @testset "Okada Line Fault" begin
        ft = fault(Val(:LineOkada), 100.0, 2.0, 33.0)
        @test ft.mesh.nξ == 50
    end

    @testset "Okada Rect Fault" begin
        ft = fault(Val(:RectOkada), 100.0, 50.0, 2.0, 2.0, 33.0)
        @test ft.mesh.nx == 50
        @test ft.mesh.nξ == 25
    end

    @testset "compose" begin
        ft = fault(Val(:LineOkada), 100.0, 2.0, 33.0)
        me = SBarbotTet4MeshEntity([rand(5) for _ in 1: 3]..., [[rand(3) for _ in 1: 5] for _ in 1: 4]..., rand(5))
        fas = compose(ft, me)
        @test fas.mv.tag |> length == 5

        ft = fault(Val(:RectOkada), 100.0, 50.0, 2.0, 2.0, 33.0)
        me = SBarbotHex8MeshEntity([rand(5) for _ in 1: 9]..., rand(), rand(5))
        fas = compose(ft, me)
        @test fas.mv.tag |> length == 5
    end
end
