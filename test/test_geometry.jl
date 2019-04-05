using Test

@testset "Simple Mesh Generator" begin

    @testset "Line" begin
        geo = SimpleLine((0., 0., 0.), (50,))
        msh = mesh(Val(:SimpleGrid), geo, 2.)
        @test msh.nξ == 25
    end

    @testset "Rect" begin
        geo = SimpleRect((0, 0, 0,), (50, 50))
        msh = mesh(Val(:SimpleGrid), geo, 3, 4)
        @test msh.nx == 16
        @test msh.nξ == 12
    end

end
