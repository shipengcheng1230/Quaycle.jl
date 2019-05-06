using Test

@testset "properties" begin

    @testset "sdof rsf properties" begin
        tmpfile = tempname()
        p = SingleDofRSFProperties(rand(9)...)
        @save_prop tmpfile p
        p2 = @read_prop tmpfile
        @test p == p2
        rm(tmpfile)
    end

    @testset "elastic rsf properties" begin
        tmpfile = tempname()
        p = ElasticRSFProperties([rand(9) for _ in 1: 4]..., rand(6)...)
        @save_prop tmpfile p
        p2 = @read_prop tmpfile
        @test p == p2
        rm(tmpfile)
    end

end
