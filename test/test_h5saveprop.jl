using Test
using HDF5

@testset "Property" begin

    @testset "sdof rsf Property" begin
        tmpfile = tempname()
        p = SingleDofRSFProperty(rand(9)...)
        @save_prop tmpfile p
        p2 = @read_prop tmpfile
        @test p == p2
        rm(tmpfile)
    end

    @testset "elastic rsf Property" begin
        tmpfile = tempname()
        p = ElasticRSFProperty([rand(9) for _ in 1: 4]..., rand(6)...)
        @save_prop tmpfile p
        p2 = @read_prop tmpfile
        @test p == p2
        rm(tmpfile)
    end

end
