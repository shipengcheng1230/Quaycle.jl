using Test
using HDF5

@testset "Property" begin
    function test_equal(p)
        tmpfile = tempname()
        @save_prop tmpfile p
        p′ = @read_prop tmpfile
        @test p == p′
    end

    ps = [
        SingleDofRSFProperty(rand(9)...),
        ElasticRSFProperty([rand(9) for _ in 1: 4]..., rand(6)...),
        DislocationCreepProperty([rand(5) for _ in 1: 6]...),
        DiffusionCreepProperty([rand(5) for _ in 1: 7]...),
        ]
    map(test_equal, ps)
end
