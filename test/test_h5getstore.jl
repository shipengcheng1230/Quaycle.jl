using Test
using HDF5

@testset "property save and read" begin
    function test_equal(p)
        tmpfile = tempname()
        @store tmpfile p
        p′ = @getprop tmpfile
        @test p == p′
        rm(tmpfile)
    end

    ps = [
        SingleDofRSFProperty(rand(9)...),
        RateStateQuasiDynamicProperty([rand(9) for _ in 1: 4]..., rand(4)...),
        DislocationCreepProperty([rand(5) for _ in 1: 10]...),
        DiffusionCreepProperty([rand(5) for _ in 1: 11]...),
        CompositePlasticDeformationProperty([rand(5) for _ in 1: 4]..., rand(3), rand(1:6, 3))
        ]
    map(test_equal, ps)
end
