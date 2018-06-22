using Base.Test
using JuEQ

@testset "Test Okada's displacement" begin
    @testset "Normal computations" begin
        u_truth = readdlm(joinpath(@__DIR__, "test_okada.dat"), ',', Float64)
        count = 1
        for x = 10.: 15., y = 20.: 30., z = -40.: -30.
            u = dc3d_okada(x, y, z, 2./3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
            @test isapprox(u_truth[count, :], u, rtol=1e-8)
            count += 1
        end
    end

    @testset "Negative depth" begin
        u = dc3d_okada(10., 20., 30., 2./3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        @test isapprox(zeros(Float64, 12), u, rtol=1e-8)
    end

    @testset "Singularity" begin
        u = dc3d_okada(-80., 0., -50., 2./3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        @test isapprox(zeros(Float64, 12), u, rtol=1e-8)
    end
end
