using DelimitedFiles

@testset "Test Okada's displacement" begin
    @testset "Normal computations" begin
        u_truth = readdlm(joinpath(@__DIR__, "test_okada.dat"), ',', Float64)
        count = 1
        for x = 10.: 15., y = 20.: 30., z = -40.: -30.
            u = dc3d_okada(x, y, z, 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
            @test isapprox(u_truth[count, :], u, rtol=1e-8)
            count += 1
        end

    end

    @testset "Negative depth" begin
        u = dc3d_okada(10., 20., 30., 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        @test isapprox(zeros(Float64, 12), u, rtol=1e-8)
    end

    @testset "Singularity" begin
        u_singular = zeros(Float64, 12)

        u = dc3d_okada(-80., 0., -50., 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        @test isapprox(u_singular, u, atol=1e-8)

        u = dc3d_okada(-70., 0., -10., 2. /3, 10., 90., [-80., 120.], [-20., 20.], [200., -150., 100.])
        @test isapprox(u_singular, u, atol=1e-8)

        u = dc3d_okada(-80., 20., -50., 2. /3, 50., 70., [-80., 120.], [-30., 25.], [200., -150., 100.])
        u_truth = [
            -32.60576975, 59.39436703, 1.61816100, -0.43682428, 1.66503564, 0.27751175,
            0.64126234, -1.05836342, 0.03630938, -0.16974972, 0.61800758, 0.22842501]
        @test isapprox(u_truth, u, rtol=1e-8)


        u = dc3d_okada(-10., 20., -50., 2. /3, 50., 90., [-80., 120.], [-30., 25.], [200., -150., 100.])
        u_truth = [
            -56.04511679, 52.38458483, 40.87786501, -0.10140811, -0.09099896, 0.02371607,
            1.60327214, -0.63369117, -1.13609685, 0.05783881, 0.88665943, 0.21946899,
        ]
        @test isapprox(u_truth, u, rtol=1e-8)

        u = dc3d_okada(-10., 20., 0., 2. /3, 0., 0., [-80., 120.], [-30., 25.], [200., -150., 100.])
        # take a look: https://github.com/JuliaLang/julia/issues/23376
        @test isapprox(u_singular, u, atol=1e-8)

        u = dc3d_okada(-100., 20., 0., 2. /3, 0., 0., [-80., -90.], [20., 25.], [200., -150., 100.])
        @test isapprox(u_singular, u, atol=1e-8)

        u = dc3d_okada(-80., 20., 0., 2. /3, 0., 0., [-80., 100.], [25., 35.], [200., -150., 100.])
        @test isapprox(u_singular, u, atol=1e-8)
    end
end
