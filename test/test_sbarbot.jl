
using Test
using DelimitedFiles
using Base.Iterators

@testset "SBarbot" begin
    epsv11 = 11e-6
    epsv12 = 5e-6
    epsv13 = 6e-6
    epsv22 = 7e-6
    epsv23 = 8e-6
    epsv33 = 13e-6

    L = 40e3
    W = 20e3
    T = 10e3
    q1 = -10e3
    q2 = 10e3
    q3 = 1e3
    theta = 40.0
    G = 1.0
    nu = 0.25

    x1s = range(-10000.0, stop=10000.0, step=1000.0)
    x2s = range(-10000.0, stop=10000.0, step=1000.0)
    x3s = range(-30.0, stop=-20.0, step=5.0)
    xxs = product(x1s, x2s, x3s)

    @testset "Displacement" begin
        u_truth = readdlm(joinpath(@__DIR__, "data/test_sbarbot_disp.dat"), ' ', Float64)
        fu = (x) -> sbarbot_disp_quad4(x..., q1, q2, q3, L, T, W, theta, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, G, nu)
        u_cal = map(fu, xxs) |> vec
        ftest = (i) -> u_cal[i] ≈ u_truth[i,:]
        @test map(ftest, 1: length(u_cal)) |> all
    end

    @testset "Stress as well as Singularity" begin
        fu = (x) -> sbarbot_stress_quad4(x..., q1, q2, q3, L, T, W, theta, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, G, nu)
        u_cal = map(fu, xxs) |> vec
        function funtest(i::Integer)
            if all(map(isnan, u_cal[i])) && all(map(isnan, u_truth[i,:]))
                return true
            else
                return u_cal[i] ≈ u_truth[i,:]
            end
        end
        @test map(funtest, 1: length(u_cal)) |> all
    end
end
