using Test
using DelimitedFiles
using Base.Iterators
using FastGaussQuadrature

# Those corresponding verifiable data are obtained from orginal matlab functions at
# https://bitbucket.org/sbarbot/bssa-2016237/src/master/
# https://bitbucket.org/sbarbot/bssa-2018058/src/default/
@testset "SBarbot Hex8" begin
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
        u_truth = readdlm(joinpath(@__DIR__, "data/test_sbarbot_disp_hex8.dat"), ' ', Float64)
        fu = (x) -> sbarbot_disp_hex8(x..., q1, q2, q3, L, T, W, theta, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, G, nu)
        u_cal = map(fu, xxs) |> vec
        ftest = (i) -> u_cal[i] ≈ u_truth[i,:]
        @test map(ftest, 1: length(u_cal)) |> all
    end

    @testset "Stress as well as Singularity" begin
        epsv11 = 11e-6
        u_truth = readdlm(joinpath(@__DIR__, "data/test_sbarbot_stress_hex8.dat"), ' ', Float64)
        fu = (x) -> sbarbot_stress_hex8(x..., q1, q2, q3, L, T, W, theta, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, G, nu)
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

@testset "SBarbot Tet4" begin
    epsv11 = 11e0
    epsv12 = 5e0
    epsv13 = 6e0
    epsv22 = 7e0
    epsv23 = 8e0
    epsv33 = 13e0
    G = 1.0
    nu = 0.25
    A = [1.0, 2.0, 3.0]
    B = [2.0, 5.0, 4.0]
    C = [11.0, 13.0, 12.0]
    D = [6.0, 5.0, 4.0]

    x1s = range(-100.0, stop=100.0, step=20.0)
    x2s = range(-100.0, stop=100.0, step=20.0)
    x3s = range(20.0, stop=30.0, step=5.0)
    xxs = product(x1s, x2s, x3s)
    qd = gausslegendre(15)

    @testset "Displacement" begin
        u_truth = readdlm(joinpath(@__DIR__, "data/test_sbarbot_disp_tet4.dat"), ' ', Float64)
        fu = (x) -> sbarbot_disp_tet4(qd, x..., A, B, C, D, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, nu)
        u_cal = map(fu, xxs) |> vec
        ftest = (i) -> u_cal[i] ≈ u_truth[i,:]
        @test map(ftest, 1: length(u_cal)) |> all
    end

    @testset "Stress" begin
        u_truth = readdlm(joinpath(@__DIR__, "data/test_sbarbot_stress_tet4.dat"), ' ', Float64)
        fu = (x) -> sbarbot_stress_tet4(qd, x..., A, B, C, D, epsv11, epsv12, epsv13, epsv22, epsv23, epsv33, G, nu)
        u_cal = map(fu, xxs) |> vec
        ftest = (i) -> u_cal[i] ≈ u_truth[i,:]
        @test map(ftest, 1: length(u_cal)) |> all
    end
end
