using Test
using DelimitedFiles
using Base.Iterators
using FastGaussQuadrature
using GmshTools
using Gmsh_SDK_jll

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
    x3s = range(-30.0, stop=-20.0, step=5.0) # in real case should not be negative
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

@testset "SBarbot Hex8 & Tet4 Convergence" begin
    f1 = tempname() * ".msh"
    rfzn = ones(Int64, 3)
    rfzh = ones(length(rfzn)) |> cumsum
    normalize!(rfzh, Inf)
    gen_gmsh_mesh(Val(:BoxHexExtrudeFromSurface), -50.0, -30.0, -10.0, 100.0, 60.0, 50.0, 5, 3, 1.0, 1.0, rfzn, rfzh; filename=f1)

    f2 = tempname() * ".msh"
    @gmsh_do begin
        gmsh.model.occ.addBox(-50.0, -30.0, -60.0, 100.0, 60.0, 50.0, 1)
        @addOption begin
            "Mesh.CharacteristicLengthMax", 20.0
            "Mesh.CharacteristicLengthMin", 20.0
        end
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.write(f2)
    end

    # Both `f1` and `f2` are two identical box shape with characteristic length of 20.
    # And we test if both displacement and stress are consistent for these two different meshes.

    mc1 = read_gmsh_mesh(Val(:SBarbotHex8), f1; phytag=-1, check=true)
    mc2 = read_gmsh_mesh(Val(:SBarbotTet4), f2; phytag=-1)

    obs = rand(3)
    u1, u2 = [zeros(Float64, 3) for _ in 1: 2]
    σ1, σ2 = [zeros(Float64, 6) for _ in 1: 2]
    uvec = Vector{Float64}(undef, 3)
    σvec = Vector{Float64}(undef, 6)
    ϵ = rand(6)
    G = 1.0
    ν = 0.25
    quadrature = gausslegendre(11)

    @inbounds @fastmath @simd for i in 1: length(mc1.tag)
        sbarbot_disp_hex8!(uvec, obs..., mc1.q1[i], mc1.q2[i], mc1.q3[i], mc1.L[i], mc1.T[i], mc1.W[i], 0.0, ϵ..., G, ν)
        sbarbot_stress_hex8!(σvec, obs..., mc1.q1[i], mc1.q2[i], mc1.q3[i], mc1.L[i], mc1.T[i], mc1.W[i], 0.0, ϵ..., G, ν)
        u1 .+= uvec
        σ1 .+= σvec
    end

    @inbounds @fastmath @simd for i in 1: length(mc2.tag)
        sbarbot_disp_tet4!(uvec, quadrature, obs..., mc2.A[i], mc2.B[i], mc2.C[i], mc2.D[i], ϵ..., ν)
        sbarbot_stress_tet4!(σvec, quadrature, obs..., mc2.A[i], mc2.B[i], mc2.C[i], mc2.D[i], ϵ..., G, ν)
        u2 .+= uvec
        σ2 .+= σvec
    end

    @test norm(u1 - u2) < 1e-3
    @test norm(σ1 - σ2) < 1e-3
    foreach(rm, [f1, f2])
end

@testset "SBarbot Quad4" begin
    epsv22 = 7e-6
    epsv23 = 8e-6
    epsv33 = 9e-6
    T = 10e3
    W = 20e3
    q2 = 15e3
    q3 = 1e-3
    phi = 0.0
    G = 1.0
    nu = 0.25

    x3 = range(-10000.0, stop=-1000.0, step=1000)
    x2 = range(-10000.0, stop=10000.0, step=2000)
    xxs = product(x2, x3)

    @testset "Stress" begin
        u_truth = readdlm(joinpath(@__DIR__, "data", "test_sbarbot_stress_quad4.dat"), ' ', Float64)
        fu = (x) -> sbarbot_stress_quad4(x..., q2, q3, T, W, phi, epsv22, epsv23, epsv33, G, nu)
        u_cal = map(fu, xxs) |> vec
        ftest = (i) -> u_cal[i] ≈ u_truth[i,:]
        @test map(ftest, 1: length(u_cal)) |> all
    end
end
