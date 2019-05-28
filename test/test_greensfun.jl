using Test
using GmshTools

@testset "unit dislocation for plane fault types" begin
    @test JuEQ.unit_dislocation(DIPPING()) == [0.0, 1.0, 0.0]
    @test JuEQ.unit_dislocation(STRIKING()) == [1.0, 0.0, 0.0]
end

@testset "Stress Green's Function between SBarbot Mesh Entity and Okada Mesh" begin
    filename = tempname() * ".msh"
    rfzn = ones(Int64, 5)
    rfzh = rand(5) |> cumsum
    normalize!(rfzh, Inf)
    gen_gmsh_mesh(Val(:BoxHexExtrudeFromSurface), -150.0, -75.0, -60.0, 300.0, 150.0, 100.0, 5, 5, 1.5, 2.0, rfzn, rfzh; filename=filename)
    ma = read_gmsh_mesh(Val(:SBarbotHex8), filename)
    λ, μ = 1.0, 1.0
    α = (λ + μ) / (λ + 2μ)
    ν = λ / 2 / (λ + μ)
    mf = gen_mesh(Val(:RectOkada), 100.0, 50.0, 20.0, 20.0, 45.0)
    s2i = LinearIndices((mf.nx, mf.nξ))

    function test_okada2sbarbot_stress(ft)
        ud = JuEQ.unit_dislocation(ft)
        st = okada_stress_gf_tensor(mf, ma, λ, μ, ft; nrept=0, buffer_ratio=0)
        for _ in 1: 5 # random check 5 position
            i, j, k = rand(1: mf.nx), rand(1: mf.nξ), rand(1: length(ma.tag))
            u = dc3d(ma.x2[k], ma.x1[k], -ma.x3[k], α, 0.0, 45.0, mf.ax[i], mf.aξ[j], ud)
            σ = JuEQ.stress_components(u, λ, μ)
            @test map(x -> st[x][k, s2i[i,j]], 1: 6) == σ
        end
    end

    function test_sbarbot2okada_traction(ft)
        st = sbarbot_stress_gf_tensor(ma, mf, λ, μ, ft, allcomp)
        for (ic, _st) in enumerate(st)
            uϵ = JuEQ.unit_strain(Val(allcomp[ic]))
            for _ in 1: 5 # random check 5 position
                i, j, k = rand(1: mf.nx), rand(1: length(ma.tag)), rand(1: mf.nξ)
                σ = sbarbot_stress_hex8(mf.y[k], mf.x[i], -mf.z[k], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
                τ = JuEQ.shear_traction_sbarbot(ft, σ, λ, μ, mf.dip)
                @test _st[s2i[i,k],j] == τ
            end
        end
    end

    function test_sbarbot_self_stress(comp)
        st = sbarbot_stress_gf_tensor(ma, λ, μ, comp)
        uϵ = JuEQ.unit_strain(Val(comp))
        indexST = Base.OneTo(6)
        for _ in 1: 5 # random check 5 position
            i, j = rand(1: length(ma.tag), 2)
            σ = sbarbot_stress_hex8(ma.x1[i], ma.x2[i], ma.x3[i], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
            JuEQ.coordinate_sbarbot2okada!(σ)
            @test σ == map(x -> st[x][i,j], indexST)
        end
    end

    allfaulttype = [STRIKING(), DIPPING()]
    allcomp = (:xx, :xy, :xz, :yz, :yy, :zz)
    foreach(test_okada2sbarbot_stress, allfaulttype)
    foreach(test_sbarbot2okada_traction, allfaulttype)
    foreach(test_sbarbot_self_stress, allcomp)
    rm(filename)
end

@testset "Shear traction consistency between 2 coordinate system" begin
    λ, μ = rand(2)
    dip = rand() * 90
    ν = λ / 2 / (λ + μ)
    rϵ = rand(6)
    # observation point is out of strain volume
    args = (10.0, 10.0, 10.0, -10.0, -10.0, 20.0, 2.0, 2.0, 2.0, 0.0, rϵ..., μ, ν)
    ϵ = sbarbot_strain_hex8(args...)
    σ = sbarbot_stress_hex8(args...)
    # construct an equivalent output from `dc3d`, ignoring the first 3 displacement
    u = [0.0, 0.0, 0.0, ϵ[4], ϵ[2], -ϵ[5], ϵ[2], ϵ[1], -ϵ[3], -ϵ[5], -ϵ[3], ϵ[6]]
    τyz1 = JuEQ.shear_traction_dc3d(DIPPING(), u, λ, μ, dip)
    τyz2 = JuEQ.shear_traction_sbarbot(DIPPING(), σ, λ, μ, dip)
    @test τyz1 ≈ τyz2
    τxy1 = JuEQ.shear_traction_dc3d(STRIKING(), u, λ, μ, dip)
    τxy2 = JuEQ.shear_traction_sbarbot(STRIKING(), σ, λ, μ, dip)
    @test τxy1 ≈ τxy2
end
