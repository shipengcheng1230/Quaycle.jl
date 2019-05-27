using Test
using GmshTools

@testset "unit dislocation for plane fault types" begin
    @test JuEQ.unit_dislocation(DIPPING()) == [0.0, 1.0, 0.0]
    @test JuEQ.unit_dislocation(STRIKING()) == [1.0, 0.0, 0.0]
end

@testset "Okada Stress Green's Function towards SBarbot Mesh Entity" begin
    filename = tempname() * ".msh"
    rfzn = ones(Int64, 10)
    rfzh = rand(10) |> cumsum
    normalize!(rfzh, Inf)
    gen_gmsh_mesh(Val(:BoxHexExtrudeFromSurface), -150.0, -75.0, -60.0, 300.0, 150.0, 100.0, 10, 10, 1.5, 2.0, rfzn, rfzh; filename=filename)
    ma = read_gmsh_mesh(Val(:SBarbotHex8), filename)
    λ, μ = 1.0, 1.0
    α = (λ + μ) / (λ + 2μ)
    mf = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 45.0)
    s2i = LinearIndices((mf.nx, mf.nξ))

    function __test__(ft)
        ud = JuEQ.unit_dislocation(ft)
        st = okada_stress_gf_tensor(mf, ma, λ, μ, ft; nrept=0, buffer_ratio=0)
        for _ in 1: 10 # random check 10 position
            i, j, k = rand(1: mf.nx), rand(1: mf.nξ), rand(1: length(ma.tag))
            u = dc3d(ma.x2[k], ma.x1[k], -ma.x3[k], α, 0.0, 45.0, mf.ax[i], mf.aξ[j], ud)
            σ = JuEQ.stress_components(u, λ, μ)
            @test map(x -> st[x][k, s2i[i,j]], 1: 6) == σ
        end
    end

    __test__(STRIKING())
    __test__(DIPPING())

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
    τyz1 = shear_traction_dc3d(DIPPING(), u, λ, μ, dip)
    τyz2 = shear_traction_sbarbot(DIPPING(), σ, λ, μ, dip)
    @test τyz1 ≈ τyz2
    τxy1 = shear_traction_dc3d(STRIKING(), u, λ, μ, dip)
    τxy2 = shear_traction_sbarbot(STRIKING(), σ, λ, μ, dip)
    @test τxy1 ≈ τxy2
end
