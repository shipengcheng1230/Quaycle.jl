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
    comp = collect(1: 6)
    λ, μ = 1.0, 1.0
    α = (λ + μ) / (λ + 2μ)
    mf = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 45.0)
    s2i = LinearIndices((mf.nx, mf.nξ))

    function __test__(ft)
        ud = JuEQ.unit_dislocation(ft)
        st = okada_stress_gf_tensor(mf, ma, λ, μ, ft, comp; nrept=0, buffer_ratio=0)
        for _ in 1: 10 # random check 10 position
            i, j, k = rand(1: mf.nx), rand(1: mf.nξ), rand(1: length(ma.tag))
            u = dc3d(ma.x2[k], ma.x1[k], -ma.x3[k], α, 0.0, 45.0, mf.ax[i], mf.aξ[j], ud)
            σ = JuEQ.stress_components(u, λ, μ)
            @test map(x -> st[x][k, s2i[i,j]], comp) == σ
        end
    end

    __test__(STRIKING())
    __test__(DIPPING())
    @test_throws AssertionError okada_stress_gf_tensor(mf, ma, λ, μ, STRIKING(), [1, 2, 3, 4, 5, 6, 1]; nrept=0, buffer_ratio=0)
    @test_throws AssertionError okada_stress_gf_tensor(mf, ma, λ, μ, STRIKING(), [-1]; nrept=0, buffer_ratio=0)
    @test_throws AssertionError okada_stress_gf_tensor(mf, ma, λ, μ, STRIKING(), [7]; nrept=0, buffer_ratio=0)

    rm(filename)
end
