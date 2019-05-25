using Test

@testset "unit dislocation for plane fault types" begin
    @test JuEQ.unit_dislocation(DIPPING()) == [0.0, 1.0, 0.0]
    @test JuEQ.unit_dislocation(STRIKING()) == [1.0, 0.0, 0.0]
end

@testset "Okada Strain Green's Function towards SBarbot Mesh Entity" begin
    filename = tempname() * ".msh"
    rfzn = ones(Int64, 10)
    rfzh = rand(10) |> cumsum
    normalize!(rfzh, Inf)
    gen_gmsh_mesh(Val(:BoxHexExtrudeFromSurface), -150.0, -75.0, -60.0, 300.0, 150.0, 100.0, 10, 10, 1.5, 2.0, rfzn, rfzh; filename=filename)
    ma = read_gmsh_mesh(Val(:SBarbotHex8), filename)
    comp = collect(1: 6)
    λ, μ = 1.0, 1.0
    α = (λ + μ) / (λ + 2μ)

    @testset "Strking Fault" begin
        ft = STRIKING()
        mf = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 90.0)
        st = okada_strain_gf_tensor(mf, ma, λ, μ, ft, comp; nrept=0, buffer_ratio=0)
        ud = JuEQ.unit_dislocation(ft)
        s2i = LinearIndices((mf.nx, mf.nξ))
        st = okada_strain_gf_tensor(mf, ma, λ, μ, ft, comp; nrept=0, buffer_ratio=0)
        for _ in 1: 10 # random check 10 position
            i, j, k = rand(1: mf.nx), rand(1: mf.nξ), rand(1: length(ma.tag))
            u = dc3d(ma.x2[k], ma.x1[k], -ma.x3[k], α, 0.0, 90.0, mf.ax[i], mf.aξ[j], ud)
            ϵ = JuEQ.strain_components(u)
            @test map(x -> st[x][s2i[i,j], k], comp) == ϵ
        end
        @test_throws AssertionError okada_strain_gf_tensor(mf, ma, λ, μ, STRIKING(), [1, 2, 3, 4, 5, 6, 1]; nrept=0, buffer_ratio=0)
        @test_throws AssertionError okada_strain_gf_tensor(mf, ma, λ, μ, STRIKING(), [-1]; nrept=0, buffer_ratio=0)
        @test_throws AssertionError okada_strain_gf_tensor(mf, ma, λ, μ, STRIKING(), [7]; nrept=0, buffer_ratio=0)
    end

    @testset "Dipping Fault" begin
        ft = DIPPING()
        mf = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 45.0)
        st = okada_strain_gf_tensor(mf, ma, λ, μ, ft, comp; nrept=0, buffer_ratio=0)
        ud = JuEQ.unit_dislocation(ft)
        s2i = LinearIndices((mf.nx, mf.nξ))
        st = okada_strain_gf_tensor(mf, ma, λ, μ, ft, comp; nrept=0, buffer_ratio=0)
        for _ in 1: 10 # random check 10 position
            i, j, k = rand(1: mf.nx), rand(1: mf.nξ), rand(1: length(ma.tag))
            u = dc3d(ma.x2[k], ma.x1[k], -ma.x3[k], α, 0.0, 45.0, mf.ax[i], mf.aξ[j], ud)
            ϵ = JuEQ.strain_components(u)
            @test map(x -> st[x][s2i[i,j], k], comp) == ϵ
        end
    end
    rm(filename)
end
