using Test
using GmshTools
using FastGaussQuadrature
using LinearAlgebra
using JuEQ:
    unit_dislocation, unit_strain,
    shear_traction_sbarbot, shear_traction_dc3d,
    stress_components, coordinate_sbarbot2okada!,
    relative_velocity!, relative_strain!,
    dτ_dt!, dσ_dt!,
    deviatoric_stress!, stress_norm!

@testset "Unit dislocation for plane fault types" begin
    @test unit_dislocation(DIPPING()) == [0.0, 1.0, 0.0]
    @test unit_dislocation(STRIKING()) == [1.0, 0.0, 0.0]
end

@testset "Unit strain for volume" begin
    @test unit_strain(Val(:xx)) == [0., 0., 0., 1., 0., 0.]
    @test unit_strain(Val(:xy)) == [0., 1., 0., 0., 0., 0.]
    @test unit_strain(Val(:xz)) == [0., 0., 0., 0., -1., 0.]
    @test unit_strain(Val(:yy)) == [1., 0., 0., 0., 0., 0.]
    @test unit_strain(Val(:yz)) == [0., 0., -1., 0., 0., 0.]
    @test unit_strain(Val(:zz)) == [0., 0., 0., 0., 0., 1.]
end

@testset "Stress green's function between SBarbotMeshEntity and OkadaMesh" begin
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
        ud = unit_dislocation(ft)
        st = okada_stress_gf_tensor(mf, ma, λ, μ, ft; nrept=0, buffer_ratio=0)
        indexST = Base.OneTo(6)
        for _ in 1: 5 # random check 5 position
            i, j, k = rand(1: mf.nx), rand(1: mf.nξ), rand(1: length(ma.tag))
            u = dc3d(ma.x2[k], ma.x1[k], -ma.x3[k], α, 0.0, 45.0, mf.ax[i], mf.aξ[j], ud)
            σ = stress_components(u, λ, μ)
            @test map(x -> st[x][k, s2i[i,j]], indexST) == σ
        end
    end

    function test_sbarbot2okada_traction(ft)
        st = sbarbot_stress_gf_tensor(ma, mf, λ, μ, ft, allcomp)
        for (ic, _st) in enumerate(st)
            uϵ = unit_strain(Val(allcomp[ic]))
            for _ in 1: 5 # random check 5 position
                i, j, k = rand(1: mf.nx), rand(1: length(ma.tag)), rand(1: mf.nξ)
                σ = sbarbot_stress_hex8(mf.y[k], mf.x[i], -mf.z[k], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
                τ = shear_traction_sbarbot(ft, σ, λ, μ, mf.dip)
                @test _st[s2i[i,k],j] == τ
            end
        end
    end

    function test_sbarbot_self_stress()
        st = sbarbot_stress_gf_tensor(ma, λ, μ, allcomp)
        indexST = Base.OneTo(6)
        for (ic, _st) in enumerate(st)
            uϵ = unit_strain(Val(allcomp[ic]))
            for _ in 1: 5 # random check 5 position
                i, j = rand(1: length(ma.tag), 2)
                σ = sbarbot_stress_hex8(ma.x1[i], ma.x2[i], ma.x3[i], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
                coordinate_sbarbot2okada!(σ)
                @test σ == map(x -> _st[x][i,j], indexST)
            end
        end
    end

    allfaulttype = [STRIKING(), DIPPING()]
    allcomp = (:xx, :xy, :xz, :yz, :yy, :zz)
    foreach(test_okada2sbarbot_stress, allfaulttype)
    foreach(test_sbarbot2okada_traction, allfaulttype)
    test_sbarbot_self_stress()
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

@testset "Consistency of shear traction of SBarbot green's function" begin
    quadrature = gausslegendre(15)
    mf = gen_mesh(Val(:RectOkada), 100.0, 20.0, 20.0, 10.0, 45.0)

    f1 = tempname() * ".msh"
    f2 = tempname() * ".msh"

    rfzn = ones(Int64, 2)
    rfzh = ones(length(rfzn)) |> cumsum
    normalize!(rfzh, Inf)
    gen_gmsh_mesh(Val(:BoxHexExtrudeFromSurface), -20.0, -20.0, -20.0, 40.0, 40.0, 40.0, 2, 2, 1.0, 1.0, rfzn, rfzh; filename=f1)

    @gmsh_do begin
        gmsh.model.occ.addBox(-20.0, -20.0, -60.0, 40.0, 40.0, 40.0, 1)
        @addOption begin
            "Mesh.CharacteristicLengthMax", 40.0
            "Mesh.CharacteristicLengthMin", 40.0
        end
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.write(f2)
    end

    mc1 = read_gmsh_mesh(Val(:SBarbotHex8), f1; phytag=-1, check=true)
    mc2 = read_gmsh_mesh(Val(:SBarbotTet4), f2; phytag=-1)

    comp = rand([:xx, :xy, :xz, :yy, :yz, :zz])
    ft = rand([DIPPING(), STRIKING()])
    λ, μ = 1.0, 1.0

    st_hex = sbarbot_stress_gf_tensor(mc1, mf, λ, μ, ft, comp)
    st_tet = sbarbot_stress_gf_tensor(mc2, mf, λ, μ, ft, comp; quadrature=quadrature)

    ϵ1 = ones(length(mc1.tag))
    σ1 = st_hex * ϵ1
    ϵ2 = ones(length(mc2.tag))
    σ2 = st_tet * ϵ2
    @test norm(σ1 - σ2) < 1e-3

    foreach(rm, [f1, f2])
end

@testset "BLAS update stress/traction rate" begin
    tmpfile = tempname() * ".msh"
    mf = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, 90.0)
    rfzn = ones(4)
    rfzh = accumulate((x, y) -> x * y, fill(1.2, size(rfzn))) |> cumsum |> x -> normalize!(x, Inf)
    mv = gen_gmsh_mesh(mf,
        Val(:BoxHexExtrudeFromSurface), -60.0, -30.0, -60.0, 120.0, 60.0, 60.0,
        4, 4, 1.0, 3.0, rfzn, rfzh; filename=tmpfile)
    me = read_gmsh_mesh(Val(:SBarbotHex8), tmpfile; phytag=1)
    λ, μ = 1.0, 1.0
    ft = rand([DIPPING(), STRIKING()])
    comp = (:xx, :xy, :xz, :yy, :yz, :zz)

    gfoo = okada_stress_gf_tensor(mf, λ, μ, ft)
    gfos = okada_stress_gf_tensor(mf, me, λ, μ, ft)
    gfso = sbarbot_stress_gf_tensor(me, mf, λ, μ, ft, comp)
    gfss = sbarbot_stress_gf_tensor(me, λ, μ, comp)
    alos = gen_alloc(mf, me, length(comp))

    ϵ = rand(length(me.tag), length(comp))
    ϵ₀ = rand(length(comp))
    v = rand(mf.nx, mf.nξ)
    vpl = rand()

    @testset "relative velocity" begin
        relative_velocity!(alos.e, vpl, v)
        @test alos.e.relvnp ≈ v .- vpl
    end

    @testset "relative strain rate" begin
        relative_strain!(alos.v, ϵ₀, ϵ)
        for i = 1: length(comp)
            @test alos.v.relϵ[i] == ϵ[:,i] .- ϵ₀[i]
        end
    end

    @testset "inelastic ⟶ elastic" begin
        fill!(alos.e.dτ_dt, 0.0)
        res = zeros(mf.nx * mf.nξ)
        for i = 1: length(comp)
            res .+= vec(alos.e.dτ_dt) + gfso[i] * alos.v.relϵ[i]
        end
        dτ_dt!(gfso, alos)
        @test res ≈ vec(alos.e.dτ_dt)
    end

    @testset "elastic ⟶ inelastic" begin
        for i = 1: 6
            res = gfos[i] * vec(alos.e.relvnp)
            dσ_dt!(gfos, alos)
            @test alos.v.dσ′_dt[i] ≈ res
        end
    end

    @testset "inelastic ⟷ inelastic" begin
        res = [zeros(length(me.tag)) for _ in 1: 6]
        for i = 1: 6
            res[i] .+= alos.v.dσ′_dt[i]
            for j = 1: length(comp)
                res[i] .+= gfss[j][i] * alos.v.relϵ[j]
            end
        end
        dσ_dt!(gfss, alos.v)
        @test map(i -> res[i] ≈ alos.v.dσ′_dt[i], 1: 6) |> all
    end

    rm(tmpfile)
end

@testset "Deviatoric stress and norm" begin
    alv = gen_alloc(10, 3, 6)
    dσ = rand(10, 6)
    for i = 1: 6
        alv.dσ′_dt[i] .= dσ[:,i]
    end
    σkk = (dσ[:,1] + dσ[:,4] + dσ[:,6]) / 3
    dσ[:,1] -= σkk
    dσ[:,4] -= σkk
    dσ[:,6] -= σkk
    deviatoric_stress!(alv)
    @test map(i -> dσ[:,i] ≈ alv.dσ′_dt[i], 1: 6) |> all
    stress_norm!(alv)
    @test alv.dς′_dt ≈ map(x -> norm(dσ[x,:]), 1: size(dσ, 1))
end
