using Test
using GmshTools
using FastGaussQuadrature
using LinearAlgebra
using JuEQ:
    unit_dislocation, unit_strain, shear_traction_td,
    shear_traction_sbarbot_on_okada, shear_traction_dc3d,
    stress_components, coordinate_sbarbot2okada!,
    ViscoelasticCompositeGreensFunction,
    relative_velocity!, relative_strain_rate!,
    dτ_dt!, dσ_dt!, deviatoric_stress!, gen_alloc

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
        st = stress_greens_func(mf, ma, λ, μ, ft; nrept=0, buffer_ratio=0)
        indexST = Base.OneTo(6)
        for _ in 1: 5 # random check 5 position
            i, j, k = rand(1: mf.nx), rand(1: mf.nξ), rand(1: length(ma.tag))
            u = dc3d(ma.x2[k], ma.x1[k], -ma.x3[k], α, 0.0, 45.0, mf.ax[i], mf.aξ[j], ud)
            σ = stress_components(u, λ, μ)
            @test map(x -> st[x][k, s2i[i,j]], indexST) == σ
        end
    end

    function test_sbarbot2okada_traction(ft)
        st = stress_greens_func(ma, mf, λ, μ, ft, allcomp)
        for (ic, _st) in enumerate(st)
            uϵ = unit_strain(Val(allcomp[ic]))
            for _ in 1: 5 # random check 5 position
                i, j, k = rand(1: mf.nx), rand(1: length(ma.tag)), rand(1: mf.nξ)
                σ = sbarbot_stress_hex8(mf.y[k], mf.x[i], -mf.z[k], ma.q1[j], ma.q2[j], ma.q3[j], ma.L[j], ma.T[j], ma.W[j], ma.θ, uϵ..., μ, ν)
                τ = shear_traction_sbarbot_on_okada(ft, σ, mf.dip)
                @test _st[s2i[i,k],j] == τ
            end
        end
    end

    function test_sbarbot_self_stress()
        st = stress_greens_func(ma, λ, μ, allcomp)
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
    τyz2 = shear_traction_sbarbot_on_okada(DIPPING(), σ, dip)
    @test τyz1 ≈ τyz2
    τxy1 = shear_traction_dc3d(STRIKING(), u, λ, μ, dip)
    τxy2 = shear_traction_sbarbot_on_okada(STRIKING(), σ, dip)
    @test τxy1 ≈ τxy2
end

@testset "Shear traction consistency between Quad4 and Tri3 dislocation" begin
    θ = rand() * 90
    mo = gen_mesh(Val(:RectOkada), 100.0, 50.0, 10.0, 10.0, θ)
    i, j = rand(1: mo.nx), rand(1: mo.nξ)
    m, n = rand(1: mo.nx), rand(1: mo.nξ)
    λ, μ = 2.0, 1.0
    α = (λ + μ) / (λ + 2μ)
    ν = λ / 2 / (λ + μ)
    disl = [1.0, 1.0, 2.0]
    ss = [1.0, 0.0, 0.0]
    ds = [0.0, cosd(θ), sind(θ)]
    ts = ss × ds
    E = Matrix{Float64}(undef, 3, 3)

    uo = dc3d(mo.x[i], mo.y[j], mo.z[j], α, mo.dep, mo.dip, mo.ax[m], mo.aξ[n], disl)
    E[1,1] = uo[4]
    E[2,2] = uo[8]
    E[3,3] = uo[12]
    E[1,2] = (uo[5] + uo[7]) / 2
    E[1,3] = (uo[6] + uo[10]) / 2
    E[2,3] = (uo[9] + uo[11]) / 2
    E[2,1] = E[1,2]
    E[3,1] = E[1,3]
    E[3,2] = E[2,3]

    cx, cy, cz = mo.x[m], mo.y[n], mo.z[n]
    p1 = [cx - mo.Δx/2, cy + mo.Δξ/2 * cosd(mo.dip), cz + mo.Δξ/2 * sind(mo.dip)]
    p2 = [cx + mo.Δx/2, cy + mo.Δξ/2 * cosd(mo.dip), cz + mo.Δξ/2 * sind(mo.dip)]
    p3 = [cx + mo.Δx/2, cy - mo.Δξ/2 * cosd(mo.dip), cz - mo.Δξ/2 * sind(mo.dip)]
    p4 = [cx - mo.Δx/2, cy - mo.Δξ/2 * cosd(mo.dip), cz - mo.Δξ/2 * sind(mo.dip)]

    σtd1 = td_stress_hs(mo.x[i], mo.y[j], mo.z[j], p1, p4, p3, disl..., λ, μ)
    σtd2 = td_stress_hs(mo.x[i], mo.y[j], mo.z[j], p3, p2, p1, disl..., λ, μ)
    σ = map(+, σtd1, σtd2) |> collect

    fdo = shear_traction_dc3d(DIPPING(), uo, λ, μ, mo.dip)
    fso = shear_traction_dc3d(STRIKING(), uo, λ, μ, mo.dip)

    fdt = shear_traction_td(DIPPING(), σ, ss, ds, ts)
    fst = shear_traction_td(STRIKING(), σ, ss, ds, ts)

    @test fdo ≈ fdt
    @test fso ≈ fst
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

    st_hex = stress_greens_func(mc1, mf, λ, μ, ft, comp)
    st_tet = stress_greens_func(mc2, mf, λ, μ, ft, comp; quadrature=quadrature)

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

    gfoo = stress_greens_func(mf, λ, μ, ft)
    gfos = stress_greens_func(mf, me, λ, μ, ft)
    gfso = stress_greens_func(me, mf, λ, μ, ft, comp)
    gfss = stress_greens_func(me, λ, μ, comp)
    gg = cat_greensfunc(gfoo, gfos, gfso, gfss)

    @testset "concatenate components of green's function" begin
        @test gg.ee == gfoo
        @test gg.ev == vcat(gfos...)
        @test gg.ve == hcat(gfso...)
        @test gg.vv == hcat([vcat(x...) for x in gfss]...)
    end

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
        relative_strain_rate!(alos.v, ϵ₀, ϵ)
        @test alos.v.reldϵ == ϵ .- ϵ₀'
    end

    @testset "inelastic ⟶ elastic" begin
        fill!(alos.e.dτ_dt, 0.0)
        res = zeros(mf.nx * mf.nξ)
        for i = 1: length(comp)
            res .+= vec(alos.e.dτ_dt) + gfso[i] * alos.v.reldϵ[:,i]
        end
        dτ_dt!(gg.ve, alos)
        @test res ≈ vec(alos.e.dτ_dt)
    end

    @testset "elastic ⟶ inelastic" begin
        dσ = Matrix{Float64}(undef, length(me.tag), 6)
        dσ_dt!(dσ, gg.ev, alos.e)
        for i = 1: 6
            res = gfos[i] * vec(alos.e.relvnp)
            @test dσ[:,i] ≈ res
        end
    end

    @testset "inelastic ⟷ inelastic" begin
        dσ = rand(length(me.tag), 6)
        res = zeros(length(me.tag), 6) + dσ
        for i = 1: 6
            for j = 1: length(comp)
                res[:,i] .+= gfss[j][i] * alos.v.reldϵ[:,j]
            end
        end
        dσ_dt!(dσ, gg.vv, alos.v)
        @test res ≈ dσ
    end

    rm(tmpfile)
end

@testset "Deviatoric stress and norm" begin
    alv = gen_alloc(10, 3, 6)
    σ = rand(10, 6)
    σkk = (σ[:,1] + σ[:,4] + σ[:,6]) / 3
    σ′ = copy(σ)
    σ′[:,1] -= σkk
    σ′[:,4] -= σkk
    σ′[:,6] -= σkk
    deviatoric_stress!(σ, alv)
    @test alv.σ′ ≈ σ′
    ς′ = map(x -> norm(σ′[x,:]), 1: alv.nume)
    @test alv.ς′ ≈ ς′
end
