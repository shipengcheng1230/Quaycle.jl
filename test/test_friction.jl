using Test
using Sundials

@testset "SDOF without damping" begin
    mp = SingleDofRSFProperty(a=0.001, b=0.0015, L=1e-4, k=50.0, vpl=1e-5, f0=0.6, v0=1e-6, η=0.0, σ=1.0)
    prob = assemble(mp, [1e-6, mp.L/1e-6], (0., 100.,); se=DieterichStateLaw(), flf=CForm())
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=5.0)
    μ = friction(CForm(), sol.u, mp)
    μ_truth = [
        0.60000000, 0.60170649, 0.60171045, 0.60096431, 0.60023676,
        0.59967431, 0.59929548, 0.59907074, 0.59895143, 0.59889360,
        0.59886748, 0.59885630, 0.59885171, 0.59884987, 0.59884915,
        0.59884888, 0.59884877, 0.59884873, 0.59884872, 0.59884871,
        0.59884871,
        ]
    @test isapprox(μ, μ_truth, rtol=1e-8)
end

@testset "SDOF with damping" begin
    mp = SingleDofRSFProperty(a=0.001, b=0.0015, L=3e-5, k=10.0, vpl=1e-5, f0=0.6, v0=1e-6, η=0.5, σ=1.0)
    prob = assemble(mp, [1e-6, mp.L/1e-6], (0., 500.,); se=DieterichStateLaw(), flf=CForm())
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=25.)
    μ = friction(CForm(), sol.u, mp)
    μ_truth = [
        0.60000000, 0.59650552, 0.59894203, 0.60106049, 0.59777486,
        0.60015600, 0.59659840, 0.59903374, 0.60109141, 0.59786773,
        0.60024044, 0.59669104, 0.59912469, 0.60109864, 0.59796052,
        0.60032384, 0.59678383, 0.59921537, 0.60106116, 0.59805324,
        0.60040602,
    ]
    @test isapprox(μ, μ_truth, rtol=1e-8)
end

@testset "SDOF with stepping loading" begin
    condition = (u, t , integrator) -> t - 10.0
    affect! = (integrator) -> integrator.p.vpl = 10.
    cb = ContinuousCallback(condition, affect!, nothing, save_positions=(false, false))

    @testset "Ruina" begin
        mp = SingleDofRSFProperty(a=0.01, b=0.005, L=10., k=1e-3, vpl=1.0, f0=0.6, v0=1., η=0., σ=1.0)
        prob = assemble(mp, [1.0, mp.L/1.0], (0., 40.,); se=RuinaStateLaw(), flf=CForm())
        # sol = solve(prob, CVODE_BDF(), reltol=1e-8, abstol=1e-8, callback=cb, saveat=1.0)
        # μ = friction(CForm(), sol.u, mp)
        # μ_truth = [
        #     0.60000000, 0.60000000, 0.60000000, 0.60000000, 0.60000000,
        #     0.60000000, 0.60000000, 0.60000000, 0.60000000, 0.60000000,
        #     0.60000000, 0.60000000, 0.60000000, 0.60000000, 0.61451424,
        #     0.61239234, 0.61152559, 0.61133979, 0.61138186, 0.61145348,
        #     0.61149781, 0.61151494, 0.61151781, 0.61151612, 0.61151421,
        #     0.61151317, 0.61151282, 0.61151279, 0.61151285, 0.61151290,
        #     0.61151292, 0.61151293, 0.61151293, 0.61151293, 0.61151293,
        #     0.61151293, 0.61151293, 0.61151293, 0.61151293, 0.61151293,
        #     0.61151293,
        # ]
        # Since Julia 1.4.0, this fails on Linux OS and Windows OS
        # It also breaks for Mac OS on Travis, but not on my personal computer.
        # @test_skip isapprox(μ, μ_truth, rtol=1e-8)
    end

    @testset "PRZ" begin
        mp = SingleDofRSFProperty(a=0.01, b=0.005, L=10., k=1e-3, vpl=1.0, f0=0.6, v0=1., η=0., σ=1.0)
        prob = assemble(mp, [1.0, mp.L/1.0], (0., 40.,); se=PrzStateLaw(), flf=CForm())
        # sol = solve(prob, CVODE_BDF(), reltol=1e-8, abstol=1e-8, callback=cb, saveat=1.0)
        # μ = friction(CForm(), sol.u, mp)
        # μ_truth = [
        #     0.60000000, 0.60001736, 0.60006450, 0.60013523, 0.60022460,
        #     0.60032861, 0.60044398, 0.60056796, 0.60069827, 0.60083298,
        #     0.60097045, 0.60958894, 0.61636558, 0.61928354, 0.61785150,
        #     0.61592716, 0.61505648, 0.61483759, 0.61485737, 0.61491944,
        #     0.61496151, 0.61497911, 0.61498283, 0.61498167, 0.61497997,
        #     0.61497897, 0.61497860, 0.61497855, 0.61497859, 0.61497863,
        #     0.61497866, 0.61497866, 0.61497866, 0.61497866, 0.61497866,
        #     0.61497866, 0.61497866, 0.61497866, 0.61497866, 0.61497866,
        #     0.61497866,
        # ]
        # @test isapprox(μ, μ_truth, rtol=1e-8)
    end
end

@testset "SDOF of variational form" begin
    mp = SingleDofRSFProperty(a=0.001, b=0.0015, L=3e-5, k=10.0, vpl=1e-5, f0=0.6, v0=1e-6, η=0.5, σ=1.0)
    prob = assemble(mp, [1e-6, mp.L/1e-6], (0., 500.,); se=DieterichStateLaw(), flf=RForm())
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=25.)
    μ = friction(RForm(), sol.u, mp)
    μ_truth = [
        0.60000000, 0.59885200, 0.59884871, 0.59884871, 0.59884871,
        0.59884871, 0.59884871, 0.59884871, 0.59884871, 0.59884871,
        0.59884871, 0.59884871, 0.59884871, 0.59884871, 0.59884871,
        0.59884871, 0.59884871, 0.59884871, 0.59884871, 0.59884871,
        0.59884871
    ]
    @test isapprox(μ, μ_truth, rtol=1e-8)
end

@testset "SDOF of info display" begin
    mp = SingleDofRSFProperty(a=0.001, b=0.0015, L=3e-5, k=10.0, vpl=1e-5, f0=0.6, v0=1e-6, η=0.0, σ=1.0)
    @test_logs (:warn, "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt.") assemble(mp, [1e-6, mp.L/1e-6], (0., 500.,); se=DieterichStateLaw(), flf=RForm())
end
