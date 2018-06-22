using Base.Test
using JuEQ

function example_dense_loading()
    tload = collect(0.: 1.: 40.)
    vload = zeros(tload)
    vload[11: end] = 10.
    DenseSDLoading(1e-3, tload, vload)
end

@testset "Dieterich" begin
    @testset "without damping" begin
        ev = DieterichLaw(0.0015, 1e-4)
        ld = ConstantSDLoading(50., 1e-5)
        rsf = UnitRSFModel(0.001, 0.6, 1e-6, 0.0, ev)
        sol = simulate(rsf, ld, (0., 100.), reltol=1e-8, abstol=1e-8, saveat=5.)
        μ = μ_equation(rsf, sol[1, :], sol[2, :])
        μ_truth = [
            0.60000000, 0.60170649, 0.60171045, 0.60096431, 0.60023676,
            0.59967431, 0.59929548, 0.59907074, 0.59895143, 0.59889360,
            0.59886748, 0.59885630, 0.59885171, 0.59884987, 0.59884915,
            0.59884888, 0.59884877, 0.59884873, 0.59884872, 0.59884871,
            0.59884871,
        ]
        @test isapprox(μ, μ_truth, rtol=1e-8)
    end

    @testset "with dampling" begin
        ev = DieterichLaw(0.0015, 3e-5)
        ld = ConstantSDLoading(10., 1e-5)
        rsf = UnitRSFModel(0.001, 0.6, 1e-6, 0.5, ev)
        sol = simulate(rsf, ld, (0., 500.), reltol=1e-8, abstol=1e-8, saveat=25.)
        μ = μ_equation(rsf, sol[1, :], sol[2, :])
        μ_truth = [
            0.60000000, 0.59650552, 0.59894203, 0.60106049, 0.59777486,
            0.60015600, 0.59659840, 0.59903374, 0.60109141, 0.59786773,
            0.60024044, 0.59669104, 0.59912469, 0.60109864, 0.59796052,
            0.60032384, 0.59678383, 0.59921537, 0.60106116, 0.59805324,
            0.60040602,
        ]
        @test isapprox(μ, μ_truth, rtol=1e-8)
    end
end

@testset "Ruina law" begin
    ld = example_dense_loading()
    ev = RuinaLaw(0.005, 10.)
    rsf = UnitRSFModel(0.01, 0.6, 1., 0., ev)
    sol = simulate(rsf, ld, (0., 40.), reltol=1e-8, abstol=1e-8, saveat=1., alg_hints=[:stiff])
    μ = μ_equation(rsf, sol[1, :], sol[2, :])
    μ_truth = [
        0.60000000, 0.59904759, 0.59818148, 0.59738987, 0.59666284,
        0.59599208, 0.59537057, 0.59479242, 0.59425263, 0.59374696,
        0.59820309, 0.60697059, 0.61398982, 0.61719027, 0.61539604,
        0.61277364, 0.61160777, 0.61132454, 0.61135336, 0.61143551,
        0.61149072, 0.61151367, 0.61151845, 0.61151688, 0.61151464,
        0.61151332, 0.61151284, 0.61151277, 0.61151283, 0.61151289,
        0.61151292, 0.61151293, 0.61151293, 0.61151293, 0.61151293,
        0.61151293, 0.61151293, 0.61151293, 0.61151293, 0.61151293,
        0.61151293,
    ]
    @test isapprox(μ, μ_truth, rtol=1e-8)
end

@testset "PRZ law" begin
    ld = example_dense_loading()
    ev = PRZLaw(0.005, 10.)
    rsf = UnitRSFModel(0.01, 0.6, 1., 0., ev)
    sol = simulate(rsf, ld, (0., 40.), reltol=1e-8, abstol=1e-8, saveat=1., alg_hints=[:stiff])
    μ = μ_equation(rsf, sol[1, :], sol[2, :])
    μ_truth = [
        0.60000000, 0.59906355, 0.59823645, 0.59749754, 0.59683103,
        0.59622488, 0.59566972, 0.59515814, 0.59468418, 0.59424300,
        0.59877105, 0.60770816, 0.61517034, 0.61930953, 0.61852021,
        0.61630011, 0.61516675, 0.61484249, 0.61483781,0.61490373,
        0.61495412, 0.61497710, 0.61498298, 0.61498223, 0.61498036,
        0.61497913, 0.61497863, 0.61497853, 0.61497857, 0.61497862,
        0.61497865, 0.61497866, 0.61497866, 0.61497866, 0.61497866,
        0.61497866, 0.61497866, 0.61497866,0.61497866, 0.61497866,
        0.61497866,
    ]
    @test isapprox(μ, μ_truth, rtol=1e-8)
end

@testset "Different ODEs settings" begin
    ld = example_dense_loading()
    ev = DieterichLaw(0.005, 10.)
    rsf = UnitRSFModel(0.01, 0.6, 1., 0., ev)

    sol_vθ = simulate(rsf, ld, (0., 40.), Val{:vθ}, reltol=1e-8, abstol=1e-8, saveat=1., alg_hints=[:stiff])
    μ = μ_equation(rsf, sol_vθ[1, :], sol_vθ[2, :])

    sol_μθ = simulate(rsf, ld, (0., 40.), Val{:μθ}, reltol=1e-8, abstol=1e-8, saveat=1., alg_hints=[:stiff])
    v = v_equation(rsf, sol_μθ[1, :], sol_μθ[2, :])

    sol_μvθ = simulate(rsf, ld, (0., 40.), Val{:μvθ}, reltol=1e-8, abstol=1e-8, saveat=1., alg_hints=[:stiff])

    @test isapprox(sol_μθ[1, :], μ, rtol=1e-6)
    @test isapprox(sol_vθ[1, :], v, rtol=1e-6)
    @test isapprox(sol_μvθ[1, :], μ, rtol=1e-6)
    @test isapprox(sol_μvθ[2, :], v, rtol=1e-6)
end
