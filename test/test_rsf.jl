using Base.Test
using JuEQ

@testset "RateStateFriction" begin
    @testset "Dieterich Law" begin

    model = UniformRSFModel(0.005, 0.010, 10.0, 0.6, 1e-3, 1.0, Dieterich)
    tload = collect(0.: 1: 40.)
    vload = ones(tload)
    vload[11: end] = 10.
    loading = DenseTimeLoading(tload, vload)
    u0 = [0.6, 1., model.L / model.vref]
    tspan = (0., 40.)
    sol = solve_model(model, loading, u0, tspan; saveat=1., reltol=1e-8, abstol=1e-8)

    μ_truth = [
        0.60000000, 0.60000000, 0.60000000, 0.60000000, 0.60000000,
        0.60000000, 0.60000000, 0.60000000, 0.60000000, 0.60000000,
        0.60418935, 0.60744265, 0.59014796, 0.58400035, 0.58748100,
        0.58986633, 0.58949737, 0.58825787, 0.58807707, 0.58845076,
        0.58862857, 0.58854378, 0.58845178, 0.58845594, 0.58849077,
        0.58849965, 0.58848946, 0.58848316, 0.58848503, 0.58848791,
        0.58848807, 0.58848707, 0.58848670, 0.58848697, 0.58848718,
        0.58848714, 0.58848706, 0.58848704, 0.58848707, 0.58848709,
        0.58848708,
        ]

    @test isapprox(sol[1, :], μ_truth, rtol=1e-8)
    end
end
