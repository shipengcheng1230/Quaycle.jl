using Test
using JuEQ: dϵ_dt

@testset "Rheology property" begin
    p1 = DislocationCreepProperty([rand(5) for _ in 1: 10]...)
    p2 = DiffusionCreepProperty([rand(5) for _ in 1: 11]...)
    pe = RateStateQuasiDynamicProperty([rand(3, 3) for _ in 1: 4]..., rand(4)...)

    pvm = compose(pe, rand(3), [:yz, :xy, :zz], p1, p2)
    @test pvm.pv.dϵind == [5, 2, 6]
    σ, τ = [rand(5) for _ in 1: 2]

    dϵ1 = @. pvm.pv.disl * τ^(pvm.pv.n - 1) * σ
    dϵ2 = @. pvm.pv.diff * σ

    dϵ1′ = dϵ_dt.(Ref(DislocationCreep()), p1.A, σ, τ, p1.n, p1.COH, p1.r, p1.α, p1.ϕ, p1.Q, p1.P, p1.Ω, p1.T)
    dϵ2′ = dϵ_dt.(Ref(DiffusionCreep()), p2.A, σ, p2.d, p2.m, p2.COH, p2.r, p2.α, p2.ϕ, p2.Q, p2.P, p2.Ω, p2.T)

    @test dϵ1 ≈ dϵ1′
    @test dϵ2 ≈ dϵ2′
end
