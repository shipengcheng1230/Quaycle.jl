@testset "rheology property" begin
    p1 = DislocationCreepProperty([rand(5) for _ in 1: 6]...)
    p2 = DiffusionCreepProperty([rand(5) for _ in 1: 7]...)
    pe = ElasticRSFProperty([rand(3, 3) for _ in 1: 4]..., rand(6)...)

    pvm = ViscoelasticMaxwellProperty(pe, p1, p2)
    σ, τ = [rand(5) for _ in 1: 2]

    dϵ1 = @. pvm.pv.disl * τ^(pvm.pv.n - 1) * σ
    dϵ2 = @. pvm.pv.diff * σ

    dϵ1′ = dϵ_dt.(Ref(DislocationCreep()), p1.A, σ, τ, p1.n, p1.fH₂0, p1.r, p1.Q, p1.T)
    dϵ2′ = dϵ_dt.(Ref(DiffusionCreep()), p2.A, σ, p2.d, p2.m, p2.fH₂0, p2.r, p2.Q, p2.T)

    @test dϵ1 ≈ dϵ1′
    @test dϵ2 ≈ dϵ2′
end
