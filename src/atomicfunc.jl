@inline function dμ_dvdθ_atomic!(i::I, ::CForm, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where {T, I}
    alloc.dμ_dθ[i] = @fastmath p.σ[i] * p.b[i] / θ[i]
    alloc.dμ_dv[i] = @fastmath p.σ[i] * p.a[i] / v[i]
end

@inline function dμ_dvdθ_atomic!(i::I, ::RForm, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where {T, I}
    ψ1 = @fastmath exp((p.f0 + p.b[i] * log(p.v0 * θ[i] / p.L[i])) / p.a[i]) / 2p.v0
    ψ2 = @fastmath p.σ[i] * ψ1 / hypot(1, v[i] * ψ1)
    alloc.dμ_dv[i] = @fastmath p.a[i] * ψ2
    alloc.dμ_dθ[i] = @fastmath p.b[i] / θ[i] * v[i] * ψ2
end

@inline function dvdθ_dt_atomic!(i::I, se::StateEvolutionLaw, dv::T, dθ::T, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where {T, I}
    dθ[i] = dθ_dt(se, v[i], θ[i], p.L[i])
    dv[i] = dv_dt(alloc.dτ_dt[i], alloc.dμ_dv[i], alloc.dμ_dθ[i], dθ[i], p.η)
end

@inline dδ_dt_atomic!(i::I, dδ::T, v::T) where {T, I} = dδ[i] = v[i]

@inline function dϵdt_reldϵ_atomic!(i::I, σ::T, dϵ::T, p::CompositePlasticDeformationProperty, alloc::StressRateAllocMatrix) where {T, I}
    σkk = zero(eltype(T))
    @fastmath @simd for j ∈ 1: alloc.numσ
        σkk += σ[i, j] * alloc.isdiag[j] # σ[i, :] ⋅ alloc.isdiag
    end
    # this is not rigorously correct if only one diagonal are considered
    # in such case, that one diagonal is considered as deviatoric component
    # otherwise, no strain rate will be on that component ever
    @fastmath σkk /= alloc.numdiag
    @fastmath @simd for j ∈ 1: alloc.numσ
        alloc.σ′[i, j] = σ[i, j] - σkk * alloc.isdiag[j] # σ′ᵢⱼ = σᵢⱼ - δᵢⱼσkk
    end
    alloc.ς′[i] = zero(eltype(T))
    @fastmath @simd for j ∈ 1: alloc.numσ
        # deviatoric stress norm: ς′ = |σ′| = √(σ′:σ′) = √(tr(σ′ᵀσ′))
        alloc.ς′[i] += alloc.σ′[i, j] ^ 2 * alloc.diagcoef[j]
    end
    @fastmath alloc.ς′[i] = √(alloc.ς′[i]) ^ p.n₋1[i] * p.disl[i] + p.diff[i]
    @fastmath @simd for j ∈ 1: alloc.numϵ
        dϵ[i, j] = alloc.ς′[i] * alloc.σ′[i, alloc.ϵ2σ[j]]
        alloc.reldϵ[i, j] = dϵ[i, j] - p.dϵref[j]
    end
end

@inline function atomic_update_vθδ!(v::T, θ::T, dv::T, dθ::T, dδ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation, flf::FrictionLawForm, se::StateEvolutionLaw) where T
    @inbounds @threads for i ∈ eachindex(v)
        dμ_dvdθ_atomic!(i, flf, v, θ, p, alloc)
        dvdθ_dt_atomic!(i, se, dv, dθ, v, θ, p, alloc)
        dδ_dt_atomic!(i, dδ, v)
    end
end

@inline function atomic_update_ϵσ!(σ, dϵ, p::CompositePlasticDeformationProperty, alloc::StressRateAllocation)
    @inbounds @threads for i ∈ size(σ, 1)
        dϵdt_reldϵ_atomic!(i, σ, dϵ, p, alloc)
    end
end

function ∂u∂t_atomic(du::ArrayPartition{T}, u::ArrayPartition{T}, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation{N}, gf::AbstractArray, flf::FrictionLawForm, se::StateEvolutionLaw,
    ) where {T, N}
    v, θ, δ = u.x # velocity, state, slip
    dv, dθ, dδ = du.x
    clamp!(θ, zero(T), Inf)
    clamp!(v, zero(T), Inf)
    relative_velocity!(alloc, p.vpl, v)
    dτ_dt!(gf, alloc)
    atomic_update_vθδ!(v, θ, dv, dθ, dδ, p, alloc, flf, se)
end

function ∂u∂t_atomic(du::ArrayPartition{T}, u::ArrayPartition{T}, p::ViscoelasticMaxwellProperty, alloc::ViscoelasticCompositeAlloc{N}, gf::ViscoelasticCompositeGreensFunction, flf::FrictionLawForm, se::StateEvolutionLaw,
    ) where {T, N}
    v, θ, ϵ, σ, δ = u.x # velocity, state, strain, stress, slip
    dv, dθ, dϵ, dσ, dδ = du.x
    clamp!(θ, zero(T), Inf)
    clamp!(v, zero(T), Inf)
    relative_velocity!(alloc.e, p.pe.vpl, v)
    atomic_update_ϵσ!(σ, dϵ, p.pv, alloc.v)
    dτ_dt!(gf.ee, alloc.e) # clear `dτ_dt`
    dτ_dt!(gf.ve, alloc) # accumulate `dτ_dt`
    dσ_dt!(dσ, gf.ev, alloc.e) # clear `dσ_dt`
    dσ_dt!(dσ, gf.vv, alloc.v) # accumulate `dσ_dt`
    atomic_update_vθδ!(v, θ, dv, dθ, dδ, p.pe, alloc.e, flf, se)
end
