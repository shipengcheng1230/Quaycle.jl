## single degree of freedom
"Just Hook's law. Notice `vpl` is the leading velocity as opposed those in [`relative_velocity!`](@ref) where `vpl` is passively."
@inline dτ_dt(K::T, v::T, vpl::T) where T = @fastmath K * (vpl - v)

"Derivative of velocity in quai-dynamic rate-and-state governing equation."
@inline dv_dt(dτdt::T, dμdv::T, dμdθ::T, dθdt::T, η::T) where T = @fastmath (dτdt - dμdθ * dθdt) / (dμdv + η)

"""
Out-place version of derivative of *velocity* (or you may call it slip rate)
and *state* (the *state* in rate-and-state friction) for single degree of freedom system.
This for now only support single *state* variable, which is most-widely used.
"""
@inline function dv_dθ_dt(::CForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, ::Vararg{T},
    ) where {T<:Number, SE<:StateEvolutionLaw}
    dμdθ = @fastmath σ * b / θ
    dμdv = @fastmath σ * a / v
    dθdt = dθ_dt(se, v, θ, L)
    dτdt = dτ_dt(k, v, vpl)
    return dv_dt(dτdt, dμdv, dμdθ, dθdt, η), dθdt
end

@inline function dv_dθ_dt(::RForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T,
    ) where {T<:Number, SE<:StateEvolutionLaw}
    ψ1 = @fastmath exp((f0 + b * log(v0 * θ / L)) / a) / 2v0
    ψ2 = @fastmath σ * ψ1 / √(1 + (v * ψ1)^2)
    dμdv = @fastmath a * ψ2
    dμdθ = @fastmath b / θ * v * ψ2
    dθdt = dθ_dt(se, v, θ, L)
    dτdt = dτ_dt(k, v, vpl)
    return dv_dt(dτdt, dμdv, dμdθ, dθdt, η), dθdt
end

## inplace derivatives
@inline function dμ_dvdθ!(::CForm, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where T
    @fastmath @inbounds @threads for i ∈ eachindex(v)
        alloc.dμ_dθ[i] = p.σ[i] * p.b[i] / θ[i]
        alloc.dμ_dv[i] = p.σ[i] * p.a[i] / v[i]
    end
end

@inline function dμ_dvdθ_atomic!(i::I, ::CForm, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where {T, I}
    alloc.dμ_dθ[i] = @fastmath p.σ[i] * p.b[i] / θ[i]
    alloc.dμ_dv[i] = @fastmath p.σ[i] * p.a[i] / v[i]
end

@inline function dμ_dvdθ!(::RForm, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where T
    @fastmath @inbounds @threads for i ∈ eachindex(v)
        ψ1 = exp((p.f0 + p.b[i] * log(p.v0 * θ[i] / p.L[i])) / p.a[i]) / 2p.v0
        ψ2 = p.σ[i] * ψ1 / hypot(1, v[i] * ψ1)
        alloc.dμ_dv[i] = p.a[i] * ψ2
        alloc.dμ_dθ[i] = p.b[i] / θ[i] * v[i] * ψ2
    end
end

@inline function dμ_dvdθ_atomic!(i::I, ::RForm, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where {T, I}
    ψ1 = @fastmath exp((p.f0 + p.b[i] * log(p.v0 * θ[i] / p.L[i])) / p.a[i]) / 2p.v0
    ψ2 = @fastmath p.σ[i] * ψ1 / hypot(1, v[i] * ψ1)
    alloc.dμ_dv[i] = @fastmath p.a[i] * ψ2
    alloc.dμ_dθ[i] = @fastmath p.b[i] / θ[i] * v[i] * ψ2
end

@inline function dvdθ_dt!(se::StateEvolutionLaw, dv::T, dθ::T, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where T
    @inbounds @threads for i ∈ eachindex(v)
        dθ[i] = dθ_dt(se, v[i], θ[i], p.L[i])
        dv[i] = dv_dt(alloc.dτ_dt[i], alloc.dμ_dv[i], alloc.dμ_dθ[i], dθ[i], p.η)
    end
end

@inline function dvdθ_dt_atomic!(i::I, se::StateEvolutionLaw, dv::T, dθ::T, v::T, θ::T, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation) where {T, I}
    dθ[i] = dθ_dt(se, v[i], θ[i], p.L[i])
    dv[i] = dv_dt(alloc.dτ_dt[i], alloc.dμ_dv[i], alloc.dμ_dθ[i], dθ[i], p.η)
end

@inline function deviatoric_stress!(σ::AbstractVecOrMat{T}, alloc::StressRateAllocation) where T
    BLAS.blascopy!(alloc.numσ * alloc.nume, σ, 1, alloc.σ′, 1)
    if alloc.numdiag > 1
        # this is not rigorously correct if only one diagonal are considered
        # in such case, that one diagonal is considered as deviatoric component
        # otherwise, no strain rate will be on that component ever
        # BLAS.gemv!('N', one(T), σ, alloc.isdiag, zero(T), alloc.ς′) # store `σkk` into `ς′`
        mul!(alloc.ς′, σ, alloc.isdiag)
        @strided @fastmath alloc.ς′ ./= alloc.numdiag
        BLAS.ger!(-one(T), alloc.ς′, alloc.isdiag, alloc.σ′)
    end

    # deviatoric stress norm: ς′ = |σ′| = √(σ′:σ′) = √(tr(σ′ᵀσ′))
    fill!(alloc.ς′, zero(T))
    @inbounds @fastmath @threads for i in eachindex(alloc.ς′)
        for j in 1: alloc.numσ
            alloc.ς′[i] += alloc.σ′[i,j] ^ 2 * alloc.diagcoef[j] # for higher precision use `hypot` or `norm`
        end
        alloc.ς′[i] = √(alloc.ς′[i])
    end
end

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

@inline function dϵ_dt!(dϵ::AbstractArray, p::CompositePlasticDeformationProperty, alloc::StressRateAllocMatrix)
    alloc.ς′ .^= p.n₋1
    @inbounds @fastmath for j = 1: alloc.numϵ
        @threads for i = 1: alloc.nume
            dϵ[i,j] = (p.disl[i] * alloc.ς′[i] + p.diff[i]) * alloc.σ′[i,alloc.ϵ2σ[j]]
        end
    end
end

@inline function dδ_dt!(dδ::T, v::T) where T
    @inbounds @fastmath @threads for i ∈ eachindex(dδ)
        dδ[i] = v[i]
    end
end

@inline dδ_dt_atomic!(i::I, dδ::T, v::T) where {T, I} = dδ[i] = v[i]

## concrete combination of derivatives
function ∂u∂t(du::ArrayPartition{T}, u::ArrayPartition{T}, p::RateStateQuasiDynamicProperty, alloc::TractionRateAllocation{N}, gf::AbstractArray, flf::FrictionLawForm, se::StateEvolutionLaw,
    ) where {T, N}
    v, θ, δ = u.x # velocity, state, slip
    dv, dθ, dδ = du.x
    clamp!(θ, zero(T), Inf)
    clamp!(v, zero(T), Inf)
    relative_velocity!(alloc, p.vpl, v)
    dτ_dt!(gf, alloc)
    atomic_update_vθδ!(v, θ, dv, dθ, dδ, p, alloc, flf, se)
end

function ∂u∂t(du::ArrayPartition{T}, u::ArrayPartition{T}, p::ViscoelasticMaxwellProperty, alloc::ViscoelasticCompositeAlloc{N}, gf::ViscoelasticCompositeGreensFunction, flf::FrictionLawForm, se::StateEvolutionLaw,
    ) where {T, N}
    v, θ, ϵ, σ, δ = u.x # velocity, state, strain, stress, slip
    dv, dθ, dϵ, dσ, dδ = du.x
    clamp!(θ, 0.0, Inf)
    clamp!(v, 0.0, Inf)
    relative_velocity!(alloc.e, p.pe.vpl, v)

    # deviatoric_stress!(σ, alloc.v)
    # dϵ_dt!(dϵ, p.pv, alloc.v)
    # relative_strain_rate!(alloc.v, p.pv.dϵref, dϵ)

    atomic_update_ϵσ!(σ, dϵ, p.pv, alloc.v)
    dτ_dt!(gf.ee, alloc.e) # clear `dτ_dt`
    dτ_dt!(gf.ve, alloc) # accumulate `dτ_dt`
    dσ_dt!(dσ, gf.ev, alloc.e) # clear `dσ_dt`
    dσ_dt!(dσ, gf.vv, alloc.v) # accumulate `dσ_dt`
    atomic_update_vθδ!(v, θ, dv, dθ, dδ, p.pe, alloc.e, flf, se)
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
