export assemble

## single degree of freedom
"Just Hook's law. Notice `vpl` is the leading velocity as opposed those in [`relative_velocity!`](@ref) where `vpl` is passively."
@inline dτ_dt(K::T, v::T, vpl::T) where {T<:Number} = K * (vpl - v)

"Derivative of velocity in quai-dynamic rate-and-state governing equation."
@inline dv_dt(dτdt::T, dμdv::T, dμdθ::T, dθdt::T, η::T) where {T<:Number} = (dτdt - dμdθ * dθdt) / (dμdv + η)

"""
Out-place version of derivative of *velocity* (or you may call it slip rate)
and *state* (the *state* in rate-and-state friction) for single degree of freedom system.
This for now only support single *state* variable, which is most-widely used.
"""
@inline function dv_dθ_dt(::CForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, ::Vararg{T},
    ) where {T<:Number, SE<:StateEvolutionLaw}
    dμdθ = σ * b / θ
    dμdv = σ * a / v
    dθdt = dθ_dt(se, v, θ, L)
    dτdt = dτ_dt(k, v, vpl)
    return dv_dt(dτdt, dμdv, dμdθ, dθdt, η), dθdt
end

@inline function dv_dθ_dt(::RForm, se::SE, v::T, θ::T, a::T, b::T, L::T, k::T, σ::T, η::T, vpl::T, f0::T, v0::T,
    ) where {T<:Number, SE<:StateEvolutionLaw}
    ψ1 = exp((f0 + b * log(v0 * clamp(θ, zero(T), Inf) / L)) / a) / 2v0
    ψ2 = σ * ψ1 / sqrt(1 + (v * ψ1)^2)
    dμdv = a * ψ2
    dμdθ = b / θ * v * ψ2
    dθdt = dθ_dt(se, v, θ, L)
    dτdt = dτ_dt(k, v, vpl)
    return dv_dt(dτdt, dμdv, dμdθ, dθdt, η), dθdt
end

"""
    assemble(p::SingleDofRSFProperty, u0::AbstractArray, tspan::NTuple;
        flf::FrictionLawForm=CForm(), se::StateEvolutionLaw=DieterichStateLaw())

Assemble the `ODEProblem` for single degree of freedom system.

## Arguments
- `p::SingleDofRSFProperty`: all system properties
- `u0::AbstractArray`: initial condition
- `tspan::NTuple`: time span for simulation
- `flf::FrictionLawForm`: form of friction law, either [`CForm`](@ref) or [`RForm`](@ref)
- `se::StateEvolutionLaw`: state evolutional law, see [`StateEvolutionLaw`](@ref)
"""
function assemble(p::SingleDofRSFProperty, u0::AbstractArray, tspan::NTuple;
    flf::FrictionLawForm=CForm(), se::StateEvolutionLaw=DieterichStateLaw()) where T<:Real
    (typeof(flf) == RForm && p.η ≈ 0) && @warn "Regularized form requires nonzero `η` to avoid `Inf` in dv/dt."
    op! = (du, u, p, t) -> du .= dv_dθ_dt(flf, se, u[1], u[2], p.a, p.b, p.L, p.k, p.σ, p.η, p.vpl, p.f0, p.v0)
    return ODEProblem(op!, u0, tspan, p)
end

## inplace derivatives
@inline function dμ_dvdθ!(::CForm, v::T, θ::T, p::ElasticRSFProperty, alloc::TractionRateAllocation) where T
    @fastmath @inbounds @threads for i = 1: prod(alloc.dims)
        alloc.dμ_dθ[i] = p.σ[i] * p.b[i] / θ[i]
        alloc.dμ_dv[i] = p.σ[i] * p.a[i] / v[i]
    end
end

@inline function dμ_dvdθ!(::RForm, v::T, θ::T, p::ElasticRSFProperty, alloc::TractionRateAllocation) where T
    @fastmath @inbounds @threads for i = 1: prod(alloc.dims)
        ψ1 = exp((p.f0 + p.b[i] * log(p.v0 * θ[i] / p.L[i])) / p.a[i]) / 2p.v0
        ψ2 = p.σ[i] * ψ1 / hypot(1, v[i] * ψ1)
        alloc.dμ_dv[i] = p.a[i] * ψ2
        alloc.dμ_dθ[i] = p.b[i] / θ[i] * v[i] * ψ2
    end
end

@inline function dvdθ_dt!(se::StateEvolutionLaw, dv::T, dθ::T, v::T, θ::T, p::ElasticRSFProperty, alloc::TractionRateAllocation) where T
    @fastmath @inbounds @threads for i = 1: prod(alloc.dims)
        dθ[i] = dθ_dt(se, v[i], θ[i], p.L[i])
        dv[i] = dv_dt(alloc.dτ_dt[i], alloc.dμ_dv[i], alloc.dμ_dθ[i], dθ[i], p.η)
    end
end

function ∂u∂t(du::ArrayPartition{T}, u::ArrayPartition{T}, p::ElasticRSFProperty, alloc::TractionRateAllocation{N}, gf::AbstractArray, flf::FrictionLawForm, se::StateEvolutionLaw,
    ) where {T, N}
    v, θ = u.x # velocity, state
    dv, dθ = du.x
    clamp!(θ, zero(T), Inf)
    clamp!(v, zero(T), Inf)
    relative_velocity!(alloc, p.vpl, v)
    dτ_dt!(gf, alloc)
    dμ_dvdθ!(flf, v, θ, p, alloc)
    dvdθ_dt!(se, dv, dθ, v, θ, p, alloc)
end

@inline function deviatoric_stress!(σ::AbstractVecOrMat, alloc::StressRateAllocation{3})
    BLAS.blascopy!(6alloc.nume, σ, 1, alloc.σ′, 1)
    @inbounds @fastmath @threads for i = 1: alloc.nume
        σkk = (σ[i,1] + σ[i,4] + σ[i,6]) / 3
        alloc.σ′[i,1] -= σkk
        alloc.σ′[i,4] -= σkk
        alloc.σ′[i,6] -= σkk
    end
    @strided alloc.ς′ .= sqrt.(vec(sum(abs2, alloc.σ′; dims=2))) # for higher precision use `hypot` or `norm`
end

@inline function dϵ_dt!(dϵ::AbstractArray, p::CompositePlasticDeformationProperty, alloc::StressRateAllocMatrix)
    @inbounds @fastmath for j = 1: alloc.numϵ
        @threads for i = 1: alloc.nume
            ςⁿ⁻¹ = alloc.ς′[i]^(p.n[i]-1)
            dϵ[i,j] = (p.disl[i] * ςⁿ⁻¹ + p.diff[i]) * alloc.σ′[i,j]
        end
    end
end

function ∂u∂t(du::ArrayPartition{T}, u::ArrayPartition{T}, p::ViscoelasticMaxwellProperty, alloc::ViscoelasticCompositeAlloc{N}, gf::ViscoelasticCompositeGreensFunction, flf::FrictionLawForm, se::StateEvolutionLaw,
    ) where {T, N}
    v, θ, ϵ, σ = u.x # velocity, state, strain, stress
    dv, dθ, dϵ, dσ = du.x
    clamp!(θ, 0.0, Inf)
    clamp!(v, 0.0, Inf)
    relative_velocity!(alloc.e, p.pe.vpl, v)
    deviatoric_stress!(σ, alloc.v)
    dϵ_dt!(dϵ, p.pv, alloc.v)
    relative_strain_rate!(alloc.v, p.pv.dϵref, dϵ)
    dτ_dt!(gf.ee, alloc.e) # clear `dτ_dt`
    dτ_dt!(gf.ve, alloc) # accumulate `dτ_dt`
    dσ_dt!(dσ, gf.ev, alloc.e) # clear `dσ_dt`
    dσ_dt!(dσ, gf.vv, alloc.v) # accumulate `dσ_dt`
    dμ_dvdθ!(flf, v, θ, p.pe, alloc.e)
    dvdθ_dt!(se, dv, dθ, v, θ, p.pe, alloc.e)
end

"""
    assemble(mf::AbstractMesh, gf::AbstractArray,
        p::ElasticRSFProperty, u0::AbstractArray, tspan::NTuple{2};
        flf::FrictionLawForm=RForm(), se::StateEvolutionLaw=DieterichStateLaw(), kwargs...)

Assemble the `ODEProblem` for elastic fault.

##  Extra Arguments
- `mf::AbstractMesh`: fault mesh, currently only support [`OkadaMesh`](@ref)
- `gf::AbstractArray`: green's function associated with `fs.mesh` and `p.λ` & `p.μ`
- `p::ElasticRSFProperty`: all system properties
- `u0::ArrayPartition`: initial condition. By rule of order in this package:
    1. *velocity*
    2. *state*


- `tspan::NTuple`: time span for simulation
- `flf::FrictionLawForm`: form of friction law, either [`CForm`](@ref) or [`RForm`](@ref)
- `se::StateEvolutionLaw`: state evolutional law, see [`StateEvolutionLaw`](@ref)
"""
function assemble(
    mf::AbstractMesh, gf::AbstractArray, p::ElasticRSFProperty, u0::ArrayPartition, tspan::NTuple{2};
    flf::FrictionLawForm=RForm(), se::StateEvolutionLaw=DieterichStateLaw(),
    )
    alloc = gen_alloc(mf)
    _assemble(p, alloc, gf, u0, tspan, flf, se)
end

"""
    assemble(fas::LithAsthSpace, gf::ViscoelasticCompositeGreensFunction,
        p::ViscoelasticMaxwellProperty, u0::ArrayPartition, tspan::NTuple{2};
        flf::FrictionLawForm=RForm(), se::StateEvolutionLaw=DieterichStateLaw(),)

Assemble the `ODEProblem` for elastic fault plus viscoelastic relaxation.

## Arguments
- `fas::LithAsthSpace`: containing fault mesh, fault type and asthenosphere mesh
- `gf::ViscoelasticCompositeGreensFunction`: green's function associated with the composite space `fas`
- `p::ViscoelasticCompositeGreensFunction`: all system properties
- `u0::ArrayPartition`: initial condition. By rule of order in this package:
    1. *velocity*
    2. *state*
    3. *strain*
    4. *stress*


- `tspan::NTuple`: time span for simulation
- `flf::FrictionLawForm`: form of friction law, either [`CForm`](@ref) or [`RForm`](@ref)
- `se::StateEvolutionLaw`: state evolutional law, see [`StateEvolutionLaw`](@ref)
"""
function assemble(
    fas::LithAsthSpace, gf::ViscoelasticCompositeGreensFunction, p::ViscoelasticMaxwellProperty, u0::ArrayPartition, tspan::NTuple{2};
    flf::FrictionLawForm=RForm(), se::StateEvolutionLaw=DieterichStateLaw(),
    )
    alloc = gen_alloc(fas.me, fas.mv, length(p.pv.dϵref))
    _assemble(p, alloc, gf, u0, tspan, flf, se)
end

function _assemble(p::AbstractProperty, alloc::AbstractAllocation, gf, u0, tspan, flf, se)
    f! = (du, u, p, t) -> ∂u∂t(du, u, p, alloc, gf, flf, se)
    return ODEProblem(f!, u0, tspan, p)
end
