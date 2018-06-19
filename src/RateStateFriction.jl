module RateStateFriction

using DifferentialEquations

export # rate-state model settings
    RSFModel,
    UniformRSFModel

export # state evolution laws
    EvolutionLaw,
    DieterichLaw, RuinaLaw, PRZLaw,
    Dieterich, Ruina, PRZ,
    dθdt

export # external loading system
    RSFTimeLoading,
    ConstantTimeLoading,
    DenseTimeLoading,
    FunctionalTimeLoading

export # ODEs derivations
    dμ_dt, dθ_dt, dv_dt

export # simulators
    solve_model, rsf_odes!

## State evolution law
abstract type EvolutionLaw end

struct DieterichLaw <: EvolutionLaw end
struct RuinaLaw <: EvolutionLaw end
struct PRZLaw <: EvolutionLaw end

const Dieterich = DieterichLaw()
const Ruina = RuinaLaw()
const PRZ = PRZLaw()

function dθdt(::DieterichLaw, v, θ, L)
    1.0 - v * θ / L
end

function dθdt(::RuinaLaw, v, θ, L)
    x = v * θ / L
    -x * log(x)
end

function dθdt(::PRZLaw, v, θ, L)
    1.0 - (v * θ / 2L) ^ 2
end

## RSF model parameters
abstract type RSFModel end

# unit system
struct UniformRSFModel{T, L} <: RSFModel where {T <: Number, L <: EvolutionLaw}
    a::T
    b::T
    L::T
    μ0::T
    k::T
    vref::T
    θlaw::L
end

## Loading types
abstract type RSFTimeLoading end

struct ConstantTimeLoading{T} <: RSFTimeLoading where {T <: Number}
    vload::T
end

struct DenseTimeLoading{A} <: RSFTimeLoading where {A <: AbstractVector}
    tload::A
    vload::A
end

struct FunctionalTimeLoading{F} <: RSFTimeLoading where {F <: Function}
    vload::F
end

## Governing equations
function μ_contrib_by_θ(rsf::RSFModel, θ)
    rsf.b * log(rsf.vref * θ / rsf.L)
end

function μ_contrib_by_v(rsf::RSFModel, v)
    rsf.a * log(v / rsf.vref)
end

function μ_equation(rsf::RSFModel, v, θ)
    μ_contrib_by_θ(rsf, state) + μ_contrib_by_v(rsf, state)
end

## Derivation of system variable
function dμdt(rsf::UniformRSFModel, ld::ConstantTimeLoading, μ, v, θ, t)
    rsf.k * (ld.vload - v)
end

function dμdt(rsf::UniformRSFModel, ld::DenseTimeLoading, μ, v, θ, t)
    nearest = indmin(abs(t - ld.tload))
    rsf.k * (ld.vload[nearest] - v)
end

function dvdt(rsf::UniformRSFModel, dμ, dθ, μ, v, θ)
    dv = (dμ - rsf.b / θ * dθ) / (rsf.a / v)
end

## Model solvor
function solve_model(
    rsf::UniformRSFModel, ld::RSFTimeLoading, u0, tspan; kwargs...
    )
    prob = ODEProblem(rsf_odes!, u0, tspan, (rsf, ld))
    sol = solve(prob; kwargs...)
end

# ODEs settings
function rsf_odes!(du, u, p, t)
    rsf, ld = p
    dμ = dμdt(rsf, ld, u[1], u[2], u[3], t)
    dθ = dθdt(rsf.θlaw, u[2], u[3], rsf.L)
    dv = dvdt(rsf, dμ, dθ, u[1], u[2], u[3])
    du[1], du[2], du[3] = dμ, dv, dθ
end

end
