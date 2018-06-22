module RateStateFriction

using DifferentialEquations

export # rate-state model settings
    RSFModel,
    UnitRSFModel

export # state evolution laws
    EvolutionLaw,
    DieterichLaw, RuinaLaw, PRZLaw

export # external loading system
    LoadingSystem,
    SingleDegreeLoading,
    ConstantSDLoading,
    DenseSDLoading

export # system derivations
    dϴdt, dμdt, dμdv, dμdθ, dvdt_dϴdt

export # simulation tools
    simulate, solve_model, build_model

export # system equations
    μ_equation, v_equation

export # odes settings
	odes_v_θ!, odes_μ_θ!, odes_μ_v_θ!

## State evolution law
abstract type StateEvolutionLaw end
abstract type OneStateEvolutionLaw <: StateEvolutionLaw end

struct DieterichLaw{T} <: OneStateEvolutionLaw where {T <: Number}
    b::T
    L::T
end

struct RuinaLaw{T} <: OneStateEvolutionLaw where {T <: Number}
    b::T
    L::T
end

struct PRZLaw{T} <: OneStateEvolutionLaw where {T <: Number}
    b::T
    L::T
end

function dθdt(ev::DieterichLaw, v, θ)
    1.0 - v * θ / ev.L
end

function dθdt(ev::RuinaLaw, v, θ)
    x = v * θ / ev.L
    -x * log(x)
end

function dθdt(ev::PRZLaw, v, θ)
    1.0 - (v * θ / 2ev.L) ^ 2
end

## RSF model parameters
abstract type RSFModel end

# unit system
struct UnitRSFModel{T, L} <: RSFModel where {T <: Number, L <: StateEvolutionLaw}
    a::T
    μref::T
    vref::T
    η::T
    θlaw::L
end

## Loading types
abstract type LoadingSystem end
abstract type SingleDegreeLoading <: LoadingSystem end

struct ConstantSDLoading{T} <: SingleDegreeLoading where {T <: Number}
    k::T
    vload::T
end

struct DenseSDLoading{T, A} <: SingleDegreeLoading where {T <: Number, A <: AbstractVector}
    k::T
    tload::A
    vload::A
end

## Derivation of system variables: μ, v, θ
function dμdt(ld::ConstantSDLoading, v::T, t::T) where {T <: Number}
    ld.k * (ld.vload - v)
end

function dμdt(ld::DenseSDLoading, v, t)
    nearest = indmin(abs(t - ld.tload))
    ld.k * (ld.vload[nearest] - v)
end

function dμdθ(ev::OneStateEvolutionLaw, θ)
    ev.b / θ
end

function dμdv(rsf::UnitRSFModel, v)
    rsf.a / v
end

function dvdt_dϴdt(rsf::UnitRSFModel, ld::SingleDegreeLoading, v, θ, t)
    dμ_dt = dμdt(ld, v, t)
    dμ_dθ = dμdθ(rsf.θlaw, θ)
    dμ_dv = dμdv(rsf, v)
    dθ_dt = dθdt(rsf.θlaw, v, θ)
    dv_dt = (dμ_dt - dμ_dθ * dθ_dt) / (dμ_dv + rsf.η)
    return dv_dt, dθ_dt
end

## Simulation procedure
function simulate(rsf::RSFModel, ld::LoadingSystem, tspan, depvar=Val{:vθ}, u0=nothing; kwargs...)
    prob = build_model(rsf, ld, tspan, depvar, u0)
    sol = solve_model(prob; kwargs...)
end

function solve_model(prob::ODEProblem; kwargs...)
    sol = solve(prob; kwargs...)
end

function build_model(rsf::UnitRSFModel, ld::SingleDegreeLoading, tspan, depvar::Type{Val{:vθ}}, u0)
    u0 == nothing && (u0 = [rsf.vref, rsf.θlaw.L / rsf.vref])
	prob = ODEProblem(odes_v_θ!, u0, tspan, (rsf, ld))
end

function build_model(rsf::UnitRSFModel, ld::SingleDegreeLoading, tspan, depvar::Type{Val{:μθ}}, u0)
	u0 == nothing && (u0 = [rsf.μref, rsf.θlaw.L / rsf.vref])
	prob = ODEProblem(odes_μ_θ!, u0, tspan, (rsf, ld))
end

function build_model(rsf::UnitRSFModel, ld::SingleDegreeLoading, tspan, depvar::Type{Val{:μvθ}}, u0)
    u0 == nothing && (u0 = [rsf.μref, rsf.vref, rsf.θlaw.L / rsf.vref])
    prob = ODEProblem(odes_μ_v_θ!, u0, tspan, (rsf, ld))
end

"""
*u* represents in order: v, θ
"""
function odes_v_θ!(du, u, p, t)
    rsf, ld = p
    du[1], du[2] = dvdt_dϴdt(rsf, ld, u[1], u[2], t)
end

"""
*u* represents in order: μ, θ
"""
function odes_μ_θ!(du, u, p, t)
	rsf, ld = p
	v = v_equation(rsf, u[1], u[2])
	du[1] = dμdt(ld, v, t)
	du[2] = dθdt(rsf.θlaw, v, u[2])
end

"""
*u* represents in order: μ, v, θ
"""
function odes_μ_v_θ!(du, u, p, t)
    rsf, ld = p
    du[1] = dμdt(ld, u[2], t)
    du[2], du[3] = dvdt_dϴdt(rsf, ld, u[2], u[3], t)
end

## Governing equations
function μ_equation(rsf::UnitRSFModel, v, θ)
    a, b, μref, vref, L = rsf.a, rsf.θlaw.b, rsf.μref, rsf.vref, rsf.θlaw.L
    μref + a * log(v / vref) + b * log(vref * θ / L)
end

function v_equation(rsf::UnitRSFModel, μ, θ)
	μ_contrib_by_θ = rsf.θlaw.b * log(rsf.vref * θ / rsf.θlaw.L)
	v = rsf.vref * exp((μ - rsf.μref - μ_contrib_by_θ) / rsf.a)
end

end
