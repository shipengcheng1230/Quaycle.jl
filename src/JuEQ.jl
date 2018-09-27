module JuEQ

using Reexport

@reexport using DifferentialEquations

using SimpleTraits
using Parameters
using FFTW
using FFTW: Plan

using Distributed
using Base.Threads
using LinearAlgebra
using SharedArrays
using FileIO
using JLD2

include("rsf.jl")
include("fault.jl")
include("bem.jl")
include("constructor.jl")
include("utils.jl")
include(joinpath(@__DIR__, "dc3d.jl"))

export NormalFault, ThrustFault, StrikeSlipFault
export PlaneFaultDomain
export fault

export DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export CForm, RForm
export friction
export MaterialProperties

export discretize, properties
export EarthquakeCycleProblem
export dc3d_okada

export max_velocity

end
