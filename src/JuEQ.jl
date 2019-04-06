module JuEQ

using Reexport

@reexport using DifferentialEquations
@reexport using Sundials

using SimpleTraits
using Parameters
using FFTW
using FFTW: Plan

using Distributed
using Base.Threads
using LinearAlgebra
using SharedArrays

include("config.jl")
include("geometry.jl")
include("greensfun.jl")

include("rsf.jl")
include("fault.jl")
include("bem.jl")
include("constructor.jl")
# include("utils.jl")
# include("dc3d.jl")


export NormalFault, ThrustFault, StrikeSlipFault
export PlaneFaultDomain
export fault

export DieterichStateLaw, RuinaStateLaw, PrzStateLaw
export CForm, RForm
export friction
export MaterialProperties

export discretize, properties, stiffness_tensor
export EarthquakeCycleProblem

export max_velocity, moment_magnitude, DECallbackSaveToFile

end
