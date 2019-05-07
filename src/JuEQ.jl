module JuEQ

using Reexport

@reexport using DifferentialEquations
@reexport using Sundials
@reexport using HDF5
@reexport using FastGaussQuadrature

using Parameters
using FFTW
using FFTW: Plan

using Distributed
using Base.Threads
using LinearAlgebra
using SharedArrays

include("config.jl")
include("properties.jl")
include("rsf.jl")
include("mesh.jl")
include("fault.jl")
include("greensfunction.jl")
include("assemble.jl")

const TOOLS = abspath(joinpath(@__DIR__, "tools"))
map(x -> include(joinpath(TOOLS, x)), filter!(x -> endswith(x, ".jl"), readdir(TOOLS)))

end # module
