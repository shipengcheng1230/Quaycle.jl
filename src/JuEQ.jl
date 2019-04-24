module JuEQ

using Reexport

@reexport using DifferentialEquations
@reexport using Sundials
@reexport using HDF5

using SimpleTraits
using Parameters
using FFTW
using FFTW: Plan

using Distributed
using Base.Threads
using LinearAlgebra
using SharedArrays

include("config.jl")
include("rsf.jl")
include("faulttype.jl")
include("geometry.jl")

include("faultspace.jl")
include("greensfun.jl")
include("gfoperator.jl")
include("properties.jl")
include("assemble.jl")

const UTILSDIR = abspath(joinpath(@__DIR__, "utils"))
map(x -> include(joinpath(UTILSDIR, x)), filter!(x -> endswith(x, ".jl"), readdir(UTILSDIR)))

end # module
