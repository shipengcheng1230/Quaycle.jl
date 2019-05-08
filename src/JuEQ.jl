module JuEQ

using Reexport

@reexport using OrdinaryDiffEq
@reexport using Sundials
@reexport using DiffEqCallbacks

using Parameters
using Requires
using FFTW
using FFTW: Plan

using Distributed
using Base.Threads
using LinearAlgebra
using SharedArrays
using ProgressMeter

include("config.jl")
include("properties.jl")
include("rsf.jl")
include("mesh.jl")
include("fault.jl")
include("greensfunction.jl")
include("assemble.jl")

include("tools/maxvelocity.jl")
include("tools/mmapsave.jl")

function __init__()
    @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" begin
        using .HDF5
        include("tools/h5savecallback.jl")
        include("tools/h5saveprop.jl")
    end
end

end # module
