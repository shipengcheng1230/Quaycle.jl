module Quaycle

using Reexport

@reexport using OrdinaryDiffEq
@reexport using DiffEqCallbacks
@reexport using RecursiveArrayTools

using Parameters
using Requires
using FFTW
using PhysicalConstants
using BlockArrays
using Strided
using WriteVTK
using Formatting
using ProgressMeter

using TensorOperations

using Distributed
using Base.Threads
using LinearAlgebra
using SharedArrays

include("index.jl")
include("friction.jl")
include("rheology.jl")
include("property.jl")
include("mesh.jl")
include("space.jl")
include("gf.jl")
include("derivative.jl")
include("assemble.jl")
include("visualize.jl")

include("tools/maxvelocity.jl")
include("tools/mmapsave.jl")

function __init__()
    @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" begin
        using .HDF5
        include("tools/h5solution.jl")
        include("tools/h5getstore.jl")
    end

    @require GmshTools="82e2f556-b1bd-5f1a-9576-f93c0da5f0ee" begin
        using .GmshTools
        include("tools/gmshtools.jl")
    end

    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
        using .PyPlot
        include("tools/pyplottools.jl")
    end

    include(joinpath(@__DIR__, "config.jl"))
end

end # module
