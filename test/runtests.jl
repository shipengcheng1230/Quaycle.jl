using Base.Filesystem
push!(LOAD_PATH, joinpath(dirname(@__DIR__), "src"))

using RateState
using Base.Test

println("Test Okada Solution")
@time include("test_okada_solution.jl")
