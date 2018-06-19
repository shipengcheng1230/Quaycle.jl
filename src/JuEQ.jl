module JuEQ

using Reexport

include("RateStateFriction.jl")

@reexport using .RateStateFriction

end
