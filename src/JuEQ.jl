module JuEQ

using Reexport

include("RateStateFriction.jl")
include("BEM.jl")

@reexport using .RateStateFriction
@reexport using .BEM

end
