## calculate stiffness matrix
dc3d_file = joinpath(dirname(@__DIR__), "src/dc3d.jl")
# need string insertion here to work around
@everywhere include($dc3d_file)
