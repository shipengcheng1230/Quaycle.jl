## static green's function

const KERNELDIR = joinpath(@__DIR__, "gfkernel")

for f in filter!(x -> endswith(x, ".jl"), readdir(KERNELDIR))
    include(abspath(joinpath(KERNELDIR, f)))
end

## okada
