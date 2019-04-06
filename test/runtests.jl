using JuEQ

const TESTDIR = @__DIR__

for f in filter!(x -> startswith(x, "test_") && endswith(x, ".jl"), readdir(TESTDIR))
    include(abspath(joinpath(TESTDIR, f)))
end
