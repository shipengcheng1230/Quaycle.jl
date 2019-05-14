using JuEQ

const TESTDIR = @__DIR__

if isempty(ARGS)
    TESTFILES = filter(x -> startswith(x, "test_") && endswith(x, ".jl"), readdir(TESTDIR))
else
    TESTFILES = filter(x -> startswith(x, "_") && endswith(x, ".jl"), readdir(TESTDIR))
end

foreach(x -> include(joinpath(TESTDIR, x)), TESTFILES)
