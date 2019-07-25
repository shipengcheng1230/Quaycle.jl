using JuEQ

# const TESTDIR = @__DIR__
# const TESTFILES = filter(x -> startswith(x, "test_") && endswith(x, ".jl"), readdir(TESTDIR))
# foreach(x -> include(joinpath(TESTDIR, x)), TESTFILES)

using GmshTools
try
    gmsh.initialize()
    @test true
catch
    @test false
finally
    gmsh.finalize()
end
