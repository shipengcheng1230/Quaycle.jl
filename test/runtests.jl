using Quaycle
using Pkg

# temporary solution to bypass General registry
# required by GmshTools v0.4.0
Pkg.add(PackageSpec(url="https://github.com/shipengcheng1230/Gmsh_SDK_jll.jl"))
using Gmsh_SDK_jll

const TESTDIR = @__DIR__
const TESTFILES = filter(x -> startswith(x, "test_") && endswith(x, ".jl"), readdir(TESTDIR))
foreach(x -> include(joinpath(TESTDIR, x)), TESTFILES)
