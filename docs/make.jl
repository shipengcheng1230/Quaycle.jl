# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, Quaycle
using Plots # to not capture precompilation output

using HDF5
using GmshTools

# ENV["GKSwstype"] = "100"

include("generate.jl")

makedocs(
    doctest=false,
    modules = [Quaycle],
    sitename = "Quaycle",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "Quasi Dynamic Simulation using BEM" => [
            "Example 1D" => "examples/generated/bp1.md",
            "Example 2D" => "examples/generated/otfsync.md",
        ],
        "Interface" => [
            "Assemble" => "interface_assemble.md",
            "Friction" => "interface_friction.md",
            "Greens Function" => "interface_greensfun.md",
            "HDF5 Utility" => "interface_HDF5.md",
            "Mesh" => "interface_mesh.md",
            "Simulation Property" => "interface_property.md",
            "Rheology" => "interface_rheology.md",
            "Visualize" => "interface_visualize.md",
        ],
        "Acknowledge" => "acknowledge.md",
    ],
)

deploydocs(
  repo = "github.com/shipengcheng1230/Quaycle.jl.git",
  # deps = Deps.pip("pymdown-extensions", "pygments", "mkdocs", "python-markdown-math", "mkdocs-material"),
  target = "build",
  # make = () -> run(`mkdocs build`),
)
