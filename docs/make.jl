# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, DocumenterMarkdown, Quaycle
using Plots # to not capture precompilation output

using HDF5
using GmshTools

ENV["GKSwstype"] = "100"

include("generate.jl")

makedocs(
    doctest=false,
    modules = [Quaycle],
    format = Markdown(),
)

deploydocs(
  repo = "github.com/shipengcheng1230/Quaycle.jl.git",
  deps = Deps.pip("pymdown-extensions", "pygments", "mkdocs", "python-markdown-math", "mkdocs-material"),
  target = "site",
  make = () -> run(`mkdocs build`),
)
