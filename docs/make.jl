# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, DocumenterMarkdown, JuEQ
using Plots # to not capture precompilation output

ENV["GKSwstype"] = "100"

include("generate.jl")

cp(joinpath(@__DIR__, "..", "LICENSE"), joinpath(@__DIR__, "src", "LICENSE.md"); force = true)

makedocs(
    doctest=false,
    modules = [JuEQ],
    sitename = "JuEQ",
    format = Markdown(),
)

deploydocs(
  repo  = "github.com/shipengcheng1230/JuEQ.jl.git",
  deps = Deps.pip("pygments", "mkdocs", "python-markdown-math", "mkdocs-windmill"),
  make = () -> run(`mkdocs build`),
  julia = "1.0",
)
