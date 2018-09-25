# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, JuEQ
using Documenter
using Literate
using Plots # to not capture precompilation output

bem_ex1 = joinpath(@__DIR__, "..", "examples", "bp1.jl")
bem_ex1_output = joinpath(@__DIR__, "src/generated")

Literate.markdown(bem_ex1, bem_ex1_output)

makedocs(
    modules = [JuEQ],
    format = :html,
    sitename = "JuEQ",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Quasi-dynamic" => [
                "Introduction" => "quasi_dynamic_intro.md"
                "Example 1: 1D fault" => "generated/bp1.md"
            ]
        ],
        "Libray" => [
            "Public" => "public_interface.md",
            "Private" => "private_interface.md",
        ],
    ],
)

deploydocs(
  repo  = "github.com/shipengcheng1230/JuEQ.jl.git",
  target = "build",
  deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
  make = nothing,
  julia = "1.0",
)
