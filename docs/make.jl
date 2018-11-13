# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, JuEQ
using Documenter
using Literate
using Plots # to not capture precompilation output

ENV["GKSwstype"] = "100"

function convert_examples_to_markdown(name::AbstractString; outdir=joinpath(@__DIR__, "src/generated"))
    ex = joinpath(@__DIR__, "..", "examples", name)
    Literate.markdown(ex, outdir)
end

examples = ["bp1.jl", "otfsync.jl"]
map(convert_examples_to_markdown, examples)

makedocs(
    modules = [JuEQ],
    format = :html,
    sitename = "JuEQ",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Quasi-dynamic" => [
                "Introduction" => "quasi_dynamic_intro.md",
                "Example 1: 1D fault" => "generated/bp1.md",
                "Example 2: 2D fault" => "generated/otfsync.md",
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
)
