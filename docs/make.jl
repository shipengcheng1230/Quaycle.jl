# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, JuEQ

makedocs(
    modules = [JuEQ],
    format = :html,
    sitename = "JuEQ",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Quasi-dynamic" => [
                "Introduction" => "quasi_dynamic_intro.md"
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
