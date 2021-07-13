using Documenter
using PYMethod

# Setup for doctests in docstrings
DocMeta.setdocmeta!(PYMethod, :DocTestSetup, recursive = true,
    quote
        using PYMethod
    end
)

makedocs(;
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [PYMethod],
    sitename = "PYMethod.jl",
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "manual/FEM.md",
        ],
        "API" => [
            "api/FEM.md",
            "api/ChangEquation.md",
        ]
    ],
    doctest = true, # :fix
)

deploydocs(
    repo = "github.com/KeitaNakamura/PYMethod.jl.git",
    devbranch = "main",
)
