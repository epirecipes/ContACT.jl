using Documenter
using Documenter.Remotes
using ContACT

makedocs(;
    modules = [ContACT],
    warnonly = true,
    authors = "Simon Frost",
    sitename = "ContACT.jl",
    repo = Remotes.GitHub("epiforecasts", "ContACT.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/getting_started.md",
            "Composition & Stratification" => "tutorials/composition.md",
            "Categorical Framework" => "tutorials/categorical.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Survey Operations" => "api/survey.md",
            "Matrix Operations" => "api/operations.md",
            "ACSet Schemas" => "api/schemas.md",
        ],
    ],
)

deploydocs(;
    repo = "github.com/epiforecasts/ContACT.jl",
    devbranch = "main",
    push_preview = true,
)
