using Documenter
using Documenter.Remotes
using ContACT

makedocs(;
    modules = [ContACT],
    authors = "Simon Frost",
    sitename = "ContACT.jl",
    repo = Remotes.GitHub("epirecipes", "ContACT.jl"),
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
        "Vignettes" => [
            "Introduction" => "vignettes/01-introduction.md",
            "Composition" => "vignettes/02-composition.md",
            "Categorical Framework" => "vignettes/03-categorical-framework.md",
            "Generalized Contact Matrices" => "vignettes/04-generalized-contact-matrices.md",
            "CoMix/SEP Reconstruction" => "vignettes/05-comix-sep-reconstruction.md",
            "POLYMOD Real-Data Reconstruction" => "vignettes/06-polymod-reconstruction.md",
            "Epidemic Bounds" => "vignettes/07-epidemic-bounds.md",
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
    repo = "github.com/epirecipes/ContACT.jl",
    devbranch = "main",
    push_preview = true,
)
