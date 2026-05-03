using Documenter
using Documenter.Remotes
using ContACT

# Copy rendered vignettes into docs/src/vignettes/ so Documenter includes them
vignette_src = joinpath(@__DIR__, "..", "vignettes")
vignette_dst = joinpath(@__DIR__, "src", "vignettes")
mkpath(vignette_dst)

for dir in readdir(vignette_src; join=true)
    isdir(dir) || continue
    name = basename(dir)
    startswith(name, "0") || continue
    dst_dir = joinpath(vignette_dst, name)
    mkpath(dst_dir)
    md_src = joinpath(dir, "index.md")
    if isfile(md_src)
        cp(md_src, joinpath(dst_dir, "index.md"); force=true)
    end
end

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
            "Introduction" => "vignettes/01-introduction/index.md",
            "Composition" => "vignettes/02-composition/index.md",
            "Categorical Framework" => "vignettes/03-categorical-framework/index.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Survey Operations" => "api/survey.md",
            "Matrix Operations" => "api/operations.md",
            "ACSet Schemas" => "api/schemas.md",
        ],
    ],
)

# Copy standalone HTML vignettes into build/vignettes/ for direct access
build_vignettes = joinpath(@__DIR__, "build", "vignettes")
mkpath(build_vignettes)
for dir in readdir(vignette_src; join=true)
    isdir(dir) || continue
    name = basename(dir)
    startswith(name, "0") || continue
    dst_dir = joinpath(build_vignettes, name)
    cp(dir, dst_dir; force=true)
end

deploydocs(;
    repo = "github.com/epirecipes/ContACT.jl",
    devbranch = "main",
    push_preview = true,
)
