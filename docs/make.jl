using Documenter

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using GaussletBases

makedocs(
    sitename = "GaussletBases.jl",
    modules = [GaussletBases],
    doctest = false,
    checkdocs = :none,
    format = Documenter.HTML(prettyurls = false, edit_link = nothing),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "First radial workflow" => "tutorials/first_radial_workflow.md",
        ],
        "How-to" => [
            "Recommended atomic setup" => "howto/recommended_atomic_setup.md",
            "Example guide" => "howto/example_guide.md",
        ],
        "Explanations" => [
            "Current atomic branch" => "explanations/current_atomic_branch.md",
            "Current ordinary branch" => "explanations/current_ordinary_branch.md",
            "Architecture" => "explanations/architecture.md",
        ],
    ],
)
