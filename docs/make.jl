using Documenter

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using GaussletBases

makedocs(
    sitename = "GaussletBases.jl",
    modules = [GaussletBases],
    doctest = true,
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
        "Reference" => [
            "Reference overview" => "reference/index.md",
            "Bases and mappings" => "reference/bases_and_mappings.md",
            "Operators and diagnostics" => "reference/operators_and_diagnostics.md",
            "Atomic and ordinary workflows" => "reference/atomic_and_ordinary.md",
            "Export layer" => "reference/export.md",
        ],
    ],
)
