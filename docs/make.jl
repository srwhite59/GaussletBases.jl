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
        "Manual" => [
            "Overview" => "manual/index.md",
            "First radial workflow" => "tutorials/first_radial_workflow.md",
            "Recommended atomic setup" => "howto/recommended_atomic_setup.md",
            "Current atomic branch" => "explanations/current_atomic_branch.md",
            "Current ordinary branch" => "explanations/current_ordinary_branch.md",
            "Example guide" => "howto/example_guide.md",
        ],
        "Reference" => [
            "Reference overview" => "reference/index.md",
            "Bases and mappings" => "reference/bases_and_mappings.md",
            "Operators and diagnostics" => "reference/operators_and_diagnostics.md",
            "Atomic and ordinary workflows" => "reference/atomic_and_ordinary.md",
            "Export layer" => "reference/export.md",
        ],
        "Developer Notes" => [
            "Overview" => "developer/index.md",
            "Architecture and current direction" => "developer/architecture.md",
            "Supporting note map" => "developer/supporting_notes.md",
        ],
    ],
)
