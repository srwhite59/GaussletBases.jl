using Documenter

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using GaussletBases

const DOCS_CI = get(ENV, "CI", "false") == "true"

makedocs(
    sitename = "GaussletBases.jl",
    modules = [GaussletBases],
    doctest = true,
    checkdocs = :none,
    format = Documenter.HTML(
        prettyurls = DOCS_CI,
        edit_link = "main",
        canonical = DOCS_CI ? "https://srwhite59.github.io/GaussletBases.jl/dev/" : nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => "manual/index.md",
        "Algorithms" => "algorithms/index.md",
        "Examples" => "examples/index.md",
        "Reference" => "reference/index.md",
        "Developer Notes" => "developer/index.md",
        hide("Qiu-White residual-Gaussian route" => "algorithms/qiu_white_residual_gaussian_route.md"),
        hide("1D distorted-gausslet PGDG refinement hierarchy" => "algorithms/distorted_gausslet_pgdg_refinement_hierarchy.md"),
        hide("Cartesian nested face construction" => "algorithms/cartesian_nested_face_construction.md"),
        hide("Cartesian nested atomic nonrecursive route" => "algorithms/cartesian_nested_atomic_nonrecursive_route.md"),
        hide("Radial interval-sampled build and extents" => "algorithms/radial_interval_sampled_build_and_extents.md"),
        hide("First radial workflow" => "tutorials/first_radial_workflow.md"),
        hide("Recommended atomic setup" => "howto/recommended_atomic_setup.md"),
        hide("Current atomic branch" => "explanations/current_atomic_branch.md"),
        hide("Current ordinary branch" => "explanations/current_ordinary_branch.md"),
        hide("Example guide" => "howto/example_guide.md"),
        hide("Bases and mappings" => "reference/bases_and_mappings.md"),
        hide("Operators and diagnostics" => "reference/operators_and_diagnostics.md"),
        hide("Atomic and ordinary workflows" => "reference/atomic_and_ordinary.md"),
        hide("Export layer" => "reference/export.md"),
        hide("Architecture and current direction" => "developer/architecture.md"),
        hide("Supporting note map" => "developer/supporting_notes.md"),
    ],
)

deploydocs(
    repo = "github.com/srwhite59/GaussletBases.jl.git",
    devbranch = "main",
)
