const _DOCS_SITE = "https://srwhite59.github.io/GaussletBases.jl"
const _DOCS_VERSION_TAG =
    r"^v[0-9]+\.[0-9]+\.[0-9]+(?:-[0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*)?$"
const _DOCS_VERSIONS = [
    "stable" => "v^",
    "v#.#",
    "v0.2.0-rc1" => "v0.2.0-rc1",
    "dev" => "dev",
]

function _documentation_target(event::String, ref::String)
    isempty(event) && isempty(ref) && return (:local, nothing)
    event == "pull_request" && return (:pull_request, "dev")
    event == "push" || error("unsupported documentation event: $(repr(event))")
    ref == "refs/heads/main" && return (:dev, "dev")
    prefix = "refs/tags/"
    startswith(ref, prefix) || error("unsupported documentation push ref: $(repr(ref))")
    tag = chop(ref; head = length(prefix), tail = 0)
    occursin('+', tag) && error("documentation tags may not contain build metadata: $(repr(tag))")
    match(_DOCS_VERSION_TAG, tag) === nothing &&
        error("documentation tag is not an exact semantic version: $(repr(tag))")
    version = tryparse(VersionNumber, tag[2:end])
    (isnothing(version) || tag != "v$(version)") &&
        error("documentation tag is not canonical: $(repr(tag))")
    return (:tag, tag)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
include(joinpath(@__DIR__, "check_manager_log.jl"))
ManagerLogPolicy.check_live_log()

using Documenter

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using GaussletBases

const DOCS_CI = get(ENV, "CI", "false") == "true"
const DOCS_DEPLOY = get(ENV, "GAUSSLETBASES_DOCS_DEPLOY", "false") == "true"
const DOCS_TARGET = _documentation_target(
    get(ENV, "GITHUB_EVENT_NAME", ""), get(ENV, "GITHUB_REF", ""))
DOCS_DEPLOY && DOCS_TARGET[1] ∉ (:dev, :tag) &&
    error("documentation deployment is not allowed for context $(DOCS_TARGET[1])")
const DOCS_CANONICAL = DOCS_CI && !isnothing(DOCS_TARGET[2]) ?
                       "$(_DOCS_SITE)/$(DOCS_TARGET[2])/" : nothing
@info "Documentation build context" context = DOCS_TARGET[1] canonical = DOCS_CANONICAL deploy = DOCS_DEPLOY

makedocs(
    sitename = "GaussletBases.jl",
    modules = [GaussletBases],
    doctest = true,
    checkdocs = :none,
    format = Documenter.HTML(
        prettyurls = DOCS_CI,
        edit_link = "main",
        canonical = DOCS_CANONICAL,
        size_threshold_ignore = [
            "developer/designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_through_pass_379.md",
            "developer/designs/cartesian_hamiltonian_producer/registry.md",
        ],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Overview" => "manual/index.md",
            "Projected q-shells (PQS)" => "manual/projected_q_shells.md",
            "Reference-density Hartree screening" => "manual/reference_density_hartree_screening.md",
        ],
        "Algorithms" => [
            "Overview" => "algorithms/index.md",
            "Cartesian PQS and IDA overview" => "algorithms/cartesian_ida_overview.md",
            "PQS shell construction" => "algorithms/pqs_shell_construction.md",
            "IDA Hamiltonian and counterpoise" => "algorithms/ida_hamiltonian_and_counterpoise.md",
        ],
        "Examples" => "examples/index.md",
        "Reference" => "reference/index.md",
        "Developer Notes" => "developer/index.md",
        hide("Qiu-White residual-Gaussian route" => "algorithms/qiu_white_residual_gaussian_route.md"),
        hide("Cartesian low-dimensional operator assembly" => "algorithms/cartesian_low_dimensional_operator_assembly.md"),
        hide("Residual-Gaussian extension" => "algorithms/residual_gaussian_extension.md"),
        hide("Atomic IDA exchange angular-sector rule" => "algorithms/atomic_ida_exchange_angular_sectors.md"),
        hide("1D distorted-gausslet PGDG refinement hierarchy" => "algorithms/distorted_gausslet_pgdg_refinement_hierarchy.md"),
        hide("Cartesian nested face construction" => "algorithms/cartesian_nested_face_construction.md"),
        hide("Cartesian nested atomic nonrecursive route" => "algorithms/cartesian_nested_atomic_nonrecursive_route.md"),
        hide("Cartesian nested diatomic box policy" => "algorithms/cartesian_nested_diatomic_box_policy.md"),
        hide("Cartesian nested endcap/panel shared-shell route" => "algorithms/cartesian_nested_endcap_panel_shared_shell.md"),
        hide("Cartesian nested diatomic coordinate distortion" => "algorithms/cartesian_nested_diatomic_coordinate_distortion.md"),
        hide("Radial interval-sampled build and extents" => "algorithms/radial_interval_sampled_build_and_extents.md"),
        hide("First radial workflow" => "tutorials/first_radial_workflow.md"),
        hide("Recommended atomic setup" => "howto/recommended_atomic_setup.md"),
        hide("Visualization utilities" => "howto/visualization.md"),
        hide("Current atomic branch" => "explanations/current_atomic_branch.md"),
        hide("Angular research track" => "explanations/angular_research_track.md"),
        hide("Current ordinary branch" => "explanations/current_ordinary_branch.md"),
        hide("Example guide" => "howto/example_guide.md"),
        hide("Bases and mappings" => "reference/bases_and_mappings.md"),
        hide("Operators and diagnostics" => "reference/operators_and_diagnostics.md"),
        hide("Atomic and ordinary workflows" => "reference/atomic_and_ordinary.md"),
        hide("Export layer" => "reference/export.md"),
        hide("Architecture and current direction" => "developer/architecture.md"),
        hide("Gausslet methods fundamentals" => "developer/architecture/gausslet_methods_fundamentals.md"),
        hide("Gausslet algorithm refresher" => "developer/architecture/gausslet_algorithm_refresher.md"),
        hide("Cartesian route migration" => "developer/cartesian/route_migration.md"),
        hide("Cartesian feature donor inventory" => "developer/cartesian/feature_donor_inventory.md"),
        hide("PQS thin route demolition history" => "developer/archive/pqs_thin_route_demolition_history.md"),
        hide("Numerical contracts" => "developer/numerical_contracts.md"),
        hide("Supporting note map" => "developer/supporting_notes.md"),
    ],
)

if DOCS_DEPLOY
    deploydocs(
        repo = "github.com/srwhite59/GaussletBases.jl.git",
        devbranch = "main",
        versions = _DOCS_VERSIONS,
    )
elseif DOCS_CI
    @info "Documentation build completed without deployment." context = DOCS_TARGET[1]
end
end
