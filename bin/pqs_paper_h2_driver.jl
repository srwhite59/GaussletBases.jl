#!/usr/bin/env julia
using GaussletBases, LinearAlgebra, SHA
const GB = GaussletBases
const INPUTS = (:system, :method, :R, :output_dir, :padding)
const FIELDS = split("""
system method R_bohr padding_bohr output_dir git_head git_dirty route route_q ns core_spacing s_factor mapping_s_standard mapping_s_effective mapping_types tail_spacing source_span
coulomb_accuracy coulomb_terms coulomb_fingerprint parent_axes parent_bounds_bohr actual_transverse_padding_bohr actual_parallel_padding_bohr parent_dimension parent_fingerprint_sha256
terminal_row_count topology direct_core_columns complete_shell_columns slab_columns compact_product_columns identity_columns final_dimension T_expectation_Ha U_left_expectation_Ha U_right_expectation_Ha H1_expectation_Ha
electronic_energy_Ha second_eigenvalue_Ha eigen_gap_Ha nuclear_repulsion_Ha total_energy_Ha eigen_residual_Ha H1_symmetry_error_Ha construction_elapsed_s construction_allocated_bytes kinetic_elapsed_s kinetic_allocated_bytes nuclear_elapsed_s nuclear_allocated_bytes solve_elapsed_s solve_allocated_bytes total_elapsed_s total_allocated_bytes peak_rss peak_rss_units tsv_path report_path""")
values = Dict{Symbol,Any}(:system => :h2plus, :method => :both, :R => 2.0, :output_dir => "/tmp/pqs_paper_h2plus_R2", :padding => 10.0)
function apply_inputs!(input)
    for (key, value) in pairs(input); name = Symbol(key)
        name in INPUTS || throw(ArgumentError("unknown driver input: $name")); values[name] = value; end
end
args = collect(ARGS)
if !isempty(args) && !occursin("=", first(args))
    loaded = Base.include(Main, abspath(popfirst!(args)))
    loaded isa NamedTuple || loaded isa AbstractDict || throw(ArgumentError("trusted config must return a NamedTuple or AbstractDict")); apply_inputs!(loaded)
end
for arg in args
    parts = split(arg, "="; limit = 2); length(parts) == 2 || throw(ArgumentError("expected name=value input"))
    apply_inputs!(Dict(Symbol(parts[1]) => Core.eval(Main, Meta.parse(parts[2])))); end
system = Symbol(values[:system]); system === :h2 && throw(ArgumentError("system=:h2 is not enabled in this H2+ gate"))
system === :h2plus || throw(ArgumentError("system must be :h2plus or :h2"))
method = Symbol(values[:method]); method in (:pqs, :wl, :both) || throw(ArgumentError("method must be :pqs, :wl, or :both"))
R, padding = Float64(values[:R]), Float64(values[:padding])
isfinite(R) && R > 0 && isfinite(padding) && padding >= 10 || throw(ArgumentError("R must be positive and padding at least 10 bohr"))
output_dir = abspath(String(values[:output_dir])); isempty(output_dir) && throw(ArgumentError("output_dir must not be empty")); mkpath(output_dir)
tsv_path, report_path = joinpath(output_dir, "h2plus.tsv"), joinpath(output_dir, "report.txt")
git_head, git_dirty = readchomp(`git rev-parse HEAD`), !isempty(readchomp(`git status --porcelain`))
nuclei = NTuple{3,Float64}[(0.0, 0.0, -R / 2), (0.0, 0.0, R / 2)]
expansion = GB._cartesian_base_resolve_coulomb_expansion(:high); coulomb = GB._cartesian_coulomb_expansion_summary(:high, expansion)
axis_data(axis) = (axis.basis.reference_center_data, axis.basis.center_data, axis.basis.integral_weight_data,
    axis.basis.coefficient_matrix, axis.centers, axis.weights, axis.overlap, axis.kinetic, axis.position, axis.x2)
function parent_fingerprint(pgdg)
    bytes = vcat((reinterpret(UInt8, vec(array)) for axis in pgdg for array in axis_data(axis))...)
    return bytes2hex(sha256(bytes))
end
function run_method(nesting)
    q = nesting === :pqs ? 5 : 3
    built = @timed begin
        system_stage = GB.cartesian_system((; atom_symbols = (:H, :H), nuclear_charges = (1.0, 1.0),
            atom_locations = Tuple(nuclei), nup = 1, ndn = 0, bond_axis = :z, bond_length = R,
            radius = R / 2 + padding, parent_axis_counts = nothing, map_backend = :pgdg_localized_experimental))
        recipe = GB.cartesian_recipe(GB._cartesian_base_route(:h2, nesting, :ordinary))
        spacing = (; q, n_s = 5, reference_spacing = 1.0, tail_spacing = 2.8,
            q_to_core_spacing_rule = :standard_pqs_ns_equals_q, core_spacing = 0.30, s_factor = 1.0,
            xmax_parallel = R / 2 + padding, xmax_transverse = padding)
        parent_inputs = (; parent_axis_bundle_backend = :pgdg_localized_experimental, parent_axis_family = :G10, coulomb_expansion = expansion)
        parent = GB.cartesian_parent(system_stage, spacing, parent_inputs, recipe)
        shells = GB.cartesian_shells(parent, spacing, recipe); units = GB.cartesian_units(parent, shells, recipe)
        transforms = GB.cartesian_transforms(units, recipe)
        setup = parent.standard_setup
        due_input = (; kind = :h2, nesting, source_span = :ordinary, symbols = ["H", "H"],
            charges = [1.0, 1.0], nup = 1, ndn = 0, locations = nuclei, radius = nothing,
            xmax_parallel = setup.xmax_parallel, xmax_transverse = setup.xmax_transverse,
            core_spacing = setup.core_spacing, reference_spacing = setup.reference_spacing, tail_spacing = setup.tail_spacing,
            q = setup.q, ns = setup.n_s)
        due = GB._cartesian_terminal_due_diligence_report(due_input, parent, transforms)
        (; parent, transforms, due)
    end
    parent, transforms, due = built.value.parent, built.value.transforms, built.value.due
    setup, basis = parent.standard_setup, transforms.terminal_basis_realization
    isnothing(basis) && error("terminal basis was not realized")
    rows = due.terminal_rows; isempty(rows) && error("complete terminal due-diligence rows are unavailable")
    required = (:region_key, :region_kind, :lowering_kind, :retained_unit_kind, :realization_kind, :realization_class, :shell_index, :index_ranges, :physical_ranges, :support_rows, :final_cols, :terminal_order, :retained_count, :final_column_range, :slab_axis, :slab_side, :slab_thickness, :slab_stack_index, :slab_stack_count, :warning_flags, :warning_summary)
    all(row -> all(name -> hasproperty(row, name), required), rows) || error("terminal due-diligence row structure is incomplete")
    hard_flags = (:missing_shell_index, :missing_physical_bounds, :slab_without_native_metadata)
    all(row -> all(flag -> flag ∉ hard_flags, row.warning_flags) && all(range -> all(isfinite, range) && first(range) <= last(range), (row.physical_ranges.x, row.physical_ranges.y, row.physical_ranges.z)), rows) || error("terminal due diligence lacks required geometry")
    all(i -> rows[i].terminal_order == i, eachindex(rows)) || error("terminal due-diligence order is not native")
    all(row -> row.final_cols == row.retained_count == length(row.final_column_range), rows) || error("terminal row dimensions disagree")
    collect(Iterators.flatten(row.final_column_range for row in rows)) == collect(1:basis.final_dimension) || error("terminal columns do not partition the basis")
    expected_family = nesting === :pqs ? :pqs_source_box : :white_lindsey_low_order; transforms.route_family === expected_family || error("staged route family differs")
    expected_realization = nesting === :pqs ? :pqs_source_modes_boundary_selection_shell_realization_contract : :white_lindsey_boundary_stratum_product_contract
    all(row -> row.region_kind !== :complete_shell || row.realization_kind === expected_realization, rows) || error("complete-shell realization differs from route")
    route = nesting === :pqs ? :pqs_source_box_first : :white_lindsey_product; setup.q == q || error("route-local q gate failed")
    expected_s = (sqrt(0.30), sqrt(0.30))
    (setup.n_s, setup.core_spacing, setup.s_factor, setup.tail_spacing, setup.mapping_s_standard, setup.mapping_s_effective) == (5, 0.30, 1.0, 2.8, expected_s, expected_s) || error("resolved parent setup differs")
    parent_basis = parent.parent_basis_object; axes = GB.CartesianParentGaussletBases.parent_axes(parent_basis)
    axis_counts = GB.CartesianParentGaussletBases.parent_axis_counts(parent_basis); parent_dimension = GB.CartesianParentGaussletBases.parent_dimension(parent_basis)
    pgdg = Tuple(GB._nested_axis_pgdg(parent.parent_axis_bundle_object, axis) for axis in (:x, :y, :z))
    bounds = (; x = extrema(pgdg[1].centers), y = extrema(pgdg[2].centers), z = extrema(pgdg[3].centers))
    mappings = (GB.mapping(axes.x), GB.mapping(axes.y), GB.mapping(axes.z))
    all(mapping -> mapping isa GB.CombinedInvsqrtMapping, mappings) || error("multicenter parent did not use combined inverse-sqrt mappings")
    all(mapping -> mapping.tail_spacing == setup.tail_spacing, mappings) || error("live mapping tail spacings disagree")
    due.geometry.parent_axis_counts == axis_counts && due.geometry.parent_physical_bounds == bounds || error("due diligence differs from live parent")
    kinetic = @timed GB._pqs_source_box_route_driver_terminal_products(basis, pgdg)
    nuclear = @timed GB._pqs_source_box_route_driver_terminal_unit_nuclear(basis, expansion, pgdg, nuclei)
    H1 = copy(kinetic.value.kinetic); foreach(unit -> H1 .+= unit, nuclear.value)
    all(isfinite, H1) || error("H1 contains nonfinite values")
    symmetry = norm(H1 - transpose(H1), Inf); symmetry <= 1.0e-10 || error("H1 symmetry gate failed")
    solved = @timed eigen(Symmetric(H1))
    electronic, second, orbital = solved.value.values[1], solved.value.values[2], solved.value.vectors[:, 1]
    residual = norm(H1 * orbital - electronic * orbital); residual <= 1.0e-9 || error("lowest-state residual gate failed")
    expectation(A) = dot(orbital, A, orbital); T_expectation = expectation(kinetic.value.kinetic)
    U_left, U_right = expectation.(nuclear.value); H1_expectation = expectation(H1)
    isapprox(H1_expectation, T_expectation + U_left + U_right; atol = 1.0e-10, rtol = 0) || error("one-body decomposition does not close")
    dimensions = due.dimensions; transverse_padding = minimum((-bounds.x[1], bounds.x[2], -bounds.y[1], bounds.y[2]))
    parallel_padding = (nuclei[1][3] - bounds.z[1], bounds.z[2] - nuclei[2][3])
    elapsed = built.time + kinetic.time + nuclear.time + solved.time; allocated = built.bytes + kinetic.bytes + nuclear.bytes + solved.bytes
    row = Dict{String,String}(
        "system" => "h2plus", "method" => String(nesting), "R_bohr" => repr(R), "padding_bohr" => repr(padding),
        "output_dir" => output_dir, "git_head" => git_head, "git_dirty" => string(git_dirty),
        "route" => String(route), "route_q" => string(setup.q),
        "ns" => string(setup.n_s), "core_spacing" => repr(setup.core_spacing), "s_factor" => repr(setup.s_factor),
        "mapping_s_standard" => repr(setup.mapping_s_standard), "mapping_s_effective" => repr(setup.mapping_s_effective),
        "mapping_types" => join(string.(nameof.(typeof.(mappings))), ","),
        "tail_spacing" => repr(setup.tail_spacing), "source_span" => String(due.geometry.source_span),
        "coulomb_accuracy" => String(coulomb.policy), "coulomb_terms" => string(coulomb.term_count),
        "coulomb_fingerprint" => string(GB._coulomb_expansion_fingerprint(expansion)), "parent_axes" => join(axis_counts, "x"),
        "parent_bounds_bohr" => repr(bounds), "actual_transverse_padding_bohr" => repr(transverse_padding),
        "actual_parallel_padding_bohr" => repr(parallel_padding), "parent_dimension" => string(parent_dimension),
        "parent_fingerprint_sha256" => parent_fingerprint(pgdg),
        "terminal_row_count" => string(length(due.terminal_rows)), "topology" => join(String.(unique(row.region_kind for row in due.terminal_rows)), ","),
        "direct_core_columns" => string(dimensions.direct_core_columns), "complete_shell_columns" => string(dimensions.complete_shell_columns),
        "slab_columns" => string(dimensions.slab_columns), "compact_product_columns" => string(dimensions.compact_product_columns), "identity_columns" => string(dimensions.identity_columns),
        "final_dimension" => string(basis.final_dimension), "T_expectation_Ha" => repr(T_expectation),
        "U_left_expectation_Ha" => repr(U_left), "U_right_expectation_Ha" => repr(U_right), "H1_expectation_Ha" => repr(H1_expectation),
        "electronic_energy_Ha" => repr(electronic), "second_eigenvalue_Ha" => repr(second),
        "eigen_gap_Ha" => repr(second - electronic), "nuclear_repulsion_Ha" => repr(inv(R)),
        "total_energy_Ha" => repr(electronic + inv(R)), "eigen_residual_Ha" => repr(residual),
        "H1_symmetry_error_Ha" => repr(symmetry), "construction_elapsed_s" => repr(built.time),
        "construction_allocated_bytes" => string(built.bytes), "kinetic_elapsed_s" => repr(kinetic.time), "kinetic_allocated_bytes" => string(kinetic.bytes),
        "nuclear_elapsed_s" => repr(nuclear.time), "nuclear_allocated_bytes" => string(nuclear.bytes), "solve_elapsed_s" => repr(solved.time),
        "solve_allocated_bytes" => string(solved.bytes), "total_elapsed_s" => repr(elapsed),
        "total_allocated_bytes" => string(allocated), "peak_rss" => string(Sys.maxrss()),
        "peak_rss_units" => "diagnostic bytes (combined process; Julia Sys.maxrss; $(Sys.KERNEL))",
        "tsv_path" => tsv_path, "report_path" => report_path)
    return (; row, due, pgdg, mappings)
end
requested = method === :both ? (:pqs, :wl) : (method,); results = [run_method(nesting) for nesting in requested]
if length(results) == 2
    results[1].row["parent_bounds_bohr"] == results[2].row["parent_bounds_bohr"] || error("independent parent bounds differ")
    results[1].row["parent_fingerprint_sha256"] == results[2].row["parent_fingerprint_sha256"] || error("parent fingerprints differ")
    all(axis_data(a) == axis_data(b) for (a, b) in zip(results[1].pgdg, results[2].pgdg)) || error("parent numerical objects differ")
end
open(tsv_path, "w") do io
    println(io, join(FIELDS, '\t'))
    foreach(result -> println(io, join((result.row[field] for field in FIELDS), '\t')), results)
end
open(report_path, "w") do io
    println(io, "Matched source-box-first PQS / White-Lindsey two-center H2+ one-body gate")
    for result in results
        println(io, "\n[", result.row["method"], "]")
        foreach(field -> println(io, field, " = ", result.row[field]), FIELDS)
        println(io, "parent_mappings = ", result.mappings, "\nterminal_due_diligence =")
        show(IOContext(io, :limit => false, :compact => false, :displaysize => (10000, 120)), MIME"text/plain"(), result.due); println(io)
    end
end
println("TSV: ", tsv_path, "\nreport: ", report_path)
