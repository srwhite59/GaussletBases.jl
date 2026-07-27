#!/usr/bin/env julia
using GaussletBases, LinearAlgebra, SHA
const GB = GaussletBases
const INPUTS = (:system, :method, :R, :output_dir, :padding, :tail_spacing)
const FIELDS = split("""
system method R_bohr padding_bohr output_dir git_head git_dirty route route_q ns core_spacing s_factor mapping_s_standard mapping_s_effective mapping_types tail_spacing source_span
coulomb_accuracy coulomb_terms coulomb_fingerprint parent_axes parent_bounds_bohr actual_transverse_padding_bohr actual_parallel_padding_bohr parent_dimension parent_fingerprint_sha256
terminal_row_count topology direct_core_columns complete_shell_columns slab_columns compact_product_columns identity_columns final_dimension T_expectation_Ha U_left_expectation_Ha U_right_expectation_Ha H1_expectation_Ha
electronic_energy_Ha second_eigenvalue_Ha eigen_gap_Ha nuclear_repulsion_Ha total_energy_Ha eigen_residual_Ha H1_symmetry_error_Ha construction_elapsed_s construction_allocated_bytes kinetic_elapsed_s kinetic_allocated_bytes nuclear_elapsed_s nuclear_allocated_bytes solve_elapsed_s solve_allocated_bytes total_elapsed_s total_allocated_bytes peak_rss peak_rss_units tsv_path report_path overlap_identity_error independent_reference_total_Ha independent_reference_error_Ha parent_resolution_error_Ha contraction_error_Ha
parent_ground_state_norm terminal_capture_fraction terminal_lost_norm capture_closure_error
supplement_fingerprint_sha256 supplement_candidate_count supplement_parent_capture_min_sv supplement_terminal_capture_min_sv residual_occupation_cutoff
residual_dimension residual_min_retained_occupation residual_max_discarded_occupation terminal_residual_orthogonality_error bare_terminal_energy_change_Ha
fixed_target_capture fixed_target_norm_error fixed_target_fingerprint_sha256 fixed_H1_expectation_Ha fixed_density_density_direct_Ha fixed_density_density_exchange_Ha
fixed_density_density_Vee_Ha fixed_electronic_energy_Ha fixed_total_energy_Ha interaction_convention interaction_symmetry_error_Ha interaction_elapsed_s interaction_allocated_bytes""")
const REFERENCE_TOTAL = -0.6026342144949465
values = Dict{Symbol,Any}(:system => :h2plus, :method => :both, :R => 2.0, :output_dir => "/tmp/pqs_paper_h2plus_R2", :padding => 10.0, :tail_spacing => 2.8)
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
system = Symbol(values[:system]); system in (:h2plus, :h2) || throw(ArgumentError("system must be :h2plus or :h2"))
method = Symbol(values[:method]); method in (:pqs, :wl, :both) || throw(ArgumentError("method must be :pqs, :wl, or :both"))
R, padding, tail_spacing = (Float64(values[key]) for key in (:R, :padding, :tail_spacing))
R == 2.0 && (padding, tail_spacing) in ((10.0, 2.8), (20.0, 2.8), (10.0, 2.0)) || throw(ArgumentError("unsupported R, padding, and tail_spacing combination"))
system === :h2 && (method, padding, tail_spacing) != (:both, 10.0, 2.8) && throw(ArgumentError("system=:h2 requires the frozen supplemented preflight"))
capture_enabled = system === :h2plus && method === :both && padding == 10.0 && tail_spacing == 2.8
output_dir = abspath(String(values[:output_dir])); isempty(output_dir) && throw(ArgumentError("output_dir must not be empty")); mkpath(output_dir)
tsv_path, report_path = joinpath(output_dir, system === :h2 ? "h2.tsv" : "h2plus.tsv"), joinpath(output_dir, "report.txt")
git_head, git_dirty = readchomp(`git rev-parse HEAD`), !isempty(readchomp(`git status --porcelain`))
nuclei = NTuple{3,Float64}[(0.0, 0.0, -R / 2), (0.0, 0.0, R / 2)]
expansion = GB._cartesian_base_resolve_coulomb_expansion(:high); coulomb = GB._cartesian_coulomb_expansion_summary(:high, expansion)
axis_data(axis) = (axis.basis.reference_center_data, axis.basis.center_data, axis.basis.integral_weight_data,
    axis.basis.coefficient_matrix, axis.centers, axis.weights, axis.overlap, axis.kinetic, axis.position, axis.x2)
parent_fingerprint(pgdg) = bytes2hex(sha256(vcat((reinterpret(UInt8, vec(array)) for axis in pgdg for array in axis_data(axis))...)))
factor_term(factors, term) = factors isa AbstractArray{<:Real,3} ? @view(factors[term, :, :]) : factors[term]
valid_parent_matrix(matrix) = all(isfinite, matrix) && norm(matrix - transpose(matrix), Inf) <= 1.0e-10
function mode_product(matrix, tensor, axis)
    order = axis == 1 ? (1, 2, 3) : axis == 2 ? (2, 1, 3) : (3, 1, 2)
    moved = PermutedDimsArray(tensor, order)
    product = matrix * reshape(moved, size(tensor, axis), :)
    return permutedims(reshape(product, size(matrix, 1), size(tensor, order[2]), size(tensor, order[3])), invperm(order))
end
function product_apply(matrices, vector, dimensions)
    tensor = reshape(vector, reverse(dimensions)); for axis in 1:3; tensor = mode_product(matrices[4 - axis], tensor, axis); end
    return vec(tensor)
end
function inverse_sqrt(overlap)
    decomposition = eigen(Symmetric(overlap))
    minimum(decomposition.values) > 0 || error("parent overlap is not positive definite")
    root = decomposition.vectors * Diagonal(inv.(sqrt.(decomposition.values))) * transpose(decomposition.vectors)
    identity_error = norm(root * overlap * root - I, Inf); identity_error <= 1.0e-10 || error("parent overlap orthogonalization failed")
    return root, identity_error
end
function fixed_charge(state)
    all(isfinite, state) && !isempty(state) || error("fixed target state is empty or nonfinite")
    q = abs2.(state); total = sum(q); isfinite(total) && total > 0 || error("fixed target charge is invalid")
    q ./= total; abs(sum(q) - 1.0) <= 1.0e-12 || error("fixed target charge does not close"); return q
end
function fixed_fingerprint(state)
    canonical = Float64.(state); pivot = argmax(abs.(canonical))
    canonical[pivot] < 0 && (canonical .*= -1)
    return bytes2hex(sha256(reinterpret(UInt8, canonical)))
end
function terminal_projection(basis, metric_parent)
    coefficients = zeros(Float64, basis.final_dimension)
    for block in basis.blocks
        source = @view metric_parent[block.support_indices]
        coefficients[block.column_range] .= isnothing(block.coefficients) ? source : transpose(block.coefficients) * source
    end
    return coefficients
end
function matrix_interaction(V, state, owner)
    size(V) == (length(state), length(state)) && all(isfinite, V) || error("$owner interaction is malformed")
    symmetry = norm(V - transpose(V), Inf); q = fixed_charge(state); J = dot(q, V, q)
    detail = "$owner dimensions=$(size(V)) fingerprint=$(GB.external_gto_overlap_fingerprint(V)) charge=$(sum(q))"
    return (; J, symmetry, detail)
end
function parent_interaction(parent, state, dimensions)
    terms = ntuple(axis -> parent.pgdg[axis].pair_factor_terms, 3); nterms = length(expansion.coefficients)
    all(axis -> size(terms[axis], 1) == nterms && all(isfinite, terms[axis]), 1:3) ||
        error("parent IDA factor resource is malformed")
    function apply_vee(vector)
        out = zeros(Float64, length(vector))
        for term in 1:nterms
            out .+= expansion.coefficients[term] .* product_apply(
                ntuple(axis -> factor_term(terms[axis], term), 3), vector, dimensions)
        end
        return out
    end
    q = fixed_charge(state); Vq = apply_vee(q)
    probe = normalize!(sin.(0.013 .* (1:length(q))) .+ 0.19 .* cos.(0.031 .* (1:length(q))))
    symmetry = abs(dot(probe, Vq) - dot(apply_vee(probe), q))
    fingerprint = bytes2hex(sha256(vcat((reinterpret(UInt8, vec(term)) for term in terms)...)))
    return (; J = dot(q, Vq), symmetry,
        detail = "parent_high135 dimensions=$(dimensions) factor_fingerprint=$fingerprint charge=$(sum(q))")
end
function record_fixed!(result, state, capture, norm_error, h1, interaction;
    elapsed = interaction.time, allocated = interaction.bytes)
    J, symmetry = interaction.value.J, interaction.value.symmetry
    all(isfinite, (capture, norm_error, h1, J, symmetry)) && -1.0e-10 <= capture <= 1.0 + 1.0e-8 &&
        norm_error <= 1.0e-10 && J >= 0 && symmetry <= 1.0e-10 || error("fixed-state numerical gate failed")
    direct, exchange = 2J, -J; vee = direct + exchange; electronic = 2h1 + vee; total = electronic + inv(R)
    isapprox(vee, J; atol = 1.0e-12, rtol = 0) &&
        isapprox(total, 2h1 + J + inv(R); atol = 1.0e-12, rtol = 0) || error("fixed-state energy accounting does not close")
    merge!(result.row, Dict("fixed_target_capture" => repr(capture), "fixed_target_norm_error" => repr(norm_error),
        "fixed_target_fingerprint_sha256" => fixed_fingerprint(state), "fixed_H1_expectation_Ha" => repr(h1),
        "fixed_density_density_direct_Ha" => repr(direct), "fixed_density_density_exchange_Ha" => repr(exchange),
        "fixed_density_density_Vee_Ha" => repr(vee), "fixed_electronic_energy_Ha" => repr(electronic),
        "fixed_total_energy_Ha" => repr(total), "interaction_convention" => "closed_shell_density_density_2J_minus_J",
        "interaction_symmetry_error_Ha" => repr(symmetry), "interaction_elapsed_s" => repr(elapsed),
        "interaction_allocated_bytes" => string(allocated)))
    push!(result.interaction_details, interaction.value.detail)
end
function run_method(nesting)
    q = nesting === :pqs ? 5 : 3
    built = @timed begin
        system_stage = GB.cartesian_system((; atom_symbols = (:H, :H), nuclear_charges = (1.0, 1.0),
            atom_locations = Tuple(nuclei), nup = 1, ndn = system === :h2 ? 1 : 0, bond_axis = :z, bond_length = R,
            radius = R / 2 + padding, parent_axis_counts = nothing, map_backend = :pgdg_localized_experimental))
        recipe = GB.cartesian_recipe(GB._cartesian_base_route(:h2, nesting, :ordinary))
        spacing = (; q, n_s = 5, reference_spacing = 1.0, tail_spacing,
            q_to_core_spacing_rule = :standard_pqs_ns_equals_q, core_spacing = 0.30, s_factor = 1.0,
            xmax_parallel = R / 2 + padding, xmax_transverse = padding)
        parent_inputs = (; parent_axis_bundle_backend = :pgdg_localized_experimental, parent_axis_family = :G10, coulomb_expansion = expansion)
        parent = GB.cartesian_parent(system_stage, spacing, parent_inputs, recipe)
        shells = GB.cartesian_shells(parent, spacing, recipe); units = GB.cartesian_units(parent, shells, recipe)
        transforms = GB.cartesian_transforms(units, recipe)
        setup = parent.standard_setup
        due_input = (; kind = :h2, nesting, source_span = :ordinary, symbols = ["H", "H"],
            charges = [1.0, 1.0], nup = 1, ndn = system === :h2 ? 1 : 0, locations = nuclei, radius = nothing,
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
    (setup.n_s, setup.core_spacing, setup.s_factor, setup.tail_spacing, setup.mapping_s_standard, setup.mapping_s_effective) == (5, 0.30, 1.0, tail_spacing, expected_s, expected_s) || error("resolved parent setup differs")
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
        "system" => String(system), "method" => String(nesting), "R_bohr" => repr(R), "padding_bohr" => repr(padding),
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
        "tsv_path" => tsv_path, "report_path" => report_path, "overlap_identity_error" => "not_applicable",
        "independent_reference_total_Ha" => repr(REFERENCE_TOTAL), "independent_reference_error_Ha" => repr(electronic + inv(R) - REFERENCE_TOTAL),
        "parent_resolution_error_Ha" => "unavailable", "contraction_error_Ha" => "unavailable",
        "parent_ground_state_norm" => "not_applicable", "terminal_capture_fraction" => "not_applicable",
        "terminal_lost_norm" => "not_applicable", "capture_closure_error" => "not_applicable")
    foreach(field -> row[field] = "not_applicable", FIELDS[(end - 22):end])
    return (; row, due, pgdg, mappings, basis, kinetic = kinetic.value.kinetic,
        unit_nuclear = nuclear.value, H1, capture_regions = String[], supplement_details = String[], interaction_details = String[])
end
requested = method === :both ? (:pqs, :wl) : (method,); results = [run_method(nesting) for nesting in requested]
if length(results) == 2
    pqs, wl = results
    pqs.row["parent_bounds_bohr"] == wl.row["parent_bounds_bohr"] || error("independent parent bounds differ")
    pqs.row["parent_fingerprint_sha256"] == wl.row["parent_fingerprint_sha256"] || error("parent fingerprints differ")
    all(axis_data(a) == axis_data(b) for (a, b) in zip(pqs.pgdg, wl.pgdg)) || error("parent numerical objects differ")
end
function run_parent(template)
    pgdg = template.pgdg; dimensions = ntuple(axis -> length(pgdg[axis].centers), 3)
    all(d -> d > 0 && isodd(d), dimensions) && prod(dimensions) == template.due.dimensions.parent_grid_size || error("invalid full-parent dimensions")
    all(axis -> all(isfinite, axis.weights) && all(>(0), axis.weights), pgdg) || error("parent weights must be finite and strictly positive")
    prepared = @timed begin
        overlaps = ntuple(axis -> pgdg[axis].overlap, 3); kinetic = ntuple(axis -> pgdg[axis].kinetic, 3)
        roots_and_errors = ntuple(axis -> inverse_sqrt(overlaps[axis]), 3); roots = ntuple(axis -> roots_and_errors[axis][1], 3)
        factors = ntuple(center -> ntuple(axis -> GB._pqs_source_box_route_driver_centered_factor_terms(
            pgdg[axis], expansion, nuclei[center][axis]), 3), 2)
        for triple in (overlaps, kinetic, roots), matrix in triple
            valid_parent_matrix(matrix) || error("parent axis operator is nonfinite or asymmetric")
        end
        for center in factors, axis in center, term in eachindex(expansion.coefficients)
            valid_parent_matrix(factor_term(axis, term)) || error("parent nuclear factor is nonfinite or asymmetric")
        end
        (; overlaps, kinetic, roots, overlap_error = maximum(last.(roots_and_errors)), factors)
    end
    overlaps, kinetic = prepared.value.overlaps, prepared.value.kinetic
    roots, factors = prepared.value.roots, prepared.value.factors
    orthogonalize(vector) = product_apply(roots, vector, dimensions)
    raw_kinetic(vector) = product_apply((kinetic[1], overlaps[2], overlaps[3]), vector, dimensions) +
        product_apply((overlaps[1], kinetic[2], overlaps[3]), vector, dimensions) +
        product_apply((overlaps[1], overlaps[2], kinetic[3]), vector, dimensions)
    function raw_nuclear(vector, center)
        out = zeros(Float64, prod(dimensions)); for term in eachindex(expansion.coefficients)
            out .-= expansion.coefficients[term] .* product_apply(
                ntuple(axis -> factor_term(factors[center][axis], term), 3), vector, dimensions)
        end
        return out
    end
    orthogonal_component(raw, vector) = orthogonalize(raw(orthogonalize(vector)))
    raw_h(vector) = raw_kinetic(vector) + raw_nuclear(vector, 1) + raw_nuclear(vector, 2)
    apply_h!(out, vector) = copyto!(out, orthogonal_component(raw_h, vector))
    solved = @timed GB._lanczos_ground_state_apply(apply_h!, prod(dimensions); tol = 1.0e-10)
    orbital, electronic = solved.value.vector, solved.value.value; reapplied = similar(orbital); apply_h!(reapplied, orbital)
    residual = norm(reapplied - electronic .* orbital); residual <= 1.0e-9 || error("recomputed full-parent residual gate failed: $residual")
    kinetic_result = @timed orthogonal_component(raw_kinetic, orbital)
    nuclear_results = @timed ntuple(center -> orthogonal_component(v -> raw_nuclear(v, center), orbital), 2)
    T_expectation = dot(orbital, kinetic_result.value); U_left, U_right = ntuple(center -> dot(orbital, nuclear_results.value[center]), 2)
    H1_expectation = dot(orbital, reapplied)
    isapprox(H1_expectation, T_expectation + U_left + U_right; atol = 1.0e-10, rtol = 0) || error("parent one-body decomposition does not close")
    isapprox(U_left, U_right; atol = 1.0e-10, rtol = 0) || error("parent nuclear expectations break inversion symmetry")
    probe_a = normalize!(sin.(0.017 .* (1:prod(dimensions))) .+ 0.31 .* cos.(0.029 .* (1:prod(dimensions))))
    probe_b = normalize!(cos.(0.023 .* (1:prod(dimensions))) .- 0.27 .* sin.(0.037 .* (1:prod(dimensions))))
    h_a, h_b = similar(probe_a), similar(probe_b); apply_h!(h_a, probe_a); apply_h!(h_b, probe_b)
    symmetry = abs(dot(probe_a, h_b) - dot(h_a, probe_b)); symmetry <= 1.0e-10 || error("full-parent matrix-free H1 symmetry gate failed")
    row = copy(template.row)
    foreach(field -> row[field] = "not_applicable", ("terminal_row_count", "direct_core_columns", "complete_shell_columns", "slab_columns", "compact_product_columns", "identity_columns"))
    merge!(row, Dict(
        "method" => "parent", "route" => "full_pgdg_parent", "route_q" => "not_applicable",
        "topology" => "full_parent", "final_dimension" => string(prod(dimensions)),
        "T_expectation_Ha" => repr(T_expectation), "U_left_expectation_Ha" => repr(U_left),
        "U_right_expectation_Ha" => repr(U_right), "H1_expectation_Ha" => repr(H1_expectation),
        "electronic_energy_Ha" => repr(electronic), "second_eigenvalue_Ha" => "unavailable",
        "eigen_gap_Ha" => "unavailable", "total_energy_Ha" => repr(electronic + inv(R)),
        "eigen_residual_Ha" => repr(residual), "H1_symmetry_error_Ha" => repr(symmetry),
        "construction_elapsed_s" => repr(prepared.time), "construction_allocated_bytes" => string(prepared.bytes),
        "kinetic_elapsed_s" => repr(kinetic_result.time), "kinetic_allocated_bytes" => string(kinetic_result.bytes),
        "nuclear_elapsed_s" => repr(nuclear_results.time), "nuclear_allocated_bytes" => string(nuclear_results.bytes),
        "solve_elapsed_s" => repr(solved.time), "solve_allocated_bytes" => string(solved.bytes),
        "total_elapsed_s" => repr(prepared.time + solved.time + kinetic_result.time + nuclear_results.time),
        "total_allocated_bytes" => string(prepared.bytes + solved.bytes + kinetic_result.bytes + nuclear_results.bytes),
        "peak_rss" => string(Sys.maxrss()), "overlap_identity_error" => repr(prepared.value.overlap_error)))
    capture_enabled && (row["parent_ground_state_norm"] = repr(dot(orbital, orbital))); return (; row, due = nothing, pgdg, mappings = template.mappings, orbital, overlaps, roots, capture_regions = String[], supplement_details = String[], interaction_details = String[])
end
if method === :both
    parent = run_parent(results[1]); results = vcat([parent], results)
    parent_total = parse(Float64, parent.row["total_energy_Ha"])
    for result in results
        total = parse(Float64, result.row["total_energy_Ha"])
        result.row["independent_reference_error_Ha"] = repr(total - REFERENCE_TOTAL)
        result.row["parent_resolution_error_Ha"] = repr(parent_total - REFERENCE_TOTAL)
        result.row["contraction_error_Ha"] = result === parent ? "not_applicable" : repr(total - parent_total)
        total + 1.0e-10 >= parent_total || error("terminal energy is below the full-parent energy")
    end
    if capture_enabled
        parent_norm = dot(parent.orbital, parent.orbital); dimensions = ntuple(axis -> length(parent.pgdg[axis].centers), 3)
        sqrt_factors = ntuple(axis -> parent.overlaps[axis] * parent.roots[axis], 3); metric_orbital = product_apply(sqrt_factors, parent.orbital, dimensions)
        block_amplitudes(block, vector, indices) = isnothing(block.coefficients) ? vector[indices] : transpose(block.coefficients) * vector[indices]
        region_key(row) = (row.region_key, row.region_kind, row.shell_index, row.slab_axis, row.slab_side, row.slab_stack_index, row.slab_stack_count)
        for result in results[2:3]
            blocks, rows = result.basis.blocks, result.due.terminal_rows
            length(blocks) == length(rows) || error("terminal blocks and due-diligence rows differ")
            sort(reduce(vcat, (block.support_indices for block in blocks))) == collect(1:prod(dimensions)) || error("terminal block supports do not partition the parent")
            all(i -> rows[i].terminal_order == i && rows[i].support_rows == length(blocks[i].support_indices) && rows[i].final_column_range == blocks[i].column_range, eachindex(blocks)) || error("terminal blocks and due-diligence native ranges differ")
            result.row["parent_fingerprint_sha256"] == parent.row["parent_fingerprint_sha256"] || error("terminal capture does not use the solved parent")
            parent_coefficients = zeros(Float64, prod(dimensions))
            regional = Dict{Tuple,Tuple{Float64,Float64,Int,Int,Set{String}}}()
            for (block, row) in zip(blocks, rows)
                indices = block.support_indices; amplitudes = block_amplitudes(block, metric_orbital, indices)
                parent_coefficients[indices] .= isnothing(block.coefficients) ? amplitudes : block.coefficients * amplitudes
                key = region_key(row); old = get(regional, key, (0.0, 0.0, 0, 0, Set{String}()))
                push!(old[5], row.warning_summary)
                regional[key] = (old[1] + sum(abs2, amplitudes), old[2],
                    old[3] + length(block.support_indices), old[4] + length(block.column_range), old[5])
            end
            residual = parent.orbital - product_apply(sqrt_factors, parent_coefficients, dimensions)
            residual_metric = product_apply(sqrt_factors, residual, dimensions)
            orthogonality = 0.0
            for (block, row) in zip(blocks, rows)
                indices = block.support_indices; amplitudes = block_amplitudes(block, residual_metric, indices)
                orthogonality += sum(abs2, amplitudes)
                key = region_key(row); old = regional[key]
                regional[key] = (old[1], old[2] + sum(abs2, residual[indices]), old[3], old[4], old[5])
            end
            capture = sum(value[1] for value in Base.values(regional)); lost = sum(abs2, residual); closure = abs(capture + lost - parent_norm)
            -1.0e-10 <= capture <= parent_norm + 1.0e-10 || error("terminal capture is outside the parent norm")
            sqrt(orthogonality) <= 1.0e-9 || error("capture residual is not terminal-orthogonal")
            closure <= 1.0e-9 || error("terminal capture norm does not close")
            abs(sum(value[2] for value in Base.values(regional)) - lost) <= 1.0e-12 || error("regional loss does not close")
            result.row["parent_ground_state_norm"] = repr(parent_norm); result.row["terminal_capture_fraction"] = repr(capture / parent_norm)
            result.row["terminal_lost_norm"] = repr(lost); result.row["capture_closure_error"] = repr(closure)
            for key in sort!(collect(keys(regional)); by = key -> -regional[key][2])
                value = regional[key]
                push!(result.capture_regions, "region=$(key) support_rows=$(value[3]) retained_cols=$(value[4]) " *
                    "captured_norm=$(value[1]) lost_norm=$(value[2]) lost_fraction=$(iszero(lost) ? 0.0 : value[2] / lost) " *
                    "warnings=$(join(sort!(collect(value[5])), ','))")
            end
        end
    end
    if system === :h2
        h2_system = (; atom_symbols = ["H", "H"], nuclear_charges = [1.0, 1.0], atom_locations = nuclei, nup = 1, ndn = 1)
        basis_input(nesting) = (; ns = 5, nesting, core_spacing = 0.30,
            xmax_parallel = R / 2 + padding, xmax_transverse = padding, s_factor = 1.0, tail_spacing,
            source_span = :ordinary, coulomb_accuracy = :high)
        workings = [@timed GB.cartesian_base_working_basis(h2_system;
            basis = basis_input(nesting), supplemented = true) for nesting in (:pqs, :wl)]
        supplement_spec = (; basis_by_center = ["cc-pVTZ", "cc-pVTZ"], lmax = 1, uncontracted = false, width_filtering = nothing)
        loaded = @timed GB.cartesian_residual_gto_supplement_basis(workings[1].value, supplement_spec)
        supplement, raw = loaded.value, loaded.value.basis; orbitals = raw.orbitals; candidate_count = length(orbitals)
        S_AA = Matrix{Float64}(GB._cartesian_supplement_cross_overlap(raw, raw))
        size(S_AA) == (candidate_count, candidate_count) && all(isfinite, S_AA) &&
            norm(S_AA - transpose(S_AA), Inf) <= 1.0e-10 || error("supplement metric is malformed")
        metric = eigen(Symmetric(S_AA)); minimum(metric.values) > 0 || error("supplement metric is not positive definite")
        Aroot = metric.vectors * Diagonal(inv.(sqrt.(metric.values))) * transpose(metric.vectors)
        ordering = GB.external_gto_ordering_fingerprint(raw); overlap_fingerprint = GB.external_gto_overlap_fingerprint(S_AA)
        fingerprint = bytes2hex(sha256(codeunits(ordering * overlap_fingerprint)))
        labels = String[orbital.label for orbital in orbitals]; centers = NTuple{3,Float64}[Tuple(Float64.(orbital.center)) for orbital in orbitals]
        powers = NTuple{3,Int}[Tuple(Int.(orbital.angular_powers)) for orbital in orbitals]; owners = [only(findall(==(center), nuclei)) for center in centers]
        all(i -> length(orbitals[i].exponents) == length(orbitals[i].coefficients) > 0 && all(isfinite, orbitals[i].exponents) &&
            all(isfinite, orbitals[i].coefficients), eachindex(orbitals)) || error("supplement primitive data are malformed")
        owner_counts = [count(==(owner), owners) for owner in 1:2]
        all(>(0), owner_counts) && sum(owner_counts) == candidate_count || error("supplement ownership does not close")
        C = GB.CartesianFinalBasisRealization; proxy = C._r3a_qw_proxy_layers(workings[1].value.parent.parent_axis_bundle_object)
        parent_cross = GB.CartesianGaussianRawBlocks.gaussian_non_nuclear_overlap_blocks(proxy, C._r3a_qw_supplement(raw)).ga.overlap
        parent_dimensions = ntuple(axis -> length(results[1].pgdg[axis].centers), 3)
        parent_solved = hcat((product_apply(parent.roots,
            product_apply(parent.roots, parent_cross[:, column], parent_dimensions), parent_dimensions) for column in axes(parent_cross, 2))...)
        function capture_singulars(K)
            values = eigvals(Symmetric(0.5 .* (K .+ transpose(K)))); minimum(values) >= -1.0e-10 &&
                maximum(values) <= 1.0 + 1.0e-10 || error("supplement capture spectrum is outside [0,1]"); return sqrt.(max.(values, 0.0))
        end
        parent_capture = capture_singulars(Aroot * transpose(parent_cross) * parent_solved * Aroot); supplemented = Any[]
        parent = results[1]; target_dimensions = ntuple(axis -> length(parent.pgdg[axis].centers), 3)
        parent_state = product_apply(parent.roots, parent.orbital, target_dimensions)
        parent_norm_error = abs(dot(parent_state, product_apply(parent.overlaps, parent_state, target_dimensions)) - 1)
        metric_target = product_apply(ntuple(axis -> parent.overlaps[axis] * parent.roots[axis], 3),
            parent.orbital, target_dimensions)
        parent_vee = @timed parent_interaction(parent, parent_state, target_dimensions)
        record_fixed!(parent, parent_state, 1.0, parent_norm_error,
            parse(Float64, parent.row["H1_expectation_Ha"]), parent_vee)
        for (bare, working_timed) in zip(results[2:3], workings)
            base = working_timed.value; base.terminal_basis.final_dimension == bare.basis.final_dimension ||
                error("supplemented working basis differs from the bare terminal basis")
            same_blocks = length(base.terminal_basis.blocks) == length(bare.basis.blocks) &&
                all(zip(base.terminal_basis.blocks, bare.basis.blocks)) do (live, staged)
                    live.support_indices == staged.support_indices && live.column_range == staged.column_range &&
                        isequal(live.coefficients, staged.coefficients)
                end
            same_blocks || error("supplemented and staged bare terminal coordinates differ")
            due = base.terminal_due_diligence
            !isnothing(due) && !isempty(due.terminal_rows) &&
                collect(Iterators.flatten(row.final_column_range for row in due.terminal_rows)) == collect(1:base.terminal_basis.final_dimension) ||
                error("supplemented terminal due diligence is incomplete")
            base_pgdg = Tuple(GB._nested_axis_pgdg(base.parent.parent_axis_bundle_object, axis) for axis in (:x, :y, :z))
            parent_fingerprint(base_pgdg) == bare.row["parent_fingerprint_sha256"] || error("supplemented and bare parents differ")
            g = terminal_projection(base.terminal_basis, metric_target); bare_capture = dot(g, g)
            bare_state = g ./ sqrt(bare_capture); base_vee = @timed GB.cartesian_base_vee(base)
            bare_ham = GB.cartesian_base_hamiltonian_assembly(
                base, (; kinetic = bare.kinetic), bare.unit_nuclear, base_vee.value)
            bare_interaction = @timed matrix_interaction(base_vee.value, bare_state, "bare_$(bare.row["method"])")
            record_fixed!(bare, bare_state, bare_capture, abs(dot(bare_state, bare_state) - 1),
                dot(bare_state, bare.H1, bare_state), bare_interaction;
                elapsed = base_vee.time + bare_interaction.time, allocated = base_vee.bytes + bare_interaction.bytes)
            augmentation = @timed GB.cartesian_residual_gto_augmentation(base, supplement); residual = augmentation.value
            residual.candidate_count == candidate_count && residual.candidate_labels == labels &&
                residual.candidate_centers == centers && residual.candidate_owner_indices == owners &&
                [count(==(owner), residual.candidate_owner_indices) for owner in 1:2] == owner_counts ||
                error("residual candidate identity differs from the loaded supplement")
            products = @timed GB.cartesian_residual_gto_augmented_products(base, supplement, residual; base_kinetic = bare.kinetic)
            unit_nuclear = @timed GB.cartesian_residual_gto_augmented_unit_nuclear(base, residual, products.value; base_unit_nuclear = bare.unit_nuclear)
            X = products.value.supplement_blocks.mixed.overlap
            products.value.supplement_blocks.self.overlap == S_AA &&
                norm(C._r3a_project_parent_ga(base.terminal_basis, parent_cross) - X, Inf) <= 1.0e-10 || error("route supplement overlap differs")
            terminal_capture = capture_singulars(Aroot * transpose(X) * X * Aroot)
            owner_spectra = [eigvals(Symmetric(S_AA[owners .== owner, owners .== owner] -
                transpose(X[:, owners .== owner]) * X[:, owners .== owner])) for owner in 1:2]
            cutoff = residual.occupation_cutoff; cutoff == 1.0e-6 || error("residual occupation cutoff differs from production policy")
            ranks = count.(Ref(>(cutoff)), owner_spectra)
            ranks == residual.owner_retained_counts || error("owner-local residual ranks differ")
            retained = vcat(([value for value in spectrum if value > cutoff] for spectrum in owner_spectra)...); discarded = vcat(([value for value in spectrum if value <= cutoff] for spectrum in owner_spectra)...)
            isapprox(sort(retained), sort(residual.residual_occupations); atol = 1.0e-10, rtol = 1.0e-8) || error("owner-local retained occupations differ")
            minimum(abs.(vcat(owner_spectra...) .- cutoff)) > 1.0e-10 || error("residual occupation is numerically marginal at the production cutoff")
            orthogonality = norm(residual.T_G + X * residual.T_A, Inf)
            RSR = GB.CartesianResidualGaussians.residual_gaussian_overlap(residual.T_G, residual.T_A, X, S_AA); identity_error = norm(RSR - I, Inf)
            orthogonality <= 1.0e-10 && identity_error <= 5.0e-8 || error("augmented residual metric gate failed")
            nG = residual.base_dimension
            norm(products.value.kinetic[1:nG, 1:nG] - bare.kinetic, Inf) <= 1.0e-10 || error("augmented kinetic G-G block differs")
            maximum(norm(unit_nuclear.value[i][1:nG, 1:nG] - bare.unit_nuclear[i], Inf) for i in 1:2) <= 1.0e-10 ||
                error("augmented nuclear G-G block differs")
            H1 = products.value.kinetic + unit_nuclear.value[1] + unit_nuclear.value[2]; matrices = (products.value.kinetic, unit_nuclear.value..., H1)
            all(matrix -> all(isfinite, matrix) && norm(matrix - transpose(matrix), Inf) <= 1.0e-10, matrices) || error("augmented one-body matrix is nonfinite or asymmetric")
            solved = @timed eigen(Symmetric(H1)); electronic, second = solved.value.values[1:2]
            orbital = solved.value.vectors[:, 1]; eigen_residual = norm(H1 * orbital - electronic * orbital)
            eigen_residual <= 1.0e-9 || error("supplemented lowest-state residual gate failed")
            electronic <= parse(Float64, bare.row["electronic_energy_Ha"]) + 1.0e-10 || error("supplemented energy is not variational")
            expectation(matrix) = dot(orbital, matrix, orbital); T = expectation(products.value.kinetic); U = expectation.(unit_nuclear.value)
            isapprox(electronic, T + sum(U); atol = 1.0e-10, rtol = 0) || error("supplemented one-body decomposition does not close")
            row = copy(bare.row); row["method"] *= "_supplemented"; row["final_dimension"] = string(size(H1, 1))
            row["T_expectation_Ha"], row["U_left_expectation_Ha"], row["U_right_expectation_Ha"] = repr(T), repr(U[1]), repr(U[2])
            row["H1_expectation_Ha"] = repr(electronic); row["electronic_energy_Ha"] = repr(electronic)
            row["second_eigenvalue_Ha"] = repr(second); row["eigen_gap_Ha"] = repr(second - electronic)
            row["total_energy_Ha"] = repr(electronic + inv(R)); row["eigen_residual_Ha"] = repr(eigen_residual)
            row["H1_symmetry_error_Ha"] = repr(norm(H1 - transpose(H1), Inf)); row["overlap_identity_error"] = repr(max(orthogonality, identity_error))
            row["construction_elapsed_s"] = repr(working_timed.time + augmentation.time + loaded.time / 2); row["construction_allocated_bytes"] = string(working_timed.bytes + augmentation.bytes + loaded.bytes ÷ 2)
            row["kinetic_elapsed_s"], row["kinetic_allocated_bytes"] = repr(products.time), string(products.bytes)
            row["nuclear_elapsed_s"], row["nuclear_allocated_bytes"] = repr(unit_nuclear.time), string(unit_nuclear.bytes)
            row["solve_elapsed_s"], row["solve_allocated_bytes"] = repr(solved.time), string(solved.bytes)
            row["total_elapsed_s"] = repr(working_timed.time + augmentation.time + loaded.time / 2 + products.time + unit_nuclear.time + solved.time); row["total_allocated_bytes"] = string(working_timed.bytes + augmentation.bytes + loaded.bytes ÷ 2 + products.bytes + unit_nuclear.bytes + solved.bytes)
            row["peak_rss"] = string(Sys.maxrss()); row["independent_reference_error_Ha"] = repr(electronic + inv(R) - REFERENCE_TOTAL)
            row["contraction_error_Ha"] = "not_applicable"; row["supplement_fingerprint_sha256"] = fingerprint
            row["supplement_candidate_count"] = string(candidate_count); row["supplement_parent_capture_min_sv"] = repr(minimum(parent_capture))
            row["supplement_terminal_capture_min_sv"] = repr(minimum(terminal_capture)); row["residual_occupation_cutoff"] = repr(cutoff)
            row["residual_dimension"] = string(residual.residual_dimension); row["residual_min_retained_occupation"] = repr(minimum(retained))
            row["residual_max_discarded_occupation"] = isempty(discarded) ? "not_applicable" : repr(maximum(discarded))
            row["terminal_residual_orthogonality_error"] = repr(orthogonality)
            row["bare_terminal_energy_change_Ha"] = repr(electronic - parse(Float64, bare.row["electronic_energy_Ha"]))
            operator_fingerprints = GB.external_gto_overlap_fingerprint.(matrices)
            details = ["supplement_identity = $((; input = supplement.input, kind = raw.supplement_kind, source_metadata = raw.metadata, labels, powers, centers, owners, owner_counts, exponents = [o.exponents for o in orbitals], coefficients = [o.coefficients for o in orbitals], normalizations = [o.primitive_normalization for o in orbitals]))",
                "supplement_metric = $((; eigenvalues = metric.values, condition = maximum(metric.values) / minimum(metric.values), fingerprint = overlap_fingerprint))",
                "parent_capture_singular_values = $(parent_capture)", "terminal_capture_singular_values = $(terminal_capture)",
                "owner_residual_spectra = $(owner_spectra)", "owner_retained_ranks = $(ranks)",
                "owner_cutoff_margins = $([minimum(abs.(spectrum .- cutoff)) for spectrum in owner_spectra])",
                "dimensions = $((; bare = nG, residual = residual.residual_dimension, augmented = size(H1, 1)))",
                "operators = $((; dimensions = size.(matrices), fingerprints = operator_fingerprints, gg_kinetic_error = norm(products.value.kinetic[1:nG, 1:nG] - bare.kinetic, Inf), gg_nuclear_error = maximum(norm(unit_nuclear.value[i][1:nG, 1:nG] - bare.unit_nuclear[i], Inf) for i in 1:2)))"]
            a = transpose(parent_cross) * parent_state; r = transpose(residual.T_G) * g + transpose(residual.T_A) * a
            augmented_capture = bare_capture + dot(r, r); augmented_capture + 1.0e-12 >= bare_capture ||
                error("supplemented target capture decreased")
            augmented_state = vcat(g, r) ./ sqrt(augmented_capture)
            augmented_vee = @timed GB.cartesian_residual_gto_augmented_vee(
                base, bare_ham, residual, products.value, unit_nuclear.value)
            gg_error = norm(augmented_vee.value[1:nG, 1:nG] - base_vee.value, Inf)
            gg_error <= 1.0e-10 || error("augmented G-G interaction differs from fresh bare IDA")
            augmented_interaction = @timed matrix_interaction(
                augmented_vee.value, augmented_state, "augmented_$(bare.row["method"])")
            result = (; row, due, pgdg = base_pgdg, mappings = bare.mappings, capture_regions = String[],
                supplement_details = details, interaction_details = String[])
            record_fixed!(result, augmented_state, augmented_capture,
                abs(dot(augmented_state, augmented_state) - 1), dot(augmented_state, H1, augmented_state),
                augmented_interaction; elapsed = augmented_vee.time + augmented_interaction.time,
                allocated = augmented_vee.bytes + augmented_interaction.bytes)
            push!(result.interaction_details, "augmented_block_parity gg_error=$gg_error categories=$((; GG = (nG, nG), GR = (nG, residual.residual_dimension), RR = (residual.residual_dimension, residual.residual_dimension)))")
            push!(supplemented, result)
        end
        results = vcat(results, supplemented)
        length(results) == 5 && [result.row["method"] for result in results] ==
            ["parent", "pqs", "wl", "pqs_supplemented", "wl_supplemented"] || error("H2 preflight row contract failed")
    end
end
open(tsv_path, "w") do io
    println(io, join(FIELDS, '\t')); foreach(result -> println(io, join((result.row[field] for field in FIELDS), '\t')), results)
end
open(report_path, "w") do io
    println(io, system === :h2 ? "Fixed-state supplemented H2 density-density interaction gate" :
        "Matched source-box-first PQS / White-Lindsey two-center H2+ one-body gate")
    for result in results
        println(io, "\n[", result.row["method"], "]")
        foreach(field -> println(io, field, " = ", result.row[field]), FIELDS)
        println(io, "parent_mappings = ", result.mappings)
        if !isnothing(result.due); println(io, "terminal_due_diligence =")
            show(IOContext(io, :limit => false,
                :compact => false, :displaysize => (10000, 120)), MIME"text/plain"(), result.due); println(io)
        end
        if !isempty(result.capture_regions); println(io, "capture_regions_by_lost_norm =")
            foreach(line -> println(io, line), result.capture_regions); end
        if !isempty(result.supplement_details); println(io, "supplement_diagnostics =")
            foreach(line -> println(io, line), result.supplement_details); end
        if !isempty(result.interaction_details); println(io, "interaction_diagnostics =")
            foreach(line -> println(io, line), result.interaction_details); end
    end
end
println("TSV: ", tsv_path, "\nreport: ", report_path)
