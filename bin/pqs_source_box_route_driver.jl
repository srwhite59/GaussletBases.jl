#!/usr/bin/env julia

using Dates, JLD2
using GaussletBases

route_kind = :be2_pqs_source_box_development_spine
atom_symbols = ("Be", "Be")
nuclear_charges = (4, 4)
atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))
radius = 15.0
parent_axis_counts = (x = 9, y = 7, z = 9)
map_backend = :pgdg_localized_experimental

q = 5
reference_spacing = 1.0
tail_spacing = 10.0
q_to_core_spacing_rule = :standard_pqs_ns_equals_q
core_spacing = nothing
route_shape = (:pqs_left, :product, :pqs_right)
product_body_rule = :centered_single_z_slab
pqs_retained_rule = :boundary_comx_product_mode_selection
product_retained_rule = :product_doside_retained_unit
terms = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)
pair_factor_normalization = :density_normalized
support_dense_direct_allowed = false
reference_only_authorities = (
    :support_row_oracle,
    :dense_parent_projection,
)
save_artifact = false
save_tsv = false
outfile = "pqs_source_box_route_driver_report.jld2"
tsvfile = "pqs_source_box_route_driver_report.tsv"

inputs = String[]
if length(ARGS) > 0
    if occursin("=", ARGS[1])
        inputs = ARGS
    else
        fullpath = isabspath(ARGS[1]) ? ARGS[1] : joinpath(pwd(), ARGS[1])
        println("including ", fullpath)
        include(fullpath)
        length(ARGS) > 1 && (inputs = ARGS[2:end])
    end
end
eval.(Meta.parse.(inputs))

function _pqs_driver_print_section(title)
    println()
    println("[", title, "]")
    return nothing
end

function _pqs_driver_print_kv(key, value)
    println(rpad(String(key), 42), "  ", value)
    return nothing
end

function _pqs_driver_write_tsv_row(io, section, key, value)
    println(io, section, '\t', key, '\t', repr(value))
    return nothing
end

pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
    ArgumentError("pair_factor_normalization must be :density_normalized or :raw_weighted"),
)

standard_setup =
    GaussletBases.CartesianContractedParentMetrics._pqs_standard_source_box_route_setup(
        ;
        nuclear_charges,
        atom_locations,
        q,
        radius,
        reference_spacing,
        tail_spacing,
        q_to_core_spacing_rule,
        core_spacing,
    )

# 1. System metadata: physical labels and the parent box description.
system_metadata = (
    atom_symbols = atom_symbols,
    nuclear_charges = standard_setup.nuclear_charges,
    atom_locations = standard_setup.atom_locations,
    radius = standard_setup.radius,
    parent_axis_counts = parent_axis_counts,
    parent_box = standard_setup.parent_box,
    parent_box_rule = standard_setup.parent_box_rule,
    map_backend = map_backend,
)

# 2. Recipe metadata: route recipe inputs, not already-derived route facts.
recipe_metadata = (
    route_kind = route_kind,
    q = q,
    n_s = standard_setup.n_s,
    core_cube_side = standard_setup.core_cube_side,
    reference_spacing = standard_setup.reference_spacing,
    tail_spacing = standard_setup.tail_spacing,
    q_to_core_spacing_rule = standard_setup.q_to_core_spacing_rule,
    core_spacing = standard_setup.core_spacing,
    q_to_core_spacing_rule_status =
        standard_setup.spacing.q_to_core_spacing_rule_status,
    route_shape = route_shape,
    product_body_rule = product_body_rule,
    pqs_source_box_rule = :mode_selected_raw_box_pqs,
    pqs_retained_rule = pqs_retained_rule,
    product_retained_rule = product_retained_rule,
    terms = terms,
    pair_factor_normalization = pair_factor_normalization,
    support_dense_direct_allowed = support_dense_direct_allowed,
    reference_only_authorities = reference_only_authorities,
)

route_skeleton =
    GaussletBases.CartesianContractedParentMetrics._pqs_pqs_product_source_box_route_skeleton(
        ;
        q,
        parent_axis_counts,
        route_shape,
        product_body_rule,
        pqs_retained_rule,
        product_retained_rule,
        pair_factor_normalization,
    )

# 3. Parent description/construction: still described, not materialized here.
parent_description = (
    status = :described_not_constructed,
    standard_setup = standard_setup,
    physical_parent_box = standard_setup.parent_box,
    physical_parent_box_rule = standard_setup.parent_box_rule,
    axis_transform_status = :pending_repo_parent_constructor,
    one_dimensional_transforms = (:x_axis_transform, :y_axis_transform, :z_axis_transform),
    parent_lattice = :raw_product_box_parent_lattice,
    parent_axis_counts = route_skeleton.parent_axis_counts,
    source_boxes = route_skeleton.source_boxes,
    pending_facts = (
        route_skeleton.pending_facts...,
        :parent_axis_counts_from_standard_parent_constructor,
    ),
)

# 4. Product-type unit split: source boxes and source dimensions from the helper.
source_boxes = route_skeleton.source_boxes
source_dimensions = route_skeleton.source_dimensions

# 5. Retained-unit definitions: retained counts and ranges from the helper.
retained_units = route_skeleton.retained_units
retained_counts = route_skeleton.retained_counts
ranges = route_skeleton.ranges
retained_dimension = route_skeleton.retained_dimension

# 6. Pair inventory: upper-triangular source-box pair plan from the helper.
pair_entries = route_skeleton.pair_entries
pair_family_counts = route_skeleton.pair_family_counts
helper_by_pair_family = route_skeleton.helper_by_pair_family

# 7. Final linear algebra plan: the intended retained operator assembly.
linear_algebra_plan = (
    retained_block_formula = "O_final[i,j] = T_i' * O_source_box_pair * T_j",
    assemble_complete_retained_matrix = true,
    dense_parent_matrix_required = false,
    product_pqs_policy = :transpose_of_pqs_product_only_after_symmetric_pair_factor_check,
    product_pqs_lower_blocks = ((:product, :pqs_left), (:product, :pqs_right)),
    finite_output_check = :required_when_operator_blocks_are_materialized,
    symmetry_error_check = :required_for_symmetric_same_route_input,
)
no_go_flags = (
    public_default_behavior = false,
    packet_fixed_block_qw_hamiltonian_adoption = false,
    mwg_ida_semantic_change = false,
    retained_weight_division = false,
    retained_pqs_weights_used = false,
    repo_side_ray_id = false,
    ecp_behavior = false,
    cr2_science_claim = false,
    shell_projection = false,
    lowdin_cleanup = false,
    support_local_shell_row_algorithm = false,
    support_coefficient_matrix = false,
)

# 8. Validation/save: this dry-run validates route metadata, not operators.
stage_table = (
    (stage = 1, name = :collect_system_metadata, status = :represented),
    (stage = 2, name = :collect_recipe_metadata, status = :represented),
    (stage = 3, name = :construct_parent_object, status = :described_not_constructed),
    (stage = 4, name = :split_parent_into_product_type_units, status = :derived_by_helper),
    (stage = 5, name = :define_each_unit, status = :derived_by_helper),
    (stage = 6, name = :loop_over_unit_pairs, status = :derived_by_helper),
    (stage = 7, name = :apply_final_linear_algebra, status = :plan_reported),
    (stage = 8, name = :validate_report_save, status = :metadata_dry_run),
)
dry_run_validation = (
    builds_real_hamiltonian = false,
    builds_route_matrices = false,
    finite_output = :not_run_metadata_only,
    symmetry_error = :not_run_metadata_only,
    reference_error = :unavailable_metadata_only,
    timing_allocation = :placeholder_only,
)
diagnostics = merge(
    route_skeleton.diagnostics,
    (
        source = :pqs_source_box_route_driver_skeleton,
        standard_setup_helper =
            :_pqs_standard_source_box_route_setup,
        standard_setup_status = standard_setup.status,
        standard_setup_diagnostics = standard_setup.diagnostics,
        route_skeleton_helper =
            :_pqs_pqs_product_source_box_route_skeleton,
        n_s = standard_setup.n_s,
        core_cube_side = standard_setup.core_cube_side,
        core_cube_side_rule = standard_setup.core_cube_side_rule,
        parent_box_rule = standard_setup.parent_box_rule,
        parent_box = standard_setup.parent_box,
        core_spacing = standard_setup.core_spacing,
        mapping_s = standard_setup.mapping_s,
        q_to_core_spacing_rule = standard_setup.q_to_core_spacing_rule,
        q_to_core_spacing_rule_status =
            standard_setup.spacing.q_to_core_spacing_rule_status,
        q_to_core_spacing_provenance = standard_setup.spacing.provenance,
        q_to_core_spacing_non_optimality_claim =
            standard_setup.spacing.non_optimality_claim,
        parent_axis_counts_status =
            :manual_fixture_pending_standard_parent_constructor,
        output_representation = :retained_two_index_density_density,
        no_go_flags = no_go_flags,
        driver_builds_real_hamiltonian = false,
        driver_builds_route_matrices = false,
    ),
)
report = (
    object_kind = :pqs_source_box_route_driver_skeleton_report,
    generated_at = string(now()),
    standard_setup = standard_setup,
    system_metadata = system_metadata,
    recipe_metadata = recipe_metadata,
    parent_description = parent_description,
    route_skeleton = route_skeleton,
    route_shape = route_skeleton.route_shape,
    retained_unit_order = route_skeleton.retained_unit_order,
    source_boxes = source_boxes,
    source_dimensions = source_dimensions,
    retained_units = retained_units,
    retained_counts = retained_counts,
    ranges = ranges,
    retained_dimension = retained_dimension,
    pair_entries = pair_entries,
    pair_family_counts = pair_family_counts,
    helper_by_pair_family = helper_by_pair_family,
    linear_algebra_plan = linear_algebra_plan,
    stage_table = stage_table,
    dry_run_validation = dry_run_validation,
    diagnostics = diagnostics,
)

println("PQS source-box route driver skeleton")
@show route_kind
@show q
@show radius
@show route_shape
@show product_body_rule
@show pair_factor_normalization
@show standard_setup.n_s
@show standard_setup.core_cube_side
@show standard_setup.spacing.q_to_core_spacing_rule_status
@show retained_counts
@show retained_dimension

_pqs_driver_print_section("system_metadata")
for field in keys(system_metadata)
    _pqs_driver_print_kv(field, getproperty(system_metadata, field))
end

_pqs_driver_print_section("recipe_metadata")
for field in keys(recipe_metadata)
    _pqs_driver_print_kv(field, getproperty(recipe_metadata, field))
end

_pqs_driver_print_section("standard_setup")
for field in (
    :status,
    :n_s,
    :core_cube_side,
    :parent_box,
    :core_spacing,
    :mapping_s,
    :mapping_s_by_atom,
    :q_to_core_spacing_rule,
)
    _pqs_driver_print_kv(field, getproperty(standard_setup, field))
end
for field in (
    :q_to_core_spacing_rule_status,
    :q_to_core_spacing_provenance,
    :q_to_core_spacing_non_optimality_claim,
)
    _pqs_driver_print_kv(field, getproperty(standard_setup.diagnostics, field))
end

_pqs_driver_print_section("parent_description")
for field in keys(parent_description)
    _pqs_driver_print_kv(field, getproperty(parent_description, field))
end

_pqs_driver_print_section("source_boxes")
for field in keys(source_boxes)
    _pqs_driver_print_kv(field, getproperty(source_boxes, field))
end

_pqs_driver_print_section("retained_units")
for unit in retained_units
    println(
        unit.unit_key,
        '\t',
        unit.unit_role,
        '\t',
        unit.retained_count,
        '\t',
        unit.retained_range,
        '\t',
        unit.source_box,
        '\t',
        unit.retained_rule_derivation,
    )
end

_pqs_driver_print_section("pair_inventory")
@show length(pair_entries)
@show pair_family_counts
for entry in pair_entries
    println(
        entry.pair_key,
        '\t',
        entry.pair_family,
        '\t',
        entry.density_density_helper,
        '\t',
        entry.transpose_policy,
    )
end

_pqs_driver_print_section("linear_algebra_plan")
for field in keys(linear_algebra_plan)
    _pqs_driver_print_kv(field, getproperty(linear_algebra_plan, field))
end

_pqs_driver_print_section("diagnostics")
for field in keys(diagnostics)
    _pqs_driver_print_kv(field, getproperty(diagnostics, field))
end

if save_artifact
    println("saving JLD2 report ", outfile)
    jldsave(outfile; report)
end

if save_tsv
    println("saving TSV report ", tsvfile)
    open(tsvfile, "w") do io
        println(io, "section\tkey\tvalue")
        for field in keys(system_metadata)
            _pqs_driver_write_tsv_row(
                io,
                "system_metadata",
                field,
                getproperty(system_metadata, field),
            )
        end
        for field in keys(recipe_metadata)
            _pqs_driver_write_tsv_row(
                io,
                "recipe_metadata",
                field,
                getproperty(recipe_metadata, field),
            )
        end
        for field in keys(standard_setup.diagnostics)
            _pqs_driver_write_tsv_row(
                io,
                "standard_setup_diagnostics",
                field,
                getproperty(standard_setup.diagnostics, field),
            )
        end
        for unit in retained_units
            _pqs_driver_write_tsv_row(io, "retained_unit", unit.unit_key, unit)
        end
        for entry in pair_entries
            _pqs_driver_write_tsv_row(io, "pair_entry", entry.pair_key, entry)
        end
        for field in keys(diagnostics)
            _pqs_driver_write_tsv_row(
                io,
                "diagnostics",
                field,
                getproperty(diagnostics, field),
            )
        end
    end
end

println("driver skeleton complete")
