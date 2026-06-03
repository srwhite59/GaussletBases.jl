#!/usr/bin/env julia

using Dates, JLD2
using GaussletBases

route_kind = :be2_pqs_source_box_development_spine
atom_symbols = ("Be", "Be")
nuclear_charges = (4, 4)
atom_locations = ((-1.15, 0.0, 0.0), (1.15, 0.0, 0.0))
parent_axis_counts = (x = 9, y = 7, z = 9)
parent_box = (x = (-8.0, 8.0), y = (-6.0, 6.0), z = (-8.0, 8.0))
map_backend = :pgdg_localized_experimental

q = 5
L = 7
product_doside = 1
pqs_retained_count = 98
product_retained_count = 25
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

function _pqs_driver_product(values)
    result = 1
    for value in values
        result *= value
    end
    return result
end

function _pqs_driver_source_box_dimension(box)
    return _pqs_driver_product(length(getproperty(box, axis)) for axis in (:x, :y, :z))
end

function _pqs_driver_range(offset::Int, count::Int)
    count > 0 || throw(ArgumentError("retained counts must be positive"))
    return (offset + 1):(offset + count)
end

function _pqs_driver_pair_family(left_kind::Symbol, right_kind::Symbol)
    if left_kind == :product_doside && right_kind == :product_doside
        return :product_product
    elseif left_kind == :pqs && right_kind == :pqs
        return :pqs_pqs
    elseif left_kind == :pqs && right_kind == :product_doside
        return :pqs_product
    elseif left_kind == :product_doside && right_kind == :pqs
        return :product_pqs
    end
    throw(ArgumentError("unsupported unit pair kinds $(left_kind), $(right_kind)"))
end

function _pqs_driver_density_density_helper(pair_family::Symbol, mode::Symbol)
    if pair_family == :pqs_pqs
        return mode == :raw_weighted ?
            :_pqs_pqs_source_box_raw_weighted_density_density_interaction_block :
            :_pqs_pqs_source_box_density_density_interaction_block
    elseif pair_family in (:pqs_product, :product_pqs)
        return mode == :raw_weighted ?
            :_pqs_product_source_box_raw_weighted_density_density_interaction_block :
            :_pqs_product_source_box_density_density_interaction_block
    elseif pair_family == :product_product
        return mode == :raw_weighted ?
            :_product_doside_source_box_raw_weighted_density_density_interaction_block :
            :_product_doside_source_box_density_density_interaction_block
    end
    throw(ArgumentError("unsupported density-density pair family $(pair_family)"))
end

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
parent_axis_counts.x >= q || throw(ArgumentError("parent_axis_counts.x must be at least q"))
parent_axis_counts.y >= q || throw(ArgumentError("parent_axis_counts.y must be at least q"))
parent_axis_counts.z >= q || throw(ArgumentError("parent_axis_counts.z must be at least q"))
parent_axis_counts.z >= product_doside ||
    throw(ArgumentError("parent_axis_counts.z must be at least product_doside"))

source_mode_dims = (q, q, q)
product_z_start = div(parent_axis_counts.z - product_doside, 2) + 1
source_boxes = (
    pqs_left = (x = 1:q, y = 1:q, z = 1:q),
    pqs_right = (
        x = 1:q,
        y = 1:q,
        z = (parent_axis_counts.z - q + 1):parent_axis_counts.z,
    ),
    product = (
        x = 1:q,
        y = 1:q,
        z = product_z_start:(product_z_start + product_doside - 1),
    ),
)

pqs_left_range = _pqs_driver_range(0, pqs_retained_count)
pqs_right_range = _pqs_driver_range(last(pqs_left_range), pqs_retained_count)
product_range = _pqs_driver_range(last(pqs_right_range), product_retained_count)
retained_dimension = last(product_range)

system_metadata = (
    atom_symbols = atom_symbols,
    nuclear_charges = nuclear_charges,
    atom_locations = atom_locations,
    parent_axis_counts = parent_axis_counts,
    parent_box = parent_box,
    map_backend = map_backend,
)
recipe_metadata = (
    route_kind = route_kind,
    q = q,
    L = L,
    pqs_source_box_rule = :mode_selected_raw_box_pqs,
    product_body_rule = :middle_product_doside_slab,
    pqs_retained_rule = pqs_retained_rule,
    product_retained_rule = product_retained_rule,
    terms = terms,
    pair_factor_normalization = pair_factor_normalization,
    support_dense_direct_allowed = support_dense_direct_allowed,
    reference_only_authorities = reference_only_authorities,
)
parent_description = (
    status = :described_not_constructed,
    axis_transform_status = :pending_repo_parent_constructor,
    one_dimensional_transforms = (:x_axis_transform, :y_axis_transform, :z_axis_transform),
    parent_lattice = :raw_product_box_parent_lattice,
    source_boxes = source_boxes,
)
retained_units = (
    (
        unit_key = :pqs_left,
        unit_role = :left_mode_selected_raw_box_pqs_unit,
        retained_unit_kind = :pqs,
        source_family = :mode_selected_raw_product_box,
        source_box = source_boxes.pqs_left,
        source_dimensions = source_mode_dims,
        source_dimension = _pqs_driver_source_box_dimension(source_boxes.pqs_left),
        retained_rule_kind = pqs_retained_rule,
        retained_range = pqs_left_range,
        retained_count = length(pqs_left_range),
        provenance_label = :pqs_left_source_modes,
        weight_semantics = :retained_columns_not_positive_quadrature_weights,
    ),
    (
        unit_key = :pqs_right,
        unit_role = :right_mode_selected_raw_box_pqs_unit,
        retained_unit_kind = :pqs,
        source_family = :mode_selected_raw_product_box,
        source_box = source_boxes.pqs_right,
        source_dimensions = source_mode_dims,
        source_dimension = _pqs_driver_source_box_dimension(source_boxes.pqs_right),
        retained_rule_kind = pqs_retained_rule,
        retained_range = pqs_right_range,
        retained_count = length(pqs_right_range),
        provenance_label = :pqs_right_source_modes,
        weight_semantics = :retained_columns_not_positive_quadrature_weights,
    ),
    (
        unit_key = :product,
        unit_role = :middle_product_doside_slab_unit,
        retained_unit_kind = :product_doside,
        source_family = :product_doside,
        source_box = source_boxes.product,
        source_dimensions = (q, q, product_doside),
        source_dimension = _pqs_driver_source_box_dimension(source_boxes.product),
        retained_rule_kind = product_retained_rule,
        retained_range = product_range,
        retained_count = length(product_range),
        provenance_label = :product_doside_source_modes,
        weight_semantics = :product_source_weights_owned_by_source_box_helpers,
    ),
)
unit_by_key = Dict(unit.unit_key => unit for unit in retained_units)
route_shape = (:pqs_left, :pqs_right, :product)

pair_entries = Any[]
for i in eachindex(route_shape), j in i:length(route_shape)
    left_key = route_shape[i]
    right_key = route_shape[j]
    left_unit = unit_by_key[left_key]
    right_unit = unit_by_key[right_key]
    pair_family = _pqs_driver_pair_family(
        left_unit.retained_unit_kind,
        right_unit.retained_unit_kind,
    )
    helper = _pqs_driver_density_density_helper(pair_family, pair_factor_normalization)
    transpose_policy =
        left_key == right_key ? :none :
        pair_family == :product_pqs ?
        :product_pqs_uses_transpose_of_pqs_product_when_pair_factors_are_symmetric :
        :lower_block_uses_transpose_when_pair_factors_are_symmetric
    push!(
        pair_entries,
        (
            pair_key = (left_key, right_key),
            pair_family = pair_family,
            density_density_helper = helper,
            source_box_algorithmic_path = true,
            fallback_oracle_path = false,
            transpose_policy = transpose_policy,
            output_representation = :retained_two_index_density_density,
        ),
    )
end
pair_entries = Tuple(pair_entries)
pair_family_counts = (
    pqs_pqs = count(entry -> entry.pair_family == :pqs_pqs, pair_entries),
    pqs_product = count(entry -> entry.pair_family == :pqs_product, pair_entries),
    product_pqs = count(entry -> entry.pair_family == :product_pqs, pair_entries),
    product_product =
        count(entry -> entry.pair_family == :product_product, pair_entries),
)
helper_by_pair_family = (
    pqs_pqs = _pqs_driver_density_density_helper(:pqs_pqs, pair_factor_normalization),
    pqs_product =
        _pqs_driver_density_density_helper(:pqs_product, pair_factor_normalization),
    product_pqs = :transpose_of_pqs_product_helper_for_lower_blocks_only,
    product_product =
        _pqs_driver_density_density_helper(:product_product, pair_factor_normalization),
)
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
stage_table = (
    (stage = 1, name = :collect_system_metadata, status = :represented),
    (stage = 2, name = :collect_recipe_metadata, status = :represented),
    (stage = 3, name = :construct_parent_object, status = :described_not_constructed),
    (stage = 4, name = :split_parent_into_product_type_units, status = :represented),
    (stage = 5, name = :define_each_unit, status = :represented),
    (stage = 6, name = :loop_over_unit_pairs, status = :represented),
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
diagnostics = (
    source = :pqs_source_box_route_driver_skeleton,
    source_box_first = true,
    source_box_algorithmic_path_true_for_every_pair =
        all(entry -> entry.source_box_algorithmic_path, pair_entries),
    retained_dimension = retained_dimension,
    retained_unit_count = length(retained_units),
    pair_count = length(pair_entries),
    pair_family_counts = pair_family_counts,
    pair_factor_normalization = pair_factor_normalization,
    helper_by_pair_family = helper_by_pair_family,
    output_representation = :retained_two_index_density_density,
    four_index_galerkin_tensor = false,
    raw_weight_division_owner =
        pair_factor_normalization == :raw_weighted ?
        :explicit_source_quadrature_weight_outer_products :
        :caller_supplied_density_normalized_pair_factors,
    retained_weight_division_allowed = false,
    product_pqs_explicit_helper_required = false,
    product_pqs_transpose_requires_symmetric_pair_factors = true,
    product_pqs_lower_block_count = 2,
    product_pqs_lower_block_helper = :transpose_of_pqs_product_helper,
    no_go_flags = no_go_flags,
)
report = (
    object_kind = :pqs_source_box_route_driver_skeleton_report,
    generated_at = string(now()),
    system_metadata = system_metadata,
    recipe_metadata = recipe_metadata,
    parent_description = parent_description,
    route_shape = route_shape,
    retained_units = retained_units,
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
@show q L product_doside
@show pair_factor_normalization
@show retained_dimension

_pqs_driver_print_section("system_metadata")
for field in keys(system_metadata)
    _pqs_driver_print_kv(field, getproperty(system_metadata, field))
end

_pqs_driver_print_section("recipe_metadata")
for field in keys(recipe_metadata)
    _pqs_driver_print_kv(field, getproperty(recipe_metadata, field))
end

_pqs_driver_print_section("parent_description")
for field in keys(parent_description)
    _pqs_driver_print_kv(field, getproperty(parent_description, field))
end

_pqs_driver_print_section("retained_units")
for unit in retained_units
    println(
        unit.unit_key,
        '\t',
        unit.unit_role,
        '\t',
        unit.retained_range,
        '\t',
        unit.source_box,
        '\t',
        unit.weight_semantics,
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
