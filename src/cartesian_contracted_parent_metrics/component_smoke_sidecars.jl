function _pqs_pqs_product_source_box_component_route_smoke(
    bundles,
    left_source_box::NTuple{3,UnitRange{Int}},
    right_source_box::NTuple{3,UnitRange{Int}},
    product_source_box::NTuple{3,UnitRange{Int}},
    metrics::NamedTuple{(:x,:y,:z)},
    ida_provenance,
    nuclear_axis_layers::NamedTuple{(:x,:y,:z)},
    nuclear_expansion::CoulombGaussianExpansion;
    source_mode_dims::NTuple{3,Int},
    centers,
    nuclear_charges,
    center_labels = nothing,
    term_coefficients::AbstractVector{<:Real},
    pair_factor_normalization::Symbol = :density_normalized,
    pair_factor_symmetry_atol::Real = 1.0e-12,
    symmetry_atol::Real = 1.0e-10,
    nuclear_symmetry_atol::Real = 1.0e-10,
    route_name::Symbol = :pqs_pqs_product_source_box_component_route_smoke,
    parent_dims = _nested_axis_lengths(bundles),
    bond_axis = nothing,
    metadata = (;),
    provenance = (;),
    dense_parent_ida_matrix = nothing,
    dense_parent_matrix_source::Symbol = :caller_supplied_dense_parent_ida_matrix,
    dense_parent_comparison_atol::Real = 1.0e-10,
    route_supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
    orthogonality_atol::Real = 1.0e-8,
)
    pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError("component route smoke requires density_normalized or raw_weighted IDA route output"),
    )
    center_values = _source_box_nuclear_attraction_center_values(centers)
    charge_values = _source_box_nuclear_attraction_charge_values(nuclear_charges)
    length(center_values) == length(charge_values) || throw(
        ArgumentError("component route smoke center and charge counts must match"),
    )
    labels = _source_box_nuclear_attraction_center_labels(
        length(center_values),
        center_labels,
    )
    component_metadata = merge(
        (
            component_route_smoke = true,
            nuclear_center_labels = labels,
            nuclear_center_count = length(center_values),
            electron_electron_source = :ida_gausslet_source_box,
        ),
        metadata,
    )
    component_provenance = merge(
        (source = :pqs_pqs_product_source_box_component_route_smoke,),
        provenance,
    )
    ida_route =
        _pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(
            bundles,
            left_source_box,
            right_source_box,
            product_source_box,
            metrics,
            ida_provenance;
            source_mode_dims,
            term_coefficients,
            pair_factor_normalization,
            pair_factor_symmetry_atol,
            symmetry_atol,
            route_name,
            parent_dims,
            bond_axis,
            metadata = component_metadata,
            provenance = component_provenance,
            route_supported_terms,
            orthogonality_atol,
        )
    nuclear_route = _pqs_pqs_product_route_shaped_nuclear_attraction_by_center(
        ida_route.route_descriptor,
        nuclear_axis_layers,
        nuclear_expansion;
        centers = center_values,
        nuclear_charges = charge_values,
        center_labels = labels,
        symmetry_atol = nuclear_symmetry_atol,
    )
    dense_parent_authority_supported =
        pair_factor_normalization == :density_normalized
    dense_parent_authority_skipped_reason =
        isnothing(dense_parent_ida_matrix) || dense_parent_authority_supported ?
        nothing : :density_normalized_authority_only
    dense_parent_authority =
        isnothing(dense_parent_ida_matrix) || !dense_parent_authority_supported ? nothing :
        _pqs_pqs_product_dense_parent_ida_authority_comparison(
            ida_route,
            dense_parent_ida_matrix;
            parent_dims,
            dense_parent_matrix_source,
            comparison_atol = dense_parent_comparison_atol,
        )
    output_finite =
        nuclear_route.output_finite &&
        ida_route.output_finite &&
        (isnothing(dense_parent_authority) || dense_parent_authority.output_finite)
    authority_max_errors = (
        electron_electron_dense_parent =
            isnothing(dense_parent_authority) ? nothing :
            dense_parent_authority.max_error,
        nuclear_total_from_center_sum = nuclear_route.total_from_center_error,
    )
    return (
        object_kind = :pqs_pqs_product_source_box_component_route_smoke,
        path = :pqs_pqs_product_source_box_component_route_smoke,
        status = :private_component_route_smoke,
        descriptor = ida_route.route_descriptor,
        route_descriptor = ida_route.route_descriptor,
        raw_box_route_producer = ida_route.raw_box_route_producer,
        route_shape = (:pqs_left, :pqs_right, :product),
        retained_units = ida_route.route_descriptor.unit_summaries,
        ranges = ida_route.ranges,
        retained_dimension = ida_route.retained_dimension,
        nuclear_attraction_by_center = nuclear_route,
        per_center_nuclear_matrices = nuclear_route.blocks_by_center,
        nuclear_total_matrix = nuclear_route.total_block,
        electron_electron_density_density = ida_route,
        electron_electron_matrix = ida_route.block,
        dense_parent_ida_authority = dense_parent_authority,
        center_labels = labels,
        centers = Tuple(center_values),
        nuclear_charges = Tuple(charge_values),
        ida_term_count = ida_route.term_count,
        pair_factor_normalization = pair_factor_normalization,
        output_finite = output_finite,
        nuclear_symmetry_error = nuclear_route.symmetry_error,
        electron_electron_symmetry_error = ida_route.symmetry_error,
        authority_max_errors = authority_max_errors,
        metadata = component_metadata,
        provenance = component_provenance,
        diagnostics = (
            source = :pqs_pqs_product_source_box_component_route_smoke,
            private_component_route_smoke = true,
            private_shadow_only = true,
            route_shape = (:pqs_left, :pqs_right, :product),
            route_name = ida_route.route_descriptor.route_name,
            route_descriptor_object_kind = ida_route.route_descriptor.object_kind,
            route_roles = ida_route.route_descriptor.roles,
            retained_ranges = ida_route.ranges,
            retained_dimension = ida_route.retained_dimension,
            retained_unit_count = ida_route.route_descriptor.retained_unit_count,
            unit_roles = map(unit -> unit.unit_key, ida_route.route_descriptor.unit_summaries),
            unit_retained_ranges =
                map(unit -> unit.retained_range, ida_route.route_descriptor.unit_summaries),
            nuclear_pair_count = nuclear_route.pair_count,
            nuclear_pair_family_counts = nuclear_route.pair_family_counts,
            electron_electron_pair_count = ida_route.pair_count,
            electron_electron_pair_family_counts = ida_route.pair_family_counts,
            helper_used_for_nuclear_pair_families =
                nuclear_route.diagnostics.helper_used_for_pair_families,
            helper_used_for_electron_electron_pair_families =
                ida_route.diagnostics.helper_used_for_pair_families,
            center_labels = labels,
            centers = Tuple(center_values),
            nuclear_charges = Tuple(charge_values),
            center_records = nuclear_route.diagnostics.center_records,
            center_count = length(center_values),
            by_center_nuclear_terms_preserved = true,
            nuclear_total_is_explicit_sum_of_center_pieces = true,
            nuclear_total_from_center_error = nuclear_route.total_from_center_error,
            ida_term_count = ida_route.term_count,
            ida_provenance_term_count = ida_provenance.term_count,
            pair_factor_normalization = pair_factor_normalization,
            density_normalized_pair_factors =
                pair_factor_normalization == :density_normalized,
            raw_weighted_pair_factors =
                pair_factor_normalization == :raw_weighted,
            density_normalized_pair_factors_generated =
                pair_factor_normalization == :raw_weighted,
            source_weight_division_owner =
                ida_route.diagnostics.source_weight_division_owner,
            source_weight_division_applied_by_helper =
                ida_route.diagnostics.source_weight_division_applied_by_helper,
            source_weight_division_shape =
                ida_route.diagnostics.source_weight_division_shape,
            source_weights_are_raw_source_weights =
                ida_route.diagnostics.source_weights_are_raw_source_weights,
            input_pair_factor_data = :ida_gausslet_source_box_provenance,
            interaction_path = :ida_gausslet_source_box,
            mwg_supplement_residual_path = false,
            real_ida_gausslet_source_box_provenance_adapted = true,
            real_mwg_ida_pair_factor_provenance_adapted = false,
            source_box_first = true,
            source_box_algorithmic_path_true_for_every_pair =
                nuclear_route.diagnostics.source_box_algorithmic_path_true_for_every_pair &&
                ida_route.diagnostics.source_box_algorithmic_path_true_for_every_pair,
            nuclear_source_box_algorithmic_pair_count =
                nuclear_route.diagnostics.source_box_algorithmic_pair_count,
            electron_electron_source_box_algorithmic_pair_count =
                ida_route.diagnostics.source_box_algorithmic_pair_count,
            output_finite = output_finite,
            nuclear_symmetry_error = nuclear_route.symmetry_error,
            electron_electron_symmetry_error = ida_route.symmetry_error,
            dense_parent_ida_authority_available = !isnothing(dense_parent_authority),
            dense_parent_ida_authority_supported_for_pair_factor_normalization =
                dense_parent_authority_supported,
            dense_parent_ida_authority_skipped_reason =
                dense_parent_authority_skipped_reason,
            dense_parent_ida_authority_max_error =
                isnothing(dense_parent_authority) ? nothing :
                dense_parent_authority.max_error,
            dense_parent_ida_authority_within_tolerance =
                isnothing(dense_parent_authority) ? nothing :
                dense_parent_authority.within_tolerance,
            dense_parent_projection_validation_only =
                !isnothing(dense_parent_authority),
            dense_parent_projection_algorithmic = false,
            output_representation = :component_matrices,
            electron_electron_output_representation = :two_index_density_density,
            electron_electron_four_index_galerkin_tensor = false,
            hamiltonian_matrix_built = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
            mwg_interaction_implemented = false,
            mwg_supplement_residual_provenance_adapted = false,
            ida_mwg_semantics_changed = false,
            mwg_ida_semantics_changed = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_division_allowed = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_oracle_used = false,
            support_local_pqs_oracle_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
            shell_row_algorithm = false,
        ),
    )
end

function _pqs_pqs_product_component_route_smoke_summary(component)
    hasproperty(component, :object_kind) &&
        component.object_kind == :pqs_pqs_product_source_box_component_route_smoke ||
        throw(
            ArgumentError("component route smoke summary requires _pqs_pqs_product_source_box_component_route_smoke output"),
        )
    hasproperty(component, :nuclear_attraction_by_center) || throw(
        ArgumentError("component route smoke summary requires nuclear_attraction_by_center"),
    )
    hasproperty(component, :electron_electron_density_density) || throw(
        ArgumentError("component route smoke summary requires electron_electron_density_density"),
    )
    hasproperty(component, :diagnostics) || throw(
        ArgumentError("component route smoke summary requires diagnostics"),
    )
    nuclear = component.nuclear_attraction_by_center
    electron_electron = component.electron_electron_density_density
    diagnostics = component.diagnostics
    no_go_diagnostics = (
        source_box_first = diagnostics.source_box_first,
        source_box_algorithmic_path_true_for_every_pair =
            diagnostics.source_box_algorithmic_path_true_for_every_pair,
        shell_projection_used = diagnostics.shell_projection_used,
        lowdin_cleanup_used = diagnostics.lowdin_cleanup_used,
        support_local_oracle_used = diagnostics.support_local_oracle_used,
        support_local_pqs_oracle_used =
            diagnostics.support_local_pqs_oracle_used,
        support_local_shell_row_algorithm =
            diagnostics.support_local_shell_row_algorithm,
        support_coefficient_matrix_used =
            diagnostics.support_coefficient_matrix_used,
        shell_row_algorithm = diagnostics.shell_row_algorithm,
        retained_pqs_weights_used = diagnostics.retained_pqs_weights_used,
        retained_pqs_weights_positive_checked =
            diagnostics.retained_pqs_weights_positive_checked,
        retained_weight_division_allowed =
            diagnostics.retained_weight_division_allowed,
        retained_pqs_weight_division_allowed =
            diagnostics.retained_pqs_weight_division_allowed,
        ida_weight_division_allowed = diagnostics.ida_weight_division_allowed,
        packet_adoption = diagnostics.packet_adoption,
        fixed_block_routing = diagnostics.fixed_block_routing,
        qwhamiltonian_consumes = diagnostics.qwhamiltonian_consumes,
        hamiltonian_matrix_built = diagnostics.hamiltonian_matrix_built,
        public_default_consumes = diagnostics.public_default_consumes,
        mwg_supplement_residual_path =
            diagnostics.mwg_supplement_residual_path,
        mwg_supplement_residual_provenance_adapted =
            diagnostics.mwg_supplement_residual_provenance_adapted,
        ecp_terms_implemented = diagnostics.ecp_terms_implemented,
        cr2_science_status_changed = diagnostics.cr2_science_status_changed,
        dense_parent_projection_algorithmic =
            diagnostics.dense_parent_projection_algorithmic,
    )
    return (
        object_kind = :pqs_pqs_product_component_route_smoke_summary,
        status = :private_component_route_smoke_summary,
        component_object_kind = component.object_kind,
        component_status = component.status,
        route_shape = component.route_shape,
        retained_units = component.retained_units,
        ranges = component.ranges,
        retained_dimension = component.retained_dimension,
        center_labels = component.center_labels,
        nuclear_charges = component.nuclear_charges,
        center_count = length(component.center_labels),
        nuclear_pair_count = nuclear.pair_count,
        nuclear_pair_family_counts = nuclear.pair_family_counts,
        electron_electron_pair_count = electron_electron.pair_count,
        electron_electron_pair_family_counts =
            electron_electron.pair_family_counts,
        helper_used_for_nuclear_pair_families =
            diagnostics.helper_used_for_nuclear_pair_families,
        helper_used_for_electron_electron_pair_families =
            diagnostics.helper_used_for_electron_electron_pair_families,
        pair_factor_normalization = component.pair_factor_normalization,
        source_weight_division_owner =
            diagnostics.source_weight_division_owner,
        source_weight_division_applied_by_helper =
            diagnostics.source_weight_division_applied_by_helper,
        source_weight_division_shape =
            diagnostics.source_weight_division_shape,
        ida_term_count = component.ida_term_count,
        finite_checks = (
            output_finite = component.output_finite,
            nuclear_output_finite = nuclear.output_finite,
            electron_electron_output_finite = electron_electron.output_finite,
        ),
        symmetry_errors = (
            nuclear = component.nuclear_symmetry_error,
            electron_electron = component.electron_electron_symmetry_error,
        ),
        nuclear_total_from_center_error =
            nuclear.total_from_center_error,
        dense_parent_ida_authority = (
            available =
                diagnostics.dense_parent_ida_authority_available,
            max_error =
                diagnostics.dense_parent_ida_authority_max_error,
            within_tolerance =
                diagnostics.dense_parent_ida_authority_within_tolerance,
            skip_reason =
                diagnostics.dense_parent_ida_authority_skipped_reason,
            validation_only =
                diagnostics.dense_parent_projection_validation_only,
        ),
        no_go_diagnostics = no_go_diagnostics,
        performance = (
            nuclear = (
                elapsed_seconds = nuclear.performance.elapsed_seconds,
                allocated_bytes = nuclear.performance.allocated_bytes,
                gc_time_seconds = nuclear.performance.gc_time_seconds,
            ),
            electron_electron = (
                elapsed_seconds =
                    electron_electron.performance.elapsed_seconds,
                allocated_bytes =
                    electron_electron.performance.allocated_bytes,
                gc_time_seconds =
                    electron_electron.performance.gc_time_seconds,
            ),
        ),
        diagnostics = (
            source = :pqs_pqs_product_component_route_smoke_summary,
            private_component_route_smoke_summary = true,
            component_object_kind = component.object_kind,
            component_status = component.status,
            route_shape = component.route_shape,
            retained_dimension = component.retained_dimension,
            center_count = length(component.center_labels),
            pair_factor_normalization = component.pair_factor_normalization,
            output_finite = component.output_finite,
            source_box_first = no_go_diagnostics.source_box_first,
            source_box_algorithmic_path_true_for_every_pair =
                no_go_diagnostics.source_box_algorithmic_path_true_for_every_pair,
            retained_pqs_weights_used =
                no_go_diagnostics.retained_pqs_weights_used,
            retained_weight_division_allowed =
                no_go_diagnostics.retained_weight_division_allowed,
            ida_weight_division_allowed =
                no_go_diagnostics.ida_weight_division_allowed,
            packet_adoption = no_go_diagnostics.packet_adoption,
            fixed_block_routing = no_go_diagnostics.fixed_block_routing,
            qwhamiltonian_consumes = no_go_diagnostics.qwhamiltonian_consumes,
            public_default_consumes =
                no_go_diagnostics.public_default_consumes,
            mwg_supplement_residual_provenance_adapted =
                no_go_diagnostics.mwg_supplement_residual_provenance_adapted,
            ecp_terms_implemented = no_go_diagnostics.ecp_terms_implemented,
            cr2_science_status_changed =
                no_go_diagnostics.cr2_science_status_changed,
        ),
    )
end

function _pqs_pqs_product_component_route_smoke_no_go_clear(no_go)
    return no_go.source_box_first &&
        no_go.source_box_algorithmic_path_true_for_every_pair &&
        !no_go.shell_projection_used &&
        !no_go.lowdin_cleanup_used &&
        !no_go.support_local_oracle_used &&
        !no_go.support_local_pqs_oracle_used &&
        !no_go.support_local_shell_row_algorithm &&
        !no_go.support_coefficient_matrix_used &&
        !no_go.shell_row_algorithm &&
        !no_go.retained_pqs_weights_used &&
        !no_go.retained_pqs_weights_positive_checked &&
        !no_go.retained_weight_division_allowed &&
        !no_go.retained_pqs_weight_division_allowed &&
        !no_go.ida_weight_division_allowed &&
        !no_go.packet_adoption &&
        !no_go.fixed_block_routing &&
        !no_go.qwhamiltonian_consumes &&
        !no_go.hamiltonian_matrix_built &&
        !no_go.public_default_consumes &&
        !no_go.mwg_supplement_residual_path &&
        !no_go.mwg_supplement_residual_provenance_adapted &&
        !no_go.ecp_terms_implemented &&
        !no_go.cr2_science_status_changed &&
        !no_go.dense_parent_projection_algorithmic
end

function _pqs_pqs_product_component_route_smoke_report_row(summary)
    hasproperty(summary, :object_kind) &&
        summary.object_kind == :pqs_pqs_product_component_route_smoke_summary ||
        throw(
            ArgumentError("component route smoke report row requires _pqs_pqs_product_component_route_smoke_summary output"),
        )
    return (
        mode = summary.pair_factor_normalization,
        route_shape = summary.route_shape,
        retained_units = hasproperty(summary, :retained_units) ?
            summary.retained_units : (),
        ranges = hasproperty(summary, :ranges) ? summary.ranges : nothing,
        retained_dimension = summary.retained_dimension,
        nuclear_pair_count = summary.nuclear_pair_count,
        nuclear_pair_family_counts =
            summary.nuclear_pair_family_counts,
        electron_electron_pair_count = summary.electron_electron_pair_count,
        electron_electron_pair_family_counts =
            summary.electron_electron_pair_family_counts,
        helper_used_for_nuclear_pair_families =
            hasproperty(summary, :helper_used_for_nuclear_pair_families) ?
            summary.helper_used_for_nuclear_pair_families : nothing,
        helper_used_for_electron_electron_pair_families =
            hasproperty(summary, :helper_used_for_electron_electron_pair_families) ?
            summary.helper_used_for_electron_electron_pair_families : nothing,
        ida_term_count = summary.ida_term_count,
        output_finite = summary.finite_checks.output_finite,
        nuclear_symmetry_error = summary.symmetry_errors.nuclear,
        electron_electron_symmetry_error =
            summary.symmetry_errors.electron_electron,
        nuclear_total_from_center_error =
            summary.nuclear_total_from_center_error,
        dense_parent_ida_authority_available =
            summary.dense_parent_ida_authority.available,
        dense_parent_ida_authority_max_error =
            summary.dense_parent_ida_authority.max_error,
        dense_parent_ida_authority_skip_reason =
            summary.dense_parent_ida_authority.skip_reason,
        source_weight_division_owner =
            summary.source_weight_division_owner,
        source_weight_division_applied_by_helper =
            summary.source_weight_division_applied_by_helper,
        source_weight_division_shape =
            hasproperty(summary, :source_weight_division_shape) ?
            summary.source_weight_division_shape : nothing,
        no_go_clear =
            _pqs_pqs_product_component_route_smoke_no_go_clear(
                summary.no_go_diagnostics,
            ),
        no_go_diagnostics = summary.no_go_diagnostics,
        performance = hasproperty(summary, :performance) ?
            summary.performance : nothing,
        mwg_supplement_residual_path =
            summary.no_go_diagnostics.mwg_supplement_residual_path,
    )
end

function _pqs_component_route_smoke_fact(facts, key::AbstractString, default = "")
    if facts isa AbstractDict
        return get(facts, String(key), default)
    end
    symbol_key = Symbol(key)
    hasproperty(facts, symbol_key) && return getproperty(facts, symbol_key)
    return default
end

function _pqs_component_route_smoke_fact_string(value)
    return value === nothing ? "" : string(value)
end

function _pqs_component_route_smoke_fact_false(facts, key::AbstractString)
    value = _pqs_component_route_smoke_fact(facts, key, nothing)
    value === false && return true
    return lowercase(_pqs_component_route_smoke_fact_string(value)) == "false"
end

function _pqs_component_route_smoke_fact_zero(value)
    value isa Real && return iszero(value)
    parsed = tryparse(Float64, _pqs_component_route_smoke_fact_string(value))
    return !isnothing(parsed) && iszero(parsed)
end

function _pqs_component_route_smoke_parse_int_tuple(value, field::Symbol)
    value isa AbstractVector && return Tuple(Int(index) for index in value)
    text = strip(_pqs_component_route_smoke_fact_string(value))
    isempty(text) && throw(
        ArgumentError("component route smoke sidecar requires $(field)"),
    )
    values = Int[
        parse(Int, match.match) for
        match in eachmatch(r"-?\d+", text)
    ]
    isempty(values) && throw(
        ArgumentError("component route smoke sidecar could not parse integer values for $(field)"),
    )
    return Tuple(values)
end

function _pqs_component_route_smoke_parse_shape(value, field::Symbol)
    values = _pqs_component_route_smoke_parse_int_tuple(value, field)
    length(values) == 2 || throw(
        ArgumentError("component route smoke sidecar requires a two-dimensional shape for $(field)"),
    )
    all(dimension -> dimension >= 0, values) || throw(
        ArgumentError("component route smoke sidecar requires nonnegative shape entries for $(field)"),
    )
    return values
end

function _pqs_component_route_smoke_unit_key(unit)
    hasproperty(unit, :unit_key) || throw(
        ArgumentError("component route smoke sidecar source-box unit requires unit_key"),
    )
    return unit.unit_key
end

function _pqs_component_route_smoke_unit_records(units)
    return Tuple(
        (
            unit_key = _pqs_component_route_smoke_unit_key(unit),
            retained_unit_kind = hasproperty(unit, :retained_unit_kind) ?
                unit.retained_unit_kind : :unknown,
            source_family = hasproperty(unit, :source_family) ?
                unit.source_family : :unknown,
            retained_rule_kind = hasproperty(unit, :retained_rule_kind) ?
                unit.retained_rule_kind : :unknown,
            retained_range = hasproperty(unit, :retained_range) ?
                unit.retained_range : nothing,
            retained_count = hasproperty(unit, :retained_count) ?
                Int(unit.retained_count) : (
                    hasproperty(unit, :retained_range) ?
                    length(unit.retained_range) : 0
                ),
            source_dimensions = hasproperty(unit, :source_dimensions) ?
                unit.source_dimensions : nothing,
            source_dimension = hasproperty(unit, :source_dimension) ?
                unit.source_dimension : nothing,
        ) for unit in units
    )
end

function _pqs_component_route_smoke_rows_units(rows)
    unit_rows = Tuple(row for row in rows if hasproperty(row, :retained_units) && !isempty(row.retained_units))
    isempty(unit_rows) && return ()
    records = _pqs_component_route_smoke_unit_records(first(unit_rows).retained_units)
    for row in unit_rows
        _pqs_component_route_smoke_unit_records(row.retained_units) == records || throw(
            ArgumentError("component route smoke report rows disagree on source-box retained-unit records"),
        )
    end
    return records
end

function _pqs_component_route_smoke_rows_ranges(rows)
    range_rows = Tuple(row for row in rows if hasproperty(row, :ranges) && !isnothing(row.ranges))
    isempty(range_rows) && return nothing
    ranges = first(range_rows).ranges
    all(row -> row.ranges == ranges, range_rows) || throw(
        ArgumentError("component route smoke report rows disagree on source-box retained ranges"),
    )
    return ranges
end

function _pqs_pqs_product_component_route_smoke_report_adapter(
    summaries,
    final_residual_mwg_facts;
    generated_at = nothing,
    source_report = nothing,
    route_shape =
        "left raw-box PQS | product/doside slab | right raw-box PQS",
    parent_dims = nothing,
    source_mode_dims = nothing,
    left_source_box = nothing,
    right_source_box = nothing,
    product_source_box = nothing,
    nuclear_centers = nothing,
    nuclear_charges = nothing,
    strict::Bool = true,
)
    rows = Tuple(
        _pqs_pqs_product_component_route_smoke_report_row(summary)
        for summary in summaries
    )
    isempty(rows) && throw(
        ArgumentError("component route smoke report adapter requires at least one source-box summary"),
    )
    source_rows_clear = all(row -> row.output_finite && row.no_go_clear, rows)
    density_rows = filter(row -> row.mode == :density_normalized, rows)
    density_authority =
        isempty(density_rows) ? nothing :
        first(density_rows).dense_parent_ida_authority_max_error
    source_retained_units = _pqs_component_route_smoke_rows_units(rows)
    source_retained_ranges = _pqs_component_route_smoke_rows_ranges(rows)

    residual_owner_vector = _pqs_component_route_smoke_fact(
        final_residual_mwg_facts,
        "residual_nucleus_indices",
        "",
    )
    owner_inferred = _pqs_component_route_smoke_fact(
        final_residual_mwg_facts,
        "owner_semantics_inferred_from_raw_to_final_support",
        "",
    )
    mwg_authority_error = _pqs_component_route_smoke_fact(
        final_residual_mwg_facts,
        "max_authority_error",
        "",
    )
    no_owner_inference =
        _pqs_component_route_smoke_fact_false(
            final_residual_mwg_facts,
            "owner_semantics_inferred_from_raw_to_final_support",
        )
    no_raw_gto_gto_mwg_blocks =
        _pqs_component_route_smoke_fact_false(
            final_residual_mwg_facts,
            "diagnostics.raw_gto_gto_mwg_interaction_blocks_used",
        )
    no_fixed_raw_gto_mwg_blocks =
        _pqs_component_route_smoke_fact_false(
            final_residual_mwg_facts,
            "diagnostics.fixed_raw_gto_mwg_interaction_blocks_used",
        )
    mwg_authority_zero =
        _pqs_component_route_smoke_fact_zero(mwg_authority_error)
    owner_vector_available =
        !isempty(_pqs_component_route_smoke_fact_string(residual_owner_vector))

    if strict
        source_rows_clear || throw(
            ArgumentError("component route smoke report adapter requires finite source-box rows with clear no-go diagnostics"),
        )
        owner_vector_available || throw(
            ArgumentError("component route smoke report adapter requires explicit residual owner metadata"),
        )
        mwg_authority_zero || throw(
            ArgumentError("component route smoke report adapter requires zero final-residual MWG authority error"),
        )
        no_owner_inference || throw(
            ArgumentError("component route smoke report adapter must not infer owners from raw_to_final support"),
        )
        no_raw_gto_gto_mwg_blocks || throw(
            ArgumentError("component route smoke report adapter must not use raw GTO/GTO MWG blocks"),
        )
        no_fixed_raw_gto_mwg_blocks || throw(
            ArgumentError("component route smoke report adapter must not use fixed/raw-GTO MWG blocks"),
        )
    end

    source_box_smoke = (
        route_shape = route_shape,
        parent_dims = parent_dims,
        source_mode_dims = source_mode_dims,
        left_source_box = left_source_box,
        right_source_box = right_source_box,
        product_source_box = product_source_box,
        nuclear_centers = nuclear_centers,
        nuclear_charges = nuclear_charges,
        pqs_source_box_fixed_side_facts_available = true,
        by_center_nuclear_attraction_available = true,
        ida_source_box_electron_electron_available = true,
        electron_electron_representation =
            "retained two-index density-density",
        mwg_supplement_residual_adapted_in_source_box_smoke = false,
        retained_units = source_retained_units,
        retained_ranges = source_retained_ranges,
        retained_unit_count = length(source_retained_units),
        source_unit_label_status =
            isempty(source_retained_units) ? :unavailable :
            :explicit_route_descriptor_unit_keys,
        source_unit_labels =
            Tuple(unit.unit_key for unit in source_retained_units),
        rows = rows,
    )
    final_residual = (
        source_report = source_report,
        route = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "route",
            "",
        ),
        residual_owner_metadata_available = owner_vector_available,
        residual_nucleus_indices = residual_owner_vector,
        residual_owner_counts = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "residual_owner_counts",
            "",
        ),
        owner_metadata_source = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "owner_metadata_source",
            "",
        ),
        owner_semantics_inferred_from_raw_to_final_support = owner_inferred,
        component_helper = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "component_helper",
            "",
        ),
        max_authority_error = mwg_authority_error,
        fixed_fixed_shape = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "fixed_fixed_shape",
            "",
        ),
        fixed_residual_shape = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "fixed_residual_shape",
            "",
        ),
        residual_residual_shape = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "residual_residual_shape",
            "",
        ),
        final_interaction_shape = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "final_interaction_shape",
            "",
        ),
        raw_gto_rows_role = _pqs_component_route_smoke_fact(
            final_residual_mwg_facts,
            "diagnostics.raw_gto_rows_role",
            "",
        ),
        raw_gto_gto_mwg_interaction_blocks_used =
            _pqs_component_route_smoke_fact(
                final_residual_mwg_facts,
                "diagnostics.raw_gto_gto_mwg_interaction_blocks_used",
                "",
            ),
        fixed_raw_gto_mwg_interaction_blocks_used =
            _pqs_component_route_smoke_fact(
                final_residual_mwg_facts,
                "diagnostics.fixed_raw_gto_mwg_interaction_blocks_used",
                "",
            ),
    )
    lane_boundaries = (
        can_share_private_smoke_report = true,
        pqs_source_box_ida_and_mwg_residual_same_algorithm = false,
        source_box_ida_lane =
            "fixed/PQS source-box retained two-index density-density",
        mwg_residual_lane =
            "ordinary final-residual MWG supplement coupling with explicit owners",
        must_remain_separate_until_reviewed = true,
        no_owner_inference_from_raw_to_final_support = no_owner_inference,
        no_raw_gto_gto_mwg_blocks = no_raw_gto_gto_mwg_blocks,
        no_fixed_raw_gto_mwg_blocks = no_fixed_raw_gto_mwg_blocks,
        no_retained_weight_or_ida_division = true,
        no_packet_fixed_block_qw_hamiltonian_public_adoption = true,
        no_ecp_scf_hf_cr2_science_claim = true,
        no_mwg_ida_semantic_change = true,
    )

    return (
        object_kind = :pqs_pqs_product_component_route_smoke_report_adapter,
        status = :private_component_route_smoke_report,
        report_status = "private_summary_only",
        title = "Be2/PQS component-route smoke summary",
        generated_at = generated_at,
        smallest_honest_smoke = "adjacent_component_facts_in_one_report",
        single_algorithmic_operator = false,
        hamiltonian_built = false,
        production_route_adoption = false,
        public_default_route = false,
        source_box_pqs_ida_component_smoke = source_box_smoke,
        final_residual_mwg_supplement_component_facts = final_residual,
        lane_boundaries = lane_boundaries,
        diagnostics = (
            source = :pqs_pqs_product_component_route_smoke_report_adapter,
            private_component_route_smoke_report_adapter = true,
            source_box_rows_all_finite_and_no_go_clear = source_rows_clear,
            density_normalized_ida_dense_parent_authority_max_error =
                density_authority,
            final_residual_mwg_authority_error_zero = mwg_authority_zero,
            final_residual_mwg_owner_vector_available =
                owner_vector_available,
            source_box_pqs_ida_and_mwg_residual_same_algorithm = false,
            no_owner_inference_from_raw_to_final_support = no_owner_inference,
            no_raw_gto_gto_mwg_blocks = no_raw_gto_gto_mwg_blocks,
            no_fixed_raw_gto_mwg_blocks = no_fixed_raw_gto_mwg_blocks,
            no_retained_weight_or_ida_division = true,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            public_default_consumes = false,
            ecp_terms_implemented = false,
            scf_hf_validation_claim = false,
            cr2_science_status_changed = false,
            science_route_adoption = false,
            mwg_ida_semantics_changed = false,
        ),
    )
end

function _pqs_component_route_smoke_source_box_sidecar_components(row)
    retained_shape = (row.retained_dimension, row.retained_dimension)
    return (
        mode = row.mode,
        pair_factor_normalization = row.mode,
        retained_matrix_shape = retained_shape,
        nuclear_attraction_by_center = (
            available = true,
            pair_count = row.nuclear_pair_count,
            pair_family_counts = row.nuclear_pair_family_counts,
            helper_by_family = row.helper_used_for_nuclear_pair_families,
            authority_error = row.nuclear_total_from_center_error,
            symmetry_error = row.nuclear_symmetry_error,
        ),
        electron_electron_density_density = (
            available = true,
            representation = :retained_two_index_density_density,
            not_four_index_galerkin_tensor = true,
            pair_count = row.electron_electron_pair_count,
            pair_family_counts = row.electron_electron_pair_family_counts,
            helper_by_family =
                row.helper_used_for_electron_electron_pair_families,
            dense_parent_authority_available =
                row.dense_parent_ida_authority_available,
            dense_parent_authority_max_error =
                row.dense_parent_ida_authority_max_error,
            dense_parent_authority_skip_reason =
                row.dense_parent_ida_authority_skip_reason,
            symmetry_error = row.electron_electron_symmetry_error,
        ),
        source_weight_division_owner = row.source_weight_division_owner,
        source_weight_division_applied_by_helper =
            row.source_weight_division_applied_by_helper,
        source_weight_division_shape = row.source_weight_division_shape,
        output_finite = row.output_finite,
        no_go_clear = row.no_go_clear,
    )
end

function _pqs_component_route_smoke_false_keys(flags)
    return Tuple(key for key in keys(flags) if !getproperty(flags, key))
end

function _pqs_component_route_smoke_true_keys(flags)
    return Tuple(key for key in keys(flags) if getproperty(flags, key))
end

function _pqs_component_route_smoke_row_performance(row)
    performance = hasproperty(row, :performance) ? row.performance : nothing
    available = !isnothing(performance)
    return (
        mode = row.mode,
        available = available,
        nuclear_elapsed_seconds =
            available ? performance.nuclear.elapsed_seconds : nothing,
        nuclear_allocated_bytes =
            available ? performance.nuclear.allocated_bytes : nothing,
        nuclear_gc_time_seconds =
            available ? performance.nuclear.gc_time_seconds : nothing,
        electron_electron_elapsed_seconds =
            available ? performance.electron_electron.elapsed_seconds : nothing,
        electron_electron_allocated_bytes =
            available ? performance.electron_electron.allocated_bytes : nothing,
        electron_electron_gc_time_seconds =
            available ? performance.electron_electron.gc_time_seconds : nothing,
    )
end

function _pqs_component_route_smoke_sidecar_inference_flags(cr2_sidecar)
    isnothing(cr2_sidecar) && return (
        available = false,
        source_metadata_sidecar_available = false,
        source_shell_mode_inventory_available = false,
        label_reconstruction_from_centers = false,
        nearest_grid_or_center_label_heuristic = false,
        source_shell_mode_inference = false,
        retained_weight_or_ida_division = false,
        route_construction_changed = false,
    )

    cr2_sidecar.object_kind ==
        :pqs_pqs_product_component_route_smoke_cr2_sidecar_schema || throw(
        ArgumentError("private route-adapter readiness requires _pqs_pqs_product_component_route_smoke_cr2_sidecar_schema output when a CR2 sidecar is supplied"),
    )

    fixed_side = cr2_sidecar.fixed_side_retained_unit_metadata
    fixed_side_label_reconstruction =
        !isnothing(fixed_side) &&
        fixed_side.labels.label_reconstruction_from_centers
    fixed_side_nearest_grid =
        !isnothing(fixed_side) &&
        fixed_side.labels.nearest_grid_or_center_label_heuristic
    fixed_side_retained_weight_division =
        !isnothing(fixed_side) &&
        fixed_side.diagnostics.retained_weight_or_ida_division

    source_shell_modes = cr2_sidecar.source_shell_mode_inventory
    source_shell_mode_inventory_available = !isnothing(source_shell_modes)
    source_shell_mode_inference =
        !isnothing(source_shell_modes) &&
        (
            source_shell_modes.diagnostics.inferred_from_centers ||
            source_shell_modes.diagnostics.inferred_from_nearest_grid ||
            source_shell_modes.diagnostics.inferred_from_support_order ||
            source_shell_modes.diagnostics.inferred_from_support_indices ||
            source_shell_modes.diagnostics.inferred_from_raw_to_final_support
        )
    source_shell_mode_retained_weight_division =
        !isnothing(source_shell_modes) &&
        source_shell_modes.diagnostics.retained_weight_or_ida_division
    source_shell_mode_route_changed =
        !isnothing(source_shell_modes) &&
        source_shell_modes.diagnostics.route_construction_changed

    return (
        available = true,
        source_metadata_sidecar_available =
            source_shell_mode_inventory_available,
        source_shell_mode_inventory_available =
            source_shell_mode_inventory_available,
        label_reconstruction_from_centers =
            cr2_sidecar.labels.label_reconstruction_from_centers ||
            fixed_side_label_reconstruction,
        nearest_grid_or_center_label_heuristic =
            cr2_sidecar.labels.nearest_grid_or_center_label_heuristic ||
            fixed_side_nearest_grid,
        source_shell_mode_inference = source_shell_mode_inference,
        retained_weight_or_ida_division =
            fixed_side_retained_weight_division ||
            source_shell_mode_retained_weight_division,
        route_construction_changed = source_shell_mode_route_changed,
    )
end

function _pqs_pqs_product_private_source_box_route_adapter_readiness_summary(
    report;
    cr2_sidecar = nothing,
    require_source_metadata_sidecar::Bool = false,
)
    report.object_kind ==
        :pqs_pqs_product_component_route_smoke_report_adapter || throw(
        ArgumentError("private route-adapter readiness requires _pqs_pqs_product_component_route_smoke_report_adapter output"),
    )
    source_box = report.source_box_pqs_ida_component_smoke
    rows = source_box.rows
    isempty(rows) && throw(
        ArgumentError("private route-adapter readiness requires at least one source-box row"),
    )
    final_residual = report.final_residual_mwg_supplement_component_facts
    sidecar_flags =
        _pqs_component_route_smoke_sidecar_inference_flags(cr2_sidecar)
    row_performance = Tuple(_pqs_component_route_smoke_row_performance(row) for row in rows)
    timing_allocation_fields_available =
        all(row -> row.available, row_performance)
    authority_available_modes = Tuple(
        row.mode for row in rows
        if row.dense_parent_ida_authority_available
    )
    authority_skip_reasons = Tuple(
        row.dense_parent_ida_authority_skip_reason for row in rows
        if !isnothing(row.dense_parent_ida_authority_skip_reason)
    )
    authority_comparison_accounted_for =
        all(
            row ->
                row.dense_parent_ida_authority_available ||
                !isnothing(row.dense_parent_ida_authority_skip_reason),
            rows,
        )
    source_box_algorithmic_rows_clear = all(
        row ->
            row.no_go_clear &&
            row.no_go_diagnostics.source_box_first &&
            row.no_go_diagnostics.source_box_algorithmic_path_true_for_every_pair,
        rows,
    )
    pair_factor_normalization_modes = Tuple(row.mode for row in rows)

    required_input_availability = (
        source_box_report = true,
        source_box_fixed_side_facts =
            source_box.pqs_source_box_fixed_side_facts_available,
        by_center_nuclear_attraction =
            source_box.by_center_nuclear_attraction_available &&
            all(row -> row.nuclear_pair_count > 0, rows),
        ida_source_box_electron_electron =
            source_box.ida_source_box_electron_electron_available &&
            all(row -> row.electron_electron_pair_count > 0, rows),
        source_box_algorithmic_path =
            source_box_algorithmic_rows_clear,
        final_residual_mwg_component_facts =
            hasproperty(report, :final_residual_mwg_supplement_component_facts),
        final_residual_mwg_owner_metadata =
            final_residual.residual_owner_metadata_available,
        final_residual_mwg_authority =
            report.diagnostics.final_residual_mwg_authority_error_zero,
        authority_comparison_accounted_for =
            authority_comparison_accounted_for,
        cr2_source_metadata_sidecar =
            !require_source_metadata_sidecar ||
            sidecar_flags.source_metadata_sidecar_available,
    )
    optional_input_availability = (
        cr2_sidecar_schema = !isnothing(cr2_sidecar),
        timing_allocation_fields = timing_allocation_fields_available,
        dense_parent_ida_authority_comparison =
            !isempty(authority_available_modes),
    )
    no_go_flags = (
        public_default_behavior = report.diagnostics.public_default_consumes,
        packet_fixed_block_qw_hamiltonian_adoption =
            report.diagnostics.packet_adoption ||
            report.diagnostics.fixed_block_routing ||
            report.diagnostics.qwhamiltonian_consumes ||
            report.diagnostics.hamiltonian_matrix_built,
        mwg_ida_semantic_change =
            report.diagnostics.mwg_ida_semantics_changed,
        retained_weight_division =
            !report.lane_boundaries.no_retained_weight_or_ida_division ||
            sidecar_flags.retained_weight_or_ida_division,
        raw_gto_gto_mwg_blocks =
            !report.lane_boundaries.no_raw_gto_gto_mwg_blocks,
        fixed_raw_gto_mwg_blocks =
            !report.lane_boundaries.no_fixed_raw_gto_mwg_blocks,
        owner_shell_ray_inference =
            !report.lane_boundaries.no_owner_inference_from_raw_to_final_support ||
            sidecar_flags.label_reconstruction_from_centers ||
            sidecar_flags.nearest_grid_or_center_label_heuristic ||
            sidecar_flags.source_shell_mode_inference,
        route_construction_changed = sidecar_flags.route_construction_changed,
        cr2_science_status_changed =
            report.diagnostics.cr2_science_status_changed ||
            report.diagnostics.scf_hf_validation_claim,
    )
    no_go_violations =
        _pqs_component_route_smoke_true_keys(no_go_flags)
    missing_required_pieces =
        _pqs_component_route_smoke_false_keys(required_input_availability)
    missing_optional_pieces =
        _pqs_component_route_smoke_false_keys(optional_input_availability)
    ready_for_next_private_adapter_pass =
        isempty(missing_required_pieces) &&
        isempty(no_go_violations) &&
        report.diagnostics.source_box_rows_all_finite_and_no_go_clear

    return (
        object_kind =
            :pqs_pqs_product_private_source_box_route_adapter_readiness_summary,
        status = ready_for_next_private_adapter_pass ?
            :private_route_adapter_inputs_ready :
            :private_route_adapter_inputs_incomplete,
        report_object_kind = report.object_kind,
        report_status = report.status,
        route_shape = source_box.route_shape,
        retained_dimension = first(rows).retained_dimension,
        retained_unit_count = source_box.retained_unit_count,
        retained_units = source_box.retained_units,
        retained_ranges = source_box.retained_ranges,
        pair_factor_normalization_modes = pair_factor_normalization_modes,
        required_input_availability = required_input_availability,
        optional_input_availability = optional_input_availability,
        missing_required_pieces = missing_required_pieces,
        missing_optional_pieces = missing_optional_pieces,
        by_center_nuclear_attraction = (
            available =
                required_input_availability.by_center_nuclear_attraction,
            pair_counts = Tuple(row.nuclear_pair_count for row in rows),
            pair_family_counts =
                Tuple(row.nuclear_pair_family_counts for row in rows),
            helper_by_family = Tuple(
                row.helper_used_for_nuclear_pair_families for row in rows
            ),
            total_from_center_errors =
                Tuple(row.nuclear_total_from_center_error for row in rows),
        ),
        ida_source_box_electron_electron = (
            available =
                required_input_availability.ida_source_box_electron_electron,
            normalization_modes = pair_factor_normalization_modes,
            representation = source_box.electron_electron_representation,
            pair_counts =
                Tuple(row.electron_electron_pair_count for row in rows),
            pair_family_counts =
                Tuple(row.electron_electron_pair_family_counts for row in rows),
            helper_by_family = Tuple(
                row.helper_used_for_electron_electron_pair_families
                for row in rows
            ),
            source_weight_division_owner =
                Tuple(row.source_weight_division_owner for row in rows),
            source_weight_division_applied_by_helper = Tuple(
                row.source_weight_division_applied_by_helper for row in rows
            ),
            authority_comparison = (
                available_modes = authority_available_modes,
                skip_reasons = authority_skip_reasons,
                accounted_for = authority_comparison_accounted_for,
            ),
        ),
        final_residual_mwg_component_facts = (
            available =
                required_input_availability.final_residual_mwg_component_facts,
            residual_owner_metadata_available =
                final_residual.residual_owner_metadata_available,
            component_helper = final_residual.component_helper,
            max_authority_error = final_residual.max_authority_error,
            raw_gto_rows_role = final_residual.raw_gto_rows_role,
            raw_gto_gto_mwg_interaction_blocks_used =
                final_residual.raw_gto_gto_mwg_interaction_blocks_used,
            fixed_raw_gto_mwg_interaction_blocks_used =
                final_residual.fixed_raw_gto_mwg_interaction_blocks_used,
        ),
        cr2_sidecar = (
            available = !isnothing(cr2_sidecar),
            source_metadata_sidecar_required =
                require_source_metadata_sidecar,
            source_metadata_sidecar_available =
                sidecar_flags.source_metadata_sidecar_available,
            source_shell_mode_inventory_available =
                sidecar_flags.source_shell_mode_inventory_available,
            inference_flags = sidecar_flags,
        ),
        timing_allocation = (
            available = timing_allocation_fields_available,
            rows = row_performance,
        ),
        no_go_flags = no_go_flags,
        no_go_violations = no_go_violations,
        ready_for_next_private_adapter_pass =
            ready_for_next_private_adapter_pass,
        diagnostics = (
            source =
                :pqs_pqs_product_private_source_box_route_adapter_readiness_summary,
            private_shadow_only = true,
            builds_route_matrices = false,
            construction_behavior_changed = false,
            source_box_ida_and_mwg_residual_same_algorithm = false,
            lanes_remain_separate = true,
            source_box_rows_all_finite_and_no_go_clear =
                report.diagnostics.source_box_rows_all_finite_and_no_go_clear,
            source_box_algorithmic_path =
                required_input_availability.source_box_algorithmic_path,
            cr2_source_metadata_sidecar_required =
                require_source_metadata_sidecar,
            cr2_source_metadata_sidecar_available =
                sidecar_flags.source_metadata_sidecar_available,
            authority_comparison_accounted_for =
                authority_comparison_accounted_for,
            timing_allocation_fields_available =
                timing_allocation_fields_available,
            public_default_consumes =
                report.diagnostics.public_default_consumes,
            packet_adoption = report.diagnostics.packet_adoption,
            fixed_block_routing = report.diagnostics.fixed_block_routing,
            qwhamiltonian_consumes = report.diagnostics.qwhamiltonian_consumes,
            hamiltonian_matrix_built =
                report.diagnostics.hamiltonian_matrix_built,
            mwg_ida_semantics_changed =
                report.diagnostics.mwg_ida_semantics_changed,
            retained_weight_or_ida_division =
                no_go_flags.retained_weight_division,
            raw_gto_gto_mwg_blocks =
                no_go_flags.raw_gto_gto_mwg_blocks,
            fixed_raw_gto_mwg_blocks =
                no_go_flags.fixed_raw_gto_mwg_blocks,
            owner_shell_ray_inference =
                no_go_flags.owner_shell_ray_inference,
            cr2_science_status_changed =
                no_go_flags.cr2_science_status_changed,
        ),
    )
end

function _pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
    report;
    fixed_side_retained_unit_metadata = nothing,
    source_shell_mode_inventory = nothing,
    provenance = (;),
    strict::Bool = true,
)
    report.object_kind == :pqs_pqs_product_component_route_smoke_report_adapter ||
        throw(
            ArgumentError("CR2 sidecar schema requires _pqs_pqs_product_component_route_smoke_report_adapter output"),
        )
    source_box = report.source_box_pqs_ida_component_smoke
    final_residual = report.final_residual_mwg_supplement_component_facts

    fixed_fixed_shape = _pqs_component_route_smoke_parse_shape(
        final_residual.fixed_fixed_shape,
        :fixed_fixed_shape,
    )
    fixed_residual_shape = _pqs_component_route_smoke_parse_shape(
        final_residual.fixed_residual_shape,
        :fixed_residual_shape,
    )
    residual_residual_shape = _pqs_component_route_smoke_parse_shape(
        final_residual.residual_residual_shape,
        :residual_residual_shape,
    )
    final_interaction_shape = _pqs_component_route_smoke_parse_shape(
        final_residual.final_interaction_shape,
        :final_interaction_shape,
    )
    fixed_dimension = fixed_fixed_shape[1]
    residual_dimension = residual_residual_shape[1]
    final_dimension = final_interaction_shape[1]
    residual_owner_indices = _pqs_component_route_smoke_parse_int_tuple(
        final_residual.residual_nucleus_indices,
        :residual_nucleus_indices,
    )
    fixed_range = fixed_dimension == 0 ? (1:0) : (1:fixed_dimension)
    residual_range = residual_dimension == 0 ?
        ((fixed_dimension + 1):fixed_dimension) :
        ((fixed_dimension + 1):(fixed_dimension + residual_dimension))
    final_range = final_dimension == 0 ? (1:0) : (1:final_dimension)
    residual_owner_rows_match =
        length(residual_owner_indices) == residual_dimension
    shape_consistent =
        fixed_fixed_shape == (fixed_dimension, fixed_dimension) &&
        fixed_residual_shape == (fixed_dimension, residual_dimension) &&
        residual_residual_shape == (residual_dimension, residual_dimension) &&
        final_interaction_shape == (final_dimension, final_dimension) &&
        final_dimension == fixed_dimension + residual_dimension

    if strict
        shape_consistent || throw(
            DimensionMismatch("CR2 sidecar schema final-residual MWG component shapes are inconsistent"),
        )
        residual_owner_rows_match || throw(
            DimensionMismatch("CR2 sidecar schema residual owner count does not match residual rows"),
        )
        report.lane_boundaries.no_owner_inference_from_raw_to_final_support ||
            throw(ArgumentError("CR2 sidecar schema cannot infer owners from raw_to_final support"))
        report.lane_boundaries.no_raw_gto_gto_mwg_blocks ||
            throw(ArgumentError("CR2 sidecar schema cannot include raw GTO/GTO MWG blocks"))
        report.lane_boundaries.no_fixed_raw_gto_mwg_blocks ||
            throw(ArgumentError("CR2 sidecar schema cannot include fixed/raw-GTO MWG blocks"))
        report.lane_boundaries.no_retained_weight_or_ida_division ||
            throw(ArgumentError("CR2 sidecar schema cannot include retained-weight or IDA division"))
    end

    source_unit_records = source_box.retained_units
    source_unit_label_status = isempty(source_unit_records) ?
        :unavailable : :explicit_route_descriptor_unit_keys
    component_rows = Tuple(
        _pqs_component_route_smoke_source_box_sidecar_components(row)
        for row in source_box.rows
    )
    fixed_side_metadata_available =
        !isnothing(fixed_side_retained_unit_metadata)
    if fixed_side_metadata_available
        fixed_side_retained_unit_metadata.object_kind ==
            :pqs_current_route_fixed_side_retained_unit_metadata || throw(
            ArgumentError("CR2 sidecar schema fixed-side metadata requires _pqs_current_route_fixed_side_retained_unit_metadata output"),
        )
        if strict
            fixed_side_retained_unit_metadata.labels.shell_label_status ==
                :unavailable || throw(
                ArgumentError("CR2 sidecar schema fixed-side metadata cannot provide shell labels"),
            )
            !fixed_side_retained_unit_metadata.labels.label_reconstruction_from_centers ||
                throw(
                    ArgumentError("CR2 sidecar schema fixed-side metadata cannot reconstruct labels from centers"),
                )
            !fixed_side_retained_unit_metadata.labels.nearest_grid_or_center_label_heuristic ||
                throw(
                    ArgumentError("CR2 sidecar schema fixed-side metadata cannot use nearest-grid label heuristics"),
                )
            !fixed_side_retained_unit_metadata.diagnostics.retained_weight_or_ida_division ||
                throw(
                    ArgumentError("CR2 sidecar schema fixed-side metadata cannot include retained-weight or IDA division"),
                )
            fixed_side_retained_unit_metadata.diagnostics.shell_realized_pqs_source_box_operator_ready_count ==
                0 || throw(
                ArgumentError("CR2 sidecar schema fixed-side metadata cannot mark shell-realized PQS fixtures as source-box-operator-ready without an explicit framework update"),
            )
            fixed_side_retained_unit_metadata.diagnostics.shell_realized_pqs_fixtures_are_metadata_oracle_only ||
                throw(
                    ArgumentError("CR2 sidecar schema fixed-side metadata requires shell-realized PQS fixtures to remain metadata/oracle-only"),
            )
        end
    end
    source_shell_mode_inventory_available =
        !isnothing(source_shell_mode_inventory)
    if source_shell_mode_inventory_available
        source_shell_mode_inventory.object_kind ==
            :pqs_current_route_source_shell_mode_inventory || throw(
            ArgumentError("CR2 sidecar schema source-shell/source-mode inventory requires _pqs_current_route_source_shell_mode_inventory output"),
        )
        if strict
            source_shell_mode_inventory.fixed_dimension == fixed_dimension ||
                throw(
                    DimensionMismatch("CR2 sidecar schema source-shell/source-mode inventory fixed dimension does not match sidecar fixed dimension"),
                )
            !source_shell_mode_inventory.diagnostics.inferred_from_centers ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot infer labels from centers"),
                )
            !source_shell_mode_inventory.diagnostics.inferred_from_nearest_grid ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot infer labels from nearest-grid searches"),
                )
            !source_shell_mode_inventory.diagnostics.inferred_from_support_order ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot infer labels from support order"),
                )
            !source_shell_mode_inventory.diagnostics.inferred_from_support_indices ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot infer labels from support indices"),
                )
            !source_shell_mode_inventory.diagnostics.inferred_from_raw_to_final_support ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot infer labels from raw_to_final support"),
                )
            !source_shell_mode_inventory.diagnostics.retained_weight_or_ida_division ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot include retained-weight or IDA division"),
                )
            !source_shell_mode_inventory.diagnostics.route_construction_changed ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot mutate route construction"),
                )
            !source_shell_mode_inventory.diagnostics.packet_adoption ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot adopt packet/fixed-block behavior"),
                )
            !source_shell_mode_inventory.diagnostics.qwhamiltonian_changed ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot change QW/Hamiltonian behavior"),
                )
            !source_shell_mode_inventory.diagnostics.hamiltonian_matrix_built ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot build Hamiltonian matrices"),
                )
            !source_shell_mode_inventory.diagnostics.public_default_consumes ||
                throw(
                    ArgumentError("CR2 sidecar schema source-shell/source-mode inventory cannot change public/default consumption"),
                )
        end
    end
    source_shell_count = source_shell_mode_inventory_available ?
        Int(source_shell_mode_inventory.source_shell_count) : 0
    source_mode_count = source_shell_mode_inventory_available ?
        Int(source_shell_mode_inventory.source_mode_count) : 0
    source_shell_mode_center_status = source_shell_mode_inventory_available ?
        source_shell_mode_inventory.center_status : :unavailable
    source_shell_mode_product_doside_only =
        source_shell_mode_inventory_available ?
        source_shell_mode_inventory.status == :product_doside_source_shell_modes_only :
        false
    source_shell_mode_support_dense_unavailable_count =
        source_shell_mode_inventory_available ?
        Int(
            source_shell_mode_inventory.diagnostics.support_dense_unavailable_column_count,
        ) : 0
    source_shell_mode_shell_realized_unavailable_count =
        source_shell_mode_inventory_available ?
        Int(
            source_shell_mode_inventory.diagnostics.shell_realized_pqs_unavailable_column_count,
        ) : 0

    return (
        object_kind =
            :pqs_pqs_product_component_route_smoke_cr2_sidecar_schema,
        status = :private_cr2_sidecar_schema,
        schema_version = :pqs_component_cr2_sidecar_private_v1,
        report_object_kind = report.object_kind,
        lanes = (
            source_box_pqs_ida_fixed_side = (
                status = :private_source_box_component_smoke,
                algorithm_lane = :source_box_pqs_ida,
                route_shape = source_box.route_shape,
                parent_dims = source_box.parent_dims,
                source_mode_dims = source_box.source_mode_dims,
                left_source_box = source_box.left_source_box,
                right_source_box = source_box.right_source_box,
                product_source_box = source_box.product_source_box,
                retained_dimension = first(source_box.rows).retained_dimension,
                retained_units = source_unit_records,
                retained_ranges = source_box.retained_ranges,
                source_unit_label_status = source_unit_label_status,
                source_unit_labels =
                    Tuple(unit.unit_key for unit in source_unit_records),
                components = component_rows,
            ),
            final_residual_mwg_supplement = (
                status = :ordinary_final_residual_component_facts,
                algorithm_lane = :final_residual_mwg_supplement,
                source_report = final_residual.source_report,
                route = final_residual.route,
                component_helper = final_residual.component_helper,
                fixed_dimension = fixed_dimension,
                residual_dimension = residual_dimension,
                final_dimension = final_dimension,
                fixed_column_range = fixed_range,
                residual_column_range = residual_range,
                final_column_range = final_range,
                residual_owner_metadata = (
                    status = :explicit,
                    residual_nucleus_indices = residual_owner_indices,
                    residual_owner_counts =
                        final_residual.residual_owner_counts,
                    owner_metadata_source =
                        final_residual.owner_metadata_source,
                    owner_semantics_inferred_from_raw_to_final_support =
                        false,
                    owner_count_matches_residual_rows =
                        residual_owner_rows_match,
                ),
                components = (
                    fixed_fixed = (
                        available = true,
                        shape = fixed_fixed_shape,
                        authority_error = final_residual.max_authority_error,
                        provenance = final_residual.component_helper,
                    ),
                    fixed_residual = (
                        available = true,
                        shape = fixed_residual_shape,
                        authority_error = final_residual.max_authority_error,
                        provenance = final_residual.component_helper,
                    ),
                    residual_residual = (
                        available = true,
                        shape = residual_residual_shape,
                        authority_error = final_residual.max_authority_error,
                        provenance = final_residual.component_helper,
                    ),
                    final_interaction = (
                        available = true,
                        shape = final_interaction_shape,
                        authority_error = final_residual.max_authority_error,
                        provenance = final_residual.component_helper,
                    ),
                ),
            ),
        ),
        fixed_side_retained_unit_metadata =
            fixed_side_retained_unit_metadata,
        source_shell_mode_inventory = source_shell_mode_inventory,
        labels = (
            source_unit_label_status = source_unit_label_status,
            source_unit_labels =
                Tuple(unit.unit_key for unit in source_unit_records),
            shell_label_status = :unavailable,
            shell_labels = (),
            label_reconstruction_from_centers = false,
            nearest_grid_or_center_label_heuristic = false,
        ),
        absences_by_contract = (
            raw_gto_gto_mwg_interaction_blocks = true,
            fixed_raw_gto_mwg_interaction_blocks = true,
            owner_inference_from_raw_to_final_support = true,
            retained_source_box_final_residual_weight_division = true,
            retained_weight_ida_division = true,
            packet_fixed_block_qw_hamiltonian_adoption = true,
            public_default_route = true,
            ecp_scf_hf_cr2_science_claim = true,
            mwg_ida_semantic_change = true,
        ),
        provenance = merge(
            (
                source = :pqs_pqs_product_component_route_smoke_cr2_sidecar_schema,
                report_adapter = report.object_kind,
                generated_at = report.generated_at,
            ),
            provenance,
        ),
        diagnostics = (
            source =
                :pqs_pqs_product_component_route_smoke_cr2_sidecar_schema,
            private_shadow_only = true,
            source_box_pqs_ida_and_mwg_residual_same_algorithm = false,
            lanes_remain_separate = true,
            source_box_unit_records_available =
                !isempty(source_unit_records),
            fixed_side_retained_unit_metadata_available =
                fixed_side_metadata_available,
            source_shell_mode_inventory_available =
                source_shell_mode_inventory_available,
            source_shell_count = source_shell_count,
            source_mode_count = source_mode_count,
            source_shell_mode_center_status =
                source_shell_mode_center_status,
            source_shell_mode_product_doside_only =
                source_shell_mode_product_doside_only,
            source_shell_mode_support_dense_unavailable_column_count =
                source_shell_mode_support_dense_unavailable_count,
            source_shell_mode_shell_realized_pqs_unavailable_column_count =
                source_shell_mode_shell_realized_unavailable_count,
            source_unit_label_status = source_unit_label_status,
            shell_label_status = :unavailable,
            label_reconstruction_from_centers = false,
            nearest_grid_or_center_label_heuristic = false,
            fixed_dimension = fixed_dimension,
            residual_dimension = residual_dimension,
            final_dimension = final_dimension,
            final_residual_shape_consistent = shape_consistent,
            residual_owner_rows_match = residual_owner_rows_match,
            no_owner_inference_from_raw_to_final_support =
                report.lane_boundaries.no_owner_inference_from_raw_to_final_support,
            no_raw_gto_gto_mwg_blocks =
                report.lane_boundaries.no_raw_gto_gto_mwg_blocks,
            no_fixed_raw_gto_mwg_blocks =
                report.lane_boundaries.no_fixed_raw_gto_mwg_blocks,
            no_retained_weight_or_ida_division =
                report.lane_boundaries.no_retained_weight_or_ida_division,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            public_default_consumes = false,
            ecp_terms_implemented = false,
            scf_hf_validation_claim = false,
            cr2_science_status_changed = false,
            mwg_ida_semantics_changed = false,
        ),
    )
end

function _write_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema_report(
    io::IO,
    sidecar,
)
    sidecar.object_kind ==
        :pqs_pqs_product_component_route_smoke_cr2_sidecar_schema || throw(
        ArgumentError("CR2 sidecar schema report writer requires _pqs_pqs_product_component_route_smoke_cr2_sidecar_schema output"),
    )
    source_lane = sidecar.lanes.source_box_pqs_ida_fixed_side
    mwg_lane = sidecar.lanes.final_residual_mwg_supplement

    println(io, "Be2/PQS CR2 sidecar schema summary")
    _pqs_component_route_smoke_print_kv(io, "schema_version", sidecar.schema_version)
    _pqs_component_route_smoke_print_kv(io, "status", sidecar.status)
    _pqs_component_route_smoke_print_kv(io, "source_report", mwg_lane.source_report)
    _pqs_component_route_smoke_print_kv(io, "single_algorithmic_operator", false)
    _pqs_component_route_smoke_print_kv(io, "production_route_adoption", false)
    _pqs_component_route_smoke_print_kv(io, "public_default_route", false)
    println(io)

    println(io, "[source_box_pqs_ida_fixed_side]")
    _pqs_component_route_smoke_print_kv(io, "algorithm_lane", source_lane.algorithm_lane)
    _pqs_component_route_smoke_print_kv(io, "route_shape", source_lane.route_shape)
    _pqs_component_route_smoke_print_kv(io, "parent_dims", source_lane.parent_dims)
    _pqs_component_route_smoke_print_kv(io, "source_mode_dims", source_lane.source_mode_dims)
    _pqs_component_route_smoke_print_kv(io, "left_source_box", source_lane.left_source_box)
    _pqs_component_route_smoke_print_kv(io, "right_source_box", source_lane.right_source_box)
    _pqs_component_route_smoke_print_kv(io, "product_source_box", source_lane.product_source_box)
    _pqs_component_route_smoke_print_kv(io, "retained_dimension", source_lane.retained_dimension)
    _pqs_component_route_smoke_print_kv(io, "retained_ranges", source_lane.retained_ranges)
    _pqs_component_route_smoke_print_kv(
        io,
        "source_unit_label_status",
        source_lane.source_unit_label_status,
    )
    _pqs_component_route_smoke_print_kv(io, "source_unit_labels", source_lane.source_unit_labels)
    for unit in source_lane.retained_units
        prefix = "unit.$(unit.unit_key)"
        _pqs_component_route_smoke_print_kv(io, "$prefix.kind", unit.retained_unit_kind)
        _pqs_component_route_smoke_print_kv(io, "$prefix.source_family", unit.source_family)
        _pqs_component_route_smoke_print_kv(io, "$prefix.retained_range", unit.retained_range)
        _pqs_component_route_smoke_print_kv(io, "$prefix.retained_count", unit.retained_count)
        _pqs_component_route_smoke_print_kv(io, "$prefix.source_dimensions", unit.source_dimensions)
    end
    for component in source_lane.components
        prefix = "component.$(component.mode)"
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.retained_matrix_shape",
            component.retained_matrix_shape,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.nuclear_available",
            component.nuclear_attraction_by_center.available,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.nuclear_pair_count",
            component.nuclear_attraction_by_center.pair_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.nuclear_pair_family_counts",
            component.nuclear_attraction_by_center.pair_family_counts,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.nuclear_authority_error",
            component.nuclear_attraction_by_center.authority_error,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.electron_electron_available",
            component.electron_electron_density_density.available,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.electron_electron_representation",
            component.electron_electron_density_density.representation,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.electron_electron_pair_count",
            component.electron_electron_density_density.pair_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.electron_electron_pair_family_counts",
            component.electron_electron_density_density.pair_family_counts,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.electron_electron_dense_parent_error",
            component.electron_electron_density_density.dense_parent_authority_max_error,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.source_weight_division_owner",
            component.source_weight_division_owner,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "$prefix.source_weight_division_applied_by_helper",
            component.source_weight_division_applied_by_helper,
        )
        _pqs_component_route_smoke_print_kv(io, "$prefix.no_go_clear", component.no_go_clear)
    end
    println(io)

    println(io, "[final_residual_mwg_supplement]")
    _pqs_component_route_smoke_print_kv(io, "algorithm_lane", mwg_lane.algorithm_lane)
    _pqs_component_route_smoke_print_kv(io, "component_helper", mwg_lane.component_helper)
    _pqs_component_route_smoke_print_kv(io, "fixed_dimension", mwg_lane.fixed_dimension)
    _pqs_component_route_smoke_print_kv(io, "residual_dimension", mwg_lane.residual_dimension)
    _pqs_component_route_smoke_print_kv(io, "final_dimension", mwg_lane.final_dimension)
    _pqs_component_route_smoke_print_kv(io, "fixed_column_range", mwg_lane.fixed_column_range)
    _pqs_component_route_smoke_print_kv(io, "residual_column_range", mwg_lane.residual_column_range)
    _pqs_component_route_smoke_print_kv(io, "final_column_range", mwg_lane.final_column_range)
    _pqs_component_route_smoke_print_kv(
        io,
        "residual_nucleus_indices",
        mwg_lane.residual_owner_metadata.residual_nucleus_indices,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "residual_owner_counts",
        mwg_lane.residual_owner_metadata.residual_owner_counts,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "owner_metadata_source",
        mwg_lane.residual_owner_metadata.owner_metadata_source,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "owner_count_matches_residual_rows",
        mwg_lane.residual_owner_metadata.owner_count_matches_residual_rows,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "owner_semantics_inferred_from_raw_to_final_support",
        mwg_lane.residual_owner_metadata.owner_semantics_inferred_from_raw_to_final_support,
    )
    for name in (:fixed_fixed, :fixed_residual, :residual_residual, :final_interaction)
        component = getproperty(mwg_lane.components, name)
        _pqs_component_route_smoke_print_kv(io, "component.$name.available", component.available)
        _pqs_component_route_smoke_print_kv(io, "component.$name.shape", component.shape)
        _pqs_component_route_smoke_print_kv(
            io,
            "component.$name.authority_error",
            component.authority_error,
        )
        _pqs_component_route_smoke_print_kv(io, "component.$name.provenance", component.provenance)
    end
    println(io)

    if !isnothing(sidecar.fixed_side_retained_unit_metadata)
        fixed_side = sidecar.fixed_side_retained_unit_metadata
        println(io, "[fixed_side_retained_unit_metadata]")
        _pqs_component_route_smoke_print_kv(io, "status", fixed_side.status)
        _pqs_component_route_smoke_print_kv(io, "schema_version", fixed_side.schema_version)
        _pqs_component_route_smoke_print_kv(io, "fixed_dimension", fixed_side.fixed_dimension)
        _pqs_component_route_smoke_print_kv(io, "unit_count", fixed_side.unit_count)
        _pqs_component_route_smoke_print_kv(
            io,
            "source_unit_label_status",
            fixed_side.labels.source_unit_label_status,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "source_unit_labels",
            fixed_side.labels.source_unit_labels,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "shell_label_status",
            fixed_side.labels.shell_label_status,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "label_reconstruction_from_centers",
            fixed_side.labels.label_reconstruction_from_centers,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "nearest_grid_or_center_label_heuristic",
            fixed_side.labels.nearest_grid_or_center_label_heuristic,
        )
        for unit in fixed_side.retained_units
            prefix = "unit.$(unit.unit_key)"
            _pqs_component_route_smoke_print_kv(io, "$prefix.category", unit.category)
            _pqs_component_route_smoke_print_kv(io, "$prefix.kind", unit.kind)
            _pqs_component_route_smoke_print_kv(
                io,
                "$prefix.retained_range",
                unit.retained_range,
            )
            _pqs_component_route_smoke_print_kv(
                io,
                "$prefix.retained_count",
                unit.retained_count,
            )
            _pqs_component_route_smoke_print_kv(
                io,
                "$prefix.support_count",
                unit.support_count,
            )
            _pqs_component_route_smoke_print_kv(
                io,
                "$prefix.source_mode_dims",
                unit.source_mode_dims,
            )
            _pqs_component_route_smoke_print_kv(
                io,
                "$prefix.representation_kind",
                unit.representation_kind,
            )
            _pqs_component_route_smoke_print_kv(
                io,
                "$prefix.shell_realized_pqs_metadata_oracle_fixture",
                unit.shell_realized_pqs_metadata_oracle_fixture,
            )
            _pqs_component_route_smoke_print_kv(
                io,
                "$prefix.shell_realized_pqs_source_box_operator_ready",
                unit.shell_realized_pqs_source_box_operator_ready,
            )
        end
        println(io)
    end

    if !isnothing(sidecar.source_shell_mode_inventory)
        source_shell_modes = sidecar.source_shell_mode_inventory
        println(io, "[source_shell_mode_inventory]")
        _pqs_component_route_smoke_print_kv(io, "status", source_shell_modes.status)
        _pqs_component_route_smoke_print_kv(
            io,
            "schema_version",
            source_shell_modes.schema_version,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "source_shell_count",
            source_shell_modes.source_shell_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "source_mode_count",
            source_shell_modes.source_mode_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "center_status",
            source_shell_modes.center_status,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "covered_unit_categories",
            source_shell_modes.covered_unit_categories,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "non_product_source_mode_status",
            source_shell_modes.non_product_source_mode_status,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "source_mode_label_status",
            source_shell_modes.source_mode_label_status,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "product_doside_only",
            source_shell_modes.status == :product_doside_source_shell_modes_only,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "product_doside_source_shell_count",
            source_shell_modes.diagnostics.product_doside_source_shell_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "product_doside_source_mode_count",
            source_shell_modes.diagnostics.product_doside_source_mode_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "support_dense_source_shell_count",
            source_shell_modes.diagnostics.support_dense_source_shell_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "support_dense_source_mode_count",
            source_shell_modes.diagnostics.support_dense_source_mode_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "support_dense_unavailable_column_count",
            source_shell_modes.diagnostics.support_dense_unavailable_column_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "shell_realized_pqs_source_shell_count",
            source_shell_modes.diagnostics.shell_realized_pqs_source_shell_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "shell_realized_pqs_source_mode_count",
            source_shell_modes.diagnostics.shell_realized_pqs_source_mode_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "shell_realized_pqs_unavailable_column_count",
            source_shell_modes.diagnostics.shell_realized_pqs_unavailable_column_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "total_unavailable_unit_count",
            source_shell_modes.diagnostics.total_unavailable_unit_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "total_unavailable_column_count",
            source_shell_modes.diagnostics.total_unavailable_column_count,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "inferred_from_centers",
            source_shell_modes.diagnostics.inferred_from_centers,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "inferred_from_nearest_grid",
            source_shell_modes.diagnostics.inferred_from_nearest_grid,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "inferred_from_support_order",
            source_shell_modes.diagnostics.inferred_from_support_order,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "inferred_from_support_indices",
            source_shell_modes.diagnostics.inferred_from_support_indices,
        )
        _pqs_component_route_smoke_print_kv(
            io,
            "inferred_from_raw_to_final_support",
            source_shell_modes.diagnostics.inferred_from_raw_to_final_support,
        )
        println(io)
    end

    println(io, "[labels]")
    _pqs_component_route_smoke_print_kv(
        io,
        "source_unit_label_status",
        sidecar.labels.source_unit_label_status,
    )
    _pqs_component_route_smoke_print_kv(io, "source_unit_labels", sidecar.labels.source_unit_labels)
    _pqs_component_route_smoke_print_kv(io, "shell_label_status", sidecar.labels.shell_label_status)
    _pqs_component_route_smoke_print_kv(io, "shell_labels", sidecar.labels.shell_labels)
    _pqs_component_route_smoke_print_kv(
        io,
        "label_reconstruction_from_centers",
        sidecar.labels.label_reconstruction_from_centers,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "nearest_grid_or_center_label_heuristic",
        sidecar.labels.nearest_grid_or_center_label_heuristic,
    )
    println(io)

    println(io, "[boundaries]")
    for key in keys(sidecar.absences_by_contract)
        _pqs_component_route_smoke_print_kv(
            io,
            string(key, "_absent_by_contract"),
            getproperty(sidecar.absences_by_contract, key),
        )
    end
    _pqs_component_route_smoke_print_kv(
        io,
        "lanes_remain_separate",
        sidecar.diagnostics.lanes_remain_separate,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "final_residual_shape_consistent",
        sidecar.diagnostics.final_residual_shape_consistent,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "residual_owner_rows_match",
        sidecar.diagnostics.residual_owner_rows_match,
    )
    _pqs_component_route_smoke_print_kv(io, "packet_adoption", sidecar.diagnostics.packet_adoption)
    _pqs_component_route_smoke_print_kv(io, "qwhamiltonian_consumes", sidecar.diagnostics.qwhamiltonian_consumes)
    _pqs_component_route_smoke_print_kv(
        io,
        "hamiltonian_matrix_built",
        sidecar.diagnostics.hamiltonian_matrix_built,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "cr2_science_status_changed",
        sidecar.diagnostics.cr2_science_status_changed,
    )
    return sidecar
end

function _write_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema_report(
    path::AbstractString,
    sidecar,
)
    open(path, "w") do io
        _write_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema_report(
            io,
            sidecar,
        )
    end
    return path
end

function _pqs_component_route_smoke_print_kv(io, key, value)
    println(io, key, "\t", value)
end

function _pqs_component_route_smoke_print_fields(io, record, fields)
    for key in fields
        _pqs_component_route_smoke_print_kv(
            io,
            string(key),
            getproperty(record, key),
        )
    end
end

function _write_pqs_pqs_product_private_source_box_route_adapter_readiness_summary(
    io::IO,
    readiness,
)
    readiness.object_kind ==
        :pqs_pqs_product_private_source_box_route_adapter_readiness_summary ||
        throw(
            ArgumentError("private route-adapter readiness writer requires _pqs_pqs_product_private_source_box_route_adapter_readiness_summary output"),
        )

    println(io, "[private_route_adapter_readiness]")
    _pqs_component_route_smoke_print_fields(
        io,
        readiness,
        (
            :status,
            :ready_for_next_private_adapter_pass,
            :missing_required_pieces,
            :missing_optional_pieces,
            :no_go_violations,
            :pair_factor_normalization_modes,
        ),
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "authority_available_modes",
        readiness.ida_source_box_electron_electron.authority_comparison.available_modes,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "authority_skip_reasons",
        readiness.ida_source_box_electron_electron.authority_comparison.skip_reasons,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "authority_comparison_accounted_for",
        readiness.ida_source_box_electron_electron.authority_comparison.accounted_for,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "timing_allocation_available",
        readiness.timing_allocation.available,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "lanes_remain_separate",
        readiness.diagnostics.lanes_remain_separate,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "source_box_ida_and_mwg_residual_same_algorithm",
        readiness.diagnostics.source_box_ida_and_mwg_residual_same_algorithm,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "public_default_consumes",
        readiness.diagnostics.public_default_consumes,
    )
    _pqs_component_route_smoke_print_kv(
        io,
        "construction_behavior_changed",
        readiness.diagnostics.construction_behavior_changed,
    )
    for key in keys(readiness.no_go_flags)
        _pqs_component_route_smoke_print_kv(
            io,
            "no_go.$key",
            getproperty(readiness.no_go_flags, key),
        )
    end
    return readiness
end

function _write_pqs_pqs_product_component_route_smoke_report(io::IO, report)
    report.object_kind == :pqs_pqs_product_component_route_smoke_report_adapter ||
        throw(
            ArgumentError("component route smoke report writer requires _pqs_pqs_product_component_route_smoke_report_adapter output"),
        )
    source_box = report.source_box_pqs_ida_component_smoke

    println(io, report.title)
    _pqs_component_route_smoke_print_kv(io, "generated_at", report.generated_at)
    _pqs_component_route_smoke_print_kv(io, "status", report.report_status)
    _pqs_component_route_smoke_print_fields(
        io,
        report,
        (
            :smallest_honest_smoke,
            :single_algorithmic_operator,
            :hamiltonian_built,
            :production_route_adoption,
            :public_default_route,
        ),
    )
    println(io)

    println(io, "[source_box_pqs_ida_component_smoke]")
    _pqs_component_route_smoke_print_fields(
        io,
        source_box,
        (
            :route_shape,
            :parent_dims,
            :source_mode_dims,
            :left_source_box,
            :right_source_box,
            :product_source_box,
            :nuclear_centers,
            :nuclear_charges,
            :pqs_source_box_fixed_side_facts_available,
            :by_center_nuclear_attraction_available,
            :ida_source_box_electron_electron_available,
            :electron_electron_representation,
            :mwg_supplement_residual_adapted_in_source_box_smoke,
            :retained_unit_count,
            :source_unit_label_status,
            :source_unit_labels,
            :retained_ranges,
        ),
    )

    for row in source_box.rows
        println(io)
        println(io, "[source_box_pqs_ida_component_smoke.", row.mode, "]")
        _pqs_component_route_smoke_print_fields(
            io,
            row,
            (
                :retained_dimension,
                :nuclear_pair_count,
                :electron_electron_pair_count,
                :ida_term_count,
                :output_finite,
                :nuclear_symmetry_error,
                :electron_electron_symmetry_error,
                :nuclear_total_from_center_error,
                :dense_parent_ida_authority_available,
                :dense_parent_ida_authority_max_error,
                :dense_parent_ida_authority_skip_reason,
                :source_weight_division_owner,
                :source_weight_division_applied_by_helper,
                :no_go_clear,
                :mwg_supplement_residual_path,
            ),
        )
    end
    println(io)

    println(io, "[final_residual_mwg_supplement_component_facts]")
    _pqs_component_route_smoke_print_fields(
        io,
        report.final_residual_mwg_supplement_component_facts,
        keys(report.final_residual_mwg_supplement_component_facts),
    )
    println(io)

    println(io, "[lane_boundaries]")
    _pqs_component_route_smoke_print_fields(
        io,
        report.lane_boundaries,
        keys(report.lane_boundaries),
    )
    println(io)

    readiness =
        _pqs_pqs_product_private_source_box_route_adapter_readiness_summary(
            report,
        )
    _write_pqs_pqs_product_private_source_box_route_adapter_readiness_summary(
        io,
        readiness,
    )
    return report
end

function _write_pqs_pqs_product_component_route_smoke_report(
    path::AbstractString,
    report,
)
    open(path, "w") do io
        _write_pqs_pqs_product_component_route_smoke_report(io, report)
    end
    return path
end
