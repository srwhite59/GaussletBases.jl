using GaussletBases

const _BE2_PQS_Q5_ACCEPTANCE_BACKEND = :pgdg_localized_experimental

mutable struct _PQSSourceMetadataAcceptanceChecks
    rows::Vector{Pair{String,Any}}
    failures::Vector{String}
end

_PQSSourceMetadataAcceptanceChecks() =
    _PQSSourceMetadataAcceptanceChecks(Pair{String,Any}[], String[])

function _pqs_source_metadata_acceptance_record!(
    checks::_PQSSourceMetadataAcceptanceChecks,
    name::AbstractString,
    value,
)
    push!(checks.rows, String(name) => value)
    return value
end

function _pqs_source_metadata_acceptance_require!(
    checks::_PQSSourceMetadataAcceptanceChecks,
    name::AbstractString,
    value::Bool,
)
    _pqs_source_metadata_acceptance_record!(checks, name, value)
    value || push!(checks.failures, String(name))
    return value
end

function _pqs_source_metadata_acceptance_read_key_value_tsv(path::AbstractString)
    isfile(path) || error("missing metadata TSV: $(path)")
    data = Dict{String,String}()
    for (line_index, line) in enumerate(readlines(path))
        text = strip(line)
        isempty(text) && continue
        startswith(text, "#") && continue
        fields = split(chomp(line), '\t')
        length(fields) >= 2 || error(
            "bad metadata row $(line_index) in $(path): $(line)",
        )
        key = strip(fields[1])
        value = strip(fields[2])
        lowercase(key) in ("key", "name") &&
            lowercase(value) in ("value", "values") &&
            continue
        data[key] = value
    end
    return data
end

function _pqs_source_metadata_acceptance_metadata_value(metadata, key)
    haskey(metadata, key) || error("metadata missing required key $(key)")
    return metadata[key]
end

function _pqs_source_metadata_acceptance_parse_symbol(text::AbstractString)
    return Symbol(replace(strip(text), ":" => ""))
end

function _pqs_source_metadata_acceptance_parse_axis_counts(text::AbstractString)
    cleaned = replace(replace(replace(text, "(" => ""), ")" => ""), "," => " ")
    values = parse.(Int, split(cleaned))
    length(values) == 3 || error(
        "parent_axis_counts must have three entries, got $(text)",
    )
    all(>(0), values) || error(
        "parent_axis_counts must be positive, got $(text)",
    )
    return (values[1], values[2], values[3])
end

function _be2_pqs_q5_source_metadata_acceptance_metadata(
    artifact_dir::AbstractString,
)
    metadata = _pqs_source_metadata_acceptance_read_key_value_tsv(
        joinpath(artifact_dir, "q_row_capture_h1_metadata.tsv"),
    )
    value(key) = _pqs_source_metadata_acceptance_metadata_value(metadata, key)
    parsed = (
        family = _pqs_source_metadata_acceptance_parse_symbol(value("family")),
        bond_length = parse(Float64, value("bond_length")),
        core_spacing = parse(Float64, value("core_spacing")),
        xmax_parallel = parse(Float64, value("xmax_parallel")),
        xmax_transverse = parse(Float64, value("xmax_transverse")),
        bond_axis =
            _pqs_source_metadata_acceptance_parse_symbol(value("bond_axis")),
        nuclear_charge = parse(Float64, value("nuclear_charge")),
        reference_spacing = parse(Float64, value("reference_spacing")),
        tail_spacing = parse(Float64, value("tail_spacing")),
        shared_q = parse(Int, value("shared_q")),
        shared_order = parse(Int, value("shared_order")),
        protected_atom_side_count =
            parse(Int, value("protected_atom_side_count")),
        q_min = parse(Int, value("q_min")),
        nside = parse(Int, value("nside")),
        parent_axis_counts =
            _pqs_source_metadata_acceptance_parse_axis_counts(
                value("parent_axis_counts"),
            ),
        parent_dimension = parse(Int, value("parent_dimension")),
    )
    parsed.family == :G10 || error("expected family G10, got $(parsed.family)")
    parsed.bond_axis == :z || error(
        "expected z bond axis, got $(parsed.bond_axis)",
    )
    parsed.nuclear_charge == 4.0 || error(
        "expected all-electron Z=4, got $(parsed.nuclear_charge)",
    )
    parsed.shared_q == 5 || error(
        "expected q=5 CR2 target, got shared_q=$(parsed.shared_q)",
    )
    parsed.shared_order == 5 || error(
        "expected shared_order=5, got $(parsed.shared_order)",
    )
    parsed.parent_axis_counts == (15, 15, 27) || error(
        "expected parent axes (15,15,27), got $(parsed.parent_axis_counts)",
    )
    parsed.parent_dimension == prod(parsed.parent_axis_counts) || error(
        "parent_dimension $(parsed.parent_dimension) does not match axes $(parsed.parent_axis_counts)",
    )
    parsed.parent_dimension == 6075 || error(
        "expected parent dimension 6075, got $(parsed.parent_dimension)",
    )
    return parsed
end

function _be2_pqs_q5_source_metadata_acceptance_construction(parsed)
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = bond_aligned_homonuclear_qw_basis(
        family = parsed.family,
        bond_length = parsed.bond_length,
        core_spacing = parsed.core_spacing,
        xmax_parallel = parsed.xmax_parallel,
        xmax_transverse = parsed.xmax_transverse,
        bond_axis = parsed.bond_axis,
        nuclear_charge = parsed.nuclear_charge,
        reference_spacing = parsed.reference_spacing,
        tail_spacing = parsed.tail_spacing,
    )
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(
        basis,
        expansion;
        gausslet_backend = _BE2_PQS_Q5_ACCEPTANCE_BACKEND,
    )
    axis_counts = GaussletBases._nested_axis_lengths(bundles)
    axis_counts == parsed.parent_axis_counts || error(
        "repo-built parent axes $(axis_counts) do not match CR2 axes $(parsed.parent_axis_counts)",
    )
    policy0 = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        basis,
        bundles;
        protected_atom_side_count = parsed.protected_atom_side_count,
        q_min = parsed.q_min,
    )
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy0.construction_plan;
        q_min = parsed.q_min,
        atom_q = 4,
        atom_order = 4,
        shared_q = parsed.shared_q,
        shared_order = parsed.shared_order,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = parsed.nside,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :support_reference,
            shared_shell_realization = :projected_q_shell,
        )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            construction,
        )
    diagnostics.metadata.shared_shell_realization == :projected_q_shell ||
        error("PQS opt-in realization was not selected")
    diagnostics.metadata.projected_q_shell_opt_in ||
        error("PQS opt-in diagnostic is false")
    diagnostics.support_coverage.coverage_ok ||
        error("PQS source did not cover the parent")
    return construction
end

function _pqs_source_metadata_acceptance_tsv_records(text::AbstractString)
    lines = split(chomp(text), '\n')
    isempty(lines) && error("empty TSV text")
    header = split(lines[1], '\t'; keepempty = true)
    records = Vector{Dict{String,String}}()
    for (line_index, line) in enumerate(lines[2:end])
        values = split(line, '\t'; keepempty = true)
        length(values) == length(header) || error(
            "TSV row $(line_index + 1) has $(length(values)) fields, expected $(length(header))",
        )
        push!(records, Dict(header .=> values))
    end
    return header, records
end

function _pqs_source_metadata_acceptance_table_records(
    writer,
    source_inventory;
    path = nothing,
)
    if isnothing(path)
        io = IOBuffer()
        writer(io, source_inventory)
        return _pqs_source_metadata_acceptance_tsv_records(
            String(take!(io)),
        )
    end
    writer(path, source_inventory)
    return _pqs_source_metadata_acceptance_tsv_records(read(path, String))
end

function _pqs_source_metadata_acceptance_all_records_equal(
    records,
    key::AbstractString,
    value::AbstractString,
)
    return all(record -> get(record, key, "") == value, records)
end

function _pqs_source_metadata_acceptance_all_records_false(
    records,
    key::AbstractString,
)
    return all(record -> lowercase(get(record, key, "")) == "false", records)
end

function _be2_pqs_q5_source_metadata_acceptance(
    artifact_dir::AbstractString;
    source_shells_table_path = nothing,
    source_modes_table_path = nothing,
)
    checks = _PQSSourceMetadataAcceptanceChecks()
    _pqs_source_metadata_acceptance_record!(
        checks,
        "source_artifact_dir",
        abspath(artifact_dir),
    )
    if !isdir(artifact_dir)
        _pqs_source_metadata_acceptance_require!(
            checks,
            "source_artifact_available",
            false,
        )
        return (
            rows = copy(checks.rows),
            row_dict = Dict(checks.rows),
            failures = copy(checks.failures),
        )
    end
    _pqs_source_metadata_acceptance_require!(
        checks,
        "source_artifact_available",
        true,
    )

    elapsed = @elapsed begin
        parsed = _be2_pqs_q5_source_metadata_acceptance_metadata(artifact_dir)
        construction =
            _be2_pqs_q5_source_metadata_acceptance_construction(parsed)
        metrics_module = GaussletBases.CartesianContractedParentMetrics
        retained_inventory =
            metrics_module._pqs_current_route_retained_unit_inventory(
                construction,
            )
        source_inventory =
            metrics_module._pqs_current_route_source_shell_mode_inventory(
                retained_inventory;
                provenance = (
                    source =
                        :be2_pqs_q5_source_metadata_real_artifact_acceptance,
                ),
            )
        relation_inventory =
            metrics_module._pqs_current_route_fixed_column_source_relation_inventory(
                retained_inventory;
                provenance = (
                    source =
                        :be2_pqs_q5_source_metadata_real_artifact_acceptance,
                ),
            )
        contract = metrics_module._pqs_source_metadata_export_contract()
        source_shell_header, source_shell_records =
            _pqs_source_metadata_acceptance_table_records(
                metrics_module._write_pqs_source_shells_table,
                source_inventory;
                path = source_shells_table_path,
            )
        source_mode_header, source_mode_records =
            _pqs_source_metadata_acceptance_table_records(
                metrics_module._write_pqs_source_modes_table,
                source_inventory;
                path = source_modes_table_path,
            )

        record = _pqs_source_metadata_acceptance_record!
        require = _pqs_source_metadata_acceptance_require!
        record(checks, "operator_assembly_run", false)
        record(checks, "qw_hamiltonian_construction_run", false)
        record(checks, "full_postprocess_export_run", false)
        require(checks, "fixed_dimension_is_1483", source_inventory.fixed_dimension == 1483)
        require(checks, "relation_fixed_dimension_is_1483", relation_inventory.fixed_dimension == 1483)
        require(checks, "source_shells_table_present", !isempty(source_shell_records))
        require(checks, "source_modes_table_present", !isempty(source_mode_records))
        require(checks, "source_shell_header_matches_contract", source_shell_header == collect(contract.source_shells_header))
        require(checks, "source_mode_header_matches_contract", source_mode_header == collect(contract.source_modes_header))
        require(checks, "source_shell_header_has_no_ray_id", !("ray_id" in source_shell_header))
        require(checks, "source_mode_header_has_no_ray_id", !("ray_id" in source_mode_header))
        require(checks, "source_shell_count_is_8", source_inventory.source_shell_count == 8)
        require(checks, "source_mode_count_is_2299", source_inventory.source_mode_count == 2299)
        require(checks, "source_shell_table_row_count_is_8", length(source_shell_records) == 8)
        require(checks, "source_mode_table_row_count_is_2299", length(source_mode_records) == 2299)
        require(checks, "product_doside_source_shells_is_3", source_inventory.diagnostics.product_doside_source_shell_count == 3)
        require(checks, "product_doside_source_modes_is_531", source_inventory.diagnostics.product_doside_source_mode_count == 531)
        require(checks, "support_dense_source_shells_is_2", source_inventory.diagnostics.support_dense_source_shell_count == 2)
        require(checks, "support_dense_source_modes_is_1458", source_inventory.diagnostics.support_dense_source_mode_count == 1458)
        require(checks, "shell_realized_pqs_source_shells_is_3", source_inventory.diagnostics.shell_realized_pqs_source_shell_count == 3)
        require(checks, "shell_realized_pqs_source_modes_is_310", source_inventory.diagnostics.shell_realized_pqs_source_mode_count == 310)
        require(checks, "unavailable_source_units_is_0", source_inventory.diagnostics.total_unavailable_unit_count == 0)
        require(checks, "unavailable_source_columns_is_0", source_inventory.diagnostics.total_unavailable_column_count == 0)
        require(checks, "source_inventory_no_center_inference", !source_inventory.diagnostics.inferred_from_centers)
        require(checks, "source_inventory_no_nearest_grid_inference", !source_inventory.diagnostics.inferred_from_nearest_grid)
        require(checks, "source_inventory_no_support_order_inference", !source_inventory.diagnostics.inferred_from_support_order)
        require(checks, "source_inventory_no_support_index_inference", !source_inventory.diagnostics.inferred_from_support_indices)
        require(checks, "source_inventory_no_raw_to_final_inference", !source_inventory.diagnostics.inferred_from_raw_to_final_support)
        for key in (
            "inferred_from_centers",
            "inferred_from_nearest_grid",
            "inferred_from_support_order",
            "inferred_from_support_indices",
            "inferred_from_raw_to_final_support",
        )
            require(
                checks,
                string("source_mode_table_", key, "_all_false"),
                _pqs_source_metadata_acceptance_all_records_false(
                    source_mode_records,
                    key,
                ),
            )
        end
        require(checks, "source_shell_shell_statuses_unavailable", all(==(:unavailable), source_inventory.source_shells.shell_label_statuses))
        require(checks, "source_shell_ray_statuses_unavailable", all(==(:unavailable), source_inventory.source_shells.ray_label_statuses))
        require(checks, "source_shell_radial_statuses_unavailable", all(==(:unavailable), source_inventory.source_shells.radial_order_statuses))
        require(checks, "source_mode_shell_statuses_unavailable", all(==(:unavailable), source_inventory.source_modes.shell_label_statuses))
        require(checks, "source_mode_ray_statuses_unavailable", all(==(:unavailable), source_inventory.source_modes.ray_label_statuses))
        require(checks, "source_mode_radial_statuses_unavailable", all(==(:unavailable), source_inventory.source_modes.radial_order_statuses))
        require(checks, "source_shell_table_shell_statuses_unavailable", _pqs_source_metadata_acceptance_all_records_equal(source_shell_records, "shell_label_status", "unavailable"))
        require(checks, "source_shell_table_ray_statuses_unavailable", _pqs_source_metadata_acceptance_all_records_equal(source_shell_records, "ray_label_status", "unavailable"))
        require(checks, "source_shell_table_radial_statuses_unavailable", _pqs_source_metadata_acceptance_all_records_equal(source_shell_records, "radial_order_status", "unavailable"))
        require(checks, "source_mode_table_shell_statuses_unavailable", _pqs_source_metadata_acceptance_all_records_equal(source_mode_records, "shell_label_status", "unavailable"))
        require(checks, "source_mode_table_ray_statuses_unavailable", _pqs_source_metadata_acceptance_all_records_equal(source_mode_records, "ray_label_status", "unavailable"))
        require(checks, "source_mode_table_radial_statuses_unavailable", _pqs_source_metadata_acceptance_all_records_equal(source_mode_records, "radial_order_status", "unavailable"))
        require(checks, "product_doside_relations_status", relation_inventory.status == :product_doside_axis_tuple_relations_only)
        require(checks, "product_doside_relation_rows_is_531", relation_inventory.row_count == 531)
        require(checks, "product_doside_relation_kinds_only", all(==(:product_axis_tuple), relation_inventory.relation_kinds))
        require(checks, "product_axis_tuples_not_rays", relation_inventory.absences_by_contract.product_axis_tuples_not_interpreted_as_ray_labels)
        require(checks, "product_axis_tuples_are_not_ray_labels", !relation_inventory.diagnostics.product_axis_tuples_are_ray_labels)
        require(checks, "fixed_column_source_relations_product_doside_only", relation_inventory.diagnostics.product_doside_row_count == 531)
        require(checks, "support_dense_relation_unavailable_columns_is_642", relation_inventory.diagnostics.support_dense_unavailable_column_count == 642)
        require(checks, "shell_realized_relation_unavailable_columns_is_310", relation_inventory.diagnostics.shell_realized_pqs_unavailable_column_count == 310)
        require(checks, "fixed_column_source_relations_no_center_inference", !relation_inventory.diagnostics.inferred_from_centers)
        require(checks, "fixed_column_source_relations_no_nearest_grid_inference", !relation_inventory.diagnostics.inferred_from_nearest_grid)
        require(checks, "fixed_column_source_relations_no_support_order_inference", !relation_inventory.diagnostics.inferred_from_support_order)
        require(checks, "fixed_column_source_relations_no_support_index_inference", !relation_inventory.diagnostics.inferred_from_support_indices)
        require(checks, "fixed_column_source_relations_no_raw_to_final_inference", !relation_inventory.diagnostics.inferred_from_raw_to_final_support)
        require(checks, "no_retained_weight_or_ida_division", !source_inventory.diagnostics.retained_weight_or_ida_division && !relation_inventory.diagnostics.retained_weight_or_ida_division)
        require(checks, "no_route_construction_change", !source_inventory.diagnostics.route_construction_changed && !relation_inventory.diagnostics.route_construction_changed)
        require(checks, "no_packet_adoption", !source_inventory.diagnostics.packet_adoption && !relation_inventory.diagnostics.packet_adoption)
        require(checks, "no_qw_hamiltonian_change", !source_inventory.diagnostics.qwhamiltonian_changed && !relation_inventory.diagnostics.qwhamiltonian_changed)
        require(checks, "no_public_default_consumes", !source_inventory.diagnostics.public_default_consumes && !relation_inventory.diagnostics.public_default_consumes)
    end
    _pqs_source_metadata_acceptance_record!(
        checks,
        "elapsed_seconds",
        elapsed,
    )
    if !isnothing(source_shells_table_path)
        _pqs_source_metadata_acceptance_record!(
            checks,
            "source_shells_table",
            source_shells_table_path,
        )
    end
    if !isnothing(source_modes_table_path)
        _pqs_source_metadata_acceptance_record!(
            checks,
            "source_modes_table",
            source_modes_table_path,
        )
    end
    return (
        rows = copy(checks.rows),
        row_dict = Dict(checks.rows),
        failures = copy(checks.failures),
    )
end
