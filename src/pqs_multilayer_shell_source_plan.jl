# Route-owned source planning for multi-layer PQS shell realizations.

function _pqs_multilayer_explicit_box_layer_specs(
    core_box::NTuple{3,UnitRange{Int}},
    outer_box::NTuple{3,UnitRange{Int}},
)
    layer_count = _pqs_multilayer_box_depth(outer_box, core_box)
    specs = Vector{NamedTuple}(undef, layer_count)
    for layer_index in 1:layer_count
        inner_box =
            layer_index == 1 ?
            core_box :
            _pqs_multilayer_core_box_at_depth(core_box, layer_index - 1)
        specs[layer_index] = (;
            layer_index,
            current_box = _pqs_multilayer_core_box_at_depth(core_box, layer_index),
            inner_box,
            provenance = :explicit_box_bridge,
        )
    end
    return specs
end

function _pqs_multilayer_region_plan_layer_specs(
    region_plan::PQSMultilayerShellRegionPlan,
)
    return [
        (;
            layer_index = layer.layer_index,
            current_box = layer.current_box,
            inner_box = layer.inner_box,
            provenance = :pqs_multilayer_shell_region_plan,
            terminal_region_key = layer.metadata.terminal_region_key,
            lowering_kind = layer.metadata.lowering_kind,
            source_cpb_count = layer.metadata.source_cpb_count,
        )
        for layer in region_plan.shell_layers
    ]
end

function _pqs_multilayer_realize_shell_source_plan(
    bundles::_CartesianNestedAxisBundles3D,
    core_box::NTuple{3,UnitRange{Int}},
    outer_box::NTuple{3,UnitRange{Int}},
    layer_specs::AbstractVector;
    bond_axis::Symbol,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}},
    source_kind::Symbol,
    metadata,
)
    metadata_tuple = NamedTuple(metadata)
    dims = _nested_axis_lengths(bundles)
    for axis in 1:3
        first(outer_box[axis]) >= 1 && last(outer_box[axis]) <= dims[axis] ||
            throw(ArgumentError("multi-layer PQS outer box must lie inside parent dimensions"))
    end
    metrics = _pqs_multilayer_axis_metrics(bundles)
    core_support_indices =
        _nested_box_support_indices(core_box[1], core_box[2], core_box[3], dims)
    core_support_states =
        [_cartesian_unflat_index(index, dims) for index in core_support_indices]

    shell_records = Vector{NamedTuple}(undef, length(layer_specs))
    for (record_index, spec) in pairs(layer_specs)
        inner_box = spec.inner_box
        current_box = spec.current_box
        q_values = length.(inner_box)
        all(q -> q == q_values[1], q_values) ||
            throw(ArgumentError("multi-layer PQS shell layers currently require cubic raw source dimensions"))
        q = q_values[1]
        layer = _nested_projected_q_shell_layer(
            bundles,
            current_box,
            inner_box;
            bond_axis,
            q,
            L = q,
            raw_source_dims = (q, q, q),
            selected_q = q,
            term_coefficients,
        )
        descriptor = _nested_projected_q_shell_staged_unit_descriptor(layer)
        shell_plan =
            CartesianContractedParentMetrics._pqs_shell_realization_plan(
                descriptor,
                metrics,
            )
        shell_final_coefficients =
            shell_plan.shell_projection_matrix * shell_plan.lowdin_cleanup
        shell_records[record_index] = (;
            layer_index = spec.layer_index,
            current_box,
            inner_box,
            raw_source_dims = (q, q, q),
            q,
            descriptor,
            shell_plan,
            shell_support_indices = descriptor.support_indices,
            shell_support_states = descriptor.support_states,
            shell_final_coefficients,
            support_count = descriptor.support_count,
            mode_count = descriptor.mode_count,
            retained_count = descriptor.retained_count,
            isometry_error = shell_plan.isometry_error,
            provenance = get(spec, :provenance, nothing),
        )
    end

    shell_support_indices =
        reduce(vcat, (record.shell_support_indices for record in shell_records); init = Int[])
    shell_support_states =
        reduce(vcat, (record.shell_support_states for record in shell_records); init = Tuple{Int,Int,Int}[])
    shell_final_coefficients =
        _pqs_multilayer_block_concatenate_shell_coefficients(shell_records)
    core_shell_duplicate_count = length(intersect(core_support_indices, shell_support_indices))
    shell_duplicate_count = _pqs_multilayer_duplicate_count(shell_support_indices)
    combined_support_indices = vcat(core_support_indices, shell_support_indices)
    intended_support_indices =
        _nested_box_support_indices(outer_box[1], outer_box[2], outer_box[3], dims)
    combined_unique = sort!(unique(combined_support_indices))
    intended_sorted = sort!(collect(intended_support_indices))
    covers_outer_box = combined_unique == intended_sorted
    blocker =
        shell_duplicate_count != 0 ? :duplicate_shell_support :
        core_shell_duplicate_count != 0 ? :core_shell_support_overlap :
        !covers_outer_box ? :combined_support_does_not_cover_outer_box :
        nothing
    status = isnothing(blocker) ?
        :available_pqs_multilayer_shell_source_plan :
        :blocked_pqs_multilayer_shell_source_plan

    summary = (;
        status,
        blocker,
        layer_count = length(layer_specs),
        core_support_count = length(core_support_indices),
        shell_support_count = length(shell_support_indices),
        shell_final_retained_count = size(shell_final_coefficients, 2),
        combined_support_count = length(combined_support_indices),
        intended_support_count = length(intended_sorted),
        shell_duplicate_count,
        core_shell_duplicate_count,
        covers_outer_box,
        final_basis_helper = :pqs_complete_core_shell_final_basis,
        collapsed_shell_sector = true,
        fixed_block_matrix_authority_used = false,
        explicit_box_bridge = get(metadata_tuple, :explicit_box_bridge, false),
        shellification_backed_geometry =
            get(metadata_tuple, :shellification_backed_geometry, false),
        h1_materialized = false,
        rhf_materialized = false,
    )

    return (;
        object_kind = :pqs_multilayer_shell_source_plan,
        status,
        blocker,
        source_kind,
        bundles,
        metrics,
        core_box,
        outer_box,
        bond_axis,
        layer_count = length(layer_specs),
        shell_records,
        core_support_indices,
        core_support_states,
        shell_support_indices,
        shell_support_states,
        shell_final_coefficients,
        summary,
        metadata = merge(
            metadata_tuple,
            (;
                source = :pqs_multilayer_shell_source_plan,
                collapsed_shell_sector = true,
                old_fixed_block_matrix_authority_used = false,
            ),
        ),
    )
end

"""
    pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)

Compatibility/probe bridge for direct explicit-box multi-layer PQS source
planning. New route work should prefer the shellification/lowering-backed
`pqs_multilayer_shell_source_plan(bundles, region_plan; ...)` entry point so
explicit boxes do not become shellification authority.

This bridge materializes each layer with the existing projected-q-shell
descriptor and shell-realization plan, then collapses the disjoint shell
supports and shell isometries into one shell sector suitable for
`CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis`.

This helper plans shell/source data only. It does not build final-basis overlap,
one-body operators, H1, IDA, RHF, driver wiring, exports, or artifacts.
"""
function pqs_multilayer_shell_source_plan(
    bundles::_CartesianNestedAxisBundles3D,
    core_box::NTuple{3,UnitRange{Int}},
    outer_box::NTuple{3,UnitRange{Int}};
    bond_axis::Symbol = :z,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    metadata = (;),
)
    return _pqs_multilayer_realize_shell_source_plan(
        bundles,
        core_box,
        outer_box,
        _pqs_multilayer_explicit_box_layer_specs(core_box, outer_box);
        bond_axis,
        term_coefficients,
        source_kind = :repeated_one_cell_projected_q_shell_layers,
        metadata = merge(
            NamedTuple(metadata),
            (;
                explicit_box_bridge = true,
                shellification_backed_geometry = false,
            ),
        ),
    )
end

function pqs_multilayer_shell_source_plan(
    bundles::_CartesianNestedAxisBundles3D,
    region_plan::PQSMultilayerShellRegionPlan;
    bond_axis::Symbol = :z,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    metadata = (;),
)
    region_plan.status === :available_pqs_multilayer_shell_region_plan ||
        throw(ArgumentError("PQS multi-layer source planning requires an available shell region plan"))
    source_plan = _pqs_multilayer_realize_shell_source_plan(
        bundles,
        region_plan.core_box,
        region_plan.outer_box,
        _pqs_multilayer_region_plan_layer_specs(region_plan);
        bond_axis,
        term_coefficients,
        source_kind = :shellification_backed_repeated_one_cell_projected_q_shell_layers,
        metadata = merge(
            (;
                shell_region_plan_source = :pqs_multilayer_shell_region_plan,
                shellification_backed_geometry = true,
            ),
            NamedTuple(metadata),
        ),
    )
    return merge(
        source_plan,
        (;
            source_kind = :shellification_backed_repeated_one_cell_projected_q_shell_layers,
            region_plan,
            region_plan_summary = region_plan.summary,
            summary = merge(
                source_plan.summary,
                (;
                    shell_region_plan_source = :pqs_multilayer_shell_region_plan,
                    shellification_backed_geometry = true,
                ),
            ),
        ),
    )
end

"""
    pqs_multilayer_complete_core_shell_final_basis(plan; metadata = (;), ...)

Build the complete core/shell final-basis realization for a route-owned
multi-layer PQS shell source plan. This helper performs only the overlap
assembly needed by `CartesianFinalBasisRealization`; it does not materialize
one-body operators, H1, IDA, density-density, RHF, driver wiring, exports, or
artifacts.
"""
function pqs_multilayer_complete_core_shell_final_basis(
    plan;
    identity_atol::Real = 1.0e-8,
    rank_atol::Real = 1.0e-10,
    metadata = (;),
)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer final basis requires a pqs_multilayer_shell_source_plan"))
    if plan.status !== :available_pqs_multilayer_shell_source_plan
        return _blocked_pqs_multilayer_complete_core_shell_final_basis(
            plan;
            blocker = isnothing(plan.blocker) ?
                :pqs_multilayer_shell_source_plan_not_available :
                plan.blocker,
            metadata,
        )
    end

    metrics = plan.metrics
    core_overlap = _pqs_multilayer_support_product_matrix(
        plan.core_support_states,
        plan.core_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    core_shell_overlap = _pqs_multilayer_support_product_matrix(
        plan.core_support_states,
        plan.shell_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    shell_overlap = _pqs_multilayer_support_product_matrix(
        plan.shell_support_states,
        plan.shell_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    merged_metadata = merge(
        NamedTuple(plan.metadata),
        NamedTuple(metadata),
        (;
            source = :pqs_multilayer_complete_core_shell_final_basis,
            input_source_plan = :pqs_multilayer_shell_source_plan,
            overlap_blocks_built_from_plan_support_states = true,
            h1_materialized = false,
            ida_data_materialized = false,
            density_density_materialized = false,
            rhf_materialized = false,
            driver_route_materialized = false,
            exports_materialized = false,
            artifacts_materialized = false,
        ),
    )
    final_basis =
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(
            core_support_indices = plan.core_support_indices,
            shell_support_indices = plan.shell_support_indices,
            core_overlap = core_overlap,
            core_shell_overlap = core_shell_overlap,
            shell_overlap = shell_overlap,
            shell_final_coefficients = plan.shell_final_coefficients,
            identity_atol = identity_atol,
            rank_atol = rank_atol,
            metadata = merged_metadata,
        )
    return merge(
        final_basis,
        (;
            source_plan_object_kind = plan.object_kind,
            source_plan_status = plan.status,
            source_plan_layer_count = plan.layer_count,
            source_plan_core_support_count = length(plan.core_support_indices),
            source_plan_shell_support_count = length(plan.shell_support_indices),
            source_plan_final_basis_helper =
                :pqs_multilayer_complete_core_shell_final_basis,
            h1_materialized = false,
            ida_data_materialized = false,
            density_density_materialized = false,
            rhf_materialized = false,
            driver_route_materialized = false,
            exports_materialized = false,
            artifacts_materialized = false,
        ),
    )
end

"""
    pqs_multilayer_complete_core_shell_h1_payload(plan; final_basis, ...)

Build the narrow complete core/shell PQS H1 assembly payload for an available
multi-layer source plan and final basis. This helper assembles only kinetic,
separated by-center electron-nuclear one-body matrices, the final H1
Hamiltonian, and the ordinary H1 solve. It does not materialize IDA,
density-density, RHF, GTO, driver wiring, exports, artifacts, or fixture-rule
policy.
"""
function pqs_multilayer_complete_core_shell_h1_payload(
    plan;
    final_basis,
    coulomb_expansion,
    center_records,
    axis_layers = nothing,
    gaussian_factor_terms_by_center = nothing,
    metadata = (;),
)
    _pqs_multilayer_property(plan, :object_kind) ===
        :pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer H1 payload requires a pqs_multilayer_shell_source_plan"))
    plan.status === :available_pqs_multilayer_shell_source_plan ||
        throw(ArgumentError("PQS multi-layer H1 payload requires an available source plan"))
    get(final_basis, :status, nothing) === :available_pqs_complete_core_shell_final_basis ||
        throw(ArgumentError("PQS multi-layer H1 payload requires an available complete core/shell final basis"))

    support_kinetic = pqs_multilayer_support_kinetic_matrix(plan)
    support_nuclear_by_center =
        pqs_multilayer_support_electron_nuclear_by_center_matrices(
            plan;
            coulomb_expansion,
            center_records,
            axis_layers,
            gaussian_factor_terms_by_center,
        )
    final_kinetic =
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix(
            final_basis,
            support_kinetic;
            term = :kinetic,
        )
    raw_base_layer_gaussian_factor_matrices_used =
        isnothing(gaussian_factor_terms_by_center) && !isnothing(axis_layers)
    nuclear_factor_source =
        isnothing(gaussian_factor_terms_by_center) ?
        :centered_axis_layer_gaussian_factor_terms :
        :pgdg_intermediate_gaussian_factor_terms
    final_nuclear_by_center = map(support_nuclear_by_center.records) do record
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix(
            final_basis,
            record.support_operator;
            term = :electron_nuclear_by_center,
            center_record = record,
            metadata = merge(
                record.metadata,
                (;
                    nuclear_factor_source,
                    support_gaussian_factor_terms_source =
                        record.gaussian_factor_terms_source,
                    raw_base_layer_gaussian_factor_matrices_used,
                ),
            ),
        )
    end
    final_hamiltonian =
        CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_electron_hamiltonian(
            final_kinetic,
            final_nuclear_by_center,
        )
    h1 = CartesianFinalBasisRealization.pqs_complete_core_shell_final_h1_solve(
        final_hamiltonian,
    )
    summary = (;
        status = :materialized_pqs_multilayer_complete_core_shell_h1_payload,
        blocker = nothing,
        final_dimension = h1.final_dimension,
        lowest_energy = h1.lowest_energy,
        solve_kind = h1.solve_kind,
        center_count = support_nuclear_by_center.center_count,
        final_basis_status = final_basis.status,
        source_plan_status = plan.status,
        support_kinetic_materialized = true,
        support_nuclear_by_center_materialized = true,
        final_one_body_transfer_materialized = true,
        hamiltonian_data_materialized = true,
        h1_solve_materialized = true,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        gto_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        fixture_rule_policy_materialized = false,
    )
    return (;
        object_kind = :pqs_multilayer_complete_core_shell_h1_payload,
        status = summary.status,
        blocker = nothing,
        final_kinetic,
        final_nuclear_by_center,
        final_hamiltonian,
        h1,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_complete_core_shell_h1_payload,
                support_space_one_body_helpers_scoped_to_h1 = true,
                ida_data_materialized = false,
                density_density_materialized = false,
                rhf_materialized = false,
                gto_materialized = false,
                driver_route_materialized = false,
                exports_materialized = false,
                artifacts_materialized = false,
                fixture_rule_policy_materialized = false,
            ),
        ),
    )
end

function _blocked_pqs_multilayer_complete_core_shell_final_basis(
    plan;
    blocker,
    metadata = (;),
)
    return (;
        object_kind = :pqs_multilayer_complete_core_shell_final_basis,
        status = :blocked_pqs_multilayer_complete_core_shell_final_basis,
        blocker,
        source_plan_object_kind = _pqs_multilayer_property(plan, :object_kind),
        source_plan_status = _pqs_multilayer_property(plan, :status),
        source_plan_blocker = _pqs_multilayer_property(plan, :blocker),
        source_plan_layer_count = _pqs_multilayer_property(plan, :layer_count, 0),
        source_plan_core_support_count =
            hasproperty(plan, :core_support_indices) ?
            length(plan.core_support_indices) : 0,
        source_plan_shell_support_count =
            hasproperty(plan, :shell_support_indices) ?
            length(plan.shell_support_indices) : 0,
        final_basis_materialized = false,
        one_body_operator_materialized = false,
        h1_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        driver_route_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
        metadata = merge(
            hasproperty(plan, :metadata) ? NamedTuple(plan.metadata) : (;),
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_complete_core_shell_final_basis,
                input_source_plan = :pqs_multilayer_shell_source_plan,
                final_basis_helper = :pqs_complete_core_shell_final_basis,
            ),
        ),
    )
end
