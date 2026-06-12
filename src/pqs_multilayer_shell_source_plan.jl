# Route-owned source planning for multi-layer PQS shell realizations.

function _pqs_multilayer_axis_metrics(bundles::_CartesianNestedAxisBundles3D)
    pgdg_x = _nested_axis_pgdg(bundles, :x)
    pgdg_y = _nested_axis_pgdg(bundles, :y)
    pgdg_z = _nested_axis_pgdg(bundles, :z)
    return (;
        x = (overlap = pgdg_x.overlap, kinetic = pgdg_x.kinetic),
        y = (overlap = pgdg_y.overlap, kinetic = pgdg_y.kinetic),
        z = (overlap = pgdg_z.overlap, kinetic = pgdg_z.kinetic),
    )
end

function _pqs_multilayer_box_depth(
    outer_box::NTuple{3,UnitRange{Int}},
    core_box::NTuple{3,UnitRange{Int}},
)
    lower_depths = ntuple(axis -> first(core_box[axis]) - first(outer_box[axis]), 3)
    upper_depths = ntuple(axis -> last(outer_box[axis]) - last(core_box[axis]), 3)
    all(depth -> depth >= 1, lower_depths) ||
        throw(ArgumentError("multi-layer PQS shell plan requires the core box to be strictly inside the outer box"))
    lower_depths == upper_depths ||
        throw(ArgumentError("multi-layer PQS shell plan currently requires symmetric shell depth on every axis"))
    all(depth -> depth == lower_depths[1], lower_depths) ||
        throw(ArgumentError("multi-layer PQS shell plan currently requires the same shell depth on every axis"))
    return lower_depths[1]
end

function _pqs_multilayer_core_box_at_depth(
    core_box::NTuple{3,UnitRange{Int}},
    depth::Int,
)
    return ntuple(axis -> (first(core_box[axis]) - depth):(last(core_box[axis]) + depth), 3)
end

function _pqs_multilayer_block_concatenate_shell_coefficients(shell_records)
    total_support = sum(record -> record.support_count, shell_records)
    total_retained = sum(record -> record.retained_count, shell_records)
    coefficients = zeros(Float64, total_support, total_retained)
    support_offset = 0
    retained_offset = 0
    for record in shell_records
        support_range = (support_offset + 1):(support_offset + record.support_count)
        retained_range = (retained_offset + 1):(retained_offset + record.retained_count)
        coefficients[support_range, retained_range] .= record.shell_final_coefficients
        support_offset += record.support_count
        retained_offset += record.retained_count
    end
    return coefficients
end

function _pqs_multilayer_duplicate_count(values)
    return length(values) - length(unique(values))
end

function _pqs_multilayer_property(plan, key::Symbol, default = nothing)
    return hasproperty(plan, key) ? getproperty(plan, key) : default
end

function _pqs_multilayer_support_product_matrix(
    left_states::AbstractVector{<:NTuple{3,Int}},
    right_states::AbstractVector{<:NTuple{3,Int}},
    operator_x::AbstractMatrix{<:Real},
    operator_y::AbstractMatrix{<:Real},
    operator_z::AbstractMatrix{<:Real},
)
    result = Matrix{Float64}(undef, length(left_states), length(right_states))
    @inbounds for (left_index, (ix, iy, iz)) in pairs(left_states)
        for (right_index, (jx, jy, jz)) in pairs(right_states)
            result[left_index, right_index] =
                Float64(operator_x[ix, jx]) *
                Float64(operator_y[iy, jy]) *
                Float64(operator_z[iz, jz])
        end
    end
    return result
end

"""
    pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)

Build a compact route-owned source plan for a direct PQS core plus repeated
one-cell surrounding shell layers. Each layer is materialized with the existing
projected-q-shell descriptor and shell-realization plan, then the disjoint shell
supports and shell isometries are collapsed into one shell sector suitable for
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
    dims = _nested_axis_lengths(bundles)
    for axis in 1:3
        first(outer_box[axis]) >= 1 && last(outer_box[axis]) <= dims[axis] ||
            throw(ArgumentError("multi-layer PQS outer box must lie inside parent dimensions"))
    end
    layer_count = _pqs_multilayer_box_depth(outer_box, core_box)
    metrics = _pqs_multilayer_axis_metrics(bundles)
    core_support_indices =
        _nested_box_support_indices(core_box[1], core_box[2], core_box[3], dims)
    core_support_states =
        [_cartesian_unflat_index(index, dims) for index in core_support_indices]

    shell_records = Vector{NamedTuple}(undef, layer_count)
    for layer_index in 1:layer_count
        inner_box =
            layer_index == 1 ?
            core_box :
            _pqs_multilayer_core_box_at_depth(core_box, layer_index - 1)
        current_box = _pqs_multilayer_core_box_at_depth(core_box, layer_index)
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
        shell_records[layer_index] = (;
            layer_index,
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
        layer_count,
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
        h1_materialized = false,
        rhf_materialized = false,
    )

    return (;
        object_kind = :pqs_multilayer_shell_source_plan,
        status,
        blocker,
        source_kind = :repeated_one_cell_projected_q_shell_layers,
        bundles,
        metrics,
        core_box,
        outer_box,
        bond_axis,
        layer_count,
        shell_records,
        core_support_indices,
        core_support_states,
        shell_support_indices,
        shell_support_states,
        shell_final_coefficients,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source = :pqs_multilayer_shell_source_plan,
                collapsed_shell_sector = true,
                old_fixed_block_matrix_authority_used = false,
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
