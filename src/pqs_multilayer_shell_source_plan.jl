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
