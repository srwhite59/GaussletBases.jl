# Private materialized White-Lindsey benchmark-route seed helpers.

function _white_lindsey_low_order_offset_ranges(
    ranges::AbstractVector{<:UnitRange{Int}},
    offset::Int,
)
    return Tuple((first(range) + offset):(last(range) + offset) for range in ranges)
end

function _white_lindsey_low_order_materialized_seed_inventory(
    sequence::_CartesianNestedShellSequence3D,
    fixed_block::_NestedFixedBlock3D,
    structure::OneCenterAtomicNestedStructureDiagnostics,
)
    length(sequence.shell_layers) == 1 || throw(
        ArgumentError("White-Lindsey materialized seed inventory expects exactly one shell layer"),
    )
    shell = only(sequence.shell_layers)
    shell isa _CartesianNestedCompleteShell3D || throw(
        ArgumentError("White-Lindsey materialized seed inventory requires a complete-shell layer"),
    )
    fixed_block.shell === sequence || throw(
        ArgumentError("White-Lindsey materialized seed inventory requires fixed_block.shell === sequence"),
    )

    shell_range = only(sequence.layer_column_ranges)
    shell_offset = first(shell_range) - 1
    face_ranges = _white_lindsey_low_order_offset_ranges(shell.face_column_ranges, shell_offset)
    edge_ranges = _white_lindsey_low_order_offset_ranges(shell.edge_column_ranges, shell_offset)
    corner_ranges = _white_lindsey_low_order_offset_ranges(shell.corner_column_ranges, shell_offset)

    core_retained_count = length(sequence.core_column_range)
    face_retained_count = sum(length, shell.face_column_ranges)
    edge_retained_count = sum(length, shell.edge_column_ranges)
    corner_retained_count = sum(length, shell.corner_column_ranges)
    shell_retained_count = length(shell_range)
    total_retained_count = size(sequence.coefficient_matrix, 2)
    fixed_block_ready =
        size(fixed_block.coefficient_matrix, 2) == total_retained_count &&
        size(fixed_block.overlap) == (total_retained_count, total_retained_count)
    overlap_ready = fixed_block_ready && all(isfinite, fixed_block.overlap)
    retained_basis_integral_weights_ready =
        length(fixed_block.weights) == total_retained_count &&
        all(isfinite, fixed_block.weights) &&
        all(>(0.0), fixed_block.weights)

    return (;
        object_kind = :white_lindsey_low_order_materialized_seed_inventory,
        route_family = :white_lindsey_low_order,
        status = :private_development_seed,
        private_development_only = true,
        parent_side_count = structure.parent_side_count,
        source_side_count = structure.working_box_side_count,
        nside = structure.nside,
        piece_counts = (;
            core = 1,
            faces = length(shell.faces),
            edges = length(shell.edges),
            corners = length(shell.corners),
        ),
        support_counts = (;
            core = length(sequence.core_indices),
            shell = length(shell.support_indices),
            total_source = length(sequence.support_indices),
        ),
        retained_counts = (;
            core = core_retained_count,
            faces = face_retained_count,
            edges = edge_retained_count,
            corners = corner_retained_count,
            shell = shell_retained_count,
            total = total_retained_count,
        ),
        retained_ranges = (;
            core = sequence.core_column_range,
            shell = shell_range,
            faces = face_ranges,
            edges = edge_ranges,
            corners = corner_ranges,
        ),
        materialized_shell_local_ranges = (;
            faces = Tuple(shell.face_column_ranges),
            edges = Tuple(shell.edge_column_ranges),
            corners = Tuple(shell.corner_column_ranges),
        ),
        fixed_block_ready,
        overlap_ready,
        retained_basis_integral_weights_ready,
        weight_semantics = :retained_basis_integral_weights,
    )
end
