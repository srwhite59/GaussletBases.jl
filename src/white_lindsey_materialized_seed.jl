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
    ;
    packet_kernel::Union{Nothing,Symbol} = nothing,
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
        packet_kernel,
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

function _white_lindsey_low_order_materialized_seed_fixture(;
    parent_side_count::Int = 7,
    nside::Int = 5,
    Z::Real = 2.0,
    d::Real = 0.2,
    tail_spacing::Real = 10.0,
    basis_family::Symbol = :G10,
    reference_spacing::Real = 1.0,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    packet_kernel::Symbol = :factorized_direct,
)
    basis = build_basis(
        MappedUniformBasisSpec(
            basis_family;
            count = parent_side_count,
            mapping = white_lindsey_atomic_mapping(; Z = Z, d = d, tail_spacing = tail_spacing),
            reference_spacing = reference_spacing,
        ),
    )
    sequence = build_one_center_atomic_full_parent_shell_sequence(
        basis;
        expansion = expansion,
        nside = nside,
        gausslet_backend = gausslet_backend,
        refinement_levels = refinement_levels,
        packet_kernel = packet_kernel,
    )
    fixed_block = _nested_fixed_block(sequence, basis, gausslet_backend)
    structure = one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = parent_side_count,
        nside = nside,
    )
    inventory = _white_lindsey_low_order_materialized_seed_inventory(
        sequence,
        fixed_block,
        structure,
        packet_kernel = packet_kernel,
    )
    return (; parent_side_count, nside, packet_kernel, basis, sequence, fixed_block, structure, inventory)
end

function _white_lindsey_low_order_materialized_seed_packet_kernel(seed, inventory = nothing)
    hasproperty(seed, :packet_kernel) && return seed.packet_kernel
    if !isnothing(inventory) && hasproperty(inventory, :packet_kernel)
        return inventory.packet_kernel
    end
    return nothing
end

function _white_lindsey_low_order_grouped_range(ranges)
    isempty(ranges) && return nothing
    return first(first(ranges)):last(last(ranges))
end

function _white_lindsey_low_order_materialized_seed_route_units(seed)
    inventory = hasproperty(seed, :inventory) ? seed.inventory : seed
    packet_kernel = _white_lindsey_low_order_materialized_seed_packet_kernel(seed, inventory)
    source_side_count = inventory.source_side_count
    retained_ranges = inventory.retained_ranges
    inner_side_count = source_side_count - 2
    source_box = ntuple(_ -> 1:source_side_count, 3)
    core_box = ntuple(_ -> 2:(source_side_count - 1), 3)

    retained_units = (
        _pqs_source_box_route_driver_unit_record(
            unit_key = :low_order_core_direct,
            unit_role = :direct_core,
            retained_unit_kind = :white_lindsey_direct_core,
            source_family = :white_lindsey_low_order_materialized_core,
            source_box = core_box,
            source_dimensions = (inner_side_count, inner_side_count, inner_side_count),
            retained_rule_kind = :direct_parent_sites,
            retained_rule_derivation = :materialized_complete_shell_direct_core,
            retained_range = retained_ranges.core,
            retained_count = inventory.retained_counts.core,
            provenance_label = :white_lindsey_low_order_direct_core,
            weight_semantics = :retained_basis_integral_weights,
        ),
        _pqs_source_box_route_driver_unit_record(
            unit_key = :low_order_face_interiors,
            unit_role = :face_interiors,
            retained_unit_kind = :white_lindsey_face_interior_2d_products,
            source_family = :white_lindsey_low_order_materialized_faces,
            source_box = source_box,
            source_dimensions = (6, inner_side_count, inner_side_count),
            retained_rule_kind = :face_interior_2d_products_of_1d_retained_side_functions,
            retained_rule_derivation = :materialized_complete_shell_face_column_ranges,
            retained_range = _white_lindsey_low_order_grouped_range(retained_ranges.faces),
            retained_count = inventory.retained_counts.faces,
            provenance_label = :white_lindsey_low_order_face_interiors,
            weight_semantics = :retained_basis_integral_weights,
        ),
        _pqs_source_box_route_driver_unit_record(
            unit_key = :low_order_edges,
            unit_role = :edges,
            retained_unit_kind = :white_lindsey_edge_1d_side_functions,
            source_family = :white_lindsey_low_order_materialized_edges,
            source_box = source_box,
            source_dimensions = (12, inner_side_count),
            retained_rule_kind = :edge_1d_retained_side_functions,
            retained_rule_derivation = :materialized_complete_shell_edge_column_ranges,
            retained_range = _white_lindsey_low_order_grouped_range(retained_ranges.edges),
            retained_count = inventory.retained_counts.edges,
            provenance_label = :white_lindsey_low_order_edges,
            weight_semantics = :retained_basis_integral_weights,
        ),
        _pqs_source_box_route_driver_unit_record(
            unit_key = :low_order_corners,
            unit_role = :corners,
            retained_unit_kind = :white_lindsey_corner_direct_single_sites,
            source_family = :white_lindsey_low_order_materialized_corners,
            source_box = source_box,
            source_dimensions = (8,),
            retained_rule_kind = :corner_direct_single_site_pieces,
            retained_rule_derivation = :materialized_complete_shell_corner_column_ranges,
            retained_range = _white_lindsey_low_order_grouped_range(retained_ranges.corners),
            retained_count = inventory.retained_counts.corners,
            provenance_label = :white_lindsey_low_order_corners,
            weight_semantics = :retained_basis_integral_weights,
        ),
    )
    unit_inventory = _pqs_source_box_route_driver_unit_inventory(retained_units)
    pair_entries = ()
    pair_family_counts = (white_lindsey_low_order = 0,)
    route_facts = (;
        source_dimensions = unit_inventory.source_dimensions,
        retained_units,
        retained_counts = unit_inventory.retained_counts,
        ranges = unit_inventory.ranges,
        retained_dimension = unit_inventory.retained_dimension,
        pair_entries,
        pair_family_counts,
    )
    return (;
        object_kind = :white_lindsey_low_order_materialized_seed_route_units,
        route_family = :white_lindsey_low_order,
        status = :private_development_seed,
        packet_kernel,
        retained_units,
        unit_inventory,
        standard_unit_inventory =
            _pqs_source_box_route_driver_standard_unit_inventory_summary(route_facts),
        retained_dimension = unit_inventory.retained_dimension,
        pair_entries,
        pair_family_counts,
        operator_pairs_materialized = false,
        weight_semantics = :retained_basis_integral_weights,
    )
end

function _white_lindsey_low_order_materialized_seed_operator_matrices(fixed_block)
    return (;
        overlap = fixed_block.overlap,
        position_x = fixed_block.position_x,
        position_y = fixed_block.position_y,
        position_z = fixed_block.position_z,
        x2_x = fixed_block.x2_x,
        x2_y = fixed_block.x2_y,
        x2_z = fixed_block.x2_z,
        kinetic = fixed_block.kinetic,
    )
end

function _white_lindsey_low_order_operator_matrix_sizes(matrices)
    return NamedTuple{keys(matrices)}(Tuple(size(matrix) for matrix in values(matrices)))
end

function _white_lindsey_low_order_operator_finite_ready(matrices)
    return NamedTuple{keys(matrices)}(Tuple(all(isfinite, matrix) for matrix in values(matrices)))
end

function _white_lindsey_low_order_operator_symmetry_errors(matrices)
    return NamedTuple{keys(matrices)}(
        Tuple(norm(matrix - transpose(matrix), Inf) for matrix in values(matrices)),
    )
end

function _white_lindsey_low_order_operator_symmetric_ready(symmetry_errors; atol::Float64)
    return NamedTuple{keys(symmetry_errors)}(
        Tuple(error <= atol for error in values(symmetry_errors)),
    )
end

function _white_lindsey_low_order_materialized_seed_operator_inventory(seed)
    fixed_block = hasproperty(seed, :fixed_block) ? seed.fixed_block : seed
    packet_kernel = _white_lindsey_low_order_materialized_seed_packet_kernel(seed)
    matrices = _white_lindsey_low_order_materialized_seed_operator_matrices(fixed_block)
    terms = keys(matrices)
    retained_dimension = size(fixed_block.overlap, 1)
    matrix_sizes = _white_lindsey_low_order_operator_matrix_sizes(matrices)
    finite_ready = _white_lindsey_low_order_operator_finite_ready(matrices)
    symmetry_errors = _white_lindsey_low_order_operator_symmetry_errors(matrices)
    symmetry_tolerance = 1.0e-10
    symmetric_ready = _white_lindsey_low_order_operator_symmetric_ready(
        symmetry_errors;
        atol = symmetry_tolerance,
    )
    overlap_identity_error = norm(
        fixed_block.overlap - Matrix{Float64}(I, retained_dimension, retained_dimension),
        Inf,
    )

    return (;
        object_kind = :white_lindsey_low_order_materialized_seed_operator_inventory,
        route_family = :white_lindsey_low_order,
        status = :private_development_seed,
        packet_kernel,
        operator_source = :nested_fixed_block,
        retained_dimension,
        terms,
        matrix_sizes,
        finite_ready,
        all_finite = all(values(finite_ready)),
        symmetry_errors,
        symmetric_ready,
        symmetry_tolerance,
        overlap_identity_error,
        overlap_identity_ready = overlap_identity_error <= symmetry_tolerance,
        operator_pairs_materialized = false,
        electron_electron_materialized = false,
    )
end
