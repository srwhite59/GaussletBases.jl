"""
    _BondAlignedHomonuclearChainSplitCandidate3D

Experimental split candidate record for one homonuclear chain subtree.

This keeps the candidate midpoint policy explicit rather than burying it in the
tree builder.
"""
struct _BondAlignedHomonuclearChainSplitCandidate3D
    split_kind::Symbol
    nucleus_ranges::Vector{UnitRange{Int}}
    midpoint_values::Vector{Float64}
    split_indices::Vector{Int}
    child_boxes::Vector{NTuple{3,UnitRange{Int}}}
    child_physical_widths::Vector{NTuple{3,Float64}}
    child_parallel_counts::Vector{Int}
    child_parallel_to_transverse_ratios::Vector{Float64}
    count_eligible::Bool
    shape_eligible::Bool
    did_split::Bool
end

"""
    _BondAlignedHomonuclearChainNodeGeometry3D

Experimental recursive split-tree geometry record for the first homonuclear
chain nesting pass.
"""
struct _BondAlignedHomonuclearChainNodeGeometry3D
    node_label::String
    parent_box::NTuple{3,UnitRange{Int}}
    working_box::NTuple{3,UnitRange{Int}}
    chain_axis::Symbol
    odd_chain_policy::Symbol
    nside::Int
    min_parallel_to_transverse_ratio::Float64
    nucleus_range::UnitRange{Int}
    chain_coordinates::Vector{Float64}
    shared_shell_count::Int
    shared_shell_dimensions::Vector{Int}
    shared_shell_provenance::Vector{_CartesianNestedShellLayerProvenance3D}
    candidate_splits::Vector{_BondAlignedHomonuclearChainSplitCandidate3D}
    accepted_candidate_index::Union{Nothing,Int}
    local_resolution_warning::Bool
    child_nodes::Vector{_BondAlignedHomonuclearChainNodeGeometry3D}
    subtree_fixed_dimension::Int
end

"""
    _CartesianNestedBondAlignedHomonuclearChainSource3D

Experimental fixed-source container for the first homonuclear chain nested
geometry line.
"""
struct _CartesianNestedBondAlignedHomonuclearChainSource3D{B}
    basis::B
    axis_bundles::_CartesianNestedAxisBundles3D
    shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    root_geometry::_BondAlignedHomonuclearChainNodeGeometry3D
    leaf_sequences::Vector{_CartesianNestedShellSequence3D}
    sequence::_CartesianNestedShellSequence3D
end

"""
    _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D

Experimental planar split candidate record for one square-lattice subtree.

The first pass keeps the candidate family explicit:

- strip splits along `x`
- strip splits along `y`
- binary for even site counts
- centered ternary for odd site counts
"""
struct _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D
    split_family::Symbol
    split_axis::Symbol
    x_coordinate_ranges::Vector{UnitRange{Int}}
    y_coordinate_ranges::Vector{UnitRange{Int}}
    split_values::Vector{Float64}
    split_indices::Vector{Int}
    child_boxes::Vector{NTuple{3,UnitRange{Int}}}
    child_physical_widths::Vector{NTuple{3,Float64}}
    child_planar_counts::Vector{NTuple{2,Int}}
    child_in_plane_aspect_ratios::Vector{Float64}
    count_eligible::Bool
    shape_eligible::Bool
    symmetry_preserving::Bool
    did_split::Bool
end

"""
    _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D

Experimental recursive split-tree geometry record for the first homonuclear
square-lattice planar nesting pass.
"""
struct _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D
    node_label::String
    parent_box::NTuple{3,UnitRange{Int}}
    working_box::NTuple{3,UnitRange{Int}}
    lattice_size::Int
    nside::Int
    min_in_plane_aspect_ratio::Float64
    x_coordinate_range::UnitRange{Int}
    y_coordinate_range::UnitRange{Int}
    x_coordinates::Vector{Float64}
    y_coordinates::Vector{Float64}
    shared_shell_count::Int
    shared_shell_dimensions::Vector{Int}
    shared_shell_provenance::Vector{_CartesianNestedShellLayerProvenance3D}
    candidate_splits::Vector{_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D}
    accepted_candidate_index::Union{Nothing,Int}
    local_resolution_warning::Bool
    child_nodes::Vector{_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D}
    subtree_fixed_dimension::Int
end

"""
    _CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D

Experimental fixed-source container for the first homonuclear square-lattice
planar geometry line.
"""
struct _CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D{B}
    basis::B
    axis_bundles::_CartesianNestedAxisBundles3D
    shell_retention_contract::CartesianNestedCompleteShellRetentionContract
    root_geometry::_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D
    leaf_sequences::Vector{_CartesianNestedShellSequence3D}
    sequence::_CartesianNestedShellSequence3D
end

function Base.show(io::IO, candidate::_BondAlignedHomonuclearChainSplitCandidate3D)
    print(
        io,
        "_BondAlignedHomonuclearChainSplitCandidate3D(kind=:",
        candidate.split_kind,
        ", nchild=",
        length(candidate.child_boxes),
        ", did_split=",
        candidate.did_split,
        ")",
    )
end

function Base.show(io::IO, node::_BondAlignedHomonuclearChainNodeGeometry3D)
    print(
        io,
        "_BondAlignedHomonuclearChainNodeGeometry3D(node=",
        node.node_label,
        ", nuclei=",
        node.nucleus_range,
        ", nchild=",
        length(node.child_nodes),
        ", did_split=",
        !isnothing(node.accepted_candidate_index),
        ", nfixed=",
        node.subtree_fixed_dimension,
        ")",
    )
end

function Base.show(io::IO, source::_CartesianNestedBondAlignedHomonuclearChainSource3D)
    print(
        io,
        "_CartesianNestedBondAlignedHomonuclearChainSource3D(nleaf=",
        length(source.leaf_sequences),
        ", nfixed=",
        size(source.sequence.coefficient_matrix, 2),
        ", nside=",
        source.shell_retention_contract.nside,
        ", did_split=",
        !isnothing(source.root_geometry.accepted_candidate_index),
        ")",
    )
end

function Base.show(io::IO, candidate::_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D)
    print(
        io,
        "_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D(family=:",
        candidate.split_family,
        ", nchild=",
        length(candidate.child_boxes),
        ", did_split=",
        candidate.did_split,
        ")",
    )
end

function Base.show(io::IO, node::_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D)
    print(
        io,
        "_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D(node=",
        node.node_label,
        ", x_range=",
        node.x_coordinate_range,
        ", y_range=",
        node.y_coordinate_range,
        ", nchild=",
        length(node.child_nodes),
        ", did_split=",
        !isnothing(node.accepted_candidate_index),
        ", nfixed=",
        node.subtree_fixed_dimension,
        ")",
    )
end

function Base.show(io::IO, source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D)
    print(
        io,
        "_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D(nleaf=",
        length(source.leaf_sequences),
        ", nfixed=",
        size(source.sequence.coefficient_matrix, 2),
        ", nside=",
        source.shell_retention_contract.nside,
        ", did_split=",
        !isnothing(source.root_geometry.accepted_candidate_index),
        ")",
    )
end

function _nested_source_contract_audit(source::_CartesianNestedBondAlignedHomonuclearChainSource3D)
    return _nested_shell_sequence_contract_audit(source.sequence, _nested_axis_lengths(source.axis_bundles))
end

function _nested_source_contract_audit(source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D)
    return _nested_shell_sequence_contract_audit(source.sequence, _nested_axis_lengths(source.axis_bundles))
end

function _nested_chain_three_child_boxes(
    box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
    left_split_index::Int,
    right_split_index::Int,
)
    axis = _nested_axis_index(chain_axis)
    interval = box[axis]
    first(interval) <= left_split_index < right_split_index < last(interval) || throw(
        ArgumentError(
            "three-child chain box construction requires two strictly interior split indices with left < right",
        ),
    )
    left_axis = first(interval):left_split_index
    middle_axis = (left_split_index + 1):right_split_index
    right_axis = (right_split_index + 1):last(interval)
    left_box =
        axis == 1 ? (left_axis, box[2], box[3]) :
        axis == 2 ? (box[1], left_axis, box[3]) :
        (box[1], box[2], left_axis)
    middle_box =
        axis == 1 ? (middle_axis, box[2], box[3]) :
        axis == 2 ? (box[1], middle_axis, box[3]) :
        (box[1], box[2], middle_axis)
    right_box =
        axis == 1 ? (right_axis, box[2], box[3]) :
        axis == 2 ? (box[1], right_axis, box[3]) :
        (box[1], box[2], right_axis)
    return [left_box, middle_box, right_box]
end

function _nested_odd_chain_policy_thresholds(
    odd_chain_policy::Symbol,
    nside::Int,
    min_parallel_to_transverse_ratio::Float64,
)
    nside >= 1 || throw(ArgumentError("odd-chain nested policy requires nside >= 1"))
    min_parallel_to_transverse_ratio > 0.0 || throw(
        ArgumentError("odd-chain nested policy requires min_parallel_to_transverse_ratio > 0"),
    )
    if odd_chain_policy == :strict_current
        return (
            label = odd_chain_policy,
            outer_parallel_count_min = nside,
            center_parallel_count_min = 1,
            total_parallel_count_min = 2 * nside + 2,
            outer_parallel_to_transverse_ratio_min = min_parallel_to_transverse_ratio,
            center_parallel_to_transverse_ratio_min = min_parallel_to_transverse_ratio,
        )
    elseif odd_chain_policy == :central_ternary_relaxed
        return (
            label = odd_chain_policy,
            outer_parallel_count_min = max(nside - 1, 3),
            center_parallel_count_min = 3,
            total_parallel_count_min = 2 * nside + 1,
            outer_parallel_to_transverse_ratio_min = min_parallel_to_transverse_ratio,
            center_parallel_to_transverse_ratio_min = min(
                min_parallel_to_transverse_ratio,
                0.35,
            ),
        )
    else
        throw(ArgumentError("unknown odd-chain nested policy $(repr(odd_chain_policy)); expected :strict_current or :central_ternary_relaxed"))
    end
end

function _nested_chain_parallel_diagnostics(
    bundles::_CartesianNestedAxisBundles3D,
    child_box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
)
    axis = _nested_axis_index(chain_axis)
    widths = _nested_box_physical_widths(bundles, child_box)
    parallel = widths[axis]
    transverse = maximum(widths[index] for index in 1:3 if index != axis)
    ratio =
        transverse > 0.0 ?
        parallel / transverse :
        (parallel > 0.0 ? Inf : 0.0)
    return (
        widths = widths,
        parallel_count = length(child_box[axis]),
        parallel_to_transverse_ratio = ratio,
    )
end

function _nested_chain_children_are_roughly_cubic(
    bundles::_CartesianNestedAxisBundles3D,
    child_boxes::AbstractVector{<:NTuple{3,UnitRange{Int}}},
    child_nucleus_ranges::AbstractVector{<:UnitRange{Int}},
    chain_axis::Symbol;
    min_parallel_to_transverse_ratio::Float64 = 0.4,
)
    min_parallel_to_transverse_ratio > 0.0 || throw(
        ArgumentError("chain anti-sliver check requires min_parallel_to_transverse_ratio > 0"),
    )
    axis = _nested_axis_index(chain_axis)
    for (child_box, nucleus_range) in zip(child_boxes, child_nucleus_ranges)
        if length(child_box[axis]) == 1 && length(nucleus_range) == 1
            continue
        end
        widths = _nested_box_physical_widths(bundles, child_box)
        parallel = widths[axis]
        transverse = maximum(widths[index] for index in 1:3 if index != axis)
        parallel > 0.0 || return false
        transverse > 0.0 || return false
        parallel >= min_parallel_to_transverse_ratio * transverse || return false
    end
    return true
end

function _nested_chain_binary_candidate(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
    nucleus_range::UnitRange{Int},
    chain_coordinates::AbstractVector{<:Real},
    split_after::Int;
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    prefer_midpoint_tie_side::Symbol = :left,
)
    midpoint = 0.5 * (Float64(chain_coordinates[split_after]) + Float64(chain_coordinates[split_after + 1]))
    axis = _nested_axis_index(chain_axis)
    parallel_centers = _nested_axis_pgdg(bundles, chain_axis).centers
    split_index = _nested_diatomic_split_plane_index(
        parallel_centers,
        working_box[axis],
        midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    left_box, right_box = _nested_diatomic_child_boxes(working_box, chain_axis, split_index)
    child_boxes = [left_box, right_box]
    child_diagnostics = [
        _nested_chain_parallel_diagnostics(bundles, box, chain_axis) for box in child_boxes
    ]
    nucleus_ranges = [first(nucleus_range):(first(nucleus_range) + split_after - 1), (first(nucleus_range) + split_after):last(nucleus_range)]
    count_eligible =
        length(working_box[axis]) > 2 * nside &&
        minimum(length(box[axis]) for box in child_boxes) >= nside
    shape_eligible =
        count_eligible &&
        _nested_chain_children_are_roughly_cubic(
            bundles,
            child_boxes,
            nucleus_ranges,
            chain_axis;
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        )
    return _BondAlignedHomonuclearChainSplitCandidate3D(
        :binary,
        collect(nucleus_ranges),
        [midpoint],
        [split_index],
        child_boxes,
        [entry.widths for entry in child_diagnostics],
        [entry.parallel_count for entry in child_diagnostics],
        [entry.parallel_to_transverse_ratio for entry in child_diagnostics],
        count_eligible,
        shape_eligible,
        count_eligible && shape_eligible,
    )
end

function _nested_chain_ternary_candidate(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
    nucleus_range::UnitRange{Int},
    chain_coordinates::AbstractVector{<:Real},
    center_index::Int;
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    prefer_midpoint_tie_side::Symbol = :left,
)
    left_midpoint = 0.5 * (Float64(chain_coordinates[center_index - 1]) + Float64(chain_coordinates[center_index]))
    right_midpoint = 0.5 * (Float64(chain_coordinates[center_index]) + Float64(chain_coordinates[center_index + 1]))
    axis = _nested_axis_index(chain_axis)
    parallel_centers = _nested_axis_pgdg(bundles, chain_axis).centers
    left_split_index = _nested_diatomic_split_plane_index(
        parallel_centers,
        working_box[axis],
        left_midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    right_split_index = _nested_diatomic_split_plane_index(
        parallel_centers,
        working_box[axis],
        right_midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    child_boxes =
        left_split_index < right_split_index ?
        _nested_chain_three_child_boxes(working_box, chain_axis, left_split_index, right_split_index) :
        NTuple{3,UnitRange{Int}}[]
    child_diagnostics = [
        _nested_chain_parallel_diagnostics(bundles, box, chain_axis) for box in child_boxes
    ]
    nucleus_ranges = [
        first(nucleus_range):(first(nucleus_range) + center_index - 2),
        (first(nucleus_range) + center_index - 1):(first(nucleus_range) + center_index - 1),
        (first(nucleus_range) + center_index):last(nucleus_range),
    ]
    thresholds = _nested_odd_chain_policy_thresholds(
        odd_chain_policy,
        nside,
        min_parallel_to_transverse_ratio,
    )
    child_parallel_counts = [entry.parallel_count for entry in child_diagnostics]
    child_parallel_to_transverse_ratios = [
        entry.parallel_to_transverse_ratio for entry in child_diagnostics
    ]
    count_eligible =
        !isempty(child_boxes) &&
        length(working_box[axis]) >= thresholds.total_parallel_count_min &&
        min(child_parallel_counts[1], child_parallel_counts[end]) >= thresholds.outer_parallel_count_min &&
        child_parallel_counts[2] >= thresholds.center_parallel_count_min
    shape_eligible =
        count_eligible &&
        child_parallel_to_transverse_ratios[1] >= thresholds.outer_parallel_to_transverse_ratio_min &&
        child_parallel_to_transverse_ratios[end] >= thresholds.outer_parallel_to_transverse_ratio_min &&
        child_parallel_to_transverse_ratios[2] >= thresholds.center_parallel_to_transverse_ratio_min
    return _BondAlignedHomonuclearChainSplitCandidate3D(
        :ternary,
        collect(nucleus_ranges),
        [left_midpoint, right_midpoint],
        [left_split_index, right_split_index],
        collect(child_boxes),
        [entry.widths for entry in child_diagnostics],
        child_parallel_counts,
        child_parallel_to_transverse_ratios,
        count_eligible,
        shape_eligible,
        count_eligible && shape_eligible,
    )
end

function _nested_chain_split_candidates(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    chain_axis::Symbol,
    nucleus_range::UnitRange{Int},
    chain_coordinates::AbstractVector{<:Real};
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    prefer_midpoint_tie_side::Symbol = :left,
)
    natoms = length(chain_coordinates)
    natoms <= 1 && return _BondAlignedHomonuclearChainSplitCandidate3D[]
    candidates = _BondAlignedHomonuclearChainSplitCandidate3D[]
    if iseven(natoms)
        for split_after in 1:(natoms - 1)
            push!(
                candidates,
                _nested_chain_binary_candidate(
                    bundles,
                    working_box,
                    chain_axis,
                    nucleus_range,
                    chain_coordinates,
                    split_after;
                    nside = nside,
                    min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
                    prefer_midpoint_tie_side = prefer_midpoint_tie_side,
                ),
            )
        end
    else
        for center_index in 2:(natoms - 1)
            push!(
                candidates,
                _nested_chain_ternary_candidate(
                    bundles,
                    working_box,
                    chain_axis,
                    nucleus_range,
                    chain_coordinates,
                    center_index;
                    nside = nside,
                    min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
                    odd_chain_policy = odd_chain_policy,
                    prefer_midpoint_tie_side = prefer_midpoint_tie_side,
                ),
            )
        end
    end
    return candidates
end

function _nested_chain_candidate_center(candidate::_BondAlignedHomonuclearChainSplitCandidate3D)
    if candidate.split_kind == :binary
        return candidate.midpoint_values[1]
    else
        return sum(candidate.midpoint_values) / length(candidate.midpoint_values)
    end
end

function _nested_choose_chain_candidate(
    candidates::AbstractVector{_BondAlignedHomonuclearChainSplitCandidate3D},
    chain_coordinates::AbstractVector{<:Real};
    atol::Float64 = 1.0e-12,
    rtol::Float64 = 1.0e-10,
)
    accepted = Int[index for index in eachindex(candidates) if candidates[index].did_split]
    isempty(accepted) && return nothing
    target = 0.5 * (Float64(first(chain_coordinates)) + Float64(last(chain_coordinates)))
    distances = Float64[abs(_nested_chain_candidate_center(candidates[index]) - target) for index in accepted]
    minimum_distance = minimum(distances)
    tied = Int[
        accepted[index] for index in eachindex(accepted) if
        isapprox(distances[index], minimum_distance; atol = atol, rtol = rtol)
    ]
    return minimum(tied)
end

function _nested_chain_node_summary(
    node::_BondAlignedHomonuclearChainNodeGeometry3D,
)
    policy_thresholds = _nested_odd_chain_policy_thresholds(
        node.odd_chain_policy,
        node.nside,
        node.min_parallel_to_transverse_ratio,
    )
    candidate_summaries = [
        (
            split_kind = candidate.split_kind,
            nucleus_ranges = candidate.nucleus_ranges,
            midpoint_values = candidate.midpoint_values,
            split_indices = candidate.split_indices,
            child_boxes = candidate.child_boxes,
            child_physical_widths = candidate.child_physical_widths,
            child_parallel_counts = candidate.child_parallel_counts,
            child_parallel_to_transverse_ratios = candidate.child_parallel_to_transverse_ratios,
            count_eligible = candidate.count_eligible,
            shape_eligible = candidate.shape_eligible,
            did_split = candidate.did_split,
            accepted = index == node.accepted_candidate_index,
        ) for (index, candidate) in pairs(node.candidate_splits)
    ]
    return (
        node_label = node.node_label,
        parent_box = node.parent_box,
        working_box = node.working_box,
        chain_axis = node.chain_axis,
        odd_chain_policy = node.odd_chain_policy,
        odd_chain_policy_thresholds = policy_thresholds,
        nucleus_range = node.nucleus_range,
        chain_coordinates = node.chain_coordinates,
        shared_shell_count = node.shared_shell_count,
        shared_shell_dimensions = node.shared_shell_dimensions,
        shared_shell_provenance = node.shared_shell_provenance,
        accepted_candidate_index = node.accepted_candidate_index,
        did_split = !isnothing(node.accepted_candidate_index),
        local_resolution_warning = node.local_resolution_warning,
        child_count = length(node.child_nodes),
        subtree_fixed_dimension = node.subtree_fixed_dimension,
        candidate_summaries = candidate_summaries,
    )
end

function _nested_chain_collect_node_summaries(
    node::_BondAlignedHomonuclearChainNodeGeometry3D,
)
    summaries = Any[_nested_chain_node_summary(node)]
    for child in node.child_nodes
        append!(summaries, _nested_chain_collect_node_summaries(child))
    end
    return summaries
end

function _nested_bond_aligned_homonuclear_chain_node(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    parent_box::NTuple{3,UnitRange{Int}},
    nucleus_range::UnitRange{Int},
    node_label::AbstractString;
    chain_axis::Symbol = :z,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    prefer_midpoint_tie_side::Symbol = :left,
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    retain_x_edge::Int = 3,
    retain_y_edge::Int = 3,
    retain_z_edge::Int = 3,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    current_box = parent_box
    shared_shell_layers = _CartesianNestedCompleteShell3D[]
    shared_shell_dimensions = Int[]
    shared_shell_provenance = _CartesianNestedShellLayerProvenance3D[]
    local_coordinates = basis.chain_coordinates[nucleus_range]

    while true
        candidates = _nested_chain_split_candidates(
            bundles,
            current_box,
            chain_axis,
            nucleus_range,
            local_coordinates;
            nside = nside,
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
            odd_chain_policy = odd_chain_policy,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        )
        accepted_candidate = _nested_choose_chain_candidate(candidates, local_coordinates)
        (accepted_candidate !== nothing || minimum(length.(current_box)) <= nside || !_nested_can_shrink_box(current_box)) && break

        inner_box = _nested_inner_box(current_box)
        shell = _nested_complete_rectangular_shell(
            bundles,
            inner_box...;
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
            x_fixed = (first(current_box[1]), last(current_box[1])),
            y_fixed = (first(current_box[2]), last(current_box[2])),
            z_fixed = (first(current_box[3]), last(current_box[3])),
            term_coefficients = term_coefficients,
        )
        push!(shared_shell_layers, shell)
        push!(shared_shell_dimensions, size(shell.coefficient_matrix, 2))
        push!(shared_shell_provenance, shell.provenance)
        current_box = inner_box
    end

    candidates = _nested_chain_split_candidates(
        bundles,
        current_box,
        chain_axis,
        nucleus_range,
        local_coordinates;
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    accepted_candidate = _nested_choose_chain_candidate(candidates, local_coordinates)

    if isnothing(accepted_candidate)
        leaf_sequence = _nested_complete_shell_sequence_for_box(
            bundles,
            current_box;
            nside = nside,
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
            term_coefficients = term_coefficients,
        )
        merged_sequence =
            isempty(shared_shell_layers) ? leaf_sequence :
            _nested_shell_sequence_from_core_block(
                bundles,
                leaf_sequence.support_indices,
                leaf_sequence.coefficient_matrix,
                shared_shell_layers,
                term_coefficients = term_coefficients,
            )
        node = _BondAlignedHomonuclearChainNodeGeometry3D(
            String(node_label),
            parent_box,
            current_box,
            chain_axis,
            odd_chain_policy,
            nside,
            min_parallel_to_transverse_ratio,
            nucleus_range,
            collect(local_coordinates),
            length(shared_shell_layers),
            shared_shell_dimensions,
            shared_shell_provenance,
            candidates,
            nothing,
            !isempty(candidates),
            _BondAlignedHomonuclearChainNodeGeometry3D[],
            size(merged_sequence.coefficient_matrix, 2),
        )
        return (
            geometry = node,
            sequence = merged_sequence,
            leaf_sequences = _CartesianNestedShellSequence3D[leaf_sequence],
        )
    end

    candidate = candidates[accepted_candidate]
    child_results = [
        _nested_bond_aligned_homonuclear_chain_node(
            basis,
            bundles,
            candidate.child_boxes[index],
            candidate.nucleus_ranges[index],
            string(node_label, "_", index);
            chain_axis = chain_axis,
            nside = nside,
            min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
            odd_chain_policy = odd_chain_policy,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
            term_coefficients = term_coefficients,
        ) for index in eachindex(candidate.child_boxes)
    ]
    child_nodes = [result.geometry for result in child_results]
    child_leaf_sequences = _CartesianNestedShellSequence3D[]
    for result in child_results
        append!(child_leaf_sequences, result.leaf_sequences)
    end
    core_support = reduce(vcat, [result.sequence.support_indices for result in child_results])
    core_coefficients = _nested_hcat_coefficient_maps(
        [result.sequence.coefficient_matrix for result in child_results],
    )
    merged_sequence =
        isempty(shared_shell_layers) ? _nested_shell_sequence_from_core_block(
            bundles,
            core_support,
            core_coefficients,
            _AbstractCartesianNestedShellLayer3D[];
            enforce_coverage = false,
            term_coefficients = term_coefficients,
        ) :
        _nested_shell_sequence_from_core_block(
            bundles,
            core_support,
            core_coefficients,
            shared_shell_layers,
            term_coefficients = term_coefficients,
        )
    node = _BondAlignedHomonuclearChainNodeGeometry3D(
        String(node_label),
        parent_box,
        current_box,
        chain_axis,
        odd_chain_policy,
        nside,
        min_parallel_to_transverse_ratio,
        nucleus_range,
        collect(local_coordinates),
        length(shared_shell_layers),
        shared_shell_dimensions,
        shared_shell_provenance,
        candidates,
        accepted_candidate,
        false,
        child_nodes,
        size(merged_sequence.coefficient_matrix, 2),
    )
    return (
        geometry = node,
        sequence = merged_sequence,
        leaf_sequences = child_leaf_sequences,
    )
end

function _nested_bond_aligned_homonuclear_chain_source(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    chain_axis::Symbol = :z,
    nside::Int = 5,
    min_parallel_to_transverse_ratio::Float64 = 0.4,
    odd_chain_policy::Symbol = :strict_current,
    prefer_midpoint_tie_side::Symbol = :left,
    retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_x_edge::Union{Nothing,Int} = nothing,
    retain_y_edge::Union{Nothing,Int} = nothing,
    retain_z_edge::Union{Nothing,Int} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    retention = _nested_resolve_complete_shell_retention(
        nside;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        retain_x_edge = retain_x_edge,
        retain_y_edge = retain_y_edge,
        retain_z_edge = retain_z_edge,
    )
    dims = _nested_axis_lengths(bundles)
    parent_box = (1:dims[1], 1:dims[2], 1:dims[3])
    result = _nested_bond_aligned_homonuclear_chain_node(
        basis,
        bundles,
        parent_box,
        1:length(basis.chain_coordinates),
        "root";
        chain_axis = chain_axis,
        nside = nside,
        min_parallel_to_transverse_ratio = min_parallel_to_transverse_ratio,
        odd_chain_policy = odd_chain_policy,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        retain_xy = retention.retain_xy,
        retain_xz = retention.retain_xz,
        retain_yz = retention.retain_yz,
        retain_x_edge = retention.retain_x_edge,
        retain_y_edge = retention.retain_y_edge,
        retain_z_edge = retention.retain_z_edge,
        term_coefficients = term_coefficients,
    )
    return _CartesianNestedBondAlignedHomonuclearChainSource3D(
        basis,
        bundles,
        retention,
        result.geometry,
        result.leaf_sequences,
        result.sequence,
    )
end

function _nested_square_lattice_planar_diagnostics(
    bundles::_CartesianNestedAxisBundles3D,
    child_box::NTuple{3,UnitRange{Int}},
)
    widths = _nested_box_physical_widths(bundles, child_box)
    planar_x = widths[1]
    planar_y = widths[2]
    longer = max(planar_x, planar_y)
    aspect =
        longer > 0.0 ?
        min(planar_x, planar_y) / longer :
        0.0
    return (
        widths = widths,
        planar_count = (length(child_box[1]), length(child_box[2])),
        in_plane_aspect_ratio = aspect,
    )
end

function _nested_square_lattice_symmetry_preserving(
    child_planar_counts::AbstractVector{<:NTuple{2,Int}},
)
    if length(child_planar_counts) == 2
        return (
            abs(child_planar_counts[1][1] - child_planar_counts[2][1]) <= 1 &&
            abs(child_planar_counts[1][2] - child_planar_counts[2][2]) <= 1
        )
    elseif length(child_planar_counts) == 3
        return (
            abs(child_planar_counts[1][1] - child_planar_counts[3][1]) <= 1 &&
            abs(child_planar_counts[1][2] - child_planar_counts[3][2]) <= 1
        )
    else
        return false
    end
end

function _nested_square_lattice_binary_candidate(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    split_axis::Symbol,
    x_coordinate_range::UnitRange{Int},
    y_coordinate_range::UnitRange{Int},
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
)
    axis = _nested_axis_index(split_axis)
    coordinate_range = split_axis == :x ? x_coordinate_range : y_coordinate_range
    coordinates = split_axis == :x ? x_coordinates : y_coordinates
    local_count = length(coordinate_range)
    local_count >= 2 || throw(ArgumentError("square-lattice binary split requires at least two coordinates on the chosen axis"))
    split_after = local_count ÷ 2
    midpoint = 0.5 * (
        Float64(coordinates[first(coordinate_range) + split_after - 1]) +
        Float64(coordinates[first(coordinate_range) + split_after])
    )
    split_index = _nested_diatomic_split_plane_index(
        _nested_axis_pgdg(bundles, split_axis).centers,
        working_box[axis],
        midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    child_boxes = collect(_nested_diatomic_child_boxes(working_box, split_axis, split_index))
    split_axis_ranges = [
        first(coordinate_range):(first(coordinate_range) + split_after - 1),
        (first(coordinate_range) + split_after):last(coordinate_range),
    ]
    x_ranges = split_axis == :x ? split_axis_ranges : [x_coordinate_range, x_coordinate_range]
    y_ranges = split_axis == :y ? split_axis_ranges : [y_coordinate_range, y_coordinate_range]
    child_diagnostics = [
        _nested_square_lattice_planar_diagnostics(bundles, box) for box in child_boxes
    ]
    child_planar_counts = [entry.planar_count for entry in child_diagnostics]
    count_eligible =
        length(working_box[axis]) > 2 * nside &&
        minimum(length(box[axis]) for box in child_boxes) >= nside
    shape_eligible =
        count_eligible &&
        all(
            entry.in_plane_aspect_ratio >= min_in_plane_aspect_ratio for
            entry in child_diagnostics
        )
    split_family = split_axis == :x ? :split_x_binary : :split_y_binary
    return _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D(
        split_family,
        split_axis,
        x_ranges,
        y_ranges,
        [midpoint],
        [split_index],
        child_boxes,
        [entry.widths for entry in child_diagnostics],
        child_planar_counts,
        [entry.in_plane_aspect_ratio for entry in child_diagnostics],
        count_eligible,
        shape_eligible,
        _nested_square_lattice_symmetry_preserving(child_planar_counts),
        count_eligible && shape_eligible,
    )
end

function _nested_square_lattice_ternary_candidate(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    split_axis::Symbol,
    x_coordinate_range::UnitRange{Int},
    y_coordinate_range::UnitRange{Int},
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
)
    axis = _nested_axis_index(split_axis)
    coordinate_range = split_axis == :x ? x_coordinate_range : y_coordinate_range
    coordinates = split_axis == :x ? x_coordinates : y_coordinates
    local_count = length(coordinate_range)
    isodd(local_count) && local_count >= 3 || throw(
        ArgumentError("square-lattice ternary split requires an odd coordinate count of at least three on the chosen axis"),
    )
    center_offset = (local_count + 1) ÷ 2
    center_index = first(coordinate_range) + center_offset - 1
    left_midpoint = 0.5 * (
        Float64(coordinates[center_index - 1]) + Float64(coordinates[center_index])
    )
    right_midpoint = 0.5 * (
        Float64(coordinates[center_index]) + Float64(coordinates[center_index + 1])
    )
    left_split_index = _nested_diatomic_split_plane_index(
        _nested_axis_pgdg(bundles, split_axis).centers,
        working_box[axis],
        left_midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    right_split_index = _nested_diatomic_split_plane_index(
        _nested_axis_pgdg(bundles, split_axis).centers,
        working_box[axis],
        right_midpoint;
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    child_boxes =
        left_split_index < right_split_index ?
        _nested_chain_three_child_boxes(
            working_box,
            split_axis,
            left_split_index,
            right_split_index,
        ) :
        NTuple{3,UnitRange{Int}}[]
    split_axis_ranges = [
        first(coordinate_range):(center_index - 1),
        center_index:center_index,
        (center_index + 1):last(coordinate_range),
    ]
    x_ranges = split_axis == :x ? split_axis_ranges : [x_coordinate_range, x_coordinate_range, x_coordinate_range]
    y_ranges = split_axis == :y ? split_axis_ranges : [y_coordinate_range, y_coordinate_range, y_coordinate_range]
    child_diagnostics = [
        _nested_square_lattice_planar_diagnostics(bundles, box) for box in child_boxes
    ]
    child_planar_counts = [entry.planar_count for entry in child_diagnostics]
    outer_count_min = max(nside - 1, 3)
    center_count_min = 3
    count_eligible =
        !isempty(child_boxes) &&
        length(working_box[axis]) >= 2 * nside + 1 &&
        min(length(child_boxes[1][axis]), length(child_boxes[end][axis])) >= outer_count_min &&
        length(child_boxes[2][axis]) >= center_count_min
    shape_eligible =
        count_eligible &&
        all(
            entry.in_plane_aspect_ratio >= min_in_plane_aspect_ratio for
            entry in child_diagnostics
        )
    split_family = split_axis == :x ? :split_x_ternary : :split_y_ternary
    return _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D(
        split_family,
        split_axis,
        x_ranges,
        y_ranges,
        [left_midpoint, right_midpoint],
        [left_split_index, right_split_index],
        collect(child_boxes),
        [entry.widths for entry in child_diagnostics],
        child_planar_counts,
        [entry.in_plane_aspect_ratio for entry in child_diagnostics],
        count_eligible,
        shape_eligible,
        _nested_square_lattice_symmetry_preserving(child_planar_counts),
        count_eligible && shape_eligible,
    )
end

function _nested_square_lattice_split_candidates(
    bundles::_CartesianNestedAxisBundles3D,
    working_box::NTuple{3,UnitRange{Int}},
    x_coordinate_range::UnitRange{Int},
    y_coordinate_range::UnitRange{Int},
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
)
    candidates = _AxisAlignedHomonuclearSquareLatticeSplitCandidate3D[]
    nx_sites = length(x_coordinate_range)
    ny_sites = length(y_coordinate_range)

    if nx_sites > 1
        push!(
            candidates,
            iseven(nx_sites) ?
            _nested_square_lattice_binary_candidate(
                bundles,
                working_box,
                :x,
                x_coordinate_range,
                y_coordinate_range,
                x_coordinates,
                y_coordinates;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            ) :
            _nested_square_lattice_ternary_candidate(
                bundles,
                working_box,
                :x,
                x_coordinate_range,
                y_coordinate_range,
                x_coordinates,
                y_coordinates;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            ),
        )
    end

    if ny_sites > 1
        push!(
            candidates,
            iseven(ny_sites) ?
            _nested_square_lattice_binary_candidate(
                bundles,
                working_box,
                :y,
                x_coordinate_range,
                y_coordinate_range,
                x_coordinates,
                y_coordinates;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            ) :
            _nested_square_lattice_ternary_candidate(
                bundles,
                working_box,
                :y,
                x_coordinate_range,
                y_coordinate_range,
                x_coordinates,
                y_coordinates;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
                prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            ),
        )
    end

    return candidates
end

function _nested_square_lattice_candidate_center(
    candidate::_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D,
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real},
)
    x_mid = 0.5 * (Float64(first(x_coordinates)) + Float64(last(x_coordinates)))
    y_mid = 0.5 * (Float64(first(y_coordinates)) + Float64(last(y_coordinates)))
    split_center = sum(candidate.split_values) / length(candidate.split_values)
    return candidate.split_axis == :x ? (split_center, y_mid) : (x_mid, split_center)
end

function _nested_choose_square_lattice_candidate(
    candidates::AbstractVector{_AxisAlignedHomonuclearSquareLatticeSplitCandidate3D},
    x_coordinates::AbstractVector{<:Real},
    y_coordinates::AbstractVector{<:Real};
    atol::Float64 = 1.0e-12,
    rtol::Float64 = 1.0e-10,
)
    accepted = Int[index for index in eachindex(candidates) if candidates[index].did_split]
    isempty(accepted) && return nothing
    target = (
        0.5 * (Float64(first(x_coordinates)) + Float64(last(x_coordinates))),
        0.5 * (Float64(first(y_coordinates)) + Float64(last(y_coordinates))),
    )
    distances = Float64[
        hypot(
            _nested_square_lattice_candidate_center(candidates[index], x_coordinates, y_coordinates)[1] - target[1],
            _nested_square_lattice_candidate_center(candidates[index], x_coordinates, y_coordinates)[2] - target[2],
        ) for index in accepted
    ]
    minimum_distance = minimum(distances)
    tied = Int[
        accepted[index] for index in eachindex(accepted) if
        isapprox(distances[index], minimum_distance; atol = atol, rtol = rtol)
    ]
    return minimum(tied)
end

function _nested_square_lattice_node_summary(
    node::_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D,
)
    candidate_summaries = [
        (
            split_family = candidate.split_family,
            split_axis = candidate.split_axis,
            x_coordinate_ranges = candidate.x_coordinate_ranges,
            y_coordinate_ranges = candidate.y_coordinate_ranges,
            split_values = candidate.split_values,
            split_indices = candidate.split_indices,
            child_boxes = candidate.child_boxes,
            child_physical_widths = candidate.child_physical_widths,
            child_planar_counts = candidate.child_planar_counts,
            child_in_plane_aspect_ratios = candidate.child_in_plane_aspect_ratios,
            count_eligible = candidate.count_eligible,
            shape_eligible = candidate.shape_eligible,
            symmetry_preserving = candidate.symmetry_preserving,
            did_split = candidate.did_split,
            accepted = index == node.accepted_candidate_index,
        ) for (index, candidate) in pairs(node.candidate_splits)
    ]
    return (
        node_label = node.node_label,
        parent_box = node.parent_box,
        working_box = node.working_box,
        lattice_size = node.lattice_size,
        min_in_plane_aspect_ratio = node.min_in_plane_aspect_ratio,
        x_coordinate_range = node.x_coordinate_range,
        y_coordinate_range = node.y_coordinate_range,
        x_coordinates = node.x_coordinates,
        y_coordinates = node.y_coordinates,
        shared_shell_count = node.shared_shell_count,
        shared_shell_dimensions = node.shared_shell_dimensions,
        shared_shell_provenance = node.shared_shell_provenance,
        accepted_candidate_index = node.accepted_candidate_index,
        did_split = !isnothing(node.accepted_candidate_index),
        local_resolution_warning = node.local_resolution_warning,
        child_count = length(node.child_nodes),
        subtree_fixed_dimension = node.subtree_fixed_dimension,
        candidate_summaries = candidate_summaries,
    )
end

function _nested_square_lattice_collect_node_summaries(
    node::_AxisAlignedHomonuclearSquareLatticeNodeGeometry3D,
)
    summaries = Any[_nested_square_lattice_node_summary(node)]
    for child in node.child_nodes
        append!(summaries, _nested_square_lattice_collect_node_summaries(child))
    end
    return summaries
end

function _nested_axis_aligned_homonuclear_square_lattice_node(
    basis,
    bundles::_CartesianNestedAxisBundles3D,
    parent_box::NTuple{3,UnitRange{Int}},
    x_coordinate_range::UnitRange{Int},
    y_coordinate_range::UnitRange{Int},
    node_label::AbstractString;
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
    retain_xy::Tuple{Int,Int} = (4, 3),
    retain_xz::Tuple{Int,Int} = (4, 3),
    retain_yz::Tuple{Int,Int} = (4, 3),
    retain_x_edge::Int = 3,
    retain_y_edge::Int = 3,
    retain_z_edge::Int = 3,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    current_box = parent_box
    shared_shell_layers = _CartesianNestedCompleteShell3D[]
    shared_shell_dimensions = Int[]
    shared_shell_provenance = _CartesianNestedShellLayerProvenance3D[]
    local_x_coordinates = basis.x_coordinates[x_coordinate_range]
    local_y_coordinates = basis.y_coordinates[y_coordinate_range]

    while true
        candidates = _nested_square_lattice_split_candidates(
            bundles,
            current_box,
            x_coordinate_range,
            y_coordinate_range,
            basis.x_coordinates,
            basis.y_coordinates;
            nside = nside,
            min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        )
        accepted_candidate = _nested_choose_square_lattice_candidate(
            candidates,
            local_x_coordinates,
            local_y_coordinates,
        )
        (accepted_candidate !== nothing || minimum(length.(current_box)) <= nside || !_nested_can_shrink_box(current_box)) && break

        inner_box = _nested_inner_box(current_box)
        shell = _nested_complete_rectangular_shell(
            bundles,
            inner_box...;
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
            x_fixed = (first(current_box[1]), last(current_box[1])),
            y_fixed = (first(current_box[2]), last(current_box[2])),
            z_fixed = (first(current_box[3]), last(current_box[3])),
            term_coefficients = term_coefficients,
        )
        push!(shared_shell_layers, shell)
        push!(shared_shell_dimensions, size(shell.coefficient_matrix, 2))
        push!(shared_shell_provenance, shell.provenance)
        current_box = inner_box
    end

    candidates = _nested_square_lattice_split_candidates(
        bundles,
        current_box,
        x_coordinate_range,
        y_coordinate_range,
        basis.x_coordinates,
        basis.y_coordinates;
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
    )
    accepted_candidate = _nested_choose_square_lattice_candidate(
        candidates,
        local_x_coordinates,
        local_y_coordinates,
    )

    if isnothing(accepted_candidate)
        leaf_sequence = _nested_complete_shell_sequence_for_box(
            bundles,
            current_box;
            nside = nside,
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
            term_coefficients = term_coefficients,
        )
        merged_sequence =
            isempty(shared_shell_layers) ? leaf_sequence :
            _nested_shell_sequence_from_core_block(
                bundles,
                leaf_sequence.support_indices,
                leaf_sequence.coefficient_matrix,
                shared_shell_layers,
                term_coefficients = term_coefficients,
            )
        node = _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D(
            String(node_label),
            parent_box,
            current_box,
            basis.lattice_size,
            nside,
            min_in_plane_aspect_ratio,
            x_coordinate_range,
            y_coordinate_range,
            collect(local_x_coordinates),
            collect(local_y_coordinates),
            length(shared_shell_layers),
            shared_shell_dimensions,
            shared_shell_provenance,
            candidates,
            nothing,
            !isempty(candidates),
            _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D[],
            size(merged_sequence.coefficient_matrix, 2),
        )
        return (
            geometry = node,
            sequence = merged_sequence,
            leaf_sequences = _CartesianNestedShellSequence3D[leaf_sequence],
        )
    end

    candidate = candidates[accepted_candidate]
    child_results = [
        _nested_axis_aligned_homonuclear_square_lattice_node(
            basis,
            bundles,
            candidate.child_boxes[index],
            candidate.x_coordinate_ranges[index],
            candidate.y_coordinate_ranges[index],
            string(node_label, "_", index);
            nside = nside,
            min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
            prefer_midpoint_tie_side = prefer_midpoint_tie_side,
            retain_xy = retain_xy,
            retain_xz = retain_xz,
            retain_yz = retain_yz,
            retain_x_edge = retain_x_edge,
            retain_y_edge = retain_y_edge,
            retain_z_edge = retain_z_edge,
            term_coefficients = term_coefficients,
        ) for index in eachindex(candidate.child_boxes)
    ]
    child_nodes = [result.geometry for result in child_results]
    child_leaf_sequences = _CartesianNestedShellSequence3D[]
    for result in child_results
        append!(child_leaf_sequences, result.leaf_sequences)
    end
    core_support = reduce(vcat, [result.sequence.support_indices for result in child_results])
    core_coefficients = _nested_hcat_coefficient_maps(
        [result.sequence.coefficient_matrix for result in child_results],
    )
    merged_sequence =
        isempty(shared_shell_layers) ? _nested_shell_sequence_from_core_block(
            bundles,
            core_support,
            core_coefficients,
            _AbstractCartesianNestedShellLayer3D[];
            enforce_coverage = false,
            term_coefficients = term_coefficients,
        ) :
        _nested_shell_sequence_from_core_block(
            bundles,
            core_support,
            core_coefficients,
            shared_shell_layers,
            term_coefficients = term_coefficients,
        )
    node = _AxisAlignedHomonuclearSquareLatticeNodeGeometry3D(
        String(node_label),
        parent_box,
        current_box,
        basis.lattice_size,
        nside,
        min_in_plane_aspect_ratio,
        x_coordinate_range,
        y_coordinate_range,
        collect(local_x_coordinates),
        collect(local_y_coordinates),
        length(shared_shell_layers),
        shared_shell_dimensions,
        shared_shell_provenance,
        candidates,
        accepted_candidate,
        false,
        child_nodes,
        size(merged_sequence.coefficient_matrix, 2),
    )
    return (
        geometry = node,
        sequence = merged_sequence,
        leaf_sequences = child_leaf_sequences,
    )
end

function _nested_axis_aligned_homonuclear_square_lattice_source(
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
    prefer_midpoint_tie_side::Symbol = :left,
    retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
    retain_x_edge::Union{Nothing,Int} = nothing,
    retain_y_edge::Union{Nothing,Int} = nothing,
    retain_z_edge::Union{Nothing,Int} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    retention = _nested_resolve_complete_shell_retention(
        nside;
        retain_xy = retain_xy,
        retain_xz = retain_xz,
        retain_yz = retain_yz,
        retain_x_edge = retain_x_edge,
        retain_y_edge = retain_y_edge,
        retain_z_edge = retain_z_edge,
    )
    dims = _nested_axis_lengths(bundles)
    parent_box = (1:dims[1], 1:dims[2], 1:dims[3])
    result = _nested_axis_aligned_homonuclear_square_lattice_node(
        basis,
        bundles,
        parent_box,
        1:length(basis.x_coordinates),
        1:length(basis.y_coordinates),
        "root";
        nside = nside,
        min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
        prefer_midpoint_tie_side = prefer_midpoint_tie_side,
        retain_xy = retention.retain_xy,
        retain_xz = retention.retain_xz,
        retain_yz = retention.retain_yz,
        retain_x_edge = retention.retain_x_edge,
        retain_y_edge = retention.retain_y_edge,
        retain_z_edge = retention.retain_z_edge,
        term_coefficients = term_coefficients,
    )
    return _CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D(
        basis,
        bundles,
        retention,
        result.geometry,
        result.leaf_sequences,
        result.sequence,
    )
end

function _nested_fixed_block(source::_CartesianNestedBondAlignedHomonuclearChainSource3D)
    return _nested_fixed_block(source.sequence, source.basis)
end

function _nested_fixed_block(source::_CartesianNestedAxisAlignedHomonuclearSquareLatticeSource3D)
    return _nested_fixed_block(source.sequence, source.basis)
end
