module CartesianCPBBlockProviders

using ..CartesianCPB
using ..CartesianParentGaussletBases

const CPB = CartesianCPB
const CPGB = CartesianParentGaussletBases
const _AXIS_ORDER = (:x, :y, :z)
const _LOCAL_ORDERING = :parent_compatible_x_slowest_z_fastest

export CPBIntervalPair3D,
       CPBOverlapAxisBlockSet,
       cpb_interval_pair,
       cpb_overlap_axis_blocks,
       summary

"""
    CPBIntervalPair3D

Metadata-only CPB pair contract for future Cartesian CPB block providers.

This object validates CPB windows against a Cartesian parent and records
axis-ordered intervals and shapes. It does not consume parent factor packets or
slice any operator matrix.
"""
struct CPBIntervalPair3D{P,L,R,M}
    parent::P
    left_cpb::L
    right_cpb::R
    metadata::M
end

"""
    CPBOverlapAxisBlockSet

Axis-level overlap blocks for one CPB interval pair. The blocks are slices of a
parent-owned overlap factor packet; no dense local 3D block is materialized.
"""
struct CPBOverlapAxisBlockSet{P,I,B,M}
    overlap_packet::P
    interval_pair::I
    axis_blocks::B
    metadata::M
end

summary(pair::CPBIntervalPair3D) = pair.metadata
summary(block_set::CPBOverlapAxisBlockSet) = block_set.metadata

function cpb_interval_pair(
    parent::CPGB.CartesianParentGaussletBasis3D,
    left_cpb::CPB.CoordinateProductBox,
    right_cpb::CPB.CoordinateProductBox,
)
    parent_intervals = _axis_named_tuple(CPGB.parent_box(parent))
    left_intervals = _axis_named_tuple(CPB.intervals(left_cpb))
    right_intervals = _axis_named_tuple(CPB.intervals(right_cpb))
    outside_parent_intervals = _outside_parent_intervals(
        parent_intervals,
        left_intervals,
        right_intervals,
    )
    status =
        isempty(outside_parent_intervals) ?
        :available_cpb_interval_pair :
        :blocked_cpb_interval_pair
    blocker =
        isempty(outside_parent_intervals) ?
        nothing :
        _outside_interval_blocker(first(outside_parent_intervals))
    return CPBIntervalPair3D(
        parent,
        left_cpb,
        right_cpb,
        (;
            object_kind = :cartesian_cpb_interval_pair_summary,
            status,
            blocker,
            parent_axis_counts = CPGB.parent_axis_counts(parent),
            parent_intervals,
            left_intervals,
            right_intervals,
            left_shape = _axis_named_tuple(CPB.shape(left_cpb)),
            right_shape = _axis_named_tuple(CPB.shape(right_cpb)),
            left_support_count = CPB.support_count(left_cpb),
            right_support_count = CPB.support_count(right_cpb),
            left_codimension = CPB.codimension(left_cpb),
            right_codimension = CPB.codimension(right_cpb),
            axis_order = _AXIS_ORDER,
            local_ordering = _LOCAL_ORDERING,
            outside_parent_intervals,
            operator_matrices_sliced = false,
            parent_factor_packet_consumed = false,
            route_driver_wiring = false,
        ),
    )
end

function _axis_named_tuple(values)
    length(values) == 3 || throw(ArgumentError("axis tuple requires three values"))
    return (x = values[1], y = values[2], z = values[3])
end

function _outside_parent_intervals(parent_intervals, left_intervals, right_intervals)
    outside = NamedTuple[]
    for side in (:left, :right)
        intervals = side === :left ? left_intervals : right_intervals
        for axis in _AXIS_ORDER
            interval = getproperty(intervals, axis)
            parent_interval = getproperty(parent_intervals, axis)
            _interval_inside_parent(interval, parent_interval) && continue
            push!(
                outside,
                (;
                    side,
                    axis,
                    interval,
                    parent_interval,
                ),
            )
        end
    end
    return Tuple(outside)
end

function _interval_inside_parent(
    interval::UnitRange{Int},
    parent_interval::UnitRange{Int},
)
    return first(parent_interval) <= first(interval) &&
           last(interval) <= last(parent_interval)
end

function _outside_interval_blocker(record)
    return Symbol("$(record.side)_$(record.axis)_interval_outside_parent")
end

function cpb_overlap_axis_blocks(
    overlap_packet::CPGB.CartesianParentAxisFactorPacket3D,
    interval_pair::CPBIntervalPair3D,
)
    packet_summary = CPGB.summary(overlap_packet)
    interval_summary = summary(interval_pair)
    blocker = _cpb_overlap_axis_blocks_blocker(
        overlap_packet,
        packet_summary,
        interval_pair,
        interval_summary,
    )
    axis_blocks =
        isnothing(blocker) ?
        _cpb_overlap_axis_block_views(overlap_packet, interval_summary) :
        nothing
    status =
        isnothing(blocker) ?
        :available_cpb_overlap_axis_blocks :
        :blocked_cpb_overlap_axis_blocks
    return CPBOverlapAxisBlockSet(
        overlap_packet,
        interval_pair,
        axis_blocks,
        _cpb_overlap_axis_block_summary(
            status,
            blocker,
            packet_summary,
            interval_summary,
            axis_blocks,
        ),
    )
end

function _cpb_overlap_axis_blocks_blocker(
    overlap_packet,
    packet_summary,
    interval_pair,
    interval_summary,
)
    packet_summary.status === :available_parent_overlap_axis_factors ||
        return isnothing(packet_summary.blocker) ?
               :unavailable_parent_overlap_axis_factors :
               packet_summary.blocker
    interval_summary.status === :available_cpb_interval_pair ||
        return isnothing(interval_summary.blocker) ?
               :unavailable_cpb_interval_pair :
               interval_summary.blocker
    overlap_packet.parent === interval_pair.parent || return :parent_mismatch
    return nothing
end

function _cpb_overlap_axis_block_views(overlap_packet, interval_summary)
    overlap_1d = overlap_packet.overlap_1d
    left = interval_summary.left_intervals
    right = interval_summary.right_intervals
    return (;
        x = view(overlap_1d.x, left.x, right.x),
        y = view(overlap_1d.y, left.y, right.y),
        z = view(overlap_1d.z, left.z, right.z),
    )
end

function _cpb_overlap_axis_block_summary(
    status::Symbol,
    blocker,
    packet_summary,
    interval_summary,
    axis_blocks,
)
    available = status === :available_cpb_overlap_axis_blocks
    return (;
        object_kind = :cartesian_cpb_overlap_axis_block_set_summary,
        status,
        blocker,
        overlap_packet_summary = packet_summary,
        interval_pair_summary = interval_summary,
        axis_blocks_available = available,
        axis_block_shapes =
            available ?
            _axis_block_shapes(axis_blocks) :
            :unavailable,
        factor_space = packet_summary.factor_space,
        factor_convention = packet_summary.factor_convention,
        normalization_convention = packet_summary.normalization_convention,
        index_domain = packet_summary.index_domain,
        index_domain_source = packet_summary.index_domain_source,
        index_domain_status = packet_summary.index_domain_status,
        axis_order = packet_summary.axis_order,
        bra_ket_order = packet_summary.bra_ket_order,
        local_ordering = interval_summary.local_ordering,
        blocks_are_views = available,
        blocks_are_copies = false,
        dense_local_block_materialized = false,
        parent_factor_packet_consumed = true,
        route_driver_wiring = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        exports_or_artifacts = false,
    )
end

function _axis_block_shapes(axis_blocks)
    return (;
        x = size(axis_blocks.x),
        y = size(axis_blocks.y),
        z = size(axis_blocks.z),
    )
end

end # module CartesianCPBBlockProviders
