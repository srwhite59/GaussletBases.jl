module CartesianCPBBlockProviders

using ..CartesianCPB
using ..CartesianParentGaussletBases

const CPB = CartesianCPB
const CPGB = CartesianParentGaussletBases
const _AXIS_ORDER = (:x, :y, :z)
const _LOCAL_ORDERING = :parent_compatible_x_slowest_z_fastest

export CPBIntervalPair3D,
       CPBOverlapAxisBlockSet,
       CPBOverlapDenseBlock,
       CPBLocalOverlapBlockRecord,
       CPBLocalOverlapBlockCollection,
       cpb_interval_pair,
       cpb_overlap_axis_blocks,
       cpb_overlap_dense_block,
       cpb_local_overlap_block_record,
       cpb_local_overlap_block_collection,
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

"""
    CPBOverlapDenseBlock

Local dense CPB product-space overlap block materialized from overlap axis
blocks. This remains local to the CPB pair and does not imply route/global
matrix placement.
"""
struct CPBOverlapDenseBlock{S,D,M}
    axis_block_set::S
    dense_block::D
    metadata::M
end

"""
    CPBLocalOverlapBlockRecord

One local CPB overlap provider output prepared for future collection and
placement layers. The source block is kept by reference; metadata remains a
compact fingerprint and does not duplicate dense numerical payloads.
"""
struct CPBLocalOverlapBlockRecord{K,S,M}
    block_key::K
    source_block::S
    metadata::M
end

"""
    CPBLocalOverlapBlockCollection

Compact collection of local CPB overlap block records. This layer does not
assign placement, infer retained transforms, or assemble a global matrix.
"""
struct CPBLocalOverlapBlockCollection{R,M}
    records::R
    metadata::M
end

"""
    CPBRetainedTransformCarry

Metadata-only retained-transform carry for one side of one local CPB overlap
record. This object validates shape/count/order metadata but does not apply the
transform or assign global placement.
"""
struct CPBRetainedTransformCarry{K,S,T,M}
    block_key::K
    source_cpb_summary::S
    transform_object::T
    metadata::M
end

"""
    CPBSourcePairPlacementRange

Metadata-only retained-column range authority for one local source-pair
record. This object validates range/dimension/count metadata but does not place
matrices or assemble route/global overlap.
"""
struct CPBSourcePairPlacementRange{K,L,R,M}
    block_key::K
    left_transform_carry::L
    right_transform_carry::R
    metadata::M
end

summary(pair::CPBIntervalPair3D) = pair.metadata
summary(block_set::CPBOverlapAxisBlockSet) = block_set.metadata
summary(block::CPBOverlapDenseBlock) = block.metadata
summary(record::CPBLocalOverlapBlockRecord) = record.metadata
summary(collection::CPBLocalOverlapBlockCollection) = collection.metadata
summary(carry::CPBRetainedTransformCarry) = carry.metadata
summary(range::CPBSourcePairPlacementRange) = range.metadata

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
    sliceability_blocker = _cpb_overlap_packet_sliceability_blocker(packet_summary)
    isnothing(sliceability_blocker) || return sliceability_blocker
    interval_summary.status === :available_cpb_interval_pair ||
        return isnothing(interval_summary.blocker) ?
               :unavailable_cpb_interval_pair :
               interval_summary.blocker
    overlap_packet.parent === interval_pair.parent || return :parent_mismatch
    return nothing
end

function _cpb_overlap_packet_sliceability_blocker(packet_summary)
    _summary_property(packet_summary, :sliceable_by_cpb) === true ||
        return :overlap_packet_not_cpb_sliceable
    _summary_property(packet_summary, :index_domain) === :parent_axis_indices ||
        return :overlap_packet_not_cpb_sliceable
    _summary_property(packet_summary, :index_domain_source) === :axis_bundle_contract ||
        return :overlap_packet_not_cpb_sliceable
    _summary_property(
        packet_summary,
        :index_domain_status,
    ) === :assumed_parent_axis_indexed_by_current_axis_bundle_contract ||
        return :overlap_packet_not_cpb_sliceable
    sliceability_source = _summary_property(packet_summary, :sliceability_source)
    isnothing(sliceability_source) ||
        sliceability_source === :index_domain_contract ||
        return :overlap_packet_not_cpb_sliceable
    sliceability_status = _summary_property(packet_summary, :sliceability_status)
    isnothing(sliceability_status) ||
        sliceability_status === :sliceable_by_cpb_parent_axis_index_contract ||
        return :overlap_packet_not_cpb_sliceable
    return nothing
end

function _summary_property(summary, name::Symbol)
    hasproperty(summary, name) && return getproperty(summary, name)
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

function cpb_overlap_dense_block(axis_block_set::CPBOverlapAxisBlockSet)
    axis_block_summary = summary(axis_block_set)
    blocker =
        axis_block_summary.status === :available_cpb_overlap_axis_blocks ?
        nothing :
        _cpb_overlap_dense_block_blocker(axis_block_summary)
    dense_block =
        isnothing(blocker) ?
        _materialize_cpb_overlap_dense_block(axis_block_set.axis_blocks) :
        nothing
    status =
        isnothing(blocker) ?
        :materialized_cpb_overlap_dense_block :
        :blocked_cpb_overlap_dense_block
    return CPBOverlapDenseBlock(
        axis_block_set,
        dense_block,
        _cpb_overlap_dense_block_summary(
            status,
            blocker,
            axis_block_summary,
            dense_block,
        ),
    )
end

function _cpb_overlap_dense_block_blocker(axis_block_summary)
    isnothing(axis_block_summary.blocker) &&
        return :unavailable_cpb_overlap_axis_blocks
    return axis_block_summary.blocker
end

function _materialize_cpb_overlap_dense_block(axis_blocks)
    nx_left, nx_right = size(axis_blocks.x)
    ny_left, ny_right = size(axis_blocks.y)
    nz_left, nz_right = size(axis_blocks.z)
    element_type = promote_type(
        eltype(axis_blocks.x),
        eltype(axis_blocks.y),
        eltype(axis_blocks.z),
    )
    dense_block = Matrix{element_type}(
        undef,
        nx_left * ny_left * nz_left,
        nx_right * ny_right * nz_right,
    )
    for ix_left in 1:nx_left, iy_left in 1:ny_left, iz_left in 1:nz_left
        left_index = _local_product_index(
            ix_left,
            iy_left,
            iz_left,
            (nx_left, ny_left, nz_left),
        )
        for ix_right in 1:nx_right, iy_right in 1:ny_right, iz_right in 1:nz_right
            right_index = _local_product_index(
                ix_right,
                iy_right,
                iz_right,
                (nx_right, ny_right, nz_right),
            )
            dense_block[left_index, right_index] =
                axis_blocks.x[ix_left, ix_right] *
                axis_blocks.y[iy_left, iy_right] *
                axis_blocks.z[iz_left, iz_right]
        end
    end
    return dense_block
end

function _local_product_index(
    ix::Integer,
    iy::Integer,
    iz::Integer,
    shape::NTuple{3,Int},
)
    _nx, ny, nz = shape
    return (Int(ix) - 1) * ny * nz + (Int(iy) - 1) * nz + Int(iz)
end

function _cpb_overlap_dense_block_summary(
    status::Symbol,
    blocker,
    axis_block_summary,
    dense_block,
)
    available = status === :materialized_cpb_overlap_dense_block
    interval_summary = axis_block_summary.interval_pair_summary
    return (;
        object_kind = :cartesian_cpb_overlap_dense_block_summary,
        status,
        blocker,
        dense_block_available = available,
        dense_block_shape =
            available ?
            size(dense_block) :
            :unavailable,
        dense_block_eltype =
            available ?
            eltype(dense_block) :
            :unavailable,
        left_shape = interval_summary.left_shape,
        right_shape = interval_summary.right_shape,
        left_support_count = interval_summary.left_support_count,
        right_support_count = interval_summary.right_support_count,
        local_ordering = axis_block_summary.local_ordering,
        source_axis_block_summary = axis_block_summary,
        factor_space = axis_block_summary.factor_space,
        factor_convention = axis_block_summary.factor_convention,
        normalization_convention = axis_block_summary.normalization_convention,
        index_domain = axis_block_summary.index_domain,
        index_domain_source = axis_block_summary.index_domain_source,
        index_domain_status = axis_block_summary.index_domain_status,
        axis_order = axis_block_summary.axis_order,
        bra_ket_order = axis_block_summary.bra_ket_order,
        dense_local_block_materialized = available,
        route_driver_wiring = false,
        global_matrix_materialized = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        exports_or_artifacts = false,
    )
end

function cpb_local_overlap_block_record(
    source_block;
    block_key = :local_overlap_block,
)
    source_summary = summary(source_block)
    record_summary = _cpb_local_overlap_block_record_summary(
        block_key,
        source_block,
        source_summary,
    )
    return CPBLocalOverlapBlockRecord(block_key, source_block, record_summary)
end

function _cpb_local_overlap_block_record_summary(
    block_key,
    source_block::CPBOverlapDenseBlock,
    dense_block_summary,
)
    axis_block_summary = dense_block_summary.source_axis_block_summary
    interval_summary = axis_block_summary.interval_pair_summary
    available =
        dense_block_summary.status === :materialized_cpb_overlap_dense_block
    status =
        available ?
        :available_cpb_local_overlap_block_record :
        :blocked_cpb_local_overlap_block_record
    blocker =
        available ?
        nothing :
        _cpb_local_overlap_source_blocker(dense_block_summary)
    return _cpb_local_overlap_block_record_summary(
        block_key,
        :cpb_overlap_dense_block,
        status,
        blocker,
        interval_summary,
        axis_block_summary,
        dense_block_summary,
        dense_block_summary.dense_block_available,
        dense_block_summary.dense_block_shape,
        dense_block_summary.local_ordering,
        dense_block_summary.factor_space,
        dense_block_summary.factor_convention,
        dense_block_summary.normalization_convention,
        dense_block_summary.index_domain,
        dense_block_summary.index_domain_source,
        dense_block_summary.index_domain_status,
    )
end

function _cpb_local_overlap_block_record_summary(
    block_key,
    source_block::CPBOverlapAxisBlockSet,
    axis_block_summary,
)
    interval_summary = axis_block_summary.interval_pair_summary
    available = axis_block_summary.status === :available_cpb_overlap_axis_blocks
    status =
        available ?
        :available_cpb_local_overlap_block_record :
        :blocked_cpb_local_overlap_block_record
    blocker =
        available ?
        nothing :
        _cpb_local_overlap_source_blocker(axis_block_summary)
    return _cpb_local_overlap_block_record_summary(
        block_key,
        :cpb_overlap_axis_block_set,
        status,
        blocker,
        interval_summary,
        axis_block_summary,
        nothing,
        false,
        :not_materialized,
        axis_block_summary.local_ordering,
        axis_block_summary.factor_space,
        axis_block_summary.factor_convention,
        axis_block_summary.normalization_convention,
        axis_block_summary.index_domain,
        axis_block_summary.index_domain_source,
        axis_block_summary.index_domain_status,
    )
end

function _cpb_local_overlap_source_blocker(source_summary)
    isnothing(source_summary.blocker) && return :unavailable_cpb_local_overlap_source_block
    return source_summary.blocker
end

function _cpb_local_overlap_block_record_summary(
    block_key,
    source_kind::Symbol,
    status::Symbol,
    blocker,
    interval_summary,
    axis_block_summary,
    dense_block_summary,
    dense_block_available::Bool,
    dense_block_shape,
    local_ordering,
    factor_space,
    factor_convention,
    normalization_convention,
    index_domain,
    index_domain_source,
    index_domain_status,
)
    return (;
        object_kind = :cartesian_cpb_local_overlap_block_record_summary,
        term = :overlap,
        block_key,
        source_kind,
        status,
        blocker,
        left_cpb_summary = (;
            intervals = interval_summary.left_intervals,
            shape = interval_summary.left_shape,
            support_count = interval_summary.left_support_count,
        ),
        right_cpb_summary = (;
            intervals = interval_summary.right_intervals,
            shape = interval_summary.right_shape,
            support_count = interval_summary.right_support_count,
        ),
        interval_pair_summary = interval_summary,
        axis_block_summary,
        dense_block_summary,
        dense_block_available,
        dense_block_shape,
        local_ordering,
        factor_space,
        factor_convention,
        normalization_convention,
        index_domain,
        index_domain_source,
        index_domain_status,
        placement_status = :unassigned,
        retained_transform_status = :unassigned,
        global_matrix_materialized = false,
        route_driver_wiring = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        exports_or_artifacts = false,
    )
end

function cpb_local_overlap_block_collection(records::Union{Tuple,AbstractVector})
    normalized_records = _cpb_local_overlap_block_records_tuple(records)
    collection_summary =
        _cpb_local_overlap_block_collection_summary(normalized_records)
    return CPBLocalOverlapBlockCollection(normalized_records, collection_summary)
end

function cpb_local_overlap_block_collection(
    record::Union{
        CPBLocalOverlapBlockRecord,
        CPBOverlapDenseBlock,
        CPBOverlapAxisBlockSet,
    },
    records...,
)
    return cpb_local_overlap_block_collection((record, records...))
end

function _cpb_local_overlap_block_records_tuple(records)
    records isa Tuple && return _cpb_local_overlap_block_records_tuple(records, Val(:tuple))
    return _cpb_local_overlap_block_records_tuple(Tuple(records), Val(:tuple))
end

function _cpb_local_overlap_block_records_tuple(records::Tuple, ::Val{:tuple})
    return Tuple(
        _cpb_local_overlap_block_record_from_source(record, index)
        for (index, record) in enumerate(records)
    )
end

function _cpb_local_overlap_block_record_from_source(
    record::CPBLocalOverlapBlockRecord,
    _index::Integer,
)
    return record
end

function _cpb_local_overlap_block_record_from_source(source_block, index::Integer)
    return cpb_local_overlap_block_record(
        source_block;
        block_key = Symbol("local_overlap_block_", index),
    )
end

function _cpb_local_overlap_block_collection_summary(records::Tuple)
    record_summaries = Tuple(summary(record) for record in records)
    blocked = filter(record -> record.status !== :available_cpb_local_overlap_block_record,
        record_summaries)
    empty = isempty(records)
    status =
        empty ?
        :blocked_cpb_local_overlap_block_collection :
        isempty(blocked) ?
        :available_cpb_local_overlap_block_collection :
        :blocked_cpb_local_overlap_block_collection
    blocker =
        empty ?
        :empty_cpb_local_overlap_block_collection :
        isempty(blocked) ?
        nothing :
        :blocked_cpb_local_overlap_block_records
    return (;
        object_kind = :cartesian_cpb_local_overlap_block_collection_summary,
        status,
        blocker,
        record_count = length(records),
        record_summaries,
        terms = Tuple(unique(record.term for record in record_summaries)),
        block_keys = Tuple(record.block_key for record in record_summaries),
        dense_block_count =
            count(record -> record.dense_block_available, record_summaries),
        blocked_record_count = length(blocked),
        blocked_record_blockers = Tuple(record.blocker for record in blocked),
        placement_status = :unassigned,
        retained_transform_status = :unassigned,
        global_matrix_materialized = false,
        route_driver_wiring = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        exports_or_artifacts = false,
    )
end

function cpb_retained_transform_carry(
    side::Symbol,
    block_key,
    source_cpb_summary,
    source_shape,
    source_ordering,
    target_retained_column_range;
    transform_object = nothing,
    transform_convention = :unavailable,
    transform_provenance = :unavailable,
)
    side in (:left, :right) ||
        throw(ArgumentError("CPB retained transform carry side must be :left or :right"))
    source_support_count = _cpb_retained_transform_source_support_count(source_shape)
    target_retained_column_count =
        _cpb_retained_transform_target_count(target_retained_column_range)
    transform_available = !isnothing(transform_object)
    transform_shape =
        transform_object isa AbstractMatrix ? size(transform_object) : :unavailable
    transform_reference_kind =
        isnothing(transform_object) ?
        :unavailable :
        transform_object isa AbstractMatrix ?
        :matrix :
        Symbol(nameof(typeof(transform_object)))
    blocker = _cpb_retained_transform_carry_blocker(
        source_support_count,
        source_ordering,
        target_retained_column_count,
        transform_object,
        transform_shape,
    )
    status =
        isnothing(blocker) ?
        :available_cpb_retained_transform_carry :
        :blocked_cpb_retained_transform_carry
    return CPBRetainedTransformCarry(
        block_key,
        source_cpb_summary,
        transform_object,
        (;
            object_kind = :cartesian_cpb_retained_transform_carry_summary,
            status,
            blocker,
            side,
            block_key,
            source_cpb_summary,
            source_shape,
            source_support_count,
            source_ordering,
            target_retained_column_range,
            target_retained_column_count,
            transform_available,
            transform_convention,
            transform_provenance,
            transform_shape,
            transform_reference_kind,
            placement_status = :unassigned,
            transform_applied = false,
            global_matrix_materialized = false,
            route_driver_wiring = false,
        ),
    )
end

function _cpb_retained_transform_source_support_count(source_shape)
    if source_shape isa NamedTuple
        all(hasproperty(source_shape, axis) for axis in _AXIS_ORDER) ||
            return :unavailable
        return prod(Int(getproperty(source_shape, axis)) for axis in _AXIS_ORDER)
    elseif source_shape isa Tuple && length(source_shape) == 3
        return prod(Int(value) for value in source_shape)
    elseif source_shape isa Integer
        return Int(source_shape)
    end
    return :unavailable
end

function _cpb_retained_transform_target_count(target_retained_column_range)
    target_retained_column_range isa AbstractUnitRange{<:Integer} ||
        return :unavailable
    return length(target_retained_column_range)
end

function _cpb_retained_transform_carry_blocker(
    source_support_count,
    source_ordering,
    target_retained_column_count,
    transform_object,
    transform_shape,
)
    isnothing(transform_object) && return :missing_retained_transform
    source_ordering === _LOCAL_ORDERING ||
        return :retained_transform_ordering_mismatch
    transform_object isa AbstractMatrix || return nothing
    source_support_count isa Integer ||
        return :retained_transform_source_shape_mismatch
    target_retained_column_count isa Integer ||
        return :retained_transform_target_count_mismatch
    transform_shape[1] == source_support_count ||
        return :retained_transform_source_shape_mismatch
    transform_shape[2] == target_retained_column_count ||
        return :retained_transform_target_count_mismatch
    return nothing
end

function cpb_source_pair_placement_range(
    block_key;
    left_column_range = nothing,
    right_column_range = nothing,
    global_dimension = nothing,
    global_dimension_source = :unavailable,
    range_source = :unavailable,
    range_provenance = :unavailable,
    left_transform_carry = nothing,
    right_transform_carry = nothing,
)
    left_column_count = _cpb_source_pair_column_count(left_column_range)
    right_column_count = _cpb_source_pair_column_count(right_column_range)
    normalized_global_dimension =
        _cpb_source_pair_global_dimension(global_dimension)
    left_transform_summary =
        _cpb_source_pair_transform_summary(left_transform_carry)
    right_transform_summary =
        _cpb_source_pair_transform_summary(right_transform_carry)
    left_transform_status =
        _cpb_source_pair_transform_status(left_transform_summary)
    right_transform_status =
        _cpb_source_pair_transform_status(right_transform_summary)
    left_transform_target_count =
        _cpb_source_pair_transform_target_count(left_transform_summary)
    right_transform_target_count =
        _cpb_source_pair_transform_target_count(right_transform_summary)
    blocker = _cpb_source_pair_placement_range_blocker(
        left_column_range,
        right_column_range,
        left_column_count,
        right_column_count,
        normalized_global_dimension,
        left_transform_status,
        right_transform_status,
        left_transform_target_count,
        right_transform_target_count,
    )
    status =
        isnothing(blocker) ?
        :available_cpb_source_pair_placement_range :
        :blocked_cpb_source_pair_placement_range
    return CPBSourcePairPlacementRange(
        block_key,
        left_transform_carry,
        right_transform_carry,
        (;
            object_kind = :cartesian_cpb_source_pair_placement_range_summary,
            status,
            blocker,
            block_key,
            left_column_range,
            right_column_range,
            left_column_count,
            right_column_count,
            global_dimension = normalized_global_dimension,
            global_dimension_source,
            range_source,
            range_provenance,
            left_transform_status,
            right_transform_status,
            left_transform_target_count,
            right_transform_target_count,
            placement_status = :unassigned,
            global_matrix_materialized = false,
            route_driver_wiring = false,
        ),
    )
end

function _cpb_source_pair_column_count(column_range)
    column_range isa AbstractUnitRange{<:Integer} || return :unavailable
    return length(column_range)
end

function _cpb_source_pair_global_dimension(global_dimension)
    isnothing(global_dimension) && return nothing
    global_dimension isa Integer || return nothing
    dimension = Int(global_dimension)
    dimension > 0 || return nothing
    return dimension
end

function _cpb_source_pair_transform_summary(transform_carry)
    isnothing(transform_carry) && return nothing
    transform_carry isa CPBRetainedTransformCarry && return summary(transform_carry)
    transform_carry isa NamedTuple && return transform_carry
    return nothing
end

function _cpb_source_pair_transform_status(transform_summary)
    isnothing(transform_summary) && return :unavailable
    return _summary_property(transform_summary, :status)
end

function _cpb_source_pair_transform_target_count(transform_summary)
    isnothing(transform_summary) && return :unavailable
    target_count = _summary_property(transform_summary, :target_retained_column_count)
    isnothing(target_count) ? :unavailable : target_count
end

function _cpb_source_pair_placement_range_blocker(
    left_column_range,
    right_column_range,
    left_column_count,
    right_column_count,
    global_dimension,
    left_transform_status,
    right_transform_status,
    left_transform_target_count,
    right_transform_target_count,
)
    isnothing(left_column_range) && return :missing_left_column_range
    isnothing(right_column_range) && return :missing_right_column_range
    isnothing(global_dimension) && return :missing_global_dimension
    _cpb_source_pair_range_inside_dimension(left_column_range, global_dimension) ||
        return :left_column_range_dimension_mismatch
    _cpb_source_pair_range_inside_dimension(right_column_range, global_dimension) ||
        return :right_column_range_dimension_mismatch
    if left_transform_status === :available_cpb_retained_transform_carry
        left_column_count == left_transform_target_count ||
            return :left_column_range_dimension_mismatch
    end
    if right_transform_status === :available_cpb_retained_transform_carry
        right_column_count == right_transform_target_count ||
            return :right_column_range_dimension_mismatch
    end
    return nothing
end

function _cpb_source_pair_range_inside_dimension(column_range, global_dimension::Integer)
    column_range isa AbstractUnitRange{<:Integer} || return false
    return first(column_range) >= 1 && last(column_range) <= global_dimension
end

end # module CartesianCPBBlockProviders
