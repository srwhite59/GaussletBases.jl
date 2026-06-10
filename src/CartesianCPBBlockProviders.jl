module CartesianCPBBlockProviders

import ..GaussletBases: _mapped_ordinary_gausslet_1d_bundle,
                         basis_representation

using ..CartesianCPB
using ..CartesianParentGaussletBases

const CPB = CartesianCPB
const CPGB = CartesianParentGaussletBases
const _AXIS_ORDER = (:x, :y, :z)
const _LOCAL_ORDERING = :parent_compatible_x_slowest_z_fastest
const _OVERLAP_PLACEMENT_SYMMETRY_POLICIES = (:explicit_blocks_only,)
const _OVERLAP_PLACEMENT_DUPLICATE_POLICIES = (:reject_duplicate_block_keys,)

export CPBIntervalPair3D,
       CPBOverlapAxisBlockSet,
       CPBAxisProductOperatorBlock,
       CPBSumOfAxisProductsOperatorBlock,
       CPBOverlapOperatorBlock,
       CPBKineticOperatorBlock,
       CPBOneBodyAxisOperatorBlock,
       CPBMixedGTOLocalOverlapBlock,
       CPBElectronElectronLocalBlock,
       CPBElectronNuclearByCenterLocalBlock,
       CPBLocalIntegralWeights,
       CPBOverlapDenseBlock,
       CPBLocalOverlapBlockRecord,
       CPBLocalOverlapBlockCollection,
       cpb_interval_pair,
       cpb_overlap_axis_blocks,
       cpb_axis_product_operator_block,
       cpb_sum_of_axis_products_operator_block,
       cpb_overlap_operator_block,
       cpb_kinetic_operator_block,
       cpb_position_operator_block,
       cpb_x2_operator_block,
       cpb_mixed_gto_overlap_block,
       cpb_electron_electron_local_block,
       cpb_electron_nuclear_by_center_local_block,
       cpb_local_integral_weights,
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
    CPBAxisProductOperatorBlock

Generic dense local CPB product-space operator block materialized from
axis-local left/right blocks. This is local to the CPB operator layer and does
not imply route/global placement or realization.
"""
struct CPBAxisProductOperatorBlock{A,D,M}
    axis_ops::A
    dense_block::D
    metadata::M
end

"""
    CPBSumOfAxisProductsOperatorBlock

Dense local CPB product-space operator block materialized from a sum of
axis-product terms. This is for simple separable one-body terms only; it does
not imply Coulomb-family Gaussian-sum kernels, realization, or route/global
placement.
"""
struct CPBSumOfAxisProductsOperatorBlock{T,D,M}
    product_terms::T
    dense_block::D
    metadata::M
end

"""
    CPBOverlapOperatorBlock

Thin overlap wrapper around the generic CPB axis-product operator block. The
overlap term supplies the local overlap axis blocks; dense product-space
materialization is performed by `cpb_axis_product_operator_block`.
"""
struct CPBOverlapOperatorBlock{P,I,S,B,M}
    overlap_packet::P
    interval_pair::I
    axis_block_set::S
    axis_product_block::B
    metadata::M
end

"""
    CPBKineticOperatorBlock

Thin CPB-local kinetic wrapper around the sum-of-axis-products primitive. It
uses parent-owned 1D kinetic factors and explicit overlap factors on inactive
directions. It does not imply WL/PQS realization or route/global placement.
"""
struct CPBKineticOperatorBlock{P,I,B,M}
    parent_axis_factor_packet::P
    interval_pair::I
    sum_axis_product_block::B
    metadata::M
end

"""
    CPBOneBodyAxisOperatorBlock

Thin CPB-local wrapper for one-axis one-body terms such as position_x and
x2_y. It uses a selected parent one-body factor on the active axis and
explicit overlap factors on inactive axes.
"""
struct CPBOneBodyAxisOperatorBlock{P,I,B,M}
    parent_axis_factor_packet::P
    interval_pair::I
    axis_product_block::B
    metadata::M
end

"""
    CPBMixedGTOLocalOverlapBlock

Provider-level CPB-local mixed gausslet/GTO overlap pilot for one supplement
orbital. This reuses the existing Cartesian/GTO axis overlap helper and does
not imply route/global placement, Hamiltonian assembly, or broader GTO
operator support.
"""
struct CPBMixedGTOLocalOverlapBlock{P,C,O,D,M}
    parent::P
    cpb::C
    orbital::O
    dense_block::D
    metadata::M
end

"""
    CPBElectronElectronLocalBlock

Provider-level CPB-local electron-electron pair-factor interaction block. This
pilot uses parent axis pair-factor terms in their existing White-Lindsey
convention and keeps the Gaussian expansion term loop inside the local
contraction. It does not apply axis integral weights, route/global placement,
Hamiltonian assembly, or WL/PQS realization.
"""
struct CPBElectronElectronLocalBlock{S,E,I,D,M}
    parent_axis_bundle_object::S
    expansion::E
    interval_pair::I
    dense_block::D
    metadata::M
end

"""
    CPBElectronNuclearByCenterLocalBlock

Provider-level CPB-local electron-nuclear Galerkin block for one center. This
pilot keeps center terms separate, keeps the Gaussian expansion term loop
inside the local contraction, and does not apply CPB integral weights or sum a
Hamiltonian.
"""
struct CPBElectronNuclearByCenterLocalBlock{S,E,C,I,D,M}
    parent_axis_bundle_object::S
    expansion::E
    center_record::C
    interval_pair::I
    dense_block::D
    metadata::M
end

"""
    CPBLocalIntegralWeights

Provider-level CPB-local product integral weights built from parent-axis 1D
integral weights. These weights are a local realization/final-density support
object; they are not applied inside raw pair-factor Coulomb blocks.
"""
struct CPBLocalIntegralWeights{S,C,W,M}
    parent_axis_bundle_object::S
    cpb::C
    local_weights::W
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
struct CPBLocalOverlapBlockCollection{R,M,I}
    records::R
    metadata::M
    identity_token::I
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

"""
    CPBOverlapPlacementFacts

Metadata-only coherence bundle for local CPB overlap placement facts. This
object records collection, transform-carry, placement-range, placement-plan,
and accumulation-rule statuses, but does not apply transforms or place
matrices.
"""
struct CPBOverlapPlacementFacts{C,T,R,M}
    collection::C
    transform_carries::T
    placement_ranges::R
    metadata::M
end

"""
    CPBReviewedOverlapPlacementPlan

Metadata-only reviewed overlap placement plan. This object records the planned
placement contract but does not apply transforms or assemble any matrix.
"""
struct CPBReviewedOverlapPlacementPlan{K,M}
    accepted_block_keys::K
    metadata::M
end

"""
    CPBPlacedOverlapBlock

Provider-level synthetic placement pilot for one local CPB overlap block. This
object may carry a dense placement matrix, but it is not route-driver adoption
and does not imply route-global overlap availability.
"""
struct CPBPlacedOverlapBlock{S,L,R,P,F,G,M}
    source_block::S
    left_transform_carry::L
    right_transform_carry::R
    placement_range::P
    placement_facts::F
    global_overlap_matrix::G
    metadata::M
end

"""
    CPBPlacedOverlapCollection

Provider-level synthetic placement pilot for a local CPB overlap block
collection. This object may carry a dense provider-level target matrix, but it
is not route-driver adoption and does not imply route-global overlap
availability.
"""
struct CPBPlacedOverlapCollection{C,T,R,F,G,M}
    collection::C
    transform_carries::T
    placement_ranges::R
    placement_facts::F
    global_overlap_matrix::G
    metadata::M
end

summary(pair::CPBIntervalPair3D) = pair.metadata
summary(block_set::CPBOverlapAxisBlockSet) = block_set.metadata
summary(block::CPBAxisProductOperatorBlock) = block.metadata
summary(block::CPBSumOfAxisProductsOperatorBlock) = block.metadata
summary(block::CPBOverlapOperatorBlock) = block.metadata
summary(block::CPBKineticOperatorBlock) = block.metadata
summary(block::CPBOneBodyAxisOperatorBlock) = block.metadata
summary(block::CPBMixedGTOLocalOverlapBlock) = block.metadata
summary(block::CPBElectronElectronLocalBlock) = block.metadata
summary(block::CPBElectronNuclearByCenterLocalBlock) = block.metadata
summary(weights::CPBLocalIntegralWeights) = weights.metadata
summary(block::CPBOverlapDenseBlock) = block.metadata
summary(record::CPBLocalOverlapBlockRecord) = record.metadata
summary(collection::CPBLocalOverlapBlockCollection) = collection.metadata
summary(carry::CPBRetainedTransformCarry) = carry.metadata
summary(range::CPBSourcePairPlacementRange) = range.metadata
summary(facts::CPBOverlapPlacementFacts) = facts.metadata
summary(plan::CPBReviewedOverlapPlacementPlan) = plan.metadata
summary(placed::CPBPlacedOverlapBlock) = placed.metadata
summary(placed::CPBPlacedOverlapCollection) = placed.metadata

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

function cpb_axis_product_operator_block(
    axis_ops;
    term = :axis_product_operator,
    source_kind = :axis_product_axis_ops,
    source_summary = :unavailable,
    factor_space = :unavailable,
    factor_convention = :unavailable,
    normalization_convention = :unavailable,
    index_domain = :unavailable,
    index_domain_source = :unavailable,
    index_domain_status = :unavailable,
    axis_order = _AXIS_ORDER,
    bra_ket_order = (:bra, :ket),
)
    blocker = _cpb_axis_product_operator_block_blocker(axis_ops)
    dense_block =
        isnothing(blocker) ?
        _materialize_cpb_axis_product_dense_block(axis_ops) :
        nothing
    status =
        isnothing(blocker) ?
        :materialized_cpb_axis_product_operator_block :
        :blocked_cpb_axis_product_operator_block
    return CPBAxisProductOperatorBlock(
        axis_ops,
        dense_block,
        _cpb_axis_product_operator_block_summary(
            status,
            blocker,
            axis_ops,
            dense_block;
            term,
            source_kind,
            source_summary,
            factor_space,
            factor_convention,
            normalization_convention,
            index_domain,
            index_domain_source,
            index_domain_status,
            axis_order,
            bra_ket_order,
        ),
    )
end

function _cpb_axis_product_operator_block_blocker(axis_ops)
    for axis in _AXIS_ORDER
        hasproperty(axis_ops, axis) || return Symbol("missing_$(axis)_axis_operator")
        axis_op = getproperty(axis_ops, axis)
        axis_op isa AbstractMatrix ||
            return Symbol("$(axis)_axis_operator_not_matrix")
        isempty(axis_op) && return Symbol("$(axis)_axis_operator_empty")
    end
    return nothing
end

function _materialize_cpb_axis_product_dense_block(axis_ops)
    nx_left, nx_right = size(axis_ops.x)
    ny_left, ny_right = size(axis_ops.y)
    nz_left, nz_right = size(axis_ops.z)
    element_type = promote_type(
        eltype(axis_ops.x),
        eltype(axis_ops.y),
        eltype(axis_ops.z),
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
                axis_ops.x[ix_left, ix_right] *
                axis_ops.y[iy_left, iy_right] *
                axis_ops.z[iz_left, iz_right]
        end
    end
    return dense_block
end

function _cpb_axis_product_operator_block_summary(
    status::Symbol,
    blocker,
    axis_ops,
    dense_block;
    term,
    source_kind,
    source_summary,
    factor_space,
    factor_convention,
    normalization_convention,
    index_domain,
    index_domain_source,
    index_domain_status,
    axis_order,
    bra_ket_order,
)
    available = status === :materialized_cpb_axis_product_operator_block
    left_shape =
        available ?
        (x = size(axis_ops.x, 1), y = size(axis_ops.y, 1), z = size(axis_ops.z, 1)) :
        :unavailable
    right_shape =
        available ?
        (x = size(axis_ops.x, 2), y = size(axis_ops.y, 2), z = size(axis_ops.z, 2)) :
        :unavailable
    return (;
        object_kind = :cartesian_cpb_axis_product_operator_block_summary,
        status,
        blocker,
        term,
        source_kind,
        source_summary,
        representation = :dense_local_cpb_product_space,
        local_ordering = _LOCAL_ORDERING,
        left_shape,
        right_shape,
        left_support_count =
            available ? left_shape.x * left_shape.y * left_shape.z : :unavailable,
        right_support_count =
            available ? right_shape.x * right_shape.y * right_shape.z : :unavailable,
        axis_operator_shapes =
            available ? _axis_block_shapes(axis_ops) : :unavailable,
        dense_block_available = available,
        dense_block_shape =
            available ?
            size(dense_block) :
            :unavailable,
        dense_block_eltype =
            available ?
            eltype(dense_block) :
            :unavailable,
        provider_level_local_matrix_materialized = available,
        dense_local_block_materialized = available,
        realization_status = :unrealized,
        route_global_status = :unassigned,
        factor_space,
        factor_convention,
        normalization_convention,
        index_domain,
        index_domain_source,
        index_domain_status,
        axis_order,
        bra_ket_order,
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function cpb_sum_of_axis_products_operator_block(
    product_terms;
    term = :sum_of_axis_products_operator,
    source_kind = :axis_product_term_list,
    source_summary = :unavailable,
    factor_space = :unavailable,
    factor_convention = :unavailable,
    normalization_convention = :unavailable,
    index_domain = :unavailable,
    index_domain_source = :unavailable,
    index_domain_status = :unavailable,
    axis_order = _AXIS_ORDER,
    bra_ket_order = (:bra, :ket),
)
    normalized_terms = Tuple(product_terms)
    product_term_summaries = _cpb_axis_product_term_summaries(normalized_terms)
    blocker = _cpb_sum_of_axis_products_operator_block_blocker(
        product_term_summaries,
    )
    dense_block =
        isnothing(blocker) ?
        _materialize_cpb_sum_of_axis_products_dense_block(normalized_terms) :
        nothing
    status =
        isnothing(blocker) ?
        :materialized_cpb_sum_of_axis_products_operator_block :
        :blocked_cpb_sum_of_axis_products_operator_block
    return CPBSumOfAxisProductsOperatorBlock(
        normalized_terms,
        dense_block,
        _cpb_sum_of_axis_products_operator_block_summary(
            status,
            blocker,
            product_term_summaries,
            dense_block;
            term,
            source_kind,
            source_summary,
            factor_space,
            factor_convention,
            normalization_convention,
            index_domain,
            index_domain_source,
            index_domain_status,
            axis_order,
            bra_ket_order,
        ),
    )
end

function _cpb_axis_product_term_summaries(product_terms::Tuple)
    return Tuple(
        _cpb_axis_product_term_summary(product_term, product_term_index)
        for (product_term_index, product_term) in enumerate(product_terms)
    )
end

function _cpb_axis_product_term_summary(product_term, product_term_index::Integer)
    label =
        hasproperty(product_term, :label) ?
        getproperty(product_term, :label) :
        Symbol("axis_product_term_", product_term_index)
    coefficient =
        hasproperty(product_term, :coefficient) ?
        getproperty(product_term, :coefficient) :
        nothing
    coefficient_available = coefficient isa Number
    axis_ops =
        hasproperty(product_term, :axis_ops) ?
        getproperty(product_term, :axis_ops) :
        nothing
    axis_product_blocker =
        isnothing(axis_ops) ?
        :missing_axis_ops :
        _cpb_axis_product_operator_block_blocker(axis_ops)
    blocker =
        coefficient_available ?
        axis_product_blocker :
        :axis_product_term_coefficient_unavailable
    available = isnothing(blocker)
    left_shape =
        available ?
        (x = size(axis_ops.x, 1), y = size(axis_ops.y, 1), z = size(axis_ops.z, 1)) :
        :unavailable
    right_shape =
        available ?
        (x = size(axis_ops.x, 2), y = size(axis_ops.y, 2), z = size(axis_ops.z, 2)) :
        :unavailable
    return (;
        object_kind = :cartesian_cpb_axis_product_term_summary,
        product_term_index = Int(product_term_index),
        label,
        status =
            available ?
            :available_cpb_axis_product_term :
            :blocked_cpb_axis_product_term,
        blocker,
        coefficient_available,
        coefficient = coefficient_available ? coefficient : :unavailable,
        axis_product_blocker,
        left_shape,
        right_shape,
        left_support_count =
            available ? left_shape.x * left_shape.y * left_shape.z : :unavailable,
        right_support_count =
            available ? right_shape.x * right_shape.y * right_shape.z : :unavailable,
        axis_operator_shapes =
            available ? _axis_block_shapes(axis_ops) : :unavailable,
        dense_block_available = available,
        dense_block_shape =
            available ?
            (
                left_shape.x * left_shape.y * left_shape.z,
                right_shape.x * right_shape.y * right_shape.z,
            ) :
            :unavailable,
        dense_block_eltype =
            available ?
            promote_type(
                typeof(coefficient),
                eltype(axis_ops.x),
                eltype(axis_ops.y),
                eltype(axis_ops.z),
            ) :
            :unavailable,
    )
end

function _cpb_sum_of_axis_products_operator_block_blocker(product_term_summaries)
    isempty(product_term_summaries) && return :empty_axis_product_term_list
    for product_term_summary in product_term_summaries
        product_term_summary.blocker ===
            :axis_product_term_coefficient_unavailable &&
            return :axis_product_term_coefficient_unavailable
    end
    for product_term_summary in product_term_summaries
        isnothing(product_term_summary.blocker) || return :axis_product_term_blocked
    end
    first_left_shape = first(product_term_summaries).left_shape
    first_right_shape = first(product_term_summaries).right_shape
    for product_term_summary in product_term_summaries
        (
            product_term_summary.left_shape == first_left_shape &&
            product_term_summary.right_shape == first_right_shape
        ) || return :axis_product_term_shape_mismatch
    end
    return nothing
end

function _materialize_cpb_sum_of_axis_products_dense_block(product_terms::Tuple)
    first_term = first(product_terms)
    first_axis_ops = getproperty(first_term, :axis_ops)
    first_dense = _materialize_cpb_axis_product_dense_block(first_axis_ops)
    dense_eltype = _cpb_sum_of_axis_products_dense_eltype(product_terms)
    dense_block = zeros(dense_eltype, size(first_dense))
    for product_term in product_terms
        coefficient = getproperty(product_term, :coefficient)
        axis_dense = _materialize_cpb_axis_product_dense_block(
            getproperty(product_term, :axis_ops),
        )
        dense_block .+= coefficient .* axis_dense
    end
    return dense_block
end

function _cpb_sum_of_axis_products_dense_eltype(product_terms::Tuple)
    element_type = Union{}
    for product_term in product_terms
        coefficient = getproperty(product_term, :coefficient)
        axis_ops = getproperty(product_term, :axis_ops)
        element_type = promote_type(
            element_type,
            typeof(coefficient),
            eltype(axis_ops.x),
            eltype(axis_ops.y),
            eltype(axis_ops.z),
        )
    end
    return element_type
end

function _cpb_sum_of_axis_products_operator_block_summary(
    status::Symbol,
    blocker,
    product_term_summaries,
    dense_block;
    term,
    source_kind,
    source_summary,
    factor_space,
    factor_convention,
    normalization_convention,
    index_domain,
    index_domain_source,
    index_domain_status,
    axis_order,
    bra_ket_order,
)
    available =
        status === :materialized_cpb_sum_of_axis_products_operator_block
    first_product_term_summary =
        isempty(product_term_summaries) ? nothing : first(product_term_summaries)
    return (;
        object_kind = :cartesian_cpb_sum_of_axis_products_operator_block_summary,
        status,
        blocker,
        term,
        source_kind,
        source_summary,
        representation = :dense_local_cpb_sum_of_axis_products,
        local_ordering = _LOCAL_ORDERING,
        product_term_count = length(product_term_summaries),
        product_term_labels = Tuple(
            product_term_summary.label
            for product_term_summary in product_term_summaries
        ),
        product_term_summaries,
        left_shape =
            available ? first_product_term_summary.left_shape : :unavailable,
        right_shape =
            available ? first_product_term_summary.right_shape : :unavailable,
        left_support_count =
            available ? first_product_term_summary.left_support_count : :unavailable,
        right_support_count =
            available ? first_product_term_summary.right_support_count : :unavailable,
        dense_block_available = available,
        dense_block_shape =
            available ?
            size(dense_block) :
            :unavailable,
        dense_block_eltype =
            available ?
            eltype(dense_block) :
            :unavailable,
        provider_level_local_matrix_materialized = available,
        dense_local_block_materialized = available,
        realization_status = :unrealized,
        route_global_status = :unassigned,
        factor_space,
        factor_convention,
        normalization_convention,
        index_domain,
        index_domain_source,
        index_domain_status,
        axis_order,
        bra_ket_order,
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function cpb_overlap_operator_block(
    overlap_packet::CPGB.CartesianParentAxisFactorPacket3D,
    interval_pair::CPBIntervalPair3D,
)
    axis_block_set = cpb_overlap_axis_blocks(overlap_packet, interval_pair)
    axis_block_summary = summary(axis_block_set)
    axis_product_block =
        axis_block_summary.status === :available_cpb_overlap_axis_blocks ?
        cpb_axis_product_operator_block(
            axis_block_set.axis_blocks;
            term = :overlap,
            source_kind = :cpb_overlap_axis_blocks,
            source_summary = axis_block_summary,
            factor_space = axis_block_summary.factor_space,
            factor_convention = axis_block_summary.factor_convention,
            normalization_convention = axis_block_summary.normalization_convention,
            index_domain = axis_block_summary.index_domain,
            index_domain_source = axis_block_summary.index_domain_source,
            index_domain_status = axis_block_summary.index_domain_status,
            axis_order = axis_block_summary.axis_order,
            bra_ket_order = axis_block_summary.bra_ket_order,
        ) :
        nothing
    axis_product_summary =
        isnothing(axis_product_block) ? nothing : summary(axis_product_block)
    blocker =
        isnothing(axis_product_summary) ?
        _cpb_overlap_dense_block_blocker(axis_block_summary) :
        axis_product_summary.blocker
    status =
        isnothing(blocker) ?
        :materialized_cpb_overlap_operator_block :
        :blocked_cpb_overlap_operator_block
    return CPBOverlapOperatorBlock(
        overlap_packet,
        interval_pair,
        axis_block_set,
        axis_product_block,
        _cpb_overlap_operator_block_summary(
            status,
            blocker,
            axis_block_summary,
            axis_product_summary,
        ),
    )
end

function _cpb_overlap_operator_block_summary(
    status::Symbol,
    blocker,
    axis_block_summary,
    axis_product_summary,
)
    available = status === :materialized_cpb_overlap_operator_block
    return (;
        object_kind = :cartesian_cpb_overlap_operator_block_summary,
        status,
        blocker,
        term = :overlap,
        source_kind = :cpb_overlap_axis_blocks,
        axis_block_summary,
        axis_product_block_summary =
            isnothing(axis_product_summary) ? :unavailable : axis_product_summary,
        representation = :dense_local_cpb_product_space,
        local_ordering = axis_block_summary.local_ordering,
        left_shape =
            available ? axis_product_summary.left_shape : :unavailable,
        right_shape =
            available ? axis_product_summary.right_shape : :unavailable,
        dense_block_available = available,
        dense_block_shape =
            available ? axis_product_summary.dense_block_shape : :unavailable,
        dense_block_eltype =
            available ? axis_product_summary.dense_block_eltype : :unavailable,
        provider_level_local_matrix_materialized = available,
        realization_status = :unrealized,
        route_global_status = :unassigned,
        factor_space = axis_block_summary.factor_space,
        factor_convention = axis_block_summary.factor_convention,
        normalization_convention = axis_block_summary.normalization_convention,
        index_domain = axis_block_summary.index_domain,
        index_domain_source = axis_block_summary.index_domain_source,
        index_domain_status = axis_block_summary.index_domain_status,
        axis_order = axis_block_summary.axis_order,
        bra_ket_order = axis_block_summary.bra_ket_order,
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function cpb_kinetic_operator_block(
    parent_axis_factor_packet::CPGB.CartesianParentAxisFactorPacket3D,
    interval_pair::CPBIntervalPair3D,
)
    packet_summary = CPGB.summary(parent_axis_factor_packet)
    interval_summary = summary(interval_pair)
    blocker = _cpb_kinetic_operator_block_blocker(
        parent_axis_factor_packet,
        packet_summary,
        interval_pair,
        interval_summary,
    )
    product_terms =
        isnothing(blocker) ?
        _cpb_kinetic_axis_product_terms(parent_axis_factor_packet, interval_summary) :
        ()
    sum_block =
        isnothing(blocker) ?
        cpb_sum_of_axis_products_operator_block(
            product_terms;
            term = :kinetic,
            source_kind = :parent_axis_factor_packet_kinetic_overlap,
            source_summary = packet_summary,
            factor_space = :parent_axis_bundle_pgdg_intermediate,
            factor_convention = :axis_bundle_one_body_kinetic_sum,
            normalization_convention =
                :not_separate_from_axis_bundle_one_body_factors,
            index_domain = :parent_axis_indices,
            index_domain_source = :axis_bundle_contract,
            index_domain_status =
                :assumed_parent_axis_indexed_by_current_axis_bundle_contract,
            axis_order = packet_summary.axis_order,
            bra_ket_order = packet_summary.bra_ket_order,
        ) :
        nothing
    sum_summary = isnothing(sum_block) ? nothing : summary(sum_block)
    blocker =
        isnothing(sum_summary) ?
        blocker :
        sum_summary.blocker
    status =
        isnothing(blocker) ?
        :materialized_cpb_kinetic_operator_block :
        :blocked_cpb_kinetic_operator_block
    return CPBKineticOperatorBlock(
        parent_axis_factor_packet,
        interval_pair,
        sum_block,
        _cpb_kinetic_operator_block_summary(
            status,
            blocker,
            packet_summary,
            interval_summary,
            sum_summary,
        ),
    )
end

function _cpb_kinetic_operator_block_blocker(
    parent_axis_factor_packet,
    packet_summary,
    interval_pair,
    interval_summary,
)
    packet_summary.status === :available_parent_overlap_axis_factors ||
        return isnothing(packet_summary.blocker) ?
               :unavailable_parent_overlap_axis_factors :
               packet_summary.blocker
    _summary_property(packet_summary, :kinetic_1d_available) === true ||
        return :missing_parent_axis_bundle_kinetic_factors
    _summary_property(packet_summary, :kinetic_sliceable_by_cpb) === true ||
        return :kinetic_packet_not_cpb_sliceable
    _summary_property(packet_summary, :kinetic_index_domain) ===
        :parent_axis_indices ||
        return :kinetic_packet_not_cpb_sliceable
    _summary_property(packet_summary, :kinetic_index_domain_source) ===
        :axis_bundle_contract ||
        return :kinetic_packet_not_cpb_sliceable
    _summary_property(
        packet_summary,
        :kinetic_index_domain_status,
    ) === :assumed_parent_axis_indexed_by_current_axis_bundle_contract ||
        return :kinetic_packet_not_cpb_sliceable
    interval_summary.status === :available_cpb_interval_pair ||
        return isnothing(interval_summary.blocker) ?
               :unavailable_cpb_interval_pair :
               interval_summary.blocker
    parent_axis_factor_packet.parent === interval_pair.parent || return :parent_mismatch
    return nothing
end

function _cpb_kinetic_axis_product_terms(parent_axis_factor_packet, interval_summary)
    overlap_1d = parent_axis_factor_packet.overlap_1d
    kinetic_1d = parent_axis_factor_packet.kinetic_1d
    left = interval_summary.left_intervals
    right = interval_summary.right_intervals
    overlap_blocks = (;
        x = view(overlap_1d.x, left.x, right.x),
        y = view(overlap_1d.y, left.y, right.y),
        z = view(overlap_1d.z, left.z, right.z),
    )
    kinetic_blocks = (;
        x = view(kinetic_1d.x, left.x, right.x),
        y = view(kinetic_1d.y, left.y, right.y),
        z = view(kinetic_1d.z, left.z, right.z),
    )
    return (
        (;
            coefficient = 1.0,
            axis_ops = (x = kinetic_blocks.x, y = overlap_blocks.y, z = overlap_blocks.z),
            label = :kinetic_x_component,
        ),
        (;
            coefficient = 1.0,
            axis_ops = (x = overlap_blocks.x, y = kinetic_blocks.y, z = overlap_blocks.z),
            label = :kinetic_y_component,
        ),
        (;
            coefficient = 1.0,
            axis_ops = (x = overlap_blocks.x, y = overlap_blocks.y, z = kinetic_blocks.z),
            label = :kinetic_z_component,
        ),
    )
end

function _cpb_kinetic_operator_block_summary(
    status::Symbol,
    blocker,
    packet_summary,
    interval_summary,
    sum_summary,
)
    available = status === :materialized_cpb_kinetic_operator_block
    return (;
        object_kind = :cartesian_cpb_kinetic_operator_block_summary,
        status,
        blocker,
        term = :kinetic,
        source_kind = :parent_axis_factor_packet_kinetic_overlap,
        parent_axis_factor_packet_summary = packet_summary,
        interval_pair_summary = interval_summary,
        sum_axis_product_block_summary =
            isnothing(sum_summary) ? :unavailable : sum_summary,
        representation = :dense_local_cpb_sum_of_axis_products,
        local_ordering = interval_summary.local_ordering,
        kinetic_factor_form = :sum_of_axis_products,
        kinetic_component_axes = (:x, :y, :z),
        kinetic_component_terms =
            (:kinetic_x_component, :kinetic_y_component, :kinetic_z_component),
        product_term_count = available ? sum_summary.product_term_count : :unavailable,
        product_term_labels =
            available ? sum_summary.product_term_labels : :unavailable,
        left_shape = available ? sum_summary.left_shape : :unavailable,
        right_shape = available ? sum_summary.right_shape : :unavailable,
        dense_block_available = available,
        dense_block_shape =
            available ? sum_summary.dense_block_shape : :unavailable,
        dense_block_eltype =
            available ? sum_summary.dense_block_eltype : :unavailable,
        provider_level_local_matrix_materialized = available,
        realization_status = :unrealized,
        route_global_status = :unassigned,
        factor_space =
            available ? sum_summary.factor_space : :unavailable,
        factor_convention =
            available ? sum_summary.factor_convention : :unavailable,
        normalization_convention =
            available ? sum_summary.normalization_convention : :unavailable,
        index_domain =
            available ? sum_summary.index_domain : :unavailable,
        index_domain_source =
            available ? sum_summary.index_domain_source : :unavailable,
        index_domain_status =
            available ? sum_summary.index_domain_status : :unavailable,
        axis_order = packet_summary.axis_order,
        bra_ket_order =
            available ? packet_summary.bra_ket_order : :unavailable,
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function cpb_position_operator_block(
    parent_axis_factor_packet::CPGB.CartesianParentAxisFactorPacket3D,
    interval_pair::CPBIntervalPair3D;
    axis::Symbol,
)
    return _cpb_one_body_axis_operator_block(
        parent_axis_factor_packet,
        interval_pair;
        factor_name = :position,
        axis,
        term = Symbol("position_", String(axis)),
        factor_convention = :axis_bundle_one_body_position,
    )
end

function cpb_x2_operator_block(
    parent_axis_factor_packet::CPGB.CartesianParentAxisFactorPacket3D,
    interval_pair::CPBIntervalPair3D;
    axis::Symbol,
)
    return _cpb_one_body_axis_operator_block(
        parent_axis_factor_packet,
        interval_pair;
        factor_name = :x2,
        axis,
        term = Symbol("x2_", String(axis)),
        factor_convention = :axis_bundle_one_body_x2,
    )
end

function _cpb_one_body_axis_operator_block(
    parent_axis_factor_packet::CPGB.CartesianParentAxisFactorPacket3D,
    interval_pair::CPBIntervalPair3D;
    factor_name::Symbol,
    axis::Symbol,
    term::Symbol,
    factor_convention::Symbol,
)
    packet_summary = CPGB.summary(parent_axis_factor_packet)
    interval_summary = summary(interval_pair)
    blocker = _cpb_one_body_axis_operator_block_blocker(
        parent_axis_factor_packet,
        packet_summary,
        interval_pair,
        interval_summary;
        factor_name,
        axis,
    )
    axis_ops =
        isnothing(blocker) ?
        _cpb_one_body_axis_ops(
            parent_axis_factor_packet,
            interval_summary;
            factor_name,
            axis,
        ) :
        nothing
    axis_product_block =
        isnothing(blocker) ?
        cpb_axis_product_operator_block(
            axis_ops;
            term,
            source_kind =
                Symbol("parent_axis_factor_packet_", String(factor_name)),
            source_summary = packet_summary,
            factor_space = :parent_axis_bundle_pgdg_intermediate,
            factor_convention,
            normalization_convention =
                :not_separate_from_axis_bundle_one_body_factors,
            index_domain = :parent_axis_indices,
            index_domain_source = :axis_bundle_contract,
            index_domain_status =
                :assumed_parent_axis_indexed_by_current_axis_bundle_contract,
            axis_order = packet_summary.axis_order,
            bra_ket_order = packet_summary.bra_ket_order,
        ) :
        nothing
    axis_product_summary =
        isnothing(axis_product_block) ? nothing : summary(axis_product_block)
    blocker =
        isnothing(axis_product_summary) ?
        blocker :
        axis_product_summary.blocker
    status =
        isnothing(blocker) ?
        :materialized_cpb_one_body_axis_operator_block :
        :blocked_cpb_one_body_axis_operator_block
    return CPBOneBodyAxisOperatorBlock(
        parent_axis_factor_packet,
        interval_pair,
        axis_product_block,
        _cpb_one_body_axis_operator_block_summary(
            status,
            blocker,
            packet_summary,
            interval_summary,
            axis_product_summary;
            factor_name,
            axis,
            term,
            factor_convention,
        ),
    )
end

function _cpb_one_body_axis_operator_block_blocker(
    parent_axis_factor_packet,
    packet_summary,
    interval_pair,
    interval_summary;
    factor_name::Symbol,
    axis::Symbol,
)
    axis in _AXIS_ORDER || return Symbol("unsupported_$(factor_name)_axis")
    packet_summary.status === :available_parent_overlap_axis_factors ||
        return isnothing(packet_summary.blocker) ?
               :unavailable_parent_overlap_axis_factors :
               packet_summary.blocker
    _summary_property(
        packet_summary,
        Symbol(String(factor_name), "_1d_available"),
    ) === true ||
        return Symbol("missing_parent_axis_bundle_$(factor_name)_factors")
    _summary_property(
        packet_summary,
        Symbol(String(factor_name), "_sliceable_by_cpb"),
    ) === true ||
        return Symbol("$(factor_name)_packet_not_cpb_sliceable")
    _summary_property(
        packet_summary,
        Symbol(String(factor_name), "_index_domain"),
    ) ===
        :parent_axis_indices ||
        return Symbol("$(factor_name)_packet_not_cpb_sliceable")
    _summary_property(
        packet_summary,
        Symbol(String(factor_name), "_index_domain_source"),
    ) ===
        :axis_bundle_contract ||
        return Symbol("$(factor_name)_packet_not_cpb_sliceable")
    _summary_property(
        packet_summary,
        Symbol(String(factor_name), "_index_domain_status"),
    ) === :assumed_parent_axis_indexed_by_current_axis_bundle_contract ||
        return Symbol("$(factor_name)_packet_not_cpb_sliceable")
    interval_summary.status === :available_cpb_interval_pair ||
        return isnothing(interval_summary.blocker) ?
               :unavailable_cpb_interval_pair :
               interval_summary.blocker
    parent_axis_factor_packet.parent === interval_pair.parent || return :parent_mismatch
    return nothing
end

function _cpb_one_body_axis_ops(
    parent_axis_factor_packet,
    interval_summary;
    factor_name::Symbol,
    axis::Symbol,
)
    overlap_1d = parent_axis_factor_packet.overlap_1d
    one_body_1d = getproperty(
        parent_axis_factor_packet,
        Symbol(String(factor_name), "_1d"),
    )
    left = interval_summary.left_intervals
    right = interval_summary.right_intervals
    overlap_blocks = (;
        x = view(overlap_1d.x, left.x, right.x),
        y = view(overlap_1d.y, left.y, right.y),
        z = view(overlap_1d.z, left.z, right.z),
    )
    one_body_blocks = (;
        x = view(one_body_1d.x, left.x, right.x),
        y = view(one_body_1d.y, left.y, right.y),
        z = view(one_body_1d.z, left.z, right.z),
    )
    return (;
        x = axis === :x ? one_body_blocks.x : overlap_blocks.x,
        y = axis === :y ? one_body_blocks.y : overlap_blocks.y,
        z = axis === :z ? one_body_blocks.z : overlap_blocks.z,
    )
end

function _cpb_one_body_axis_operator_block_summary(
    status::Symbol,
    blocker,
    packet_summary,
    interval_summary,
    axis_product_summary;
    factor_name::Symbol,
    axis::Symbol,
    term::Symbol,
    factor_convention::Symbol,
)
    available = status === :materialized_cpb_one_body_axis_operator_block
    return (;
        object_kind = :cartesian_cpb_one_body_axis_operator_block_summary,
        status,
        blocker,
        term,
        source_kind = Symbol("parent_axis_factor_packet_", String(factor_name)),
        parent_axis_factor_packet_summary = packet_summary,
        interval_pair_summary = interval_summary,
        axis_product_block_summary =
            isnothing(axis_product_summary) ? :unavailable : axis_product_summary,
        representation = :dense_local_cpb_product_space,
        local_ordering = interval_summary.local_ordering,
        one_body_factor_name = factor_name,
        active_axis = axis,
        operator_factor_form =
            Symbol(
                String(factor_name),
                "_on_selected_axis_overlap_on_inactive_axes",
            ),
        left_shape = available ? axis_product_summary.left_shape : :unavailable,
        right_shape = available ? axis_product_summary.right_shape : :unavailable,
        dense_block_available = available,
        dense_block_shape =
            available ? axis_product_summary.dense_block_shape : :unavailable,
        dense_block_eltype =
            available ? axis_product_summary.dense_block_eltype : :unavailable,
        provider_level_local_matrix_materialized = available,
        realization_status = :unrealized,
        route_global_status = :unassigned,
        factor_space =
            available ? axis_product_summary.factor_space : :unavailable,
        factor_convention =
            available ? factor_convention : :unavailable,
        normalization_convention =
            available ? axis_product_summary.normalization_convention : :unavailable,
        index_domain =
            available ? axis_product_summary.index_domain : :unavailable,
        index_domain_source =
            available ? axis_product_summary.index_domain_source : :unavailable,
        index_domain_status =
            available ? axis_product_summary.index_domain_status : :unavailable,
        axis_order = packet_summary.axis_order,
        bra_ket_order =
            available ? packet_summary.bra_ket_order : :unavailable,
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function cpb_mixed_gto_overlap_block(
    parent::CPGB.CartesianParentGaussletBasis3D,
    cpb::CPB.CoordinateProductBox,
    orbital,
)
    blocker = _cpb_mixed_gto_overlap_block_blocker(parent, cpb, orbital)
    dense_block =
        isnothing(blocker) ?
        _materialize_cpb_mixed_gto_overlap_block(parent, cpb, orbital) :
        nothing
    status =
        isnothing(blocker) ?
        :materialized_cpb_mixed_gto_local_overlap_block :
        :blocked_cpb_mixed_gto_local_overlap_block
    return CPBMixedGTOLocalOverlapBlock(
        parent,
        cpb,
        orbital,
        dense_block,
        _cpb_mixed_gto_overlap_block_summary(
            status,
            blocker,
            parent,
            cpb,
            orbital,
            dense_block,
        ),
    )
end

function _cpb_mixed_gto_overlap_block_blocker(
    parent::CPGB.CartesianParentGaussletBasis3D,
    cpb::CPB.CoordinateProductBox,
    orbital,
)
    _is_cartesian_gaussian_shell_orbital(orbital) ||
        return :unsupported_supplement_orbital_record
    length(orbital.angular_powers) == 3 ||
        return :invalid_supplement_angular_powers
    length(orbital.center) == 3 || return :invalid_supplement_center
    outside_interval = _outside_cpb_parent_interval(parent, cpb)
    isnothing(outside_interval) ||
        return Symbol("cpb_$(outside_interval.axis)_interval_outside_parent")
    orbital.primitive_normalization === :axiswise_normalized_cartesian_gaussian ||
        return :unsupported_supplement_primitive_normalization
    all(power -> power >= 0, orbital.angular_powers) ||
        return :invalid_supplement_angular_powers
    isempty(orbital.exponents) && return :missing_supplement_primitive_data
    isempty(orbital.coefficients) && return :missing_supplement_primitive_data
    length(orbital.exponents) == length(orbital.coefficients) ||
        return :supplement_primitive_count_mismatch
    all(exponent -> isfinite(exponent) && exponent > 0.0, orbital.exponents) ||
        return :invalid_supplement_exponents
    all(isfinite, orbital.coefficients) ||
        return :invalid_supplement_coefficients
    return nothing
end

function _outside_cpb_parent_interval(
    parent::CPGB.CartesianParentGaussletBasis3D,
    cpb::CPB.CoordinateProductBox,
)
    parent_intervals = _axis_named_tuple(CPGB.parent_box(parent))
    cpb_intervals = _axis_named_tuple(CPB.intervals(cpb))
    for axis in _AXIS_ORDER
        interval = getproperty(cpb_intervals, axis)
        parent_interval = getproperty(parent_intervals, axis)
        _interval_inside_parent(interval, parent_interval) && continue
        return (;
            axis,
            interval,
            parent_interval,
        )
    end
    return nothing
end

function _materialize_cpb_mixed_gto_overlap_block(
    parent::CPGB.CartesianParentGaussletBasis3D,
    cpb::CPB.CoordinateProductBox,
    orbital,
)
    axis_cross = _cpb_mixed_gto_axis_overlap_tables(parent, orbital)
    intervals = _axis_named_tuple(CPB.intervals(cpb))
    local_shape = CPB.shape(cpb)
    dense_block = zeros(Float64, CPB.support_count(cpb), 1)
    for (local_ix, ix) in enumerate(intervals.x),
            (local_iy, iy) in enumerate(intervals.y),
            (local_iz, iz) in enumerate(intervals.z)
        row = _local_product_index(local_ix, local_iy, local_iz, local_shape)
        value = 0.0
        for primitive in eachindex(orbital.coefficients)
            value +=
                Float64(orbital.coefficients[primitive]) *
                axis_cross.x[ix, primitive] *
                axis_cross.y[iy, primitive] *
                axis_cross.z[iz, primitive]
        end
        dense_block[row, 1] = value
    end
    return dense_block
end

function _cpb_mixed_gto_axis_overlap_tables(
    parent::CPGB.CartesianParentGaussletBasis3D,
    orbital,
)
    axes = CPGB.parent_axes(parent)
    axis_representations = (;
        x = basis_representation(axes.x; operators = (:overlap,)),
        y = basis_representation(axes.y; operators = (:overlap,)),
        z = basis_representation(axes.z; operators = (:overlap,)),
    )
    return (;
        x = _cpb_mixed_gto_axis_overlap_table(
            axis_representations.x,
            orbital,
            :x,
        ),
        y = _cpb_mixed_gto_axis_overlap_table(
            axis_representations.y,
            orbital,
            :y,
        ),
        z = _cpb_mixed_gto_axis_overlap_table(
            axis_representations.z,
            orbital,
            :z,
        ),
    )
end

function _is_cartesian_gaussian_shell_orbital(orbital)
    return hasproperty(orbital, :label) &&
           hasproperty(orbital, :angular_powers) &&
           hasproperty(orbital, :center) &&
           hasproperty(orbital, :exponents) &&
           hasproperty(orbital, :coefficients) &&
           hasproperty(orbital, :primitive_normalization)
end

function _cpb_mixed_gto_axis_overlap_table(axis_representation, orbital, axis::Symbol)
    helper = getproperty(
        parentmodule(@__MODULE__),
        :_cartesian_basis_supplement_axis_primitive_cross,
    )
    return helper(axis_representation, orbital, axis)
end

function _cpb_mixed_gto_overlap_block_summary(
    status::Symbol,
    blocker,
    parent::CPGB.CartesianParentGaussletBasis3D,
    cpb::CPB.CoordinateProductBox,
    orbital,
    dense_block,
)
    available = status === :materialized_cpb_mixed_gto_local_overlap_block
    local_shape = _axis_named_tuple(CPB.shape(cpb))
    primitive_count =
        hasproperty(orbital, :exponents) ? length(orbital.exponents) : :unavailable
    return (;
        object_kind = :cartesian_cpb_mixed_gto_local_overlap_block_summary,
        status,
        blocker,
        term = :mixed_gto_overlap,
        source_kind = :mixed_gausslet_gto_supplement_overlap,
        supplement_representation_kind = :cartesian_gaussian_shell_orbital,
        orbital_label = _mixed_gto_property(orbital, :label),
        angular_powers = _mixed_gto_property(orbital, :angular_powers),
        center = _mixed_gto_property(orbital, :center),
        primitive_count,
        primitive_normalization =
            _mixed_gto_property(orbital, :primitive_normalization),
        formula_source = :GaussianAnalyticIntegrals_polynomial_gaussian,
        axis_kernel_source = :existing_cartesian_basis_supplement_axis_cross,
        parent_axis_counts = CPGB.parent_axis_counts(parent),
        parent_axis_counts_source = :parent_object_parent_axis_counts,
        cpb_summary = (;
            role = CPB.role(cpb),
            intervals = _axis_named_tuple(CPB.intervals(cpb)),
            shape = local_shape,
            support_count = CPB.support_count(cpb),
            codimension = CPB.codimension(cpb),
        ),
        local_shape,
        left_shape = local_shape,
        right_shape = (gto = 1,),
        dense_block_available = available,
        dense_block_shape = available ? size(dense_block) : :unavailable,
        dense_block_eltype = available ? eltype(dense_block) : :unavailable,
        local_ordering = _LOCAL_ORDERING,
        representation = :dense_local_cpb_product_space_by_gto_orbital,
        galerkin_operator = true,
        provider_level_local_matrix_materialized = available,
        realization_status = :unrealized,
        route_global_status = :unassigned,
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function _mixed_gto_property(orbital, name::Symbol)
    hasproperty(orbital, name) && return getproperty(orbital, name)
    return :unavailable
end

function cpb_electron_electron_local_block(
    parent_axis_bundle_object,
    expansion,
    interval_pair::CPBIntervalPair3D,
)
    interval_summary = summary(interval_pair)
    source_summary = CPGB.summary(
        CPGB.parent_coulomb_axis_source_summary(
            interval_pair.parent,
            parent_axis_bundle_object,
            expansion,
        ),
    )
    pair_factor_terms = _cpb_electron_electron_pair_factor_terms(
        parent_axis_bundle_object,
    )
    blocker = _cpb_electron_electron_local_block_blocker(
        source_summary,
        interval_summary,
        pair_factor_terms,
        expansion,
    )
    dense_block =
        isnothing(blocker) ?
        _materialize_cpb_electron_electron_local_dense_block(
            pair_factor_terms,
            expansion,
            interval_summary,
        ) :
        nothing
    status =
        isnothing(blocker) ?
        :materialized_cpb_electron_electron_local_block :
        :blocked_cpb_electron_electron_local_block
    return CPBElectronElectronLocalBlock(
        parent_axis_bundle_object,
        expansion,
        interval_pair,
        dense_block,
        _cpb_electron_electron_local_block_summary(
            status,
            blocker,
            source_summary,
            interval_summary,
            pair_factor_terms,
            expansion,
            dense_block,
        ),
    )
end

function _cpb_electron_electron_pair_factor_terms(parent_axis_bundle_object)
    return (;
        x = _cpb_axis_pair_factor_terms(parent_axis_bundle_object, :x),
        y = _cpb_axis_pair_factor_terms(parent_axis_bundle_object, :y),
        z = _cpb_axis_pair_factor_terms(parent_axis_bundle_object, :z),
    )
end

function _cpb_axis_pair_factor_terms(parent_axis_bundle_object, axis::Symbol)
    axis_bundle = _cpb_axis_bundle(parent_axis_bundle_object, axis)
    pgdg_intermediate = _summary_property(axis_bundle, :pgdg_intermediate)
    return _summary_property(pgdg_intermediate, :pair_factor_terms)
end

function _cpb_axis_bundle(parent_axis_bundle_object, axis::Symbol)
    axis === :x && hasproperty(parent_axis_bundle_object, :x) &&
        return getproperty(parent_axis_bundle_object, :x)
    axis === :y && hasproperty(parent_axis_bundle_object, :y) &&
        return getproperty(parent_axis_bundle_object, :y)
    axis === :z && hasproperty(parent_axis_bundle_object, :z) &&
        return getproperty(parent_axis_bundle_object, :z)
    bundle_field = Symbol("bundle_", String(axis))
    hasproperty(parent_axis_bundle_object, bundle_field) &&
        return getproperty(parent_axis_bundle_object, bundle_field)
    return nothing
end

function _cpb_electron_electron_local_block_blocker(
    source_summary,
    interval_summary,
    pair_factor_terms,
    expansion,
)
    interval_summary.status === :available_cpb_interval_pair ||
        return isnothing(interval_summary.blocker) ?
               :unavailable_cpb_interval_pair :
               interval_summary.blocker
    source_summary.axis_bundle_available ||
        return :missing_parent_axis_bundle_object
    source_summary.expansion_available ||
        return :missing_coulomb_gaussian_expansion
    source_summary.expansion_coefficients_available ||
        return :missing_coulomb_expansion_coefficients
    source_summary.expansion_exponents_available ||
        return :missing_coulomb_expansion_exponents
    coefficients = _summary_property(expansion, :coefficients)
    coefficients isa AbstractVector ||
        return :missing_coulomb_expansion_coefficients
    for axis in _AXIS_ORDER
        terms = getproperty(pair_factor_terms, axis)
        terms isa AbstractArray{<:Real,3} ||
            return Symbol("missing_$(axis)_pair_factor_terms")
        isempty(terms) && return Symbol("$(axis)_pair_factor_terms_empty")
        size(terms, 1) == length(coefficients) ||
            return :pair_factor_term_count_mismatch
        parent_axis_count = _cpb_axis_value(
            interval_summary.parent_axis_counts,
            axis,
        )
        (
            size(terms, 2) == parent_axis_count &&
            size(terms, 3) == parent_axis_count
        ) || return Symbol("$(axis)_pair_factor_terms_size_mismatch")
    end
    return nothing
end

function _materialize_cpb_electron_electron_local_dense_block(
    pair_factor_terms,
    expansion,
    interval_summary,
)
    coefficients = Float64[Float64(value) for value in expansion.coefficients]
    left = interval_summary.left_intervals
    right = interval_summary.right_intervals
    left_shape = interval_summary.left_shape
    right_shape = interval_summary.right_shape
    dense_block = zeros(
        Float64,
        interval_summary.left_support_count,
        interval_summary.right_support_count,
    )

    # Pair-factor terms already match the existing WL interaction oracle. Keep
    # alpha as the inner, stride-1 first index of the term-first factor arrays.
    # Do not apply/divide axis weights here; weights belong at a later
    # realization or final density interpretation boundary.
    for ix_left in left.x, iy_left in left.y, iz_left in left.z
        row = _local_product_index(
            ix_left - first(left.x) + 1,
            iy_left - first(left.y) + 1,
            iz_left - first(left.z) + 1,
            (left_shape.x, left_shape.y, left_shape.z),
        )
        for ix_right in right.x, iy_right in right.y, iz_right in right.z
            column = _local_product_index(
                ix_right - first(right.x) + 1,
                iy_right - first(right.y) + 1,
                iz_right - first(right.z) + 1,
                (right_shape.x, right_shape.y, right_shape.z),
            )
            value = 0.0
            @inbounds for term in eachindex(coefficients)
                value +=
                    coefficients[term] *
                    pair_factor_terms.x[term, ix_left, ix_right] *
                    pair_factor_terms.y[term, iy_left, iy_right] *
                    pair_factor_terms.z[term, iz_left, iz_right]
            end
            dense_block[row, column] = value
        end
    end
    return dense_block
end

function _cpb_electron_electron_local_block_summary(
    status::Symbol,
    blocker,
    source_summary,
    interval_summary,
    pair_factor_terms,
    expansion,
    dense_block,
)
    available = status === :materialized_cpb_electron_electron_local_block
    coefficients = _summary_property(expansion, :coefficients)
    coefficient_count = coefficients isa AbstractVector ? length(coefficients) : 0
    return (;
        object_kind = :cartesian_cpb_electron_electron_local_block_summary,
        status,
        blocker,
        term = :electron_electron_pair_factor_interaction,
        source_kind = :parent_axis_bundle_pair_factor_terms,
        parent_coulomb_source_summary = source_summary,
        interval_pair_summary = interval_summary,
        representation = :dense_local_cpb_electron_electron_pair_factor_interaction,
        pair_factor_source = :axis_pgdg_intermediate_pair_factor_terms,
        pair_factor_weighting = :wl_pair_factor_terms_existing_convention,
        axis_integral_weights_applied = false,
        axis_integral_weights_deferred = true,
        weight_application_stage = :realization_or_final_density_interpretation,
        retained_pqs_weights = false,
        factor_source_path = :axis_pgdg_intermediate_pair_factor_terms,
        gaussian_expansion_loop = :inner_local_contraction,
        gaussian_term_count = coefficient_count,
        coefficient_count,
        pair_factor_term_shapes =
            _cpb_pair_factor_term_shapes(pair_factor_terms),
        local_ordering = interval_summary.local_ordering,
        left_shape = available ? interval_summary.left_shape : :unavailable,
        right_shape = available ? interval_summary.right_shape : :unavailable,
        left_support_count =
            available ? interval_summary.left_support_count : :unavailable,
        right_support_count =
            available ? interval_summary.right_support_count : :unavailable,
        dense_block_available = available,
        dense_block_shape = available ? size(dense_block) : :unavailable,
        dense_block_eltype = available ? eltype(dense_block) : :unavailable,
        provider_level_local_matrix_materialized = available,
        provider_level_coulomb_block_materialized = available,
        provider_level_electron_electron_block_materialized = available,
        cpb_local_coulomb_kernel_implemented = available,
        cpb_local_electron_electron_kernel_implemented = available,
        realization_status = :unrealized,
        route_global_status = :unassigned,
        factor_space =
            available ? :parent_axis_bundle_pgdg_intermediate : :unavailable,
        factor_convention =
            available ? :axis_bundle_electron_electron_pair_factor_terms : :unavailable,
        index_domain =
            available ? :parent_axis_indices : :unavailable,
        index_domain_source =
            available ? :axis_bundle_contract : :unavailable,
        index_domain_status =
            available ?
            :assumed_parent_axis_indexed_by_current_axis_bundle_contract :
            :unavailable,
        axis_order = _AXIS_ORDER,
        bra_ket_order = (:density_left, :density_right),
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_assembly = false,
        hamiltonian_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function _cpb_pair_factor_term_shapes(pair_factor_terms)
    return (;
        x = _cpb_factor_term_shape(pair_factor_terms.x),
        y = _cpb_factor_term_shape(pair_factor_terms.y),
        z = _cpb_factor_term_shape(pair_factor_terms.z),
    )
end

function _cpb_factor_term_shape(terms)
    terms isa AbstractArray ? size(terms) : :unavailable
end

function _cpb_axis_value(values, axis::Symbol)
    axis === :x && return values[1]
    axis === :y && return values[2]
    axis === :z && return values[3]
    throw(ArgumentError("unsupported Cartesian axis $(axis)"))
end

function cpb_electron_nuclear_by_center_local_block(
    parent_axis_bundle_object,
    expansion,
    center_record,
    interval_pair::CPBIntervalPair3D,
)
    interval_summary = summary(interval_pair)
    center_summary = _cpb_electron_nuclear_center_summary(center_record)
    axis_terms = _cpb_electron_nuclear_axis_terms(
        parent_axis_bundle_object,
        expansion,
        center_summary,
    )
    blocker = _cpb_electron_nuclear_by_center_local_block_blocker(
        interval_summary,
        center_summary,
        axis_terms,
        expansion,
    )
    dense_block =
        isnothing(blocker) ?
        _materialize_cpb_electron_nuclear_by_center_dense_block(
            axis_terms,
            expansion,
            center_summary,
            interval_summary,
        ) :
        nothing
    status =
        isnothing(blocker) ?
        :materialized_cpb_electron_nuclear_by_center_local_block :
        :blocked_cpb_electron_nuclear_by_center_local_block
    return CPBElectronNuclearByCenterLocalBlock(
        parent_axis_bundle_object,
        expansion,
        center_record,
        interval_pair,
        dense_block,
        _cpb_electron_nuclear_by_center_local_block_summary(
            status,
            blocker,
            center_summary,
            interval_summary,
            axis_terms,
            expansion,
            dense_block,
        ),
    )
end

function _cpb_electron_nuclear_center_summary(center_record)
    isnothing(center_record) && return (;
        status = :blocked_electron_nuclear_center_record,
        blocker = :missing_electron_nuclear_center_record,
        center_key = :unavailable,
        center_index = :unavailable,
        charge = :unavailable,
        location = :unavailable,
        x = :unavailable,
        y = :unavailable,
        z = :unavailable,
    )
    charge = _summary_property(center_record, :charge)
    isnothing(charge) && (charge = _summary_property(center_record, :nuclear_charge))
    isnothing(charge) && (charge = _summary_property(center_record, :Z))
    location = _summary_property(center_record, :location)
    if isnothing(location)
        x = _summary_property(center_record, :x)
        y = _summary_property(center_record, :y)
        z = _summary_property(center_record, :z)
        if !isnothing(x) && !isnothing(y) && !isnothing(z)
            location = (x, y, z)
        end
    end
    location_tuple =
        _cpb_electron_nuclear_location_tuple(location)
    blocker =
        isnothing(charge) ?
        :missing_electron_nuclear_center_charge :
        (
            isnothing(location_tuple) ?
            :missing_electron_nuclear_center_location :
            nothing
        )
    return (;
        status =
            isnothing(blocker) ?
            :available_electron_nuclear_center_record :
            :blocked_electron_nuclear_center_record,
        blocker,
        center_key =
            something(_summary_property(center_record, :center_key), :unavailable),
        center_index =
            something(_summary_property(center_record, :center_index), :unavailable),
        charge =
            isnothing(charge) ? :unavailable : Float64(charge),
        location =
            isnothing(location_tuple) ? :unavailable : location_tuple,
        x = isnothing(location_tuple) ? :unavailable : location_tuple[1],
        y = isnothing(location_tuple) ? :unavailable : location_tuple[2],
        z = isnothing(location_tuple) ? :unavailable : location_tuple[3],
    )
end

function _cpb_electron_nuclear_location_tuple(location)
    isnothing(location) && return nothing
    length(location) == 3 || return nothing
    return (
        Float64(location[1]),
        Float64(location[2]),
        Float64(location[3]),
    )
end

function _cpb_electron_nuclear_axis_terms(
    parent_axis_bundle_object,
    expansion,
    center_summary,
)
    exponents = _summary_property(expansion, :exponents)
    center_summary.status === :available_electron_nuclear_center_record ||
        return (x = nothing, y = nothing, z = nothing)
    exponents isa AbstractVector || return (x = nothing, y = nothing, z = nothing)
    return (;
        x = _cpb_electron_nuclear_axis_terms(
            parent_axis_bundle_object,
            :x,
            exponents,
            center_summary.x,
        ),
        y = _cpb_electron_nuclear_axis_terms(
            parent_axis_bundle_object,
            :y,
            exponents,
            center_summary.y,
        ),
        z = _cpb_electron_nuclear_axis_terms(
            parent_axis_bundle_object,
            :z,
            exponents,
            center_summary.z,
        ),
    )
end

function _cpb_electron_nuclear_axis_terms(
    parent_axis_bundle_object,
    axis::Symbol,
    exponents,
    center_coordinate,
)
    axis_bundle = _cpb_axis_bundle(parent_axis_bundle_object, axis)
    basis = _summary_property(axis_bundle, :basis)
    backend = _summary_property(axis_bundle, :backend)
    isnothing(basis) && return nothing
    isnothing(backend) && return nothing
    centered_bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents,
        center = center_coordinate,
        backend,
    )
    pgdg_intermediate = _summary_property(centered_bundle, :pgdg_intermediate)
    return _summary_property(pgdg_intermediate, :gaussian_factor_terms)
end

function _cpb_electron_nuclear_by_center_local_block_blocker(
    interval_summary,
    center_summary,
    axis_terms,
    expansion,
)
    interval_summary.status === :available_cpb_interval_pair ||
        return isnothing(interval_summary.blocker) ?
               :unavailable_cpb_interval_pair :
               interval_summary.blocker
    center_summary.status === :available_electron_nuclear_center_record ||
        return center_summary.blocker
    isnothing(expansion) && return :missing_coulomb_gaussian_expansion
    coefficients = _summary_property(expansion, :coefficients)
    exponents = _summary_property(expansion, :exponents)
    coefficients isa AbstractVector ||
        return :missing_coulomb_expansion_coefficients
    exponents isa AbstractVector ||
        return :missing_coulomb_expansion_exponents
    for axis in _AXIS_ORDER
        terms = getproperty(axis_terms, axis)
        terms isa AbstractArray{<:Real,3} ||
            return Symbol("missing_$(axis)_electron_nuclear_axis_factor_terms")
        isempty(terms) &&
            return Symbol("$(axis)_electron_nuclear_axis_factor_terms_empty")
        size(terms, 1) == length(coefficients) ||
            return :electron_nuclear_axis_factor_term_count_mismatch
        parent_axis_count = _cpb_axis_value(
            interval_summary.parent_axis_counts,
            axis,
        )
        (
            size(terms, 2) == parent_axis_count &&
            size(terms, 3) == parent_axis_count
        ) || return Symbol("$(axis)_electron_nuclear_axis_factor_terms_size_mismatch")
    end
    return nothing
end

function _materialize_cpb_electron_nuclear_by_center_dense_block(
    axis_terms,
    expansion,
    center_summary,
    interval_summary,
)
    coefficients = Float64[-center_summary.charge * Float64(value) for value in expansion.coefficients]
    left = interval_summary.left_intervals
    right = interval_summary.right_intervals
    left_shape = interval_summary.left_shape
    right_shape = interval_summary.right_shape
    dense_block = zeros(
        Float64,
        interval_summary.left_support_count,
        interval_summary.right_support_count,
    )

    # Keep alpha inside the by-center Galerkin contraction as the stride-1
    # term-first index. Do not apply CPB local integral weights and do not sum
    # over centers in this provider block.
    for ix_left in left.x, iy_left in left.y, iz_left in left.z
        row = _local_product_index(
            ix_left - first(left.x) + 1,
            iy_left - first(left.y) + 1,
            iz_left - first(left.z) + 1,
            (left_shape.x, left_shape.y, left_shape.z),
        )
        for ix_right in right.x, iy_right in right.y, iz_right in right.z
            column = _local_product_index(
                ix_right - first(right.x) + 1,
                iy_right - first(right.y) + 1,
                iz_right - first(right.z) + 1,
                (right_shape.x, right_shape.y, right_shape.z),
            )
            value = 0.0
            @inbounds for term in eachindex(coefficients)
                value +=
                    coefficients[term] *
                    axis_terms.x[term, ix_left, ix_right] *
                    axis_terms.y[term, iy_left, iy_right] *
                    axis_terms.z[term, iz_left, iz_right]
            end
            dense_block[row, column] = value
        end
    end
    return dense_block
end

function _cpb_electron_nuclear_by_center_local_block_summary(
    status::Symbol,
    blocker,
    center_summary,
    interval_summary,
    axis_terms,
    expansion,
    dense_block,
)
    available = status === :materialized_cpb_electron_nuclear_by_center_local_block
    coefficients = _summary_property(expansion, :coefficients)
    coefficient_count = coefficients isa AbstractVector ? length(coefficients) : 0
    return (;
        object_kind = :cartesian_cpb_electron_nuclear_by_center_local_block_summary,
        status,
        blocker,
        term = :electron_nuclear_by_center_galerkin,
        source_kind = :parent_axis_bundle_per_center_gaussian_factor_terms,
        center_key = center_summary.center_key,
        center_index = center_summary.center_index,
        charge = center_summary.charge,
        center_location = center_summary.location,
        by_center = true,
        centers_summed = false,
        cpb_integral_weights_applied = false,
        galerkin_operator = true,
        ida_mwg_semantics = false,
        representation = :dense_local_cpb_electron_nuclear_by_center_galerkin,
        factor_source_path = :axis_pgdg_intermediate_gaussian_factor_terms,
        gaussian_expansion_loop = :inner_local_contraction,
        gaussian_term_count = coefficient_count,
        coefficient_count,
        axis_factor_term_shapes = _cpb_pair_factor_term_shapes(axis_terms),
        local_ordering = interval_summary.local_ordering,
        left_shape = available ? interval_summary.left_shape : :unavailable,
        right_shape = available ? interval_summary.right_shape : :unavailable,
        left_support_count =
            available ? interval_summary.left_support_count : :unavailable,
        right_support_count =
            available ? interval_summary.right_support_count : :unavailable,
        dense_block_available = available,
        dense_block_shape = available ? size(dense_block) : :unavailable,
        dense_block_eltype = available ? eltype(dense_block) : :unavailable,
        provider_level_local_matrix_materialized = available,
        provider_level_coulomb_block_materialized = available,
        provider_level_electron_nuclear_block_materialized = available,
        cpb_local_coulomb_kernel_implemented = available,
        cpb_local_electron_nuclear_kernel_implemented = available,
        realization_status = :unrealized,
        route_global_status = :unassigned,
        factor_space =
            available ? :parent_axis_bundle_pgdg_intermediate : :unavailable,
        factor_convention =
            available ?
            :axis_bundle_electron_nuclear_per_center_gaussian_factor_terms :
            :unavailable,
        index_domain =
            available ? :parent_axis_indices : :unavailable,
        index_domain_source =
            available ? :axis_bundle_contract : :unavailable,
        index_domain_status =
            available ?
            :assumed_parent_axis_indexed_by_current_axis_bundle_contract :
            :unavailable,
        axis_order = _AXIS_ORDER,
        bra_ket_order = (:bra, :ket),
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_assembly = false,
        hamiltonian_data_materialized = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function cpb_local_integral_weights(parent_axis_bundle_object, cpb::CPB.CoordinateProductBox)
    axis_weights = _cpb_axis_integral_weights(parent_axis_bundle_object)
    blocker = _cpb_local_integral_weights_blocker(axis_weights, cpb)
    local_weights =
        isnothing(blocker) ?
        _materialize_cpb_local_integral_weights(axis_weights, cpb) :
        nothing
    status =
        isnothing(blocker) ?
        :available_cpb_local_integral_weights :
        :blocked_cpb_local_integral_weights
    return CPBLocalIntegralWeights(
        parent_axis_bundle_object,
        cpb,
        local_weights,
        _cpb_local_integral_weights_summary(
            status,
            blocker,
            axis_weights,
            cpb,
            local_weights,
        ),
    )
end

function cpb_local_integral_weights(axis_bundles; cpb::CPB.CoordinateProductBox)
    return cpb_local_integral_weights(axis_bundles, cpb)
end

function _cpb_axis_integral_weights(parent_axis_bundle_object)
    return (;
        x = _cpb_axis_integral_weights(parent_axis_bundle_object, :x),
        y = _cpb_axis_integral_weights(parent_axis_bundle_object, :y),
        z = _cpb_axis_integral_weights(parent_axis_bundle_object, :z),
    )
end

function _cpb_axis_integral_weights(parent_axis_bundle_object, axis::Symbol)
    axis_bundle = _cpb_axis_bundle(parent_axis_bundle_object, axis)
    pgdg_intermediate = _summary_property(axis_bundle, :pgdg_intermediate)
    return _summary_property(pgdg_intermediate, :weights)
end

function _cpb_local_integral_weights_blocker(axis_weights, cpb::CPB.CoordinateProductBox)
    intervals = CPB.intervals(cpb)
    for (axis_index, axis) in enumerate(_AXIS_ORDER)
        weights = getproperty(axis_weights, axis)
        weights isa AbstractVector{<:Real} ||
            return Symbol("missing_$(axis)_axis_integral_weights")
        isempty(weights) && return Symbol("$(axis)_axis_integral_weights_empty")
        last(intervals[axis_index]) <= length(weights) ||
            return Symbol("$(axis)_axis_integral_weight_length_mismatch")
    end
    return nothing
end

function _materialize_cpb_local_integral_weights(axis_weights, cpb::CPB.CoordinateProductBox)
    intervals = CPB.intervals(cpb)
    shape = CPB.shape(cpb)
    local_weights = Vector{Float64}(undef, CPB.support_count(cpb))
    for ix in intervals[1], iy in intervals[2], iz in intervals[3]
        local_index = _local_product_index(
            ix - first(intervals[1]) + 1,
            iy - first(intervals[2]) + 1,
            iz - first(intervals[3]) + 1,
            shape,
        )
        local_weights[local_index] =
            Float64(axis_weights.x[ix]) *
            Float64(axis_weights.y[iy]) *
            Float64(axis_weights.z[iz])
    end
    return local_weights
end

function _cpb_local_integral_weights_summary(
    status::Symbol,
    blocker,
    axis_weights,
    cpb::CPB.CoordinateProductBox,
    local_weights,
)
    available = status === :available_cpb_local_integral_weights
    return (;
        object_kind = :cartesian_cpb_local_integral_weights_summary,
        status,
        blocker,
        weight_source = :axis_pgdg_intermediate_weights,
        weight_kind = :basis_function_integral_weights,
        squared_self_integral = false,
        retained_pqs_weights = false,
        ida_mwg_semantics = false,
        cpb,
        cpb_role = CPB.role(cpb),
        local_shape = _axis_named_tuple(CPB.shape(cpb)),
        local_ordering = _LOCAL_ORDERING,
        local_weight_count = available ? length(local_weights) : 0,
        axis_weight_counts = _cpb_axis_weight_counts(axis_weights),
        local_weights_available = available,
        local_weights_eltype = available ? eltype(local_weights) : :unavailable,
        route_driver_wiring = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        hamiltonian_assembly = false,
        hamiltonian_data_materialized = false,
        exports_or_artifacts = false,
    )
end

function _cpb_axis_weight_counts(axis_weights)
    return (;
        x = _cpb_axis_weight_count(axis_weights.x),
        y = _cpb_axis_weight_count(axis_weights.y),
        z = _cpb_axis_weight_count(axis_weights.z),
    )
end

function _cpb_axis_weight_count(weights)
    weights isa AbstractVector ? length(weights) : :unavailable
end

function cpb_overlap_dense_block(axis_block_set::CPBOverlapAxisBlockSet)
    axis_block_summary = summary(axis_block_set)
    blocker =
        axis_block_summary.status === :available_cpb_overlap_axis_blocks ?
        nothing :
        _cpb_overlap_dense_block_blocker(axis_block_summary)
    axis_product_block =
        isnothing(blocker) ?
        cpb_axis_product_operator_block(
            axis_block_set.axis_blocks;
            term = :overlap,
            source_kind = :cpb_overlap_axis_blocks,
            source_summary = axis_block_summary,
            factor_space = axis_block_summary.factor_space,
            factor_convention = axis_block_summary.factor_convention,
            normalization_convention = axis_block_summary.normalization_convention,
            index_domain = axis_block_summary.index_domain,
            index_domain_source = axis_block_summary.index_domain_source,
            index_domain_status = axis_block_summary.index_domain_status,
            axis_order = axis_block_summary.axis_order,
            bra_ket_order = axis_block_summary.bra_ket_order,
        ) :
        nothing
    axis_product_summary =
        isnothing(axis_product_block) ? nothing : summary(axis_product_block)
    blocker =
        isnothing(axis_product_summary) ? blocker : axis_product_summary.blocker
    dense_block =
        isnothing(axis_product_block) ? nothing : axis_product_block.dense_block
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
    return CPBLocalOverlapBlockCollection(
        normalized_records,
        collection_summary,
        Ref(nothing),
    )
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
    transform_available = transform_object isa AbstractMatrix
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
    transform_object isa AbstractMatrix ||
        return :unsupported_retained_transform_reference
    source_ordering === _LOCAL_ORDERING ||
        return :retained_transform_ordering_mismatch
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

function cpb_reviewed_overlap_placement_plan(;
    placement_plan_kind = :reviewed_overlap_placement_plan,
    accumulation_rule = nothing,
    symmetry_policy = :explicit_blocks_only,
    duplicate_record_policy = :reject_duplicate_block_keys,
    accepted_block_keys = (),
    required_global_dimension_source = :unavailable,
    local_ordering_contract = _LOCAL_ORDERING,
)
    normalized_block_keys =
        _cpb_reviewed_overlap_placement_plan_block_keys(accepted_block_keys)
    accepted_record_count =
        normalized_block_keys === :unavailable ? 0 : length(normalized_block_keys)
    blocker = _cpb_reviewed_overlap_placement_plan_blocker(
        accumulation_rule,
        symmetry_policy,
        duplicate_record_policy,
        normalized_block_keys,
        required_global_dimension_source,
        local_ordering_contract,
    )
    status =
        isnothing(blocker) ?
        :available_cpb_reviewed_overlap_placement_plan :
        :blocked_cpb_reviewed_overlap_placement_plan
    metadata = (;
        object_kind = :cartesian_cpb_reviewed_overlap_placement_plan_summary,
        status,
        blocker,
        placement_plan_kind,
        accumulation_rule =
            isnothing(accumulation_rule) ? :unavailable : accumulation_rule,
        symmetry_policy,
        duplicate_record_policy,
        accepted_block_keys =
            normalized_block_keys === :unavailable ? () : normalized_block_keys,
        accepted_record_count,
        required_global_dimension_source,
        local_ordering_contract,
        transform_application_implemented = false,
        numerical_placement_implemented = false,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        route_driver_wiring = false,
    )
    return CPBReviewedOverlapPlacementPlan(
        metadata.accepted_block_keys,
        metadata,
    )
end

function _cpb_reviewed_overlap_placement_plan_block_keys(accepted_block_keys)
    isnothing(accepted_block_keys) && return :unavailable
    accepted_block_keys isa Tuple && return accepted_block_keys
    accepted_block_keys isa AbstractVector && return Tuple(accepted_block_keys)
    return (accepted_block_keys,)
end

function _cpb_reviewed_overlap_placement_plan_blocker(
    accumulation_rule,
    symmetry_policy,
    duplicate_record_policy,
    accepted_block_keys,
    required_global_dimension_source,
    local_ordering_contract,
)
    isnothing(accumulation_rule) && return :missing_accumulation_rule
    accepted_block_keys === :unavailable && return :missing_accepted_record_inventory
    isempty(accepted_block_keys) && return :missing_accepted_record_inventory
    symmetry_policy in _OVERLAP_PLACEMENT_SYMMETRY_POLICIES ||
        return :unsupported_overlap_symmetry_policy
    duplicate_record_policy in _OVERLAP_PLACEMENT_DUPLICATE_POLICIES ||
        return :unsupported_overlap_duplicate_record_policy
    required_global_dimension_source === :unavailable &&
        return :missing_required_global_dimension_source
    isnothing(required_global_dimension_source) &&
        return :missing_required_global_dimension_source
    local_ordering_contract === _LOCAL_ORDERING ||
        return :unsupported_overlap_local_ordering_contract
    return nothing
end

function cpb_overlap_placement_facts(
    collection::CPBLocalOverlapBlockCollection;
    transform_carries = (),
    placement_ranges = (),
    placement_plan = nothing,
    accumulation_rule = nothing,
)
    collection_summary = summary(collection)
    collection_available =
        collection_summary.status === :available_cpb_local_overlap_block_collection
    normalized_transform_carries =
        _cpb_overlap_placement_items_tuple(transform_carries)
    normalized_placement_ranges =
        _cpb_overlap_placement_items_tuple(placement_ranges)
    transform_lookup =
        _cpb_overlap_placement_transform_lookup(normalized_transform_carries)
    range_lookup =
        _cpb_overlap_placement_range_lookup(normalized_placement_ranges)
    record_fact_summaries = Tuple(
        _cpb_overlap_placement_record_fact_summary(
            record_summary,
            transform_lookup,
            range_lookup,
        )
        for record_summary in collection_summary.record_summaries
    )
    effective_accumulation_rule =
        _cpb_overlap_placement_effective_accumulation_rule(
            placement_plan,
            accumulation_rule,
        )
    placement_record_inventory =
        _cpb_overlap_placement_record_inventory_summary(
            collection_available,
            collection_summary,
            placement_plan,
        )
    local_ordering_contract =
        _cpb_overlap_placement_local_ordering_contract_summary(
            collection_available,
            record_fact_summaries,
            placement_plan,
        )
    global_dimension_source_contract =
        _cpb_overlap_placement_global_dimension_source_contract_summary(
            collection_available,
            record_fact_summaries,
            placement_plan,
        )
    missing_requirements =
        _cpb_overlap_placement_missing_requirements(
            collection_available,
            record_fact_summaries,
            placement_plan,
            effective_accumulation_rule,
        )
    available_requirements =
        _cpb_overlap_placement_available_requirements(
            collection_available,
            record_fact_summaries,
            placement_plan,
            effective_accumulation_rule,
        )
    propagated_blocker =
        _cpb_overlap_placement_specific_record_blocker(record_fact_summaries)
    placement_plan_status = _cpb_overlap_placement_plan_status(placement_plan)
    placement_plan_blocker = _cpb_overlap_placement_plan_blocker(placement_plan)
    blocker =
        collection_available ?
        (
            isnothing(propagated_blocker) ?
            (
                placement_record_inventory.status ===
                :blocked_cpb_overlap_placement_record_inventory ?
                placement_record_inventory.blocker :
                local_ordering_contract.status ===
                :blocked_cpb_overlap_local_ordering_contract ?
                local_ordering_contract.blocker :
                global_dimension_source_contract.status ===
                :blocked_cpb_overlap_global_dimension_source_contract ?
                global_dimension_source_contract.blocker :
                (
                    isempty(missing_requirements) ?
                    :placement_not_implemented :
                    :missing_placement_or_retained_transform
                )
            ) :
            propagated_blocker
        ) :
        _cpb_overlap_placement_collection_blocker(collection_summary)
    accumulation_rule_status =
        isnothing(effective_accumulation_rule) ?
        :missing_accumulation_rule :
        :available_accumulation_rule
    metadata = (;
        object_kind = :cartesian_cpb_overlap_placement_facts_summary,
        status = :blocked_cpb_overlap_placement_facts,
        blocker,
        collection_available,
        record_count = collection_summary.record_count,
        block_keys = collection_summary.block_keys,
        record_fact_summaries,
        placement_plan_status,
        placement_plan_kind = _cpb_overlap_placement_plan_kind(placement_plan),
        symmetry_policy = _cpb_overlap_placement_plan_symmetry_policy(placement_plan),
        placement_plan_review_status =
            _cpb_overlap_placement_plan_review_status(placement_plan),
        placement_plan_blocker,
        placement_plan_is_reviewed =
            _cpb_overlap_placement_plan_is_reviewed(placement_plan),
        placement_record_inventory_status = placement_record_inventory.status,
        placement_record_inventory_blocker = placement_record_inventory.blocker,
        accepted_block_keys = placement_record_inventory.accepted_block_keys,
        provided_block_keys = placement_record_inventory.provided_block_keys,
        rejected_block_keys = placement_record_inventory.rejected_block_keys,
        duplicate_block_keys = placement_record_inventory.duplicate_block_keys,
        duplicate_record_policy =
            placement_record_inventory.duplicate_record_policy,
        local_ordering_contract_status = local_ordering_contract.status,
        local_ordering_contract_blocker = local_ordering_contract.blocker,
        local_ordering_contract = local_ordering_contract.local_ordering_contract,
        provided_local_orderings = local_ordering_contract.provided_local_orderings,
        mismatched_local_ordering_block_keys =
            local_ordering_contract.mismatched_local_ordering_block_keys,
        global_dimension_source_contract_status =
            global_dimension_source_contract.status,
        global_dimension_source_contract_blocker =
            global_dimension_source_contract.blocker,
        required_global_dimension_source =
            global_dimension_source_contract.required_global_dimension_source,
        provided_global_dimension_sources =
            global_dimension_source_contract.provided_global_dimension_sources,
        mismatched_global_dimension_source_block_keys =
            global_dimension_source_contract.mismatched_global_dimension_source_block_keys,
        accumulation_rule_status,
        accumulation_rule =
            isnothing(effective_accumulation_rule) ?
            :unavailable :
            effective_accumulation_rule,
        available_requirements,
        missing_requirements,
        global_overlap_status = :blocked,
        global_overlap_blocker = blocker,
        placement_engine_implemented = false,
        transform_application_implemented = false,
        global_matrix_materialized = false,
        route_driver_wiring = false,
    )
    return CPBOverlapPlacementFacts(
        collection,
        normalized_transform_carries,
        normalized_placement_ranges,
        metadata,
    )
end

function _cpb_overlap_placement_items_tuple(items)
    isnothing(items) && return ()
    items isa Tuple && return items
    items isa AbstractVector && return Tuple(items)
    return (items,)
end

function _cpb_overlap_placement_transform_lookup(transform_carries::Tuple)
    lookup = Dict{Any,Any}()
    for transform_carry in transform_carries
        transform_summary =
            _cpb_overlap_placement_transform_summary(transform_carry)
        isnothing(transform_summary) && continue
        block_key = _summary_property(transform_summary, :block_key)
        side = _summary_property(transform_summary, :side)
        isnothing(block_key) && continue
        side in (:left, :right) || continue
        lookup[(block_key, side)] = transform_summary
    end
    return lookup
end

function _cpb_overlap_placement_range_lookup(placement_ranges::Tuple)
    lookup = Dict{Any,Any}()
    for placement_range in placement_ranges
        range_summary = _cpb_overlap_placement_range_summary(placement_range)
        isnothing(range_summary) && continue
        block_key = _summary_property(range_summary, :block_key)
        isnothing(block_key) && continue
        lookup[block_key] = range_summary
    end
    return lookup
end

function _cpb_overlap_placement_transform_summary(transform_carry)
    transform_carry isa CPBRetainedTransformCarry && return summary(transform_carry)
    transform_carry isa NamedTuple && return transform_carry
    return nothing
end

function _cpb_overlap_placement_range_summary(placement_range)
    placement_range isa CPBSourcePairPlacementRange && return summary(placement_range)
    placement_range isa NamedTuple && return placement_range
    return nothing
end

function _cpb_overlap_placement_record_fact_summary(
    record_summary,
    transform_lookup,
    range_lookup,
)
    block_key = record_summary.block_key
    left_transform_summary = get(transform_lookup, (block_key, :left), nothing)
    right_transform_summary = get(transform_lookup, (block_key, :right), nothing)
    placement_range_summary = get(range_lookup, block_key, nothing)
    left_transform_status =
        _cpb_overlap_placement_transform_fact_status(left_transform_summary)
    right_transform_status =
        _cpb_overlap_placement_transform_fact_status(right_transform_summary)
    placement_range_status =
        _cpb_overlap_placement_range_fact_status(placement_range_summary)
    return (;
        block_key,
        record_status = record_summary.status,
        record_blocker = record_summary.blocker,
        left_transform_status,
        left_transform_blocker =
            _cpb_overlap_placement_fact_blocker(left_transform_summary),
        right_transform_status,
        right_transform_blocker =
            _cpb_overlap_placement_fact_blocker(right_transform_summary),
        placement_range_status,
        placement_range_blocker =
            _cpb_overlap_placement_fact_blocker(placement_range_summary),
        left_column_range =
            _cpb_overlap_placement_range_property(
                placement_range_summary,
                :left_column_range,
            ),
        right_column_range =
            _cpb_overlap_placement_range_property(
                placement_range_summary,
                :right_column_range,
            ),
        global_dimension =
            _cpb_overlap_placement_range_property(
                placement_range_summary,
                :global_dimension,
            ),
        global_dimension_source =
            _cpb_overlap_placement_range_property(
                placement_range_summary,
                :global_dimension_source,
                :unavailable,
            ),
        dense_block_available = record_summary.dense_block_available,
        dense_block_shape = record_summary.dense_block_shape,
        left_cpb_summary =
            _cpb_overlap_placement_record_property(
                record_summary,
                :left_cpb_summary,
                :unavailable,
            ),
        right_cpb_summary =
            _cpb_overlap_placement_record_property(
                record_summary,
                :right_cpb_summary,
                :unavailable,
            ),
        local_ordering = record_summary.local_ordering,
    )
end

function _cpb_overlap_placement_transform_fact_status(transform_summary)
    isnothing(transform_summary) && return :missing_retained_transform
    status = _summary_property(transform_summary, :status)
    isnothing(status) ? :unavailable_retained_transform : status
end

function _cpb_overlap_placement_range_fact_status(range_summary)
    isnothing(range_summary) && return :missing_source_pair_placement_range
    status = _summary_property(range_summary, :status)
    isnothing(status) ? :unavailable_source_pair_placement_range : status
end

function _cpb_overlap_placement_fact_blocker(fact_summary)
    isnothing(fact_summary) && return nothing
    return _summary_property(fact_summary, :blocker)
end

function _cpb_overlap_placement_range_property(
    range_summary,
    property::Symbol,
    default = nothing,
)
    isnothing(range_summary) && return default
    value = _summary_property(range_summary, property)
    isnothing(value) ? default : value
end

function _cpb_overlap_placement_record_property(
    record_summary,
    property::Symbol,
    default = nothing,
)
    value = _summary_property(record_summary, property)
    isnothing(value) ? default : value
end

function _cpb_overlap_placement_collection_blocker(collection_summary)
    blocker = _summary_property(collection_summary, :blocker)
    isnothing(blocker) ? :missing_local_overlap_collection : blocker
end

function _cpb_overlap_placement_record_inventory_summary(
    collection_available::Bool,
    collection_summary,
    placement_plan,
)
    provided_block_keys =
        collection_available ? collection_summary.block_keys : ()
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    plan_is_reviewed =
        _cpb_overlap_placement_plan_is_reviewed(placement_plan, plan_summary)
    plan_status =
        isnothing(plan_summary) ?
        :unavailable :
        _summary_property(plan_summary, :status)
    duplicate_record_policy =
        isnothing(plan_summary) ?
        :unavailable :
        _summary_property(
            plan_summary,
            :duplicate_record_policy,
        )
    accepted_block_keys =
        plan_is_reviewed && !isnothing(plan_summary) ?
        _summary_property(plan_summary, :accepted_block_keys) :
        ()
    isnothing(accepted_block_keys) && (accepted_block_keys = ())
    duplicate_block_keys =
        _cpb_overlap_placement_duplicate_block_keys(provided_block_keys)
    rejected_block_keys =
        plan_is_reviewed &&
        plan_status === :available_cpb_reviewed_overlap_placement_plan ?
        _cpb_overlap_placement_rejected_block_keys(
            provided_block_keys,
            accepted_block_keys,
        ) :
        ()
    status, blocker =
        !collection_available ?
        (:unavailable_cpb_overlap_placement_record_inventory,
            :missing_local_overlap_collection) :
        !plan_is_reviewed ?
        (:not_checked_cpb_overlap_placement_record_inventory,
            nothing) :
        plan_status !== :available_cpb_reviewed_overlap_placement_plan ?
        (:blocked_cpb_overlap_placement_record_inventory,
            _summary_property(plan_summary, :blocker)) :
        !isempty(duplicate_block_keys) &&
        duplicate_record_policy === :reject_duplicate_block_keys ?
        (:blocked_cpb_overlap_placement_record_inventory,
            :duplicate_overlap_placement_record) :
        !isempty(rejected_block_keys) ?
        (:blocked_cpb_overlap_placement_record_inventory,
            :unaccepted_overlap_placement_record) :
        (:available_cpb_overlap_placement_record_inventory,
            nothing)
    return (;
        status,
        blocker,
        accepted_block_keys,
        provided_block_keys,
        rejected_block_keys,
        duplicate_block_keys,
        duplicate_record_policy,
    )
end

function _cpb_overlap_placement_rejected_block_keys(
    provided_block_keys::Tuple,
    accepted_block_keys::Tuple,
)
    return Tuple(
        block_key for block_key in unique(provided_block_keys)
        if !(block_key in accepted_block_keys)
    )
end

function _cpb_overlap_placement_duplicate_block_keys(block_keys::Tuple)
    duplicate_block_keys = Any[]
    seen_block_keys = Set{Any}()
    pushed_block_keys = Set{Any}()
    for block_key in block_keys
        if block_key in seen_block_keys
            if !(block_key in pushed_block_keys)
                push!(duplicate_block_keys, block_key)
                push!(pushed_block_keys, block_key)
            end
        else
            push!(seen_block_keys, block_key)
        end
    end
    return Tuple(duplicate_block_keys)
end

function _cpb_overlap_placement_local_ordering_contract_summary(
    collection_available::Bool,
    record_fact_summaries::Tuple,
    placement_plan,
)
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    plan_is_reviewed =
        _cpb_overlap_placement_plan_is_reviewed(placement_plan, plan_summary)
    plan_status =
        isnothing(plan_summary) ?
        :unavailable :
        _summary_property(plan_summary, :status)
    local_ordering_contract =
        plan_is_reviewed && !isnothing(plan_summary) ?
        _summary_property(plan_summary, :local_ordering_contract) :
        :unavailable
    provided_local_orderings =
        collection_available ?
        Tuple(unique(record.local_ordering for record in record_fact_summaries)) :
        ()
    mismatched_block_keys =
        plan_is_reviewed &&
        plan_status === :available_cpb_reviewed_overlap_placement_plan ?
        Tuple(
            record.block_key for record in record_fact_summaries
            if record.record_status === :available_cpb_local_overlap_block_record &&
               record.local_ordering !== local_ordering_contract
        ) :
        ()
    status, blocker =
        !collection_available ?
        (:unavailable_cpb_overlap_local_ordering_contract,
            :missing_local_overlap_collection) :
        !plan_is_reviewed ?
        (:not_checked_cpb_overlap_local_ordering_contract, nothing) :
        plan_status !== :available_cpb_reviewed_overlap_placement_plan ?
        (:blocked_cpb_overlap_local_ordering_contract,
            _summary_property(plan_summary, :blocker)) :
        !isempty(mismatched_block_keys) ?
        (:blocked_cpb_overlap_local_ordering_contract,
            :overlap_local_ordering_contract_mismatch) :
        (:available_cpb_overlap_local_ordering_contract, nothing)
    return (;
        status,
        blocker,
        local_ordering_contract,
        provided_local_orderings,
        mismatched_local_ordering_block_keys = mismatched_block_keys,
    )
end

function _cpb_overlap_placement_global_dimension_source_contract_summary(
    collection_available::Bool,
    record_fact_summaries::Tuple,
    placement_plan,
)
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    plan_is_reviewed =
        _cpb_overlap_placement_plan_is_reviewed(placement_plan, plan_summary)
    plan_status =
        isnothing(plan_summary) ?
        :unavailable :
        _summary_property(plan_summary, :status)
    required_global_dimension_source =
        plan_is_reviewed && !isnothing(plan_summary) ?
        _summary_property(plan_summary, :required_global_dimension_source) :
        :unavailable
    available_range_records = Tuple(
        record for record in record_fact_summaries
        if record.placement_range_status ===
           :available_cpb_source_pair_placement_range
    )
    provided_global_dimension_sources =
        collection_available ?
        Tuple(unique(record.global_dimension_source for record in available_range_records)) :
        ()
    mismatched_block_keys =
        plan_is_reviewed &&
        plan_status === :available_cpb_reviewed_overlap_placement_plan ?
        Tuple(
            record.block_key for record in available_range_records
            if record.global_dimension_source !== required_global_dimension_source
        ) :
        ()
    status, blocker =
        !collection_available ?
        (:unavailable_cpb_overlap_global_dimension_source_contract,
            :missing_local_overlap_collection) :
        !plan_is_reviewed ?
        (:not_checked_cpb_overlap_global_dimension_source_contract, nothing) :
        plan_status !== :available_cpb_reviewed_overlap_placement_plan ?
        (:blocked_cpb_overlap_global_dimension_source_contract,
            _summary_property(plan_summary, :blocker)) :
        isempty(available_range_records) ?
        (:not_checked_cpb_overlap_global_dimension_source_contract, nothing) :
        !isempty(mismatched_block_keys) ?
        (:blocked_cpb_overlap_global_dimension_source_contract,
            :global_dimension_source_mismatch) :
        (:available_cpb_overlap_global_dimension_source_contract, nothing)
    return (;
        status,
        blocker,
        required_global_dimension_source,
        provided_global_dimension_sources,
        mismatched_global_dimension_source_block_keys = mismatched_block_keys,
    )
end

function _cpb_overlap_placement_missing_requirements(
    collection_available::Bool,
    record_fact_summaries::Tuple,
    placement_plan,
    accumulation_rule,
)
    missing = Symbol[]
    if !collection_available
        push!(missing, :missing_local_overlap_collection)
    else
        if any(_cpb_overlap_placement_missing_transform, record_fact_summaries)
            push!(missing, :missing_retained_transform)
        end
        if any(_cpb_overlap_placement_missing_left_range, record_fact_summaries)
            push!(missing, :missing_left_column_range)
        end
        if any(_cpb_overlap_placement_missing_right_range, record_fact_summaries)
            push!(missing, :missing_right_column_range)
        end
        if any(_cpb_overlap_placement_missing_global_dimension, record_fact_summaries)
            push!(missing, :missing_global_dimension)
        end
    end
    isnothing(placement_plan) && push!(missing, :missing_placement_plan)
    isnothing(accumulation_rule) && push!(missing, :missing_accumulation_rule)
    return Tuple(unique(missing))
end

function _cpb_overlap_placement_available_requirements(
    collection_available::Bool,
    record_fact_summaries::Tuple,
    placement_plan,
    accumulation_rule,
)
    available = Symbol[]
    collection_available && push!(available, :local_cpb_overlap_collection)
    if collection_available && !isempty(record_fact_summaries)
        all(_cpb_overlap_placement_available_transforms, record_fact_summaries) &&
            push!(available, :retained_transform)
        all(_cpb_overlap_placement_available_left_range, record_fact_summaries) &&
            push!(available, :left_column_range)
        all(_cpb_overlap_placement_available_right_range, record_fact_summaries) &&
            push!(available, :right_column_range)
        all(_cpb_overlap_placement_available_global_dimension, record_fact_summaries) &&
            push!(available, :global_dimension)
    end
    isnothing(placement_plan) || push!(available, :placement_plan)
    isnothing(accumulation_rule) || push!(available, :accumulation_rule)
    return Tuple(available)
end

function _cpb_overlap_placement_missing_transform(record_fact_summary)
    return record_fact_summary.left_transform_status === :missing_retained_transform ||
           record_fact_summary.right_transform_status === :missing_retained_transform ||
           record_fact_summary.left_transform_blocker === :missing_retained_transform ||
           record_fact_summary.right_transform_blocker === :missing_retained_transform
end

function _cpb_overlap_placement_missing_left_range(record_fact_summary)
    return record_fact_summary.placement_range_status ===
           :missing_source_pair_placement_range ||
           record_fact_summary.placement_range_blocker ===
           :missing_left_column_range
end

function _cpb_overlap_placement_missing_right_range(record_fact_summary)
    return record_fact_summary.placement_range_status ===
           :missing_source_pair_placement_range ||
           record_fact_summary.placement_range_blocker ===
           :missing_right_column_range
end

function _cpb_overlap_placement_missing_global_dimension(record_fact_summary)
    return record_fact_summary.placement_range_status ===
           :missing_source_pair_placement_range ||
           record_fact_summary.placement_range_blocker === :missing_global_dimension
end

function _cpb_overlap_placement_available_transforms(record_fact_summary)
    return record_fact_summary.left_transform_status ===
           :available_cpb_retained_transform_carry &&
           record_fact_summary.right_transform_status ===
           :available_cpb_retained_transform_carry
end

function _cpb_overlap_placement_available_left_range(record_fact_summary)
    return record_fact_summary.placement_range_status ===
           :available_cpb_source_pair_placement_range &&
           !isnothing(record_fact_summary.left_column_range)
end

function _cpb_overlap_placement_available_right_range(record_fact_summary)
    return record_fact_summary.placement_range_status ===
           :available_cpb_source_pair_placement_range &&
           !isnothing(record_fact_summary.right_column_range)
end

function _cpb_overlap_placement_available_global_dimension(record_fact_summary)
    return record_fact_summary.placement_range_status ===
           :available_cpb_source_pair_placement_range &&
           record_fact_summary.global_dimension isa Integer
end

function _cpb_overlap_placement_specific_record_blocker(record_fact_summaries::Tuple)
    for record_fact_summary in record_fact_summaries
        record_fact_summary.record_status ===
        :available_cpb_local_overlap_block_record ||
            return something(
                record_fact_summary.record_blocker,
                :blocked_cpb_local_overlap_block_record,
            )
        for blocker in (
            record_fact_summary.left_transform_blocker,
            record_fact_summary.right_transform_blocker,
            record_fact_summary.placement_range_blocker,
        )
            isnothing(blocker) && continue
            _cpb_overlap_placement_missing_blocker(blocker) && continue
            return blocker
        end
    end
    return nothing
end

function _cpb_overlap_placement_missing_blocker(blocker::Symbol)
    return blocker in (
        :missing_retained_transform,
        :missing_left_column_range,
        :missing_right_column_range,
        :missing_global_dimension,
        :missing_source_pair_placement_range,
    )
end

function _cpb_overlap_placement_plan_summary(placement_plan)
    placement_plan isa CPBReviewedOverlapPlacementPlan && return summary(placement_plan)
    placement_plan isa NamedTuple && return placement_plan
    return nothing
end

function _cpb_overlap_placement_plan_is_reviewed(
    placement_plan,
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan),
)
    placement_plan isa CPBReviewedOverlapPlacementPlan && return true
    isnothing(plan_summary) && return false
    _summary_property(
        plan_summary,
        :object_kind,
    ) === :cartesian_cpb_reviewed_overlap_placement_plan_summary && return true
    status = _summary_property(plan_summary, :status)
    return status in (
        :available_cpb_reviewed_overlap_placement_plan,
        :blocked_cpb_reviewed_overlap_placement_plan,
    )
end

function _cpb_overlap_placement_plan_status(placement_plan)
    isnothing(placement_plan) && return :missing_placement_plan
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    if !isnothing(plan_summary)
        status = _summary_property(plan_summary, :status)
        status === :available_cpb_reviewed_overlap_placement_plan &&
            return :available_placement_plan
        status === :blocked_cpb_reviewed_overlap_placement_plan &&
            return :blocked_placement_plan
    end
    return :available_placement_plan
end

function _cpb_overlap_placement_plan_review_status(placement_plan)
    isnothing(placement_plan) && return :missing_placement_plan
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    _cpb_overlap_placement_plan_is_reviewed(placement_plan, plan_summary) ||
        return :placeholder_placement_plan_compatibility
    plan_status =
        isnothing(plan_summary) ?
        :unavailable :
        _summary_property(plan_summary, :status)
    plan_status === :available_cpb_reviewed_overlap_placement_plan &&
        return :reviewed_placement_plan
    plan_status === :blocked_cpb_reviewed_overlap_placement_plan &&
        return :blocked_reviewed_placement_plan
    return :blocked_reviewed_placement_plan
end

function _cpb_overlap_placement_plan_blocker(placement_plan)
    isnothing(placement_plan) && return :missing_placement_plan
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    _cpb_overlap_placement_plan_is_reviewed(placement_plan, plan_summary) ||
        return nothing
    blocker =
        isnothing(plan_summary) ?
        nothing :
        _summary_property(plan_summary, :blocker)
    return blocker
end

function _cpb_overlap_placement_effective_accumulation_rule(
    placement_plan,
    accumulation_rule,
)
    isnothing(accumulation_rule) || return accumulation_rule
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    isnothing(plan_summary) && return nothing
    plan_accumulation_rule =
        _summary_property(plan_summary, :accumulation_rule)
    plan_accumulation_rule === :unavailable && return nothing
    return plan_accumulation_rule
end

function _cpb_overlap_placement_plan_kind(placement_plan)
    isnothing(placement_plan) && return :unavailable
    placement_plan isa Symbol && return placement_plan
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    if !isnothing(plan_summary)
        placement_plan_kind =
            _summary_property(plan_summary, :placement_plan_kind)
        isnothing(placement_plan_kind) || return placement_plan_kind
    end
    if placement_plan isa NamedTuple
        kind = _summary_property(placement_plan, :kind)
        isnothing(kind) || return kind
    end
    return Symbol(nameof(typeof(placement_plan)))
end

function _cpb_overlap_placement_plan_symmetry_policy(placement_plan)
    isnothing(placement_plan) && return :unavailable
    plan_summary = _cpb_overlap_placement_plan_summary(placement_plan)
    isnothing(plan_summary) && return :unavailable
    symmetry_policy = _summary_property(plan_summary, :symmetry_policy)
    isnothing(symmetry_policy) ? :unavailable : symmetry_policy
end

function cpb_place_overlap_block(
    source_block,
    left_transform_carry,
    right_transform_carry,
    placement_range,
    placement_facts,
)
    source_summary = _cpb_overlap_pilot_summary(source_block)
    left_transform_summary = _cpb_overlap_placement_transform_summary(
        left_transform_carry,
    )
    right_transform_summary = _cpb_overlap_placement_transform_summary(
        right_transform_carry,
    )
    range_summary = _cpb_overlap_placement_range_summary(placement_range)
    facts_summary = _cpb_overlap_pilot_facts_summary(placement_facts)
    blocker = _cpb_place_overlap_block_blocker(
        source_block,
        source_summary,
        left_transform_carry,
        left_transform_summary,
        right_transform_carry,
        right_transform_summary,
        placement_range,
        range_summary,
        placement_facts,
        facts_summary,
    )
    global_overlap_matrix =
        isnothing(blocker) ?
        _cpb_place_overlap_block_matrix(
            source_block,
            left_transform_carry,
            right_transform_carry,
            range_summary,
        ) :
        nothing
    status =
        isnothing(blocker) ?
        :materialized_cpb_overlap_placement_pilot :
        :blocked_cpb_overlap_placement_pilot
    return CPBPlacedOverlapBlock(
        source_block,
        left_transform_carry,
        right_transform_carry,
        placement_range,
        placement_facts,
        global_overlap_matrix,
        _cpb_place_overlap_block_summary(
            status,
            blocker,
            source_block,
            source_summary,
            left_transform_summary,
            right_transform_summary,
            range_summary,
            facts_summary,
            global_overlap_matrix,
        ),
    )
end

function _cpb_overlap_pilot_summary(object)
    (
        object isa CPBOverlapDenseBlock ||
        object isa CPBOverlapAxisBlockSet ||
        object isa CPBLocalOverlapBlockRecord ||
        object isa CPBLocalOverlapBlockCollection
    ) && return summary(object)
    return nothing
end

function _cpb_overlap_pilot_facts_summary(placement_facts)
    placement_facts isa CPBOverlapPlacementFacts && return summary(placement_facts)
    return nothing
end

function _cpb_place_overlap_block_blocker(
    source_block,
    _source_summary,
    left_transform_carry,
    left_transform_summary,
    right_transform_carry,
    right_transform_summary,
    placement_range,
    range_summary,
    placement_facts,
    facts_summary,
)
    placement_facts isa CPBOverlapPlacementFacts ||
        return :missing_reviewed_overlap_placement_facts
    source_block isa CPBOverlapDenseBlock ||
        return :local_overlap_source_not_dense
    facts_summary.placement_plan_review_status === :reviewed_placement_plan ||
        return :placement_facts_not_reviewed
    isempty(facts_summary.missing_requirements) ||
        return :placement_requirements_missing
    facts_summary.blocker === :placement_not_implemented ||
        return facts_summary.blocker
    isnothing(left_transform_summary) && return :retained_transform_unavailable
    isnothing(right_transform_summary) && return :retained_transform_unavailable
    isnothing(range_summary) && return :placement_range_unavailable
    left_transform_summary.status === :available_cpb_retained_transform_carry ||
        return :retained_transform_unavailable
    right_transform_summary.status === :available_cpb_retained_transform_carry ||
        return :retained_transform_unavailable
    left_transform_carry.transform_object isa AbstractMatrix ||
        return :retained_transform_unavailable
    right_transform_carry.transform_object isa AbstractMatrix ||
        return :retained_transform_unavailable
    range_summary.status === :available_cpb_source_pair_placement_range ||
        return :placement_range_unavailable
    facts_summary.accumulation_rule === :add_explicit_blocks_into_ranges ||
        return :unsupported_overlap_accumulation_rule
    facts_summary.symmetry_policy === :explicit_blocks_only ||
        return :unsupported_overlap_symmetry_policy
    facts_summary.duplicate_record_policy === :reject_duplicate_block_keys ||
        return :unsupported_overlap_duplicate_record_policy
    dense_shape = size(source_block.dense_block)
    size(left_transform_carry.transform_object, 1) == dense_shape[1] ||
        return :retained_transform_shape_mismatch
    size(right_transform_carry.transform_object, 1) == dense_shape[2] ||
        return :retained_transform_shape_mismatch
    size(left_transform_carry.transform_object, 2) == range_summary.left_column_count ||
        return :retained_transform_shape_mismatch
    size(right_transform_carry.transform_object, 2) == range_summary.right_column_count ||
        return :retained_transform_shape_mismatch
    return nothing
end

function _cpb_place_overlap_block_matrix(
    source_block::CPBOverlapDenseBlock,
    left_transform_carry::CPBRetainedTransformCarry,
    right_transform_carry::CPBRetainedTransformCarry,
    range_summary,
)
    retained_block =
        left_transform_carry.transform_object' *
        source_block.dense_block *
        right_transform_carry.transform_object
    element_type = promote_type(
        eltype(retained_block),
        eltype(source_block.dense_block),
        eltype(left_transform_carry.transform_object),
        eltype(right_transform_carry.transform_object),
    )
    global_overlap_matrix = zeros(
        element_type,
        range_summary.global_dimension,
        range_summary.global_dimension,
    )
    global_overlap_matrix[
        range_summary.left_column_range,
        range_summary.right_column_range,
    ] .+= retained_block
    return global_overlap_matrix
end

function _cpb_place_overlap_block_summary(
    status::Symbol,
    blocker,
    source_block,
    source_summary,
    left_transform_summary,
    right_transform_summary,
    range_summary,
    facts_summary,
    global_overlap_matrix,
)
    available = status === :materialized_cpb_overlap_placement_pilot
    retained_block_shape =
        available ?
        (
            left_transform_summary.target_retained_column_count,
            right_transform_summary.target_retained_column_count,
        ) :
        :not_materialized
    return (;
        object_kind = :cartesian_cpb_placed_overlap_block_summary,
        status,
        blocker,
        source_kind =
            source_block isa CPBOverlapDenseBlock ?
            :cpb_overlap_dense_block :
            :unavailable,
        source_block_status =
            isnothing(source_summary) ? :unavailable : source_summary.status,
        source_dense_block_shape =
            source_block isa CPBOverlapDenseBlock ?
            size(source_block.dense_block) :
            :unavailable,
        left_transform_status =
            isnothing(left_transform_summary) ?
            :unavailable :
            left_transform_summary.status,
        right_transform_status =
            isnothing(right_transform_summary) ?
            :unavailable :
            right_transform_summary.status,
        placement_range_status =
            isnothing(range_summary) ? :unavailable : range_summary.status,
        placement_facts_status =
            isnothing(facts_summary) ? :unavailable : facts_summary.status,
        placement_facts_blocker =
            isnothing(facts_summary) ? :unavailable : facts_summary.blocker,
        placement_plan_review_status =
            isnothing(facts_summary) ?
            :unavailable :
            facts_summary.placement_plan_review_status,
        accumulation_rule =
            isnothing(facts_summary) ? :unavailable : facts_summary.accumulation_rule,
        symmetry_policy =
            isnothing(facts_summary) ? :unavailable : facts_summary.symmetry_policy,
        duplicate_record_policy =
            isnothing(facts_summary) ? :unavailable : facts_summary.duplicate_record_policy,
        block_key =
            isnothing(range_summary) ? :unavailable : range_summary.block_key,
        left_column_range =
            isnothing(range_summary) ? :unavailable : range_summary.left_column_range,
        right_column_range =
            isnothing(range_summary) ? :unavailable : range_summary.right_column_range,
        global_dimension =
            isnothing(range_summary) ? :unavailable : range_summary.global_dimension,
        global_dimension_source =
            isnothing(range_summary) ? :unavailable : range_summary.global_dimension_source,
        retained_block_shape,
        provider_level_global_overlap_matrix_shape =
            available ? size(global_overlap_matrix) : :not_materialized,
        global_overlap_matrix_shape = :route_global_matrix_not_materialized,
        provider_level_matrix_materialized = available,
        provider_level_overlap_matrix_materialized = available,
        provider_level_pilot = true,
        synthetic_fixture_only = true,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        route_driver_wiring = false,
        route_global_overlap_stage_source = false,
        route_global_overlap_available = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

function cpb_place_overlap_collection(
    collection::CPBLocalOverlapBlockCollection,
    transform_carries,
    placement_ranges,
    placement_facts,
)
    return cpb_place_overlap_collection(
        collection;
        transform_carries,
        placement_ranges,
        placement_facts,
    )
end

function cpb_place_overlap_collection(
    collection::CPBLocalOverlapBlockCollection;
    transform_carries = (),
    placement_ranges = (),
    placement_facts = nothing,
)
    normalized_transform_carries =
        _cpb_overlap_placement_items_tuple(transform_carries)
    normalized_placement_ranges =
        _cpb_overlap_placement_items_tuple(placement_ranges)
    transform_lookup =
        _cpb_overlap_placement_transform_object_lookup(normalized_transform_carries)
    range_lookup =
        _cpb_overlap_placement_range_object_lookup(normalized_placement_ranges)
    facts_summary = _cpb_overlap_pilot_facts_summary(placement_facts)
    collection_summary = summary(collection)
    record_placement_summaries = Tuple(
        _cpb_place_overlap_collection_record_summary(
            record,
            transform_lookup,
            range_lookup,
            placement_facts,
            facts_summary,
        )
        for record in collection.records
    )
    blocker = _cpb_place_overlap_collection_blocker(
        collection,
        collection_summary,
        placement_facts,
        facts_summary,
        record_placement_summaries,
    )
    global_overlap_matrix =
        isnothing(blocker) ?
        _cpb_place_overlap_collection_matrix(
            collection,
            transform_lookup,
            range_lookup,
        ) :
        nothing
    status =
        isnothing(blocker) ?
        :materialized_cpb_overlap_collection_placement_pilot :
        :blocked_cpb_overlap_collection_placement_pilot
    return CPBPlacedOverlapCollection(
        collection,
        normalized_transform_carries,
        normalized_placement_ranges,
        placement_facts,
        global_overlap_matrix,
        _cpb_place_overlap_collection_summary(
            status,
            blocker,
            collection_summary,
            facts_summary,
            record_placement_summaries,
            global_overlap_matrix,
        ),
    )
end

function _cpb_overlap_placement_transform_object_lookup(transform_carries::Tuple)
    lookup = Dict{Any,Any}()
    for transform_carry in transform_carries
        transform_summary =
            _cpb_overlap_placement_transform_summary(transform_carry)
        isnothing(transform_summary) && continue
        block_key = _summary_property(transform_summary, :block_key)
        side = _summary_property(transform_summary, :side)
        isnothing(block_key) && continue
        side in (:left, :right) || continue
        lookup[(block_key, side)] = transform_carry
    end
    return lookup
end

function _cpb_overlap_placement_range_object_lookup(placement_ranges::Tuple)
    lookup = Dict{Any,Any}()
    for placement_range in placement_ranges
        range_summary = _cpb_overlap_placement_range_summary(placement_range)
        isnothing(range_summary) && continue
        block_key = _summary_property(range_summary, :block_key)
        isnothing(block_key) && continue
        lookup[block_key] = placement_range
    end
    return lookup
end

function _cpb_place_overlap_collection_record_summary(
    record::CPBLocalOverlapBlockRecord,
    transform_lookup,
    range_lookup,
    placement_facts,
    facts_summary,
)
    record_summary = summary(record)
    block_key = record.block_key
    left_transform_carry = get(transform_lookup, (block_key, :left), nothing)
    right_transform_carry = get(transform_lookup, (block_key, :right), nothing)
    placement_range = get(range_lookup, block_key, nothing)
    left_transform_summary =
        _cpb_overlap_placement_transform_summary(left_transform_carry)
    right_transform_summary =
        _cpb_overlap_placement_transform_summary(right_transform_carry)
    range_summary = _cpb_overlap_placement_range_summary(placement_range)
    source_summary = _cpb_overlap_pilot_summary(record.source_block)
    blocker = _cpb_place_overlap_block_blocker(
        record.source_block,
        source_summary,
        left_transform_carry,
        left_transform_summary,
        right_transform_carry,
        right_transform_summary,
        placement_range,
        range_summary,
        placement_facts,
        facts_summary,
    )
    status =
        isnothing(blocker) ?
        :materialized_cpb_overlap_collection_record_pilot :
        :blocked_cpb_overlap_collection_record_pilot
    return (;
        block_key,
        status,
        blocker,
        source_kind =
            record.source_block isa CPBOverlapDenseBlock ?
            :cpb_overlap_dense_block :
            :unavailable,
        source_block_status =
            isnothing(source_summary) ? :unavailable : source_summary.status,
        source_dense_block_shape =
            record.source_block isa CPBOverlapDenseBlock ?
            size(record.source_block.dense_block) :
            :unavailable,
        left_transform_status =
            isnothing(left_transform_summary) ?
            :unavailable :
            left_transform_summary.status,
        right_transform_status =
            isnothing(right_transform_summary) ?
            :unavailable :
            right_transform_summary.status,
        placement_range_status =
            isnothing(range_summary) ? :unavailable : range_summary.status,
        left_column_range =
            isnothing(range_summary) ? :unavailable : range_summary.left_column_range,
        right_column_range =
            isnothing(range_summary) ? :unavailable : range_summary.right_column_range,
        global_dimension =
            isnothing(range_summary) ? :unavailable : range_summary.global_dimension,
        local_ordering =
            _summary_property(record_summary, :local_ordering),
    )
end

function _cpb_place_overlap_collection_blocker(
    collection::CPBLocalOverlapBlockCollection,
    collection_summary,
    placement_facts,
    facts_summary,
    record_placement_summaries::Tuple,
)
    placement_facts isa CPBOverlapPlacementFacts ||
        return :missing_reviewed_overlap_placement_facts
    collection_summary.status === :available_cpb_local_overlap_block_collection ||
        return _cpb_overlap_placement_collection_blocker(collection_summary)
    _cpb_place_overlap_collection_facts_aligned(
        collection,
        collection_summary,
        placement_facts,
        facts_summary,
    ) || return :placement_facts_collection_mismatch
    facts_summary.placement_plan_review_status === :reviewed_placement_plan ||
        return :placement_facts_not_reviewed
    coverage = _cpb_place_overlap_collection_record_coverage(
        record_placement_summaries,
    )
    (
        isempty(coverage.missing_left_transform_block_keys) &&
        isempty(coverage.missing_right_transform_block_keys)
    ) || return :missing_retained_transform
    isempty(coverage.missing_placement_range_block_keys) ||
        return :missing_placement_range
    isempty(facts_summary.missing_requirements) ||
        return :placement_requirements_missing
    facts_summary.blocker === :placement_not_implemented ||
        return facts_summary.blocker
    facts_summary.accumulation_rule === :add_explicit_blocks_into_ranges ||
        return :unsupported_overlap_accumulation_rule
    facts_summary.symmetry_policy === :explicit_blocks_only ||
        return :unsupported_overlap_symmetry_policy
    facts_summary.duplicate_record_policy === :reject_duplicate_block_keys ||
        return :unsupported_overlap_duplicate_record_policy
    dimension_blocker =
        _cpb_place_overlap_collection_global_dimension_blocker(
            record_placement_summaries,
        )
    isnothing(dimension_blocker) || return dimension_blocker
    for record_summary in record_placement_summaries
        record_summary.status === :materialized_cpb_overlap_collection_record_pilot ||
            return record_summary.blocker
    end
    return nothing
end

function _cpb_place_overlap_collection_facts_aligned(
    collection::CPBLocalOverlapBlockCollection,
    collection_summary,
    placement_facts::CPBOverlapPlacementFacts,
    facts_summary,
)
    placement_facts.collection.identity_token === collection.identity_token ||
        return false
    facts_block_keys = _summary_property(facts_summary, :block_keys)
    isnothing(facts_block_keys) && return false
    return facts_block_keys == collection_summary.block_keys
end

function _cpb_place_overlap_collection_record_coverage(
    record_placement_summaries::Tuple,
)
    return (;
        missing_left_transform_block_keys = Tuple(
            record.block_key for record in record_placement_summaries
            if record.left_transform_status === :unavailable
        ),
        missing_right_transform_block_keys = Tuple(
            record.block_key for record in record_placement_summaries
            if record.right_transform_status === :unavailable
        ),
        missing_placement_range_block_keys = Tuple(
            record.block_key for record in record_placement_summaries
            if record.placement_range_status === :unavailable
        ),
    )
end

function _cpb_place_overlap_collection_global_dimension_blocker(
    record_placement_summaries::Tuple,
)
    global_dimensions = Tuple(
        record.global_dimension for record in record_placement_summaries
        if record.global_dimension isa Integer
    )
    isempty(global_dimensions) && return nothing
    all(==(first(global_dimensions)), global_dimensions) ||
        return :placement_range_global_dimension_mismatch
    return nothing
end

function _cpb_place_overlap_collection_matrix(
    collection::CPBLocalOverlapBlockCollection,
    transform_lookup,
    range_lookup,
)
    retained_blocks = Tuple(
        _cpb_place_overlap_collection_retained_block(
            record,
            transform_lookup,
            range_lookup,
        )
        for record in collection.records
    )
    element_type = promote_type(
        Tuple(eltype(retained.retained_block) for retained in retained_blocks)...,
    )
    global_dimension = first(retained_blocks).range_summary.global_dimension
    global_overlap_matrix = zeros(element_type, global_dimension, global_dimension)
    for retained in retained_blocks
        global_overlap_matrix[
            retained.range_summary.left_column_range,
            retained.range_summary.right_column_range,
        ] .+= retained.retained_block
    end
    return global_overlap_matrix
end

function _cpb_place_overlap_collection_retained_block(
    record::CPBLocalOverlapBlockRecord,
    transform_lookup,
    range_lookup,
)
    block_key = record.block_key
    left_transform_carry = transform_lookup[(block_key, :left)]
    right_transform_carry = transform_lookup[(block_key, :right)]
    placement_range = range_lookup[block_key]
    retained_block =
        left_transform_carry.transform_object' *
        record.source_block.dense_block *
        right_transform_carry.transform_object
    return (;
        range_summary = summary(placement_range),
        retained_block,
    )
end

function _cpb_place_overlap_collection_summary(
    status::Symbol,
    blocker,
    collection_summary,
    facts_summary,
    record_placement_summaries::Tuple,
    global_overlap_matrix,
)
    available = status === :materialized_cpb_overlap_collection_placement_pilot
    coverage = _cpb_place_overlap_collection_record_coverage(
        record_placement_summaries,
    )
    return (;
        object_kind = :cartesian_cpb_placed_overlap_collection_summary,
        status,
        blocker,
        record_count = collection_summary.record_count,
        placed_record_count =
            count(
                record ->
                    record.status ===
                    :materialized_cpb_overlap_collection_record_pilot,
                record_placement_summaries,
            ),
        blocked_record_count =
            count(
                record ->
                    record.status !==
                    :materialized_cpb_overlap_collection_record_pilot,
                record_placement_summaries,
            ),
        block_keys = collection_summary.block_keys,
        record_placement_summaries,
        missing_left_transform_block_keys =
            coverage.missing_left_transform_block_keys,
        missing_right_transform_block_keys =
            coverage.missing_right_transform_block_keys,
        missing_placement_range_block_keys =
            coverage.missing_placement_range_block_keys,
        placement_facts_status =
            isnothing(facts_summary) ? :unavailable : facts_summary.status,
        placement_facts_blocker =
            isnothing(facts_summary) ? :unavailable : facts_summary.blocker,
        placement_plan_review_status =
            isnothing(facts_summary) ?
            :unavailable :
            facts_summary.placement_plan_review_status,
        accumulation_rule =
            isnothing(facts_summary) ? :unavailable : facts_summary.accumulation_rule,
        symmetry_policy =
            isnothing(facts_summary) ? :unavailable : facts_summary.symmetry_policy,
        duplicate_record_policy =
            isnothing(facts_summary) ? :unavailable : facts_summary.duplicate_record_policy,
        provider_level_global_overlap_matrix_shape =
            available ? size(global_overlap_matrix) : :not_materialized,
        global_overlap_matrix_shape = :route_global_matrix_not_materialized,
        provider_level_matrix_materialized = available,
        provider_level_overlap_matrix_materialized = available,
        provider_level_pilot = true,
        synthetic_fixture_only = true,
        global_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        route_global_matrix_materialized = false,
        route_global_overlap_matrix_materialized = false,
        route_driver_wiring = false,
        route_global_overlap_stage_source = false,
        route_global_overlap_available = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        exports_or_artifacts = false,
    )
end

end # module CartesianCPBBlockProviders
