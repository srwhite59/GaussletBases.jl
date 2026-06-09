module CartesianParentAxisFactors

using ..CartesianParentGaussletBases

const CPGB = CartesianParentGaussletBases
const _AXIS_ORDER = (:x, :y, :z)
const _BRA_KET_ORDER = (:bra, :ket)

export CartesianParentAxisFactorPacket3D,
       parent_overlap_axis_factor_packet,
       summary

"""
    CartesianParentAxisFactorPacket3D

Internal parent-owned axis factor packet. This first contract is overlap-only:
it carries parent-axis overlap factors and reports all other factor categories
as unavailable or not requested.
"""
struct CartesianParentAxisFactorPacket3D{P,O,M}
    parent::P
    overlap_1d::O
    metadata::M
end

summary(packet::CartesianParentAxisFactorPacket3D) = packet.metadata

function parent_overlap_axis_factor_packet(
    parent::CPGB.CartesianParentGaussletBasis3D,
    parent_axis_bundle_object,
)
    parent_axis_counts = CPGB.parent_axis_counts(parent)
    overlap_1d = _overlap_1d_from_axis_bundle(parent_axis_bundle_object)
    blocker = _overlap_1d_blocker(overlap_1d, parent_axis_counts)
    status =
        isnothing(blocker) ?
        :available_parent_overlap_axis_factors :
        :blocked_parent_overlap_axis_factors
    overlap_value = isnothing(blocker) ? overlap_1d : nothing
    return CartesianParentAxisFactorPacket3D(
        parent,
        overlap_value,
        _parent_axis_factor_packet_summary(
            status,
            blocker,
            parent_axis_counts,
            overlap_value,
        ),
    )
end

function parent_overlap_axis_factor_packet(
    parent::CPGB.CartesianParentGaussletBasis3D;
    parent_axis_bundle_object = nothing,
)
    return parent_overlap_axis_factor_packet(parent, parent_axis_bundle_object)
end

function _parent_axis_factor_packet_summary(
    status::Symbol,
    blocker,
    parent_axis_counts::NTuple{3,Int},
    overlap_1d,
)
    overlap_available = !isnothing(overlap_1d)
    return (;
        object_kind = :cartesian_parent_axis_factor_packet_summary,
        packet_kind = :overlap_only_parent_axis_factors,
        status,
        blocker,
        parent_axis_counts,
        parent_axis_counts_source = :parent_object_parent_axis_counts,
        overlap_status =
            overlap_available ?
            :available_parent_overlap_axis_factors :
            :missing_parent_axis_bundle_overlap_factors,
        overlap_1d_available = overlap_available,
        factor_space =
            overlap_available ?
            :parent_axis_bundle_pgdg_intermediate :
            :unavailable,
        factor_convention =
            overlap_available ?
            :axis_bundle_one_body_overlap :
            :unavailable,
        normalization_convention =
            overlap_available ?
            :not_separate_from_axis_bundle_one_body_overlap :
            :unavailable,
        index_domain =
            overlap_available ?
            :parent_axis_indices :
            :unavailable,
        axis_order = _AXIS_ORDER,
        bra_ket_order =
            overlap_available ?
            _BRA_KET_ORDER :
            :unavailable,
        sliceable_by_cpb = overlap_available,
        category_availability = (;
            overlap =
                overlap_available ?
                :available_parent_overlap_axis_factors :
                :missing_parent_axis_bundle_overlap_factors,
            kinetic = :not_requested_parent_kinetic_axis_factors,
            position = :not_requested_parent_position_axis_factors,
            x2 = :not_requested_parent_x2_axis_factors,
            coulomb = :not_requested_parent_coulomb_axis_factors,
        ),
        full_3d_parent_matrices = false,
        cpb_slicing_implemented = false,
        route_driver_wiring = false,
        hamiltonian_data_materialized = false,
        coulomb_data_materialized = false,
        ida_mwg_semantics = false,
        exports_or_artifacts = false,
    )
end

function _overlap_1d_blocker(overlap_1d, parent_axis_counts::NTuple{3,Int})
    isnothing(overlap_1d) && return :missing_parent_axis_bundle_overlap_factors
    for (axis, count) in zip(_AXIS_ORDER, parent_axis_counts)
        matrix = getproperty(overlap_1d, axis)
        matrix isa AbstractMatrix ||
            return Symbol("$(axis)_overlap_axis_factor_not_matrix")
        size(matrix) == (count, count) ||
            return Symbol("$(axis)_overlap_axis_factor_size_mismatch")
    end
    return nothing
end

function _overlap_1d_from_axis_bundle(bundle)
    isnothing(bundle) && return nothing
    bundle_x = _axis_bundle_property(bundle, :bundle_x, :x)
    bundle_y = _axis_bundle_property(bundle, :bundle_y, :y)
    bundle_z = _axis_bundle_property(bundle, :bundle_z, :z)
    overlap_x = _axis_overlap_from_bundle_axis(bundle_x)
    overlap_y = _axis_overlap_from_bundle_axis(bundle_y)
    overlap_z = _axis_overlap_from_bundle_axis(bundle_z)
    any(isnothing, (overlap_x, overlap_y, overlap_z)) && return nothing
    return (x = overlap_x, y = overlap_y, z = overlap_z)
end

function _axis_bundle_property(bundle, primary::Symbol, fallback::Symbol)
    hasproperty(bundle, primary) && return getproperty(bundle, primary)
    hasproperty(bundle, fallback) && return getproperty(bundle, fallback)
    return nothing
end

function _property(object, name::Symbol)
    isnothing(object) && return nothing
    hasproperty(object, name) && return getproperty(object, name)
    return nothing
end

function _axis_overlap_from_bundle_axis(axis)
    isnothing(axis) && return nothing
    axis isa AbstractMatrix && return axis
    pgdg_intermediate = _property(axis, :pgdg_intermediate)
    overlap = _property(pgdg_intermediate, :overlap)
    !isnothing(overlap) && return overlap
    return _property(axis, :overlap)
end

end # module CartesianParentAxisFactors
