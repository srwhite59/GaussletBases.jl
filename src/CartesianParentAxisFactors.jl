const _AXIS_ORDER = (:x, :y, :z)
const _BRA_KET_ORDER = (:bra, :ket)

export CartesianParentAxisFactorPacket3D,
       parent_overlap_axis_factor_packet,
       summary

"""
    CartesianParentAxisFactorPacket3D

Internal parent-owned axis factor packet. The initial production contract
requires overlap factors and may also carry same-source kinetic factors when
the parent axis bundle already owns them. Position and x2 factors are also
carried when present in the same parent-axis bundle.
"""
struct CartesianParentAxisFactorPacket3D{P,O,K,X,Q,M}
    parent::P
    overlap_1d::O
    kinetic_1d::K
    position_1d::X
    x2_1d::Q
    metadata::M
end

CartesianParentAxisFactorPacket3D(parent, overlap_1d, metadata) =
    CartesianParentAxisFactorPacket3D(
        parent,
        overlap_1d,
        nothing,
        nothing,
        nothing,
        metadata,
    )

CartesianParentAxisFactorPacket3D(parent, overlap_1d, kinetic_1d, metadata) =
    CartesianParentAxisFactorPacket3D(
        parent,
        overlap_1d,
        kinetic_1d,
        nothing,
        nothing,
        metadata,
    )

summary(packet::CartesianParentAxisFactorPacket3D) = packet.metadata

function parent_overlap_axis_factor_packet(
    parent::CartesianParentGaussletBasis3D,
    parent_axis_bundle_object,
)
    axis_counts = parent_axis_counts(parent)
    overlap_1d = _overlap_1d_from_axis_bundle(parent_axis_bundle_object)
    kinetic_1d = _kinetic_1d_from_axis_bundle(parent_axis_bundle_object)
    position_1d = _position_1d_from_axis_bundle(parent_axis_bundle_object)
    x2_1d = _x2_1d_from_axis_bundle(parent_axis_bundle_object)
    blocker = _overlap_1d_blocker(overlap_1d, axis_counts)
    kinetic_blocker = _kinetic_1d_blocker(kinetic_1d, axis_counts)
    position_blocker = _position_1d_blocker(position_1d, axis_counts)
    x2_blocker = _x2_1d_blocker(x2_1d, axis_counts)
    status =
        isnothing(blocker) ?
        :available_parent_overlap_axis_factors :
        :blocked_parent_overlap_axis_factors
    overlap_value = isnothing(blocker) ? overlap_1d : nothing
    kinetic_value = isnothing(kinetic_blocker) ? kinetic_1d : nothing
    position_value = isnothing(position_blocker) ? position_1d : nothing
    x2_value = isnothing(x2_blocker) ? x2_1d : nothing
    return CartesianParentAxisFactorPacket3D(
        parent,
        overlap_value,
        kinetic_value,
        position_value,
        x2_value,
        _parent_axis_factor_packet_summary(
            status,
            blocker,
            axis_counts,
            overlap_value,
            kinetic_value,
            kinetic_blocker,
            position_value,
            position_blocker,
            x2_value,
            x2_blocker,
        ),
    )
end

function parent_overlap_axis_factor_packet(
    parent::CartesianParentGaussletBasis3D;
    parent_axis_bundle_object = nothing,
)
    return parent_overlap_axis_factor_packet(parent, parent_axis_bundle_object)
end

function _parent_axis_factor_packet_summary(
    status::Symbol,
    blocker,
    parent_axis_counts::NTuple{3,Int},
    overlap_1d,
    kinetic_1d,
    kinetic_blocker,
    position_1d,
    position_blocker,
    x2_1d,
    x2_blocker,
)
    overlap_available = !isnothing(overlap_1d)
    kinetic_available = !isnothing(kinetic_1d)
    position_available = !isnothing(position_1d)
    x2_available = !isnothing(x2_1d)
    return (;
        object_kind = :cartesian_parent_axis_factor_packet_summary,
        packet_kind =
            (kinetic_available || position_available || x2_available) ?
            :overlap_one_body_parent_axis_factors :
            :overlap_only_parent_axis_factors,
        status,
        blocker,
        parent_axis_counts,
        parent_axis_counts_source = :parent_object_parent_axis_counts,
        overlap_status =
            overlap_available ?
            :available_parent_overlap_axis_factors :
            :missing_parent_axis_bundle_overlap_factors,
        overlap_1d_available = overlap_available,
        kinetic_status =
            kinetic_available ?
            :available_parent_kinetic_axis_factors :
            (
                isnothing(kinetic_blocker) ?
                :missing_parent_axis_bundle_kinetic_factors :
                kinetic_blocker
            ),
        kinetic_1d_available = kinetic_available,
        kinetic_factor_space =
            kinetic_available ?
            :parent_axis_bundle_pgdg_intermediate :
            :unavailable,
        kinetic_factor_convention =
            kinetic_available ?
            :axis_bundle_one_body_kinetic :
            :unavailable,
        kinetic_index_domain =
            kinetic_available ?
            :parent_axis_indices :
            :unavailable,
        kinetic_index_domain_source =
            kinetic_available ?
            :axis_bundle_contract :
            :unavailable,
        kinetic_index_domain_status =
            kinetic_available ?
            :assumed_parent_axis_indexed_by_current_axis_bundle_contract :
            :unavailable,
        kinetic_sliceable_by_cpb = kinetic_available,
        kinetic_sliceability_source =
            kinetic_available ?
            :index_domain_contract :
            :unavailable,
        kinetic_sliceability_status =
            kinetic_available ?
            :sliceable_by_cpb_parent_axis_index_contract :
            :unavailable,
        position_status =
            position_available ?
            :available_parent_position_axis_factors :
            (
                isnothing(position_blocker) ?
                :missing_parent_axis_bundle_position_factors :
                position_blocker
            ),
        position_1d_available = position_available,
        position_factor_space =
            position_available ?
            :parent_axis_bundle_pgdg_intermediate :
            :unavailable,
        position_factor_convention =
            position_available ?
            :axis_bundle_one_body_position :
            :unavailable,
        position_index_domain =
            position_available ?
            :parent_axis_indices :
            :unavailable,
        position_index_domain_source =
            position_available ?
            :axis_bundle_contract :
            :unavailable,
        position_index_domain_status =
            position_available ?
            :assumed_parent_axis_indexed_by_current_axis_bundle_contract :
            :unavailable,
        position_sliceable_by_cpb = position_available,
        position_sliceability_source =
            position_available ?
            :index_domain_contract :
            :unavailable,
        position_sliceability_status =
            position_available ?
            :sliceable_by_cpb_parent_axis_index_contract :
            :unavailable,
        x2_status =
            x2_available ?
            :available_parent_x2_axis_factors :
            (
                isnothing(x2_blocker) ?
                :missing_parent_axis_bundle_x2_factors :
                x2_blocker
            ),
        x2_1d_available = x2_available,
        x2_factor_space =
            x2_available ?
            :parent_axis_bundle_pgdg_intermediate :
            :unavailable,
        x2_factor_convention =
            x2_available ?
            :axis_bundle_one_body_x2 :
            :unavailable,
        x2_index_domain =
            x2_available ?
            :parent_axis_indices :
            :unavailable,
        x2_index_domain_source =
            x2_available ?
            :axis_bundle_contract :
            :unavailable,
        x2_index_domain_status =
            x2_available ?
            :assumed_parent_axis_indexed_by_current_axis_bundle_contract :
            :unavailable,
        x2_sliceable_by_cpb = x2_available,
        x2_sliceability_source =
            x2_available ?
            :index_domain_contract :
            :unavailable,
        x2_sliceability_status =
            x2_available ?
            :sliceable_by_cpb_parent_axis_index_contract :
            :unavailable,
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
        index_domain_source =
            overlap_available ?
            :axis_bundle_contract :
            :unavailable,
        index_domain_status =
            overlap_available ?
            :assumed_parent_axis_indexed_by_current_axis_bundle_contract :
            :unavailable,
        axis_order = _AXIS_ORDER,
        bra_ket_order =
            overlap_available ?
            _BRA_KET_ORDER :
            :unavailable,
        sliceable_by_cpb = overlap_available,
        sliceability_source =
            overlap_available ?
            :index_domain_contract :
            :unavailable,
        sliceability_status =
            overlap_available ?
            :sliceable_by_cpb_parent_axis_index_contract :
            :unavailable,
        category_availability = (;
            overlap =
                overlap_available ?
                :available_parent_overlap_axis_factors :
                :missing_parent_axis_bundle_overlap_factors,
            kinetic =
                kinetic_available ?
                :available_parent_kinetic_axis_factors :
                (
                    isnothing(kinetic_blocker) ?
                    :missing_parent_axis_bundle_kinetic_factors :
                    kinetic_blocker
                ),
            position =
                position_available ?
                :available_parent_position_axis_factors :
                (
                    isnothing(position_blocker) ?
                    :missing_parent_axis_bundle_position_factors :
                    position_blocker
                ),
            x2 =
                x2_available ?
                :available_parent_x2_axis_factors :
                (
                    isnothing(x2_blocker) ?
                    :missing_parent_axis_bundle_x2_factors :
                    x2_blocker
                ),
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

function _kinetic_1d_blocker(kinetic_1d, parent_axis_counts::NTuple{3,Int})
    return _one_body_1d_blocker(
        kinetic_1d,
        parent_axis_counts,
        :kinetic,
        :missing_parent_axis_bundle_kinetic_factors,
    )
end

function _position_1d_blocker(position_1d, parent_axis_counts::NTuple{3,Int})
    return _one_body_1d_blocker(
        position_1d,
        parent_axis_counts,
        :position,
        :missing_parent_axis_bundle_position_factors,
    )
end

function _x2_1d_blocker(x2_1d, parent_axis_counts::NTuple{3,Int})
    return _one_body_1d_blocker(
        x2_1d,
        parent_axis_counts,
        :x2,
        :missing_parent_axis_bundle_x2_factors,
    )
end

function _one_body_1d_blocker(
    one_body_1d,
    parent_axis_counts::NTuple{3,Int},
    factor_name::Symbol,
    missing_blocker::Symbol,
)
    isnothing(one_body_1d) && return missing_blocker
    for (axis, count) in zip(_AXIS_ORDER, parent_axis_counts)
        matrix = getproperty(one_body_1d, axis)
        matrix isa AbstractMatrix ||
            return Symbol("$(axis)_$(factor_name)_axis_factor_not_matrix")
        size(matrix) == (count, count) ||
            return Symbol("$(axis)_$(factor_name)_axis_factor_size_mismatch")
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

function _kinetic_1d_from_axis_bundle(bundle)
    return _one_body_1d_from_axis_bundle(bundle, :kinetic)
end

function _position_1d_from_axis_bundle(bundle)
    return _one_body_1d_from_axis_bundle(bundle, :position)
end

function _x2_1d_from_axis_bundle(bundle)
    return _one_body_1d_from_axis_bundle(bundle, :x2)
end

function _one_body_1d_from_axis_bundle(bundle, factor_name::Symbol)
    isnothing(bundle) && return nothing
    bundle_x = _axis_bundle_property(bundle, :bundle_x, :x)
    bundle_y = _axis_bundle_property(bundle, :bundle_y, :y)
    bundle_z = _axis_bundle_property(bundle, :bundle_z, :z)
    factor_x = _axis_one_body_factor_from_bundle_axis(bundle_x, factor_name)
    factor_y = _axis_one_body_factor_from_bundle_axis(bundle_y, factor_name)
    factor_z = _axis_one_body_factor_from_bundle_axis(bundle_z, factor_name)
    any(isnothing, (factor_x, factor_y, factor_z)) && return nothing
    return (x = factor_x, y = factor_y, z = factor_z)
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

function _axis_kinetic_from_bundle_axis(axis)
    return _axis_one_body_factor_from_bundle_axis(axis, :kinetic)
end

function _axis_one_body_factor_from_bundle_axis(axis, factor_name::Symbol)
    isnothing(axis) && return nothing
    pgdg_intermediate = _property(axis, :pgdg_intermediate)
    pgdg_factor = _property(pgdg_intermediate, factor_name)
    !isnothing(pgdg_factor) && return pgdg_factor
    return _property(axis, factor_name)
end
