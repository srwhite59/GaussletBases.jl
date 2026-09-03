const _AXIS_ORDER = (:x, :y, :z)
const _BRA_KET_ORDER = (:bra, :ket)

export CartesianParentAxisFactorPacket3D,
       CartesianParentCoulombSourceSummary,
       parent_overlap_axis_factor_packet,
       parent_coulomb_axis_source_summary,
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

struct CartesianParentCoulombSourceSummary{M}
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
summary(summary_object::CartesianParentCoulombSourceSummary) =
    summary_object.metadata

function parent_coulomb_axis_source_summary(
    parent::CartesianParentGaussletBasis3D,
    parent_axis_bundle_object,
    expansion;
    nuclear_charges = nothing,
    atom_locations = nothing,
    center_table = nothing,
    parent_qw_basis_object = nothing,
)
    parent_axis_counts_value = parent_axis_counts(parent)
    axis_bundle_available = !isnothing(parent_axis_bundle_object)
    expansion_available = !isnothing(expansion)
    expansion_coefficients =
        _property(expansion, :coefficients)
    expansion_exponents =
        _property(expansion, :exponents)
    expansion_coefficients_available =
        expansion_coefficients isa AbstractVector && !isempty(expansion_coefficients)
    expansion_exponents_available =
        expansion_exponents isa AbstractVector && !isempty(expansion_exponents)
    gaussian_factor_terms_available =
        _axis_coulomb_source_available(parent_axis_bundle_object, :gaussian_factor_terms)
    gaussian_factors_available =
        _axis_coulomb_source_available(parent_axis_bundle_object, :gaussian_factors)
    pair_factor_terms_available =
        _axis_coulomb_source_available(parent_axis_bundle_object, :pair_factor_terms)
    pair_factors_available =
        _axis_coulomb_source_available(parent_axis_bundle_object, :pair_factors)
    pair_factor_terms_raw_available =
        _axis_coulomb_source_available(parent_axis_bundle_object, :pair_factor_terms_raw)
    axis_exponents_available =
        _axis_coulomb_source_available(parent_axis_bundle_object, :exponents)
    center_metadata =
        _parent_coulomb_center_metadata(
            nuclear_charges,
            atom_locations,
            center_table,
            parent_qw_basis_object,
            parent,
        )
    electron_nuclear_center_metadata_available = center_metadata.available
    electron_nuclear_axis_terms_available = false
    electron_electron_pair_axis_terms_available =
        pair_factor_terms_available || pair_factors_available
    missing_sources = Symbol[]
    axis_bundle_available || push!(missing_sources, :missing_parent_axis_bundle_object)
    expansion_available || push!(missing_sources, :missing_coulomb_gaussian_expansion)
    expansion_coefficients_available ||
        push!(missing_sources, :missing_coulomb_expansion_coefficients)
    expansion_exponents_available ||
        push!(missing_sources, :missing_coulomb_expansion_exponents)
    electron_electron_pair_axis_terms_available ||
        push!(missing_sources, :missing_electron_electron_pair_axis_terms)
    electron_nuclear_center_metadata_available ||
        push!(missing_sources, :missing_electron_nuclear_center_metadata)
    missing_sources_tuple = Tuple(unique(missing_sources))
    status =
        axis_bundle_available &&
        expansion_available &&
        expansion_coefficients_available &&
        expansion_exponents_available &&
        electron_electron_pair_axis_terms_available ?
        :available_parent_coulomb_axis_source_summary :
        :blocked_parent_coulomb_axis_source_summary
    blocker =
        status === :available_parent_coulomb_axis_source_summary ||
        isempty(missing_sources_tuple) ?
        nothing :
        first(missing_sources_tuple)

    return CartesianParentCoulombSourceSummary((;
        object_kind = :cartesian_parent_coulomb_axis_source_summary,
        status,
        blocker,
        parent_axis_counts = parent_axis_counts_value,
        parent_axis_counts_source = :parent_object_parent_axis_counts,
        axis_bundle_available,
        expansion_available,
        expansion_coefficients_available,
        expansion_exponents_available,
        gaussian_expansion_source =
            expansion_available ? :coulomb_gaussian_expansion_object : :unavailable,
        electron_nuclear_center_metadata_available,
        electron_nuclear_center_metadata_source = center_metadata.source,
        electron_nuclear_axis_terms_available,
        electron_nuclear_axis_terms_status =
            electron_nuclear_axis_terms_available ?
            :available_parent_electron_nuclear_axis_terms :
            :missing_per_center_electron_nuclear_axis_term_tables,
        electron_electron_pair_axis_terms_available,
        electron_electron_pair_axis_terms_status =
            electron_electron_pair_axis_terms_available ?
            :available_parent_electron_electron_pair_axis_terms :
            :missing_electron_electron_pair_axis_terms,
        electron_electron_cpb_pair_pair_source_record_available = false,
        gaussian_factor_terms_available,
        gaussian_factors_available,
        pair_factor_terms_available,
        pair_factors_available,
        pair_factor_terms_raw_available,
        axis_exponents_available,
        missing_sources = missing_sources_tuple,
        source_paths = (;
            electron_nuclear_gaussian_factor_terms =
                :axis_pgdg_intermediate_gaussian_factor_terms,
            electron_nuclear_gaussian_factors =
                :axis_pgdg_intermediate_gaussian_factors,
            electron_electron_pair_factor_terms =
                :axis_pgdg_intermediate_pair_factor_terms,
            electron_electron_pair_factors =
                :axis_pgdg_intermediate_pair_factors,
            electron_electron_pair_factor_terms_raw =
                :axis_pgdg_intermediate_pair_factor_terms_raw,
            axis_exponents = :axis_pgdg_intermediate_exponents,
            expansion_coefficients = :expansion_coefficients,
            expansion_exponents = :expansion_exponents,
            nuclear_charges = center_metadata.nuclear_charges_source,
            atom_locations = center_metadata.atom_locations_source,
            center_table = center_metadata.center_table_source,
        ),
        numerical_coulomb_blocks_materialized = false,
        cpb_local_coulomb_kernel_implemented = false,
        wl_pqs_realization = false,
        route_global_placement = false,
        route_driver_wiring = false,
        hamiltonian_assembly = false,
        ida_mwg_pqs_semantics = false,
        exports_or_artifacts = false,
    ))
end

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

function _axis_coulomb_source_available(bundle, factor_name::Symbol)
    isnothing(bundle) && return false
    bundle_x = _axis_bundle_property(bundle, :bundle_x, :x)
    bundle_y = _axis_bundle_property(bundle, :bundle_y, :y)
    bundle_z = _axis_bundle_property(bundle, :bundle_z, :z)
    return all(
        axis -> !isnothing(_axis_coulomb_factor_from_bundle_axis(axis, factor_name)),
        (bundle_x, bundle_y, bundle_z),
    )
end

function _axis_coulomb_factor_from_bundle_axis(axis, factor_name::Symbol)
    isnothing(axis) && return nothing
    pgdg_intermediate = _property(axis, :pgdg_intermediate)
    pgdg_factor = _property(pgdg_intermediate, factor_name)
    !isnothing(pgdg_factor) && return pgdg_factor
    return _property(axis, factor_name)
end

function _parent_coulomb_center_metadata(
    nuclear_charges,
    atom_locations,
    center_table,
    parent_qw_basis_object,
    parent::CartesianParentGaussletBasis3D,
)
    if !isnothing(center_table)
        return (;
            available = true,
            source = :center_table,
            nuclear_charges_source =
                _center_table_has_property(center_table, :nuclear_charge) ?
                :center_table_nuclear_charge :
                :unavailable,
            atom_locations_source =
                _center_table_has_property(center_table, :location) ?
                :center_table_location :
                :unavailable,
            center_table_source = :center_table,
        )
    elseif !isnothing(nuclear_charges) && !isnothing(atom_locations)
        return (;
            available = true,
            source = :nuclear_charges_and_atom_locations,
            nuclear_charges_source = :report_nuclear_charges,
            atom_locations_source = :report_atom_locations,
            center_table_source = :unavailable,
        )
    elseif !isnothing(parent_qw_basis_object) &&
           hasproperty(parent_qw_basis_object, :nuclei) &&
           hasproperty(parent_qw_basis_object, :nuclear_charges)
        return (;
            available = true,
            source = :parent_qw_basis_object_nuclei_and_charges,
            nuclear_charges_source = :parent_qw_basis_object_nuclear_charges,
            atom_locations_source = :parent_qw_basis_object_nuclei,
            center_table_source = :unavailable,
        )
    elseif hasproperty(parent.metadata, :nuclei) &&
           hasproperty(parent.metadata, :nuclear_charges)
        return (;
            available = true,
            source = :parent_metadata_nuclei_and_charges,
            nuclear_charges_source = :parent_metadata_nuclear_charges,
            atom_locations_source = :parent_metadata_nuclei,
            center_table_source = :unavailable,
        )
    end
    return (;
        available = false,
        source = :unavailable,
        nuclear_charges_source = :unavailable,
        atom_locations_source = :unavailable,
        center_table_source = :unavailable,
    )
end

function _center_table_has_property(center_table, property_name::Symbol)
    for center in center_table
        hasproperty(center, property_name) || return false
    end
    return true
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
