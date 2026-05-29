module CartesianContractedParents

import ..GaussletBases: _CartesianCoefficientMap,
                         _BondAlignedDiatomicAtomGrowthConstructionRegion3D,
                         _BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D,
                         _BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
                         _CartesianNestedEndcapPanelOwnedUnits3D,
                         _CartesianNestedEndcapPanelShellLayer3D,
                         _CartesianNestedProductStagedByCenterUnit3D,
                         _CartesianNestedProjectedQShellStagedUnitDescriptor3D,
                         _NestedFixedBlock3D,
                         _cartesian_coefficient_map_storage,
                         _nested_staged_by_center_sidecar
import ..GaussletBases.CartesianParentGaussletBases:
    CartesianParentGaussletBasis3D,
    cartesian_parent_gausslet_basis,
    parent_dimension

export CartesianContractionUnit3D,
       CartesianContractionRule3D,
       CartesianContractionRuleInventory3D,
       CartesianContractedParent3D,
       CartesianContractedParentStructuralAudit,
       CartesianShellRegionInventory3D,
       CartesianShellRegionRetention3D,
       CartesianShellRegion3D,
       cartesian_contraction_rule,
       cartesian_contraction_rule_inventory,
       cartesian_contraction_unit_from_rule,
       cartesian_contracted_parent,
       cartesian_shell_region,
       cartesian_shell_region_inventory,
       contracted_parent_contraction_rules,
       contracted_parent_rule_inventory,
       contracted_parent_basis,
       contracted_parent_coefficients,
       contracted_parent_units,
       contracted_parent_metadata,
       contracted_parent_parent_dimension,
       contracted_parent_dimension,
       contraction_unit_role,
       contraction_unit_support_indices,
       contraction_unit_column_range,
       contraction_unit_metadata,
       contraction_unit_rule,
       contraction_rule_family,
       contraction_rule_kind,
       contraction_rule_support_summary,
       contraction_rule_column_range,
       contraction_rule_source_dimension,
       contraction_rule_retained_dimension,
       contraction_rule_transform_rule,
       contraction_rule_cleanup_rule,
       contraction_rule_metric_capability,
       contracted_parent_unit_column_ranges,
       contracted_parent_unit_support_indices,
       contracted_parent_support_indices,
       contracted_parent_structural_audit

"""
    CartesianShellRegionRetention3D

Metadata-only retention/contraction facts for a shell region. These fields are
kept separate from region geometry so a region can be audited without implying
that a particular contraction rule already drives construction.
"""
struct CartesianShellRegionRetention3D{D}
    retention_rule::Symbol
    cleanup_rule::Symbol
    preferred_contraction_rule::Symbol
    expected_unit_family::Symbol
    metric_capability::Symbol
    required_payload_fields::Tuple{Vararg{Symbol}}
    missing_payload_fields::Tuple{Vararg{Symbol}}
    diagnostics::D
end

"""
    CartesianShellRegion3D

Metadata-only descriptor for a finite Cartesian shell/region between a full
parent lattice and final contraction units. Region facts are deliberately
separate from retention/contraction facts; this object does not build
coefficient maps, metric packets, fixed blocks, QW operators, or Hamiltonians.
"""
struct CartesianShellRegion3D{S,G,R,D,P}
    region_family::Symbol
    role::Union{Nothing,Symbol}
    status::Symbol
    box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    inner_exclusion_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    support_indices::Union{Nothing,Vector{Int}}
    support_summary::S
    ownership_coverage_contract::Symbol
    geometry::G
    retention::R
    current_route_consumes::Bool
    descriptor_drives_builder::Bool
    descriptor_only::Bool
    diagnostics::D
    provenance::P
end

"""
    CartesianShellRegionInventory3D

Ordered metadata-only inventory of shell regions for one construction route or
diagnostic source.
"""
struct CartesianShellRegionInventory3D{R,S,C,D,P}
    route::Symbol
    regions::R
    region_count::Int
    region_order::Vector{Symbol}
    status_counts::C
    support_summary::S
    current_route_consumes_count::Int
    descriptor_only_count::Int
    diagnostics::D
    provenance::P
end

"""
    CartesianContractionRule3D

Internal metadata record describing how a contracted-parent unit is meant to
be constructed. This is deliberately descriptive only: coefficient matrices,
fixed-block builders, metric packets, QW operators, and Hamiltonians remain
owned by their existing route-specific implementations.
"""
struct CartesianContractionRule3D{S,L,D,P}
    rule_family::Symbol
    kind::Symbol
    role::Union{Nothing,Symbol}
    support_indices::Vector{Int}
    support_summary::S
    local_geometry::L
    column_range::Union{Nothing,UnitRange{Int}}
    source_dimension::Int
    retained_dimension::Int
    transform_rule::Symbol
    cleanup_rule::Symbol
    metric_capability::Symbol
    diagnostics::D
    provenance::P
end

"""
    CartesianContractionRuleInventory3D

Metadata-only parent/rule inventory for contracted-parent construction rules.
It summarizes rule families, support coverage, retained dimensions, and metric
capabilities without building coefficient maps, metric packets, or operators.
"""
struct CartesianContractionRuleInventory3D{R,F,C,S,D,P}
    parent_dimension::Int
    contracted_dimension::Union{Nothing,Int}
    unit_count::Int
    rule_count::Int
    rules::R
    rule_family_counts::F
    metric_capabilities::C
    total_source_dimension::Int
    total_retained_dimension::Int
    support_summary::S
    rule_support_summaries::Vector{Any}
    every_unit_has_rule_metadata::Bool
    every_unit_rule_derivable::Bool
    metadata_only_rule_count::Int
    prototype_rule_count::Int
    any_metadata_only_rule::Bool
    any_prototype_rule::Bool
    diagnostics::D
    provenance::P
end

function _contraction_rule_support_summary(
    support_indices::AbstractVector{<:Integer};
    parent_dimension::Union{Nothing,Int} = nothing,
)
    values = Int[Int(index) for index in support_indices]
    unique_values = unique(values)
    duplicate_count = length(values) - length(unique_values)
    if isnothing(parent_dimension)
        return (
            parent_dimension = nothing,
            entry_count = length(values),
            unique_count = length(unique_values),
            duplicate_count,
            outside_count = nothing,
            missing_count = nothing,
            support_complete = nothing,
            coverage_checked = false,
        )
    end
    valid_range = 1:Int(parent_dimension)
    inside = Int[value for value in values if value in valid_range]
    unique_inside = unique(inside)
    return (
        parent_dimension = Int(parent_dimension),
        entry_count = length(values),
        unique_count = length(unique_inside),
        duplicate_count = length(inside) - length(unique_inside),
        outside_count = length(values) - length(inside),
        missing_count = Int(parent_dimension) - length(unique_inside),
        support_complete = length(unique_inside) == Int(parent_dimension),
        coverage_checked = true,
    )
end

function _shell_region_support_summary(
    support_indices::AbstractVector{<:Integer};
    parent_dimension::Union{Nothing,Int} = nothing,
)
    return merge(
        _contraction_rule_support_summary(support_indices; parent_dimension),
        (support_indices_available = true,),
    )
end

function _shell_region_support_summary(
    support_indices::Nothing;
    parent_dimension::Union{Nothing,Int} = nothing,
    support_count::Union{Nothing,Int} = nothing,
)
    complete =
        !isnothing(parent_dimension) && !isnothing(support_count) ?
        Int(parent_dimension) == Int(support_count) : nothing
    return (
        parent_dimension = isnothing(parent_dimension) ? nothing : Int(parent_dimension),
        entry_count = isnothing(support_count) ? nothing : Int(support_count),
        unique_count = isnothing(support_count) ? nothing : Int(support_count),
        duplicate_count = nothing,
        outside_count = nothing,
        missing_count = complete === true ? 0 : nothing,
        support_complete = complete,
        coverage_checked = false,
        support_indices_available = false,
    )
end

function _shell_region_named_property(metadata, name::Symbol, default)
    hasproperty(metadata, name) && return getproperty(metadata, name)
    return default
end

function _shell_region_retention(;
    retention_rule::Symbol,
    cleanup_rule::Symbol,
    preferred_contraction_rule::Symbol,
    expected_unit_family::Symbol,
    metric_capability::Symbol,
    required_payload_fields::Tuple{Vararg{Symbol}},
    missing_payload_fields::Tuple{Vararg{Symbol}} = (),
    diagnostics = (;),
)
    return CartesianShellRegionRetention3D(
        retention_rule,
        cleanup_rule,
        preferred_contraction_rule,
        expected_unit_family,
        metric_capability,
        required_payload_fields,
        missing_payload_fields,
        diagnostics,
    )
end

function _atom_growth_region_family(role::Symbol)
    role in (:left_atom_box, :right_atom_box) && return :atom_core_cube
    role == :contact_cap && return :shared_midpoint_slab_cap
    role == :regular_shared_molecular_shell && return :rectangular_molecular_shell
    role == :outer_mismatch_shared_molecular_shell && return :outer_boundary_shell
    return :support_dense_leftover_debug_region
end

function _atom_growth_region_status(role::Symbol)
    role in (
        :left_atom_box,
        :right_atom_box,
        :regular_shared_molecular_shell,
        :contact_cap,
        :outer_mismatch_shared_molecular_shell,
    ) && return :clean
    return :debug_oracle
end

function _atom_growth_region_retention_rule(role::Symbol)
    role in (:left_atom_box, :right_atom_box) && return :protected_atom_cubic_shell
    role == :contact_cap && return :shared_contact_cap
    role == :regular_shared_molecular_shell && return :policy_selected_shared_exterior
    role == :outer_mismatch_shared_molecular_shell &&
        return :outermost_mismatch_shared_molecular_shell
    return :explicit_support_dense_coefficients
end

function _atom_growth_region_contraction_rule(role::Symbol)
    role in (:left_atom_box, :right_atom_box) && return :complete_shell_sequence
    role == :contact_cap && return :contact_cap_owned_slab
    role == :regular_shared_molecular_shell && return :policy_selected_shared_exterior
    role == :outer_mismatch_shared_molecular_shell &&
        return :outer_mismatch_boundary_slab_set
    return :support_dense_fallback
end

function _atom_growth_region_expected_unit_family(role::Symbol)
    role == :regular_shared_molecular_shell && return :policy_dependent_shell_region
    return :support_dense_fallback
end

function _atom_growth_region_metric_capability(role::Symbol)
    role == :regular_shared_molecular_shell && return :metadata_only_policy_dependent
    return :support_local_product
end

function _atom_growth_region_payload_fields(role::Symbol)
    role in (:left_atom_box, :right_atom_box) && return (
        required = (:box, :support_indices, :coefficient_matrix, :column_range),
        missing = (:coefficient_matrix, :column_range),
    )
    role == :regular_shared_molecular_shell && return (
        required = (:box, :inner_exclusion_box, :policy_selected_layer, :column_range),
        missing = (:policy_selected_layer, :column_range),
    )
    role == :contact_cap && return (
        required = (:box, :support_indices, :coefficient_matrix, :column_range),
        missing = (:coefficient_matrix, :column_range),
    )
    role == :outer_mismatch_shared_molecular_shell && return (
        required = (:box, :inner_exclusion_box, :piece_descriptors, :column_range),
        missing = (:piece_descriptors, :column_range),
    )
    return (
        required = (:support_indices, :coefficient_matrix, :column_range),
        missing = (:coefficient_matrix, :column_range),
    )
end

function _source_build_region_family(build::_BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D)
    build.region_category == :atom_local && return :atom_core_cube
    build.region_category == :contact_cap && return :shared_midpoint_slab_cap
    build.region_category == :outer_mismatch && return :outer_boundary_shell
    if build.region_category == :shared_exterior
        build.primitive_family == :shared_endcap_panel_shell_layer &&
            return :endcap_panel_shared_exterior
        build.primitive_family == :projected_q_shell &&
            return :projected_q_shell_boundary_modes
        return :rectangular_molecular_shell
    end
    return :support_dense_leftover_debug_region
end

function _source_build_region_status(build::_BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D)
    build.primitive_family == :shared_endcap_panel_shell_layer && return :transitional
    build.primitive_family == :projected_q_shell && return :prototype
    build.unsupported_reason !== nothing && return :debug_oracle
    return :clean
end

function _source_build_coverage_contract(build::_BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D)
    build.primitive_family in (:shared_endcap_panel_shell_layer, :projected_q_shell) &&
        return :boundary_only
    return :disjoint_partition_piece
end

function _source_build_retention_spec(build::_BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D)
    family = build.primitive_family
    if family == :shared_endcap_panel_shell_layer
        return (
            retention_rule = :old_endcap_panel_product_split,
            cleanup_rule = :locally_orthonormal_product_doside,
            preferred_contraction_rule = :old_endcap_panel_product_split,
            expected_unit_family = :product_owned_unit,
            metric_capability = :product_staged_metric_contraction,
            required_payload_fields = (
                :owned_units,
                :unit_column_ranges,
                :support_indices,
                :support_states,
                :staged_axes,
                :axis_function_indices,
                :coefficient_matrix,
                :column_range,
            ),
            missing_payload_fields = (),
            payload_ready = true,
        )
    elseif family == :projected_q_shell
        return (
            retention_rule = :boundary_comx_product_modes_raw_boundary_projection,
            cleanup_rule = :full_rank_symmetric_lowdin,
            preferred_contraction_rule = :pqs_boundary_projection_from_filled_box,
            expected_unit_family = :projected_q_shell_descriptor,
            metric_capability = :pqs_low_order_product_metric_prototype,
            required_payload_fields = (
                :support_indices,
                :support_states,
                :axis_local_coefficients,
                :boundary_mode_indices,
                :cleanup_transform,
                :column_range,
                :fixed_block_sidecar_payload,
            ),
            missing_payload_fields = (
                :fixed_block_sidecar_payload,
                :product_staged_metric_payload,
            ),
            payload_ready = false,
        )
    elseif family == :atom_local_complete_shell_sequence
        return (
            retention_rule = :protected_atom_cubic_shell,
            cleanup_rule = :complete_shell_sequence_cleanup,
            preferred_contraction_rule = :complete_shell_sequence,
            expected_unit_family = :support_dense_fallback,
            metric_capability = :support_local_product,
            required_payload_fields = (:box_or_sequence, :coefficient_matrix, :column_range),
            missing_payload_fields = (),
            payload_ready = true,
        )
    elseif family == :contact_cap_owned_slab
        return (
            retention_rule = :shared_contact_cap,
            cleanup_rule = :explicit_direct_slab,
            preferred_contraction_rule = :contact_cap_owned_slab,
            expected_unit_family = :support_dense_fallback,
            metric_capability = :support_local_product,
            required_payload_fields = (:box, :support_indices, :coefficient_matrix, :column_range),
            missing_payload_fields = (),
            payload_ready = true,
        )
    elseif family == :outer_mismatch_boundary_slab_set
        return (
            retention_rule = :outermost_mismatch_shared_molecular_shell,
            cleanup_rule = :explicit_direct_slab_set,
            preferred_contraction_rule = :outer_mismatch_boundary_slab_set,
            expected_unit_family = :support_dense_fallback,
            metric_capability = :support_local_product,
            required_payload_fields = (
                :piece_descriptors,
                :support_indices,
                :coefficient_matrix,
                :column_range,
            ),
            missing_payload_fields = (),
            payload_ready = true,
        )
    end
    return (
        retention_rule = :explicit_support_dense_coefficients,
        cleanup_rule = :external_or_already_cleaned,
        preferred_contraction_rule = :support_dense_fallback,
        expected_unit_family = :support_dense_fallback,
        metric_capability = :support_local_product,
        required_payload_fields = (:support_indices, :coefficient_matrix, :column_range),
        missing_payload_fields = build.unsupported_reason === nothing ? () : (:coefficient_matrix,),
        payload_ready = build.unsupported_reason === nothing,
    )
end

function cartesian_shell_region(
    build::_BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    spec = _source_build_retention_spec(build)
    retention = _shell_region_retention(
        retention_rule = spec.retention_rule,
        cleanup_rule = spec.cleanup_rule,
        preferred_contraction_rule = spec.preferred_contraction_rule,
        expected_unit_family = spec.expected_unit_family,
        metric_capability = spec.metric_capability,
        required_payload_fields = spec.required_payload_fields,
        missing_payload_fields = spec.missing_payload_fields,
        diagnostics = (
            source = :bond_aligned_diatomic_high_order_recipe_region_source_build,
            payload_ready_for_current_metric_execution = spec.payload_ready,
        ),
    )
    geometry = (
        source = :bond_aligned_diatomic_high_order_recipe_region_source_build,
        region_order_index = build.region_order_index,
        region_category = build.region_category,
        recipe_family = build.recipe_family,
        q = build.q,
        order = build.order,
        mapped_primitive = build.mapped_primitive,
        primitive_family = build.primitive_family,
        region_support_count = build.region_support_count,
        built_support_count = build.built_support_count,
        retained_count = build.retained_count,
        column_range = build.column_range,
    )
    diagnostics = (
        metadata_only = true,
        built = build.built,
        unsupported_reason = build.unsupported_reason,
        support_coverage_ok = build.support_coverage.coverage_ok,
        coefficient_map_generated = false,
        metric_packet_generated = false,
        fixed_block_sidecar_installed = false,
    )
    return CartesianShellRegion3D(
        _source_build_region_family(build),
        build.region_role,
        _source_build_region_status(build),
        nothing,
        nothing,
        nothing,
        _shell_region_support_summary(
            nothing;
            parent_dimension,
            support_count = build.built_support_count,
        ),
        _source_build_coverage_contract(build),
        geometry,
        retention,
        build.active_builder_consumes,
        false,
        !build.active_builder_consumes,
        diagnostics,
        (
            source = :bond_aligned_diatomic_high_order_recipe_region_source_build,
            original_metadata = build.metadata,
        ),
    )
end

function cartesian_shell_region(
    region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    payload_fields = _atom_growth_region_payload_fields(region.role)
    retention = _shell_region_retention(
        retention_rule = _atom_growth_region_retention_rule(region.role),
        cleanup_rule = :region_rule_dependent_or_external,
        preferred_contraction_rule = _atom_growth_region_contraction_rule(region.role),
        expected_unit_family = _atom_growth_region_expected_unit_family(region.role),
        metric_capability = _atom_growth_region_metric_capability(region.role),
        required_payload_fields = payload_fields.required,
        missing_payload_fields = payload_fields.missing,
        diagnostics = (
            source = :bond_aligned_diatomic_atom_growth_construction_region,
            payload_ready_for_current_metric_execution = false,
        ),
    )
    geometry = (
        source = :bond_aligned_diatomic_atom_growth_construction_region,
        order_index = region.order_index,
        box = region.box,
        inner_exclusion_box = region.inner_exclusion_box,
        atom_side = _shell_region_named_property(region.metadata, :atom_side, nothing),
        atom_axis_index = _shell_region_named_property(region.metadata, :atom_axis_index, nothing),
        shell_offset = _shell_region_named_property(region.metadata, :shell_offset, nothing),
        contact_policy = _shell_region_named_property(region.metadata, :contact_policy, nothing),
        contact_gap_count = _shell_region_named_property(region.metadata, :contact_gap_count, nothing),
        mismatch_low_counts = _shell_region_named_property(region.metadata, :low_counts, nothing),
        mismatch_high_counts = _shell_region_named_property(region.metadata, :high_counts, nothing),
    )
    diagnostics = (
        metadata_only = true,
        source_region_role = region.role,
        support_count = length(region.support_indices),
        coefficient_map_generated = false,
        metric_packet_generated = false,
        fixed_block_sidecar_installed = false,
    )
    return CartesianShellRegion3D(
        _atom_growth_region_family(region.role),
        region.role,
        _atom_growth_region_status(region.role),
        region.box,
        region.inner_exclusion_box,
        copy(region.support_indices),
        _shell_region_support_summary(region.support_indices; parent_dimension),
        :disjoint_partition_piece,
        geometry,
        retention,
        false,
        false,
        true,
        diagnostics,
        (
            source = :bond_aligned_diatomic_atom_growth_construction_region,
            original_metadata = region.metadata,
        ),
    )
end

function cartesian_shell_region(
    owned_units::_CartesianNestedEndcapPanelOwnedUnits3D;
    parent_dimension::Union{Nothing,Int} = nothing,
    active_builder_consumes::Bool = true,
)
    retention = _shell_region_retention(
        retention_rule = :old_endcap_panel_product_split,
        cleanup_rule = :locally_orthonormal_product_doside,
        preferred_contraction_rule = :old_endcap_panel_product_split,
        expected_unit_family = :product_owned_unit,
        metric_capability = :product_staged_metric_contraction,
        required_payload_fields = (
            :owned_units,
            :unit_column_ranges,
            :support_indices,
            :support_states,
            :staged_axes,
            :axis_function_indices,
            :coefficient_matrix,
            :column_range,
        ),
        missing_payload_fields = (),
        diagnostics = (
            source = :nested_endcap_panel_owned_units,
            payload_ready_for_current_metric_execution = true,
        ),
    )
    geometry = (
        source = :nested_endcap_panel_owned_units,
        current_box = owned_units.current_box,
        inner_box = owned_units.inner_box,
        bond_axis = owned_units.bond_axis,
        q = owned_units.q,
        L = owned_units.L,
        unit_count = length(owned_units.units),
        support_contract = owned_units.support_contract,
        coefficient_contract = owned_units.coefficient_contract,
        retained_count = owned_units.audit.retained_count,
    )
    diagnostics = (
        metadata_only = true,
        transitional_current_active_implementation = true,
        coverage_ok = owned_units.audit.coverage_ok,
        support_count = length(owned_units.expected_support_indices),
        owned_unit_count = length(owned_units.units),
        coefficient_map_generated = false,
        metric_packet_generated = false,
        fixed_block_sidecar_installed = false,
    )
    return CartesianShellRegion3D(
        :endcap_panel_shared_exterior,
        :shared_endcap_panel_shell_layer,
        :transitional,
        owned_units.current_box,
        owned_units.inner_box,
        copy(owned_units.expected_support_indices),
        _shell_region_support_summary(owned_units.expected_support_indices; parent_dimension),
        :boundary_only,
        geometry,
        retention,
        active_builder_consumes,
        false,
        !active_builder_consumes,
        diagnostics,
        (
            source = :nested_endcap_panel_owned_units,
            support_contract = owned_units.support_contract,
            coefficient_contract = owned_units.coefficient_contract,
        ),
    )
end

function cartesian_shell_region(
    layer::_CartesianNestedEndcapPanelShellLayer3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    region = cartesian_shell_region(
        layer.owned_units;
        parent_dimension,
        active_builder_consumes = true,
    )
    return CartesianShellRegion3D(
        region.region_family,
        region.role,
        region.status,
        region.box,
        region.inner_exclusion_box,
        copy(layer.support_indices),
        _shell_region_support_summary(layer.support_indices; parent_dimension),
        region.ownership_coverage_contract,
        merge(
            region.geometry,
            (
                source = :nested_endcap_panel_shell_layer,
                retained_count = size(layer.coefficient_matrix, 2),
                packet_available = true,
            ),
        ),
        region.retention,
        region.current_route_consumes,
        region.descriptor_drives_builder,
        region.descriptor_only,
        merge(
            region.diagnostics,
            (
                source_layer_retained_count = size(layer.coefficient_matrix, 2),
                packet_available = true,
            ),
        ),
        (
            source = :nested_endcap_panel_shell_layer,
            layer_provenance = layer.provenance,
        ),
    )
end

function cartesian_shell_region(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    retention = _shell_region_retention(
        retention_rule = :boundary_comx_product_modes_raw_boundary_projection,
        cleanup_rule = :full_rank_symmetric_lowdin,
        preferred_contraction_rule = :pqs_boundary_projection_from_filled_box,
        expected_unit_family = :projected_q_shell_descriptor,
        metric_capability = :pqs_low_order_product_metric_prototype,
        required_payload_fields = (
            :support_indices,
            :support_states,
            :axis_local_coefficients,
            :boundary_mode_indices,
            :cleanup_transform,
            :column_range,
            :fixed_block_sidecar_payload,
        ),
        missing_payload_fields = (
            :column_range,
            :fixed_block_sidecar_payload,
            :product_staged_metric_payload,
        ),
        diagnostics = (
            source = :projected_q_shell_staged_unit_descriptor,
            payload_ready_for_current_metric_execution = false,
        ),
    )
    geometry = (
        source = :projected_q_shell_staged_unit_descriptor,
        current_box = descriptor.current_box,
        inner_box = descriptor.inner_box,
        bond_axis = descriptor.bond_axis,
        q = descriptor.q,
        L = descriptor.L,
        selection_rule = descriptor.selection_rule,
        axis_intervals = descriptor.axis_intervals,
        boundary_mode_count = descriptor.mode_count,
        boundary_column_count = length(descriptor.boundary_column_indices),
        cleanup_rank_count = descriptor.cleanup_rank_count,
        cleanup_rank_drop_count = descriptor.cleanup_rank_drop_count,
    )
    diagnostics = (
        metadata_only = true,
        prototype_only = true,
        support_count = descriptor.support_count,
        retained_count = descriptor.retained_count,
        mode_count = descriptor.mode_count,
        coefficient_map_generated = false,
        metric_packet_generated = false,
        fixed_block_sidecar_installed =
            descriptor.active_consumption.fixed_block_sidecar_installed,
        product_doside_unit = false,
        active_consumption = descriptor.active_consumption,
    )
    return CartesianShellRegion3D(
        :projected_q_shell_boundary_modes,
        descriptor.role,
        :prototype,
        descriptor.current_box,
        descriptor.inner_box,
        copy(descriptor.support_indices),
        _shell_region_support_summary(descriptor.support_indices; parent_dimension),
        :boundary_only,
        geometry,
        retention,
        false,
        false,
        true,
        diagnostics,
        (
            source = :projected_q_shell_staged_unit_descriptor,
            original_diagnostics = descriptor.diagnostics,
        ),
    )
end

function _shell_region_status_counts(regions)
    counts = Dict{Symbol,Int}()
    for region in regions
        counts[region.status] = get(counts, region.status, 0) + 1
    end
    return (; (key => counts[key] for key in sort(collect(keys(counts))))...)
end

function _shell_region_inventory_support_summary(
    regions;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    counts = [
        region.support_summary.entry_count
        for region in regions
        if !isnothing(region.support_summary.entry_count)
    ]
    total_count = isempty(counts) ? nothing : sum(counts)
    complete =
        !isnothing(parent_dimension) && !isnothing(total_count) ?
        Int(parent_dimension) == Int(total_count) : nothing
    all_count_summaries = all(!region.support_summary.support_indices_available for region in regions)
    return (
        parent_dimension = isnothing(parent_dimension) ? nothing : Int(parent_dimension),
        region_support_entry_count = total_count,
        support_complete_by_region_counts = complete,
        support_indices_available_for_all_regions =
            all(region.support_summary.support_indices_available for region in regions),
        count_only_summaries_for_all_regions = all_count_summaries,
    )
end

function cartesian_shell_region_inventory(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    effective_parent_dimension =
        isnothing(parent_dimension) ? length(construction.sequence.support_indices) : parent_dimension
    regions = [
        cartesian_shell_region(build; parent_dimension = effective_parent_dimension)
        for build in construction.region_builds
    ]
    return CartesianShellRegionInventory3D(
        :bond_aligned_diatomic_high_order_recipe_source_construction,
        Tuple(regions),
        length(regions),
        Symbol[region.role for region in regions],
        _shell_region_status_counts(regions),
        _shell_region_inventory_support_summary(
            regions;
            parent_dimension = effective_parent_dimension,
        ),
        count(region -> region.current_route_consumes, regions),
        count(region -> region.descriptor_only, regions),
        (
            metadata_only = true,
            route_uses_existing_builder = true,
            descriptor_driven_builder_count =
                count(region -> region.descriptor_drives_builder, regions),
            coefficient_maps_changed = false,
            metric_execution_changed = false,
        ),
        (
            source = :bond_aligned_diatomic_high_order_recipe_source_construction,
            construction_metadata = construction.metadata,
        ),
    )
end

function _staged_axis_rule_summary(axis)
    return (
        kind = axis.kind,
        fixed_index = axis.fixed_index,
        interval = axis.interval,
        coefficient_shape = size(axis.coefficient_matrix),
    )
end

function _product_staged_rule_family(kind::Symbol)
    kind == :product_doside && return :product_owned_unit
    kind == :support_dense && return :support_dense_fallback
    return :staged_unit
end

function _product_staged_transform_rule(kind::Symbol)
    kind == :product_doside && return :two_active_axis_product_doside
    kind == :support_dense && return :explicit_support_dense_coefficients
    return :staged_unit_coefficients
end

function _product_staged_cleanup_rule(kind::Symbol)
    kind == :product_doside && return :locally_orthonormal_product_doside
    kind == :support_dense && return :external_or_already_cleaned
    return :unspecified
end

function _product_staged_metric_capability(kind::Symbol)
    kind == :product_doside && return :product_staged_metric_contraction
    kind == :support_dense && return :support_local_product
    return :support_local_product
end

function _symbol_count_pairs(values)
    counts = Dict{Symbol,Int}()
    for value in values
        counts[value] = get(counts, value, 0) + 1
    end
    return sort!(collect(pairs(counts)); by = pair -> string(first(pair)))
end

_sorted_unique_symbols(values) = sort!(collect(Set(values)); by = string)

function _diagnostic_bool(diagnostics, name::Symbol)
    hasproperty(diagnostics, name) || return false
    return getproperty(diagnostics, name) === true
end

function _rule_metadata_only(rule::CartesianContractionRule3D)
    _diagnostic_bool(rule.diagnostics, :metadata_only) && return true
    hasproperty(rule.diagnostics, :original_diagnostics) || return false
    return _diagnostic_bool(rule.diagnostics.original_diagnostics, :metadata_only)
end

function _rule_prototype_only(rule::CartesianContractionRule3D)
    _diagnostic_bool(rule.diagnostics, :prototype_only) && return true
    rule.metric_capability == :pqs_low_order_product_metric_prototype && return true
    hasproperty(rule.diagnostics, :original_diagnostics) || return false
    return _diagnostic_bool(rule.diagnostics.original_diagnostics, :prototype_only)
end

function cartesian_contraction_rule(
    unit::_CartesianNestedProductStagedByCenterUnit3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    diagnostics = merge(
        (
            source = :nested_product_staged_by_center_unit,
            coefficient_shape = size(unit.coefficient_matrix),
            axis_function_count = length(unit.axis_function_indices),
            metadata_only = false,
            prototype_only = false,
            coefficient_contract = hasproperty(unit.provenance, :coefficient_contract) ?
                                   unit.provenance.coefficient_contract :
                                   nothing,
        ),
        unit.diagnostics,
    )
    local_geometry = (
        axes = map(_staged_axis_rule_summary, unit.axes),
        axis_function_index_count = length(unit.axis_function_indices),
    )
    return CartesianContractionRule3D(
        _product_staged_rule_family(unit.kind),
        unit.kind,
        unit.role,
        copy(unit.support_indices),
        _contraction_rule_support_summary(
            unit.support_indices;
            parent_dimension,
        ),
        local_geometry,
        unit.column_range,
        size(unit.coefficient_matrix, 1),
        length(unit.column_range),
        _product_staged_transform_rule(unit.kind),
        _product_staged_cleanup_rule(unit.kind),
        _product_staged_metric_capability(unit.kind),
        diagnostics,
        (
            source = :nested_product_staged_by_center_sidecar,
            staged_unit = unit,
            original_provenance = unit.provenance,
        ),
    )
end

function cartesian_contraction_rule(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    full_block_dimension = prod(length.(descriptor.current_box))
    diagnostics = (
        source = :projected_q_shell_staged_unit_descriptor,
        support_count = descriptor.support_count,
        boundary_mode_count = descriptor.mode_count,
        retained_count = descriptor.retained_count,
        boundary_column_count = length(descriptor.boundary_column_indices),
        cleanup_rank_count = descriptor.cleanup_rank_count,
        cleanup_rank_drop_count = descriptor.cleanup_rank_drop_count,
        cleanup_cutoff = descriptor.cleanup_cutoff,
        metadata_only = true,
        prototype_only = true,
        contracted_parent_unit_installed = false,
        non_contracts = descriptor.non_contracts,
        active_consumption = descriptor.active_consumption,
        original_diagnostics = descriptor.diagnostics,
    )
    local_geometry = (
        current_box = descriptor.current_box,
        inner_box = descriptor.inner_box,
        bond_axis = descriptor.bond_axis,
        q = descriptor.q,
        L = descriptor.L,
        axis_intervals = descriptor.axis_intervals,
        axis_local_coefficient_shapes = map(size, descriptor.axis_local_coefficients),
        cleanup_matrix_size = descriptor.cleanup_matrix_size,
    )
    return CartesianContractionRule3D(
        :projected_q_shell_boundary_modes,
        descriptor.kind,
        descriptor.role,
        copy(descriptor.support_indices),
        _contraction_rule_support_summary(
            descriptor.support_indices;
            parent_dimension,
        ),
        local_geometry,
        nothing,
        full_block_dimension,
        descriptor.retained_count,
        :boundary_comx_product_modes_raw_boundary_projection,
        :full_rank_symmetric_lowdin,
        :pqs_low_order_product_metric_prototype,
        diagnostics,
        (
            source = :projected_q_shell_staged_unit_descriptor,
            descriptor,
        ),
    )
end

function cartesian_contraction_rule_inventory(
    rules::AbstractVector{<:CartesianContractionRule3D};
    parent_dimension::Integer,
    contracted_dimension::Union{Nothing,Integer} = nothing,
    unit_count::Integer = length(rules),
    every_unit_has_rule_metadata::Bool = false,
    every_unit_rule_derivable::Bool = false,
    provenance = (; source = :cartesian_contraction_rule_collection),
)
    rule_values = collect(rules)
    parent_dim = Int(parent_dimension)
    contracted_dim = isnothing(contracted_dimension) ? nothing : Int(contracted_dimension)
    support_indices = Int[]
    for rule in rule_values
        append!(support_indices, rule.support_indices)
    end
    family_counts = _symbol_count_pairs(rule.rule_family for rule in rule_values)
    metric_capabilities = _sorted_unique_symbols(rule.metric_capability for rule in rule_values)
    metadata_only_count = count(_rule_metadata_only, rule_values)
    prototype_count = count(_rule_prototype_only, rule_values)
    all_rules_have_column_ranges = all(rule -> !isnothing(rule.column_range), rule_values)
    diagnostics = (
        source = :cartesian_contraction_rule_inventory,
        parent_level_unit_inventory =
            every_unit_rule_derivable &&
            Int(unit_count) == length(rule_values) &&
            all_rules_have_column_ranges,
        all_rules_have_column_ranges,
        rule_family_counts = family_counts,
        metric_capabilities,
        metadata_only_rule_count = metadata_only_count,
        prototype_rule_count = prototype_count,
        any_metadata_only_rule = metadata_only_count > 0,
        any_prototype_rule = prototype_count > 0,
        q_shell_rule_present = any(
            rule -> rule.rule_family == :projected_q_shell_boundary_modes,
            rule_values,
        ),
        q_shell_installed_as_contracted_parent_unit = any(
            rule ->
                rule.rule_family == :projected_q_shell_boundary_modes &&
                !isnothing(rule.column_range),
            rule_values,
        ),
    )
    return CartesianContractionRuleInventory3D(
        parent_dim,
        contracted_dim,
        Int(unit_count),
        length(rule_values),
        rule_values,
        family_counts,
        metric_capabilities,
        sum(rule.source_dimension for rule in rule_values),
        sum(rule.retained_dimension for rule in rule_values),
        _contraction_rule_support_summary(support_indices; parent_dimension = parent_dim),
        Any[rule.support_summary for rule in rule_values],
        every_unit_has_rule_metadata,
        every_unit_rule_derivable,
        metadata_only_count,
        prototype_count,
        metadata_only_count > 0,
        prototype_count > 0,
        diagnostics,
        provenance,
    )
end

"""
    CartesianContractionUnit3D

Internal provenance record for a group of contracted-parent columns.

The global coefficient matrix on `CartesianContractedParent3D` remains the
source of truth. Unit support and column ranges are structural metadata only:
support may overlap, may be incomplete, and does not imply orthonormality.
"""
struct CartesianContractionUnit3D{M}
    role::Symbol
    support_indices::Vector{Int}
    column_range::UnitRange{Int}
    metadata::M
end

function CartesianContractionUnit3D(
    role::Symbol,
    support_indices::AbstractVector{<:Integer},
    column_range::UnitRange{<:Integer};
    metadata = (;),
)
    isempty(column_range) && throw(
        ArgumentError("Cartesian contraction units require a non-empty column range"),
    )
    return CartesianContractionUnit3D(
        role,
        Int[Int(index) for index in support_indices],
        Int(first(column_range)):Int(last(column_range)),
        metadata,
    )
end

function cartesian_contraction_unit_from_rule(
    rule::CartesianContractionRule3D,
    payload::_CartesianNestedProductStagedByCenterUnit3D;
    metadata = (;),
)
    rule.kind == payload.kind || throw(
        ArgumentError("contraction rule kind $(rule.kind) does not match staged payload kind $(payload.kind)"),
    )
    rule.role == payload.role || throw(
        ArgumentError("contraction rule role $(rule.role) does not match staged payload role $(payload.role)"),
    )
    rule.support_indices == payload.support_indices || throw(
        ArgumentError("contraction rule support does not match staged payload support"),
    )
    rule.column_range == payload.column_range || throw(
        ArgumentError("contraction rule column range does not match staged payload column range"),
    )
    rule.retained_dimension == length(payload.column_range) || throw(
        ArgumentError("contraction rule retained dimension does not match staged payload column range"),
    )
    return CartesianContractionUnit3D(
        payload.role,
        payload.support_indices,
        payload.column_range;
        metadata = merge(
            (
                source = :nested_product_staged_by_center_sidecar,
                rule_driven_unit_creation = true,
                rule_family = rule.rule_family,
                rule_kind = rule.kind,
                rule_metric_capability = rule.metric_capability,
                staged_by_center_kind = payload.kind,
                staged_by_center_unit = payload,
                contraction_rule = rule,
            ),
            metadata,
        ),
    )
end

function cartesian_contraction_rule(
    unit::CartesianContractionUnit3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    if hasproperty(unit.metadata, :contraction_rule)
        return unit.metadata.contraction_rule
    elseif hasproperty(unit.metadata, :staged_by_center_unit)
        return cartesian_contraction_rule(
            unit.metadata.staged_by_center_unit;
            parent_dimension,
        )
    end
    diagnostics = (
        source = :cartesian_contraction_unit,
        metadata = unit.metadata,
        metadata_only = false,
        prototype_only = false,
    )
    return CartesianContractionRule3D(
        :support_dense_fallback,
        :support_dense,
        unit.role,
        copy(unit.support_indices),
        _contraction_rule_support_summary(
            unit.support_indices;
            parent_dimension,
        ),
        (;),
        unit.column_range,
        length(unit.support_indices),
        length(unit.column_range),
        :explicit_support_dense_coefficients,
        :external_or_already_cleaned,
        :support_local_product,
        diagnostics,
        (source = :cartesian_contraction_unit,),
    )
end

"""
    CartesianContractedParent3D

Internal identity object for columns formed as linear combinations of a full
Cartesian parent gausslet lattice.

This object deliberately stores no backend state, overlap/H/V matrices,
Gaussian supplements, residual columns, or operator packets. The coefficient
matrix is the source of truth; units are provenance only.
"""
struct CartesianContractedParent3D{P<:CartesianParentGaussletBasis3D,C<:_CartesianCoefficientMap,U<:AbstractVector,M}
    parent::P
    coefficient_matrix::C
    units::U
    metadata::M
end

function _validate_unit_column_ranges(
    units::AbstractVector{<:CartesianContractionUnit3D},
    contracted_dimension::Int,
)
    for unit in units
        first(unit.column_range) >= 1 || throw(
            ArgumentError("contraction unit column ranges must start inside contracted columns"),
        )
        last(unit.column_range) <= contracted_dimension || throw(
            ArgumentError("contraction unit column ranges must lie inside contracted columns"),
        )
    end
    return units
end

function CartesianContractedParent3D(
    parent::CartesianParentGaussletBasis3D,
    coefficients::AbstractMatrix{<:Real};
    units::AbstractVector{<:CartesianContractionUnit3D} = CartesianContractionUnit3D[],
    metadata = (;),
)
    coefficient_matrix = _cartesian_coefficient_map_storage(coefficients)
    size(coefficient_matrix, 1) == parent_dimension(parent) || throw(
        DimensionMismatch("contracted parent coefficient rows must match parent dimension"),
    )
    unit_values = CartesianContractionUnit3D[units...]
    _validate_unit_column_ranges(unit_values, size(coefficient_matrix, 2))
    return CartesianContractedParent3D(
        parent,
        coefficient_matrix,
        unit_values,
        metadata,
    )
end

function _contracted_parent_units_from_staged_sidecar(sidecar)
    if hasproperty(sidecar, :units)
        parent_dim = hasproperty(sidecar, :dims) ? prod(sidecar.dims) : nothing
        return CartesianContractionUnit3D[
            cartesian_contraction_unit_from_rule(
                cartesian_contraction_rule(unit; parent_dimension = parent_dim),
                unit,
            ) for unit in sidecar.units
        ]
    elseif hasproperty(sidecar, :block_column_ranges)
        return CartesianContractionUnit3D[
            CartesianContractionUnit3D(
                Symbol(:staged_block_, block_index),
                sidecar.block_support_indices[block_index],
                sidecar.block_column_ranges[block_index];
                metadata = (
                    source = :nested_staged_by_center_sidecar,
                    staged_by_center_kind = :support_dense,
                    staged_by_center_block_index = block_index,
                    staged_by_center_coefficients = sidecar.block_coefficients[block_index],
                    staged_by_center_support_states = sidecar.block_support_states[block_index],
                ),
            ) for block_index in eachindex(sidecar.block_column_ranges)
        ]
    end
    throw(ArgumentError("unsupported staged by-center sidecar for Cartesian contracted parent adapter"))
end

function CartesianContractedParent3D(
    fixed_block::_NestedFixedBlock3D;
    metadata = (source = :nested_fixed_block,),
)
    parent = cartesian_parent_gausslet_basis(fixed_block)
    coefficients = fixed_block.coefficient_matrix
    staged_sidecar = _nested_staged_by_center_sidecar(fixed_block)
    if !isnothing(staged_sidecar)
        path = hasproperty(staged_sidecar, :units) ? :product_staged_factorized : :staged_factorized
        return CartesianContractedParent3D(
            parent,
            coefficients;
            units = _contracted_parent_units_from_staged_sidecar(staged_sidecar),
            metadata = merge(
                metadata,
                (
                    staged_by_center_path = path,
                    staged_by_center_sidecar = staged_sidecar,
                ),
            ),
        )
    end
    unit = CartesianContractionUnit3D(
        :nested_fixed_block,
        fixed_block.support_indices,
        1:size(coefficients, 2);
        metadata = (source = :nested_fixed_block,),
    )
    return CartesianContractedParent3D(
        parent,
        coefficients;
        units = [unit],
        metadata,
    )
end

cartesian_contracted_parent(
    parent::CartesianParentGaussletBasis3D,
    coefficients::AbstractMatrix{<:Real};
    kwargs...,
) = CartesianContractedParent3D(parent, coefficients; kwargs...)

cartesian_contracted_parent(fixed_block::_NestedFixedBlock3D; kwargs...) =
    CartesianContractedParent3D(fixed_block; kwargs...)

contracted_parent_basis(parent::CartesianContractedParent3D) = parent.parent
contracted_parent_coefficients(parent::CartesianContractedParent3D) = parent.coefficient_matrix
contracted_parent_units(parent::CartesianContractedParent3D) = parent.units
contracted_parent_metadata(parent::CartesianContractedParent3D) = parent.metadata
contracted_parent_parent_dimension(parent::CartesianContractedParent3D) =
    parent_dimension(parent.parent)
contracted_parent_dimension(parent::CartesianContractedParent3D) =
    size(parent.coefficient_matrix, 2)

function contracted_parent_contraction_rules(parent::CartesianContractedParent3D)
    parent_dim = contracted_parent_parent_dimension(parent)
    return [
        contraction_unit_rule(unit; parent_dimension = parent_dim) for unit in parent.units
    ]
end

function contracted_parent_rule_inventory(parent::CartesianContractedParent3D)
    units = contracted_parent_units(parent)
    return cartesian_contraction_rule_inventory(
        contracted_parent_contraction_rules(parent);
        parent_dimension = contracted_parent_parent_dimension(parent),
        contracted_dimension = contracted_parent_dimension(parent),
        unit_count = length(units),
        every_unit_has_rule_metadata = all(
            unit -> hasproperty(unit.metadata, :contraction_rule),
            units,
        ),
        every_unit_rule_derivable = true,
        provenance = (
            source = :cartesian_contracted_parent,
            contracted_parent_metadata = parent.metadata,
        ),
    )
end

contraction_unit_role(unit::CartesianContractionUnit3D) = unit.role
contraction_unit_support_indices(unit::CartesianContractionUnit3D) = unit.support_indices
contraction_unit_column_range(unit::CartesianContractionUnit3D) = unit.column_range
contraction_unit_metadata(unit::CartesianContractionUnit3D) = unit.metadata
contraction_unit_rule(
    unit::CartesianContractionUnit3D;
    parent_dimension::Union{Nothing,Int} = nothing,
) = cartesian_contraction_rule(unit; parent_dimension)

contraction_rule_family(rule::CartesianContractionRule3D) = rule.rule_family
contraction_rule_kind(rule::CartesianContractionRule3D) = rule.kind
contraction_rule_support_summary(rule::CartesianContractionRule3D) = rule.support_summary
contraction_rule_column_range(rule::CartesianContractionRule3D) = rule.column_range
contraction_rule_source_dimension(rule::CartesianContractionRule3D) = rule.source_dimension
contraction_rule_retained_dimension(rule::CartesianContractionRule3D) = rule.retained_dimension
contraction_rule_transform_rule(rule::CartesianContractionRule3D) = rule.transform_rule
contraction_rule_cleanup_rule(rule::CartesianContractionRule3D) = rule.cleanup_rule
contraction_rule_metric_capability(rule::CartesianContractionRule3D) = rule.metric_capability

contracted_parent_unit_column_ranges(parent::CartesianContractedParent3D) =
    UnitRange{Int}[unit.column_range for unit in parent.units]

contracted_parent_unit_support_indices(parent::CartesianContractedParent3D) =
    Vector{Int}[copy(unit.support_indices) for unit in parent.units]

function contracted_parent_support_indices(parent::CartesianContractedParent3D)
    indices = Int[]
    for unit in parent.units
        append!(indices, unit.support_indices)
    end
    return indices
end

struct CartesianContractedParentStructuralAudit
    parent_dimension::Int
    contracted_dimension::Int
    unit_count::Int
    support_entry_count::Int
    unique_support_count::Int
    duplicate_support_count::Int
    missing_support_count::Int
    outside_support_count::Int
    column_entry_count::Int
    unique_column_count::Int
    duplicate_column_count::Int
    missing_column_count::Int
    outside_column_count::Int
    support_complete::Bool
    column_ranges_cover_contract::Bool
    structural_ok::Bool
end

function _entry_counts(values::AbstractVector{Int}, valid_range::UnitRange{Int})
    inside = Int[value for value in values if value in valid_range]
    outside_count = length(values) - length(inside)
    unique_inside = unique(inside)
    duplicate_count = length(inside) - length(unique_inside)
    missing_count = length(valid_range) - length(unique_inside)
    return (
        entry_count = length(values),
        unique_count = length(unique_inside),
        duplicate_count = duplicate_count,
        missing_count = missing_count,
        outside_count = outside_count,
    )
end

function contracted_parent_structural_audit(parent::CartesianContractedParent3D)
    parent_dim = contracted_parent_parent_dimension(parent)
    contracted_dim = contracted_parent_dimension(parent)
    support_counts = _entry_counts(
        contracted_parent_support_indices(parent),
        1:parent_dim,
    )
    columns = Int[]
    for unit in parent.units
        append!(columns, collect(unit.column_range))
    end
    column_counts = _entry_counts(columns, 1:contracted_dim)
    support_complete = support_counts.missing_count == 0 && support_counts.outside_count == 0
    column_ranges_cover_contract =
        column_counts.missing_count == 0 &&
        column_counts.outside_count == 0 &&
        column_counts.duplicate_count == 0
    structural_ok = support_counts.outside_count == 0 && column_ranges_cover_contract
    return CartesianContractedParentStructuralAudit(
        parent_dim,
        contracted_dim,
        length(parent.units),
        support_counts.entry_count,
        support_counts.unique_count,
        support_counts.duplicate_count,
        support_counts.missing_count,
        support_counts.outside_count,
        column_counts.entry_count,
        column_counts.unique_count,
        column_counts.duplicate_count,
        column_counts.missing_count,
        column_counts.outside_count,
        support_complete,
        column_ranges_cover_contract,
        structural_ok,
    )
end

end
