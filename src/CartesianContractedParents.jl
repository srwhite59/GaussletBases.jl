module CartesianContractedParents

using LinearAlgebra: I, norm

import ..GaussletBases: _CartesianCoefficientMap,
                         _BondAlignedDiatomicAtomGrowthConstructionRegion3D,
                         _BondAlignedDiatomicHighOrderRecipeRegionSourceBuild3D,
                         _BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
                         _CartesianNestedEndcapPanelOwnedUnits3D,
                         _CartesianNestedEndcapPanelShellLayer3D,
                         _CartesianNestedProductStagedByCenterSidecar3D,
                         _CartesianNestedProductStagedByCenterUnit3D,
                         _CartesianNestedProjectedQShellLayer3D,
                         _CartesianNestedProjectedQShellStagedUnitDescriptor3D,
                         _NestedFixedBlock3D,
                         _cartesian_coefficient_map_storage,
                         _cartesian_unflat_index,
                         _nested_parent_axis_counts,
                         _nested_projected_q_shell_descriptor_seed_coefficients,
                         _nested_projected_q_shell_staged_unit_descriptor,
                         _nested_product_staged_generic_unit,
                         _nested_product_staged_unit_from_owned_unit,
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

"""
    _CartesianResolvedContractionPayload3D

Private execution-facing normalization record for an already-existing
contraction payload. Region and rule records validate intent; this object is
only a thin description of a concrete payload that an existing metric kernel can
consume. It does not generate coefficient maps or install new sidecars.
"""
struct _CartesianResolvedContractionPayload3D{S,H,D,P}
    metric_path::Symbol
    ready_for_metric_execution::Bool
    payload_kind::Symbol
    column_range::Union{Nothing,UnitRange{Int}}
    support_indices::Union{Nothing,Vector{Int}}
    support_states::S
    payload::H
    missing_fields::Tuple{Vararg{Symbol}}
    diagnostics::D
    provenance::P
end

"""
    _CartesianPacketBuildSource3D

Private metadata-only description of the resolved payloads a future Cartesian
packet builder would consume. It records dimensions, coverage, packet-field
candidate fields, and axis requirements without checking operator data,
building packet matrices, or changing the existing `_nested_shell_packet(...)`
authority.
"""
struct _CartesianPacketBuildSource3D{R,C,S,K,A,M,D,P}
    parent_dimension::Int
    contracted_dimension::Int
    resolved_payloads::R
    column_ranges::Vector{UnitRange{Int}}
    column_coverage::C
    support_union_summary::S
    payload_kind_counts::K
    candidate_packet_fields::Tuple{Vararg{Symbol}}
    missing_packet_fields::Tuple{Vararg{Symbol}}
    axis_operator_requirements::A
    unit_descriptors::M
    diagnostics::D
    provenance::P
end

"""
    _CartesianPacketBuildPlan3D

Private metadata-only construction plan around a packet build source. This is
an audit/planning object only: it must not build overlap, kinetic, local,
Gaussian, MWG, interaction, QW, or Hamiltonian data.
"""
struct _CartesianPacketBuildPlan3D{S,D,P}
    source::S
    current_builder_authority::Symbol
    descriptor_drives_builder::Bool
    numerical_packet_matrices_built::Bool
    diagnostics::D
    provenance::P
end

"""
    _CartesianExecutableProjectedQShellPayload3D

Private executable fixture for a projected q-shell descriptor. This is not a
fixed-block sidecar and is not consumed by default routes; it is only the first
resolved-payload shape needed to prove PQS metric dispatch contracts. Its
coefficients are support-local boundary-row coefficients, not a full
parent-dimension coefficient map.
"""
struct _CartesianExecutableProjectedQShellPayload3D{C,D,P}
    kind::Symbol
    role::Union{Nothing,Symbol}
    column_range::UnitRange{Int}
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    support_coefficient_matrix::_CartesianCoefficientMap
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D
    current_box::NTuple{3,UnitRange{Int}}
    inner_box::NTuple{3,UnitRange{Int}}
    q::Int
    L::Int
    bond_axis::Symbol
    axis_intervals::NTuple{3,UnitRange{Int}}
    boundary_mode_indices::Vector{NTuple{3,Int}}
    boundary_column_indices::Vector{Int}
    cleanup_transform::C
    cleanup_diagnostics::D
    diagnostics::P
end

"""
    _CartesianProjectedQShellSidecarFixture3D

Private sidecar-shaped fixture for executable PQS payload plumbing. This is a
real container for PQS payloads, but it is deliberately not installed into
fixed blocks or consumed by default metric/QW/Hamiltonian routes.
"""
struct _CartesianProjectedQShellSidecarFixture3D{U,P,D}
    dims::NTuple{3,Int}
    payloads::U
    provenance::P
    diagnostics::D
end

"""
    _CartesianRawProductSource3D

Private metadata-only description of a local raw product source space. This is
the "source" half of the planned raw-product-source -> retained-transform ->
retained-unit contract; it does not build operator matrices or change metric
execution.
"""
struct _CartesianRawProductSource3D{S,P,D}
    source_id::Symbol
    role::Symbol
    parent_dims::NTuple{3,Int}
    local_box::NTuple{3,UnitRange{Int}}
    axis_intervals::NTuple{3,UnitRange{Int}}
    source_dimension::Int
    support_indices::Union{Nothing,Vector{Int}}
    support_summary::S
    raw_source_weight_role::Symbol
    provenance::P
    diagnostics::D
end

"""
    _CartesianRetainedTransform3D

Private metadata-only retained transform description for one raw product
source. `transform_matrix === nothing` means the full raw-to-retained
transform is not materialized; callers must inspect `transform_kind`,
`transform_stages`, and diagnostics such as
`full_raw_to_retained_matrix_materialized`. Non-identity transforms must carry
an explicit weight role so IDA-like weight division is never inferred silently.
"""
struct _CartesianRetainedTransform3D{M,C,P,D}
    source_id::Symbol
    retained_dimension::Int
    transform_kind::Symbol
    transform_matrix::M
    transform_stages::Tuple{Vararg{Symbol}}
    cleanup_diagnostics::C
    retained_column_weight_role::Symbol
    provenance::P
    diagnostics::D
end

"""
    _CartesianRawProductSourcePairOperatorPacket3D

Placeholder for future source_i/source_j raw operator blocks. This pass records
shape, symmetry, term capability, and provenance only; `operator_matrices` must
remain `nothing`.
"""
struct _CartesianRawProductSourcePairOperatorPacket3D{P,D}
    left_source_id::Symbol
    right_source_id::Symbol
    operator_kind::Symbol
    supported_terms::Tuple{Vararg{Symbol}}
    symmetry_status::Symbol
    backend::Symbol
    operator_matrices::Nothing
    provenance::P
    diagnostics::D
end

"""
    _CartesianRawProductSourcePairPlan3D

Private all-pairs planning record for raw product source operator packets. The
plan owns source/transform lookup tables and upper-triangular placeholder pair
packets, but it does not build raw or retained operator matrices.
"""
struct _CartesianRawProductSourcePairPlan3D{S,T,P,D}
    operator_kind::Symbol
    supported_terms::Tuple{Vararg{Symbol}}
    symmetry_status::Symbol
    source_ids::Vector{Symbol}
    raw_sources::S
    retained_transforms::T
    pair_packets::P
    pair_keys::Vector{NTuple{2,Symbol}}
    diagnostics::D
end

"""
    _CartesianResolvedRawProductSourcePair3D

Private dereferenced view of one raw-source pair-packet inside a pair plan. It
links the packet ids back to the plan-owned raw source and retained-transform
records and records only planning/audit diagnostics.
"""
struct _CartesianResolvedRawProductSourcePair3D{L,R,LT,RT,P,D}
    pair_key::NTuple{2,Symbol}
    pair_index::Int
    left_raw_source::L
    right_raw_source::R
    left_retained_transform::LT
    right_retained_transform::RT
    pair_packet::P
    diagnostics::D
end

"""
    _CartesianRawProductSourcePairPlanAudit3D

Private metadata-only audit of a raw-source pair plan. This confirms lookup,
upper-triangle, placeholder, and weight-role contracts without building raw or
retained operator matrices.
"""
struct _CartesianRawProductSourcePairPlanAudit3D{R,D}
    resolved_pairs::R
    diagnostics::D
end

"""
    _CartesianRawProductSourceLowOrderOperatorPacket3D

Private fixture/reference packet for one raw low-order operator block. This is
not installed in pair plans and does not apply retained transforms; it only
records a raw source-space matrix for a deliberately scoped helper path.
"""
struct _CartesianRawProductSourceLowOrderOperatorPacket3D{M,P,D}
    left_source_id::Symbol
    right_source_id::Symbol
    operator_kind::Symbol
    term::Symbol
    source_dimensions::NTuple{2,Int}
    symmetry_status::Symbol
    backend::Symbol
    raw_operator_matrix::M
    provenance::P
    diagnostics::D
end

"""
    _CartesianRetainedLowOrderOperatorBlock3D

Private retained-space block built from one raw low-order packet and explicit
materialized retained transforms. This is fixture/reference-only and is not a
metric execution or operator assembly path.
"""
struct _CartesianRetainedLowOrderOperatorBlock3D{M,P,D}
    left_source_id::Symbol
    right_source_id::Symbol
    operator_kind::Symbol
    term::Symbol
    retained_dimensions::NTuple{2,Int}
    symmetry_status::Symbol
    retained_operator_matrix::M
    provenance::P
    diagnostics::D
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

function _resolved_payload_path(rule::CartesianContractionRule3D)
    if rule.rule_family == :product_owned_unit &&
       rule.kind == :product_doside &&
       rule.metric_capability == :product_staged_metric_contraction
        return (
            metric_path = :product_staged_metric_contraction,
            linear_vector_path = :product_staged_axis_projection,
            block_role = :product,
            unsupported = false,
            prototype = false,
        )
    elseif rule.rule_family == :support_dense_fallback &&
           rule.kind == :support_dense &&
           rule.metric_capability == :support_local_product
        return (
            metric_path = :support_local_product,
            linear_vector_path = :support_local_fallback,
            block_role = :fallback,
            unsupported = false,
            prototype = false,
        )
    end
    prototype = _rule_prototype_only(rule)
    return (
        metric_path = prototype ? :unsupported_prototype : :unsupported,
        linear_vector_path = :unsupported,
        block_role = prototype ? :unsupported_prototype : :unsupported,
        unsupported = true,
        prototype = prototype,
    )
end

function _validate_resolved_payload_rule(
    rule::CartesianContractionRule3D,
    payload::_CartesianNestedProductStagedByCenterUnit3D,
)
    rule.kind == payload.kind || throw(
        ArgumentError("resolved payload kind $(payload.kind) does not match contraction rule kind $(rule.kind)"),
    )
    rule.role == payload.role || throw(
        ArgumentError("resolved payload role $(payload.role) does not match contraction rule role $(rule.role)"),
    )
    rule.support_indices == payload.support_indices || throw(
        ArgumentError("resolved payload support does not match contraction rule support"),
    )
    rule.column_range == payload.column_range || throw(
        ArgumentError("resolved payload column range does not match contraction rule column range"),
    )
    rule.retained_dimension == length(payload.column_range) || throw(
        ArgumentError("resolved payload retained dimension does not match its column range"),
    )
    return nothing
end

function _cartesian_resolved_contraction_payload(
    payload::_CartesianNestedProductStagedByCenterUnit3D;
    parent_dimension::Union{Nothing,Int} = nothing,
    rule::Union{Nothing,CartesianContractionRule3D} = nothing,
)
    resolved_rule = isnothing(rule) ?
                    cartesian_contraction_rule(payload; parent_dimension) :
                    rule
    _validate_resolved_payload_rule(resolved_rule, payload)
    path = _resolved_payload_path(resolved_rule)
    missing_fields = path.unsupported ? (:supported_product_or_support_dense_payload,) : ()
    return _CartesianResolvedContractionPayload3D(
        path.metric_path,
        !path.unsupported,
        payload.kind,
        payload.column_range,
        copy(payload.support_indices),
        copy(payload.support_states),
        payload,
        missing_fields,
        (
            source = :nested_product_staged_by_center_unit,
            rule_family = resolved_rule.rule_family,
            rule_kind = resolved_rule.kind,
            metric_capability = resolved_rule.metric_capability,
            linear_vector_path = path.linear_vector_path,
            block_role = path.block_role,
            unsupported = path.unsupported,
            prototype = path.prototype,
            coefficient_shape = size(payload.coefficient_matrix),
            payload_ready_for_current_metric_execution = !path.unsupported,
        ),
        (
            source = :nested_product_staged_by_center_unit,
            contraction_rule = resolved_rule,
        ),
    )
end

function _cartesian_resolved_contraction_payload(rule::CartesianContractionRule3D)
    path = _resolved_payload_path(rule)
    missing_fields =
        path.prototype ? (:installed_executable_payload, :fixed_block_sidecar_payload) :
        (:executable_payload,)
    return _CartesianResolvedContractionPayload3D(
        path.metric_path,
        false,
        rule.kind,
        rule.column_range,
        copy(rule.support_indices),
        nothing,
        nothing,
        missing_fields,
        (
            source = :cartesian_contraction_rule,
            rule_family = rule.rule_family,
            rule_kind = rule.kind,
            metric_capability = rule.metric_capability,
            linear_vector_path = path.linear_vector_path,
            block_role = path.block_role,
            unsupported = true,
            prototype = path.prototype,
            payload_ready_for_current_metric_execution = false,
        ),
        (
            source = :cartesian_contraction_rule,
            contraction_rule = rule,
        ),
    )
end

function _cartesian_resolved_contraction_payload(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    rule = cartesian_contraction_rule(descriptor; parent_dimension)
    return _CartesianResolvedContractionPayload3D(
        :unsupported_prototype,
        false,
        descriptor.kind,
        nothing,
        copy(descriptor.support_indices),
        copy(descriptor.support_states),
        descriptor,
        (
            :column_range,
            :fixed_block_sidecar_payload,
            :product_staged_metric_payload,
        ),
        (
            source = :projected_q_shell_staged_unit_descriptor,
            rule_family = rule.rule_family,
            rule_kind = rule.kind,
            metric_capability = rule.metric_capability,
            linear_vector_path = :unsupported,
            block_role = :unsupported_prototype,
            unsupported = true,
            prototype = true,
            payload_ready_for_current_metric_execution = false,
            prototype_metric_capability = :pqs_low_order_product_metric_prototype,
        ),
        (
            source = :projected_q_shell_staged_unit_descriptor,
            contraction_rule = rule,
        ),
    )
end

function _cartesian_executable_projected_q_shell_payload_fixture(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D;
    column_range::UnitRange{Int},
    parent_dimension::Union{Nothing,Int} = nothing,
)
    descriptor.kind == :projected_q_shell || throw(
        ArgumentError("projected q-shell executable payload fixture requires a projected_q_shell descriptor"),
    )
    length(column_range) == descriptor.retained_count || throw(
        DimensionMismatch("projected q-shell executable payload column range must match retained count"),
    )
    first(column_range) >= 1 || throw(
        ArgumentError("projected q-shell executable payload column range must start at a positive index"),
    )
    descriptor.cleanup_rank_drop_count == 0 || throw(
        ArgumentError("projected q-shell executable payload requires full-rank symmetric Lowdin cleanup"),
    )
    size(descriptor.cleanup_transform) == descriptor.cleanup_matrix_size || throw(
        DimensionMismatch("projected q-shell executable payload cleanup transform has inconsistent dimensions"),
    )
    descriptor.cleanup_matrix_size == (descriptor.mode_count, descriptor.retained_count) ||
        throw(
            DimensionMismatch("projected q-shell executable payload cleanup dimensions must match mode/retained counts"),
        )
    support_summary = _contraction_rule_support_summary(
        descriptor.support_indices;
        parent_dimension,
    )
    support_summary.outside_count in (nothing, 0) || throw(
        ArgumentError("projected q-shell executable payload support lies outside the parent dimension"),
    )
    seed = _nested_projected_q_shell_descriptor_seed_coefficients(descriptor)
    size(seed) == (descriptor.support_count, descriptor.mode_count) || throw(
        DimensionMismatch("projected q-shell seed coefficient shape does not match descriptor metadata"),
    )
    support_coefficients = Matrix{Float64}(seed * descriptor.cleanup_transform)
    size(support_coefficients) == (descriptor.support_count, descriptor.retained_count) ||
        throw(
            DimensionMismatch("projected q-shell executable payload coefficient shape does not match descriptor metadata"),
        )
    all(isfinite, support_coefficients) || throw(
        ArgumentError("projected q-shell executable payload coefficients must be finite"),
    )
    diagnostics = (
        source = :projected_q_shell_executable_payload_fixture,
        fixture_only = true,
        production_supported = false,
        coefficient_scope = :support_local_boundary_rows,
        parent_dimension_coefficient_map = false,
        fixed_block_sidecar_installed = false,
        default_builder_consumes = false,
        metric_packet_consumes = false,
        support_summary = support_summary,
        support_count = descriptor.support_count,
        support_coefficient_shape = size(support_coefficients),
        mode_count = descriptor.mode_count,
        retained_count = descriptor.retained_count,
        cleanup_method = descriptor.cleanup_method,
        cleanup_rank_count = descriptor.cleanup_rank_count,
        cleanup_rank_drop_count = descriptor.cleanup_rank_drop_count,
        cleanup_cutoff = descriptor.cleanup_cutoff,
        metric_capability = :pqs_low_order_support_local_reference,
        retained_column_weight_role = :debug_reference_only,
        retained_weight_semantics = :debug_reference_only,
        retained_weights_used_for_ida_division = false,
        ida_weight_division_allowed = false,
        quadrature_weight_semantics_claimed = false,
        active_interaction_path = :none_fixture_only,
        supported_metric_terms = (:overlap, :weights, :position_x, :position_y, :position_z),
        unsupported_metric_terms = (
            :kinetic,
            :x2,
            :nuclear_one_body,
            :gaussian_sum,
            :pair_sum,
            :interaction,
        ),
        pqs_product_optimized_path_ready = false,
    )
    return _CartesianExecutableProjectedQShellPayload3D(
        :projected_q_shell,
        descriptor.role,
        column_range,
        copy(descriptor.support_indices),
        copy(descriptor.support_states),
        _cartesian_coefficient_map_storage(support_coefficients),
        descriptor,
        descriptor.current_box,
        descriptor.inner_box,
        descriptor.q,
        descriptor.L,
        descriptor.bond_axis,
        descriptor.axis_intervals,
        copy(descriptor.boundary_mode_indices),
        copy(descriptor.boundary_column_indices),
        Matrix{Float64}(descriptor.cleanup_transform),
        (
            method = descriptor.cleanup_method,
            matrix_size = descriptor.cleanup_matrix_size,
            eigenvalues = Float64[value for value in descriptor.cleanup_eigenvalues],
            rank_count = descriptor.cleanup_rank_count,
            rank_drop_count = descriptor.cleanup_rank_drop_count,
            cutoff = descriptor.cleanup_cutoff,
        ),
        diagnostics,
    )
end

function _cartesian_resolved_contraction_payload(
    payload::_CartesianExecutableProjectedQShellPayload3D;
    parent_dimension::Union{Nothing,Int} = nothing,
    fixed_block_sidecar_installed::Bool = false,
)
    support_summary = _contraction_rule_support_summary(
        payload.support_indices;
        parent_dimension,
    )
    support_summary.outside_count in (nothing, 0) || throw(
        ArgumentError("projected q-shell resolved payload support lies outside the parent dimension"),
    )
    return _CartesianResolvedContractionPayload3D(
        :pqs_low_order_support_local_reference,
        true,
        payload.kind,
        payload.column_range,
        copy(payload.support_indices),
        copy(payload.support_states),
        payload,
        (),
        (
            source = :projected_q_shell_executable_payload_fixture,
            rule_family = :projected_q_shell_boundary_modes,
            rule_kind = payload.kind,
            metric_capability = :pqs_low_order_support_local_reference,
            linear_vector_path = :pqs_support_local_reference,
            block_role = :pqs,
            unsupported = false,
            prototype = false,
            fixture_only = true,
            production_supported = false,
            coefficient_scope = :support_local_boundary_rows,
            parent_dimension_coefficient_map = false,
            retained_column_weight_role = :debug_reference_only,
            retained_weight_semantics = :debug_reference_only,
            retained_weights_used_for_ida_division = false,
            ida_weight_division_allowed = false,
            quadrature_weight_semantics_claimed = false,
            active_interaction_path = :none_fixture_only,
            payload_ready_for_current_metric_execution = true,
            default_builder_consumes = false,
            fixed_block_sidecar_installed = fixed_block_sidecar_installed,
            support_summary = support_summary,
        ),
        (
            source = :projected_q_shell_executable_payload_fixture,
            descriptor = payload.descriptor,
        ),
    )
end

function _cartesian_projected_q_shell_sidecar_fixture(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D;
    column_range::UnitRange{Int},
    dims::NTuple{3,Int},
)
    parent_dimension = prod(dims)
    payload = _cartesian_executable_projected_q_shell_payload_fixture(
        descriptor;
        column_range,
        parent_dimension,
    )
    diagnostics = (
        source = :projected_q_shell_sidecar_fixture,
        fixture_only = true,
        production_supported = false,
        fixed_block_sidecar_installed = false,
        default_builder_consumes = false,
        metric_packet_consumes = false,
        qw_consumes = false,
        hamiltonian_consumes = false,
        payload_count = 1,
        metric_capability = :pqs_low_order_support_local_reference,
        support_local_reference_only = true,
        pqs_product_optimized_path_ready = false,
    )
    return _CartesianProjectedQShellSidecarFixture3D(
        dims,
        [payload],
        (
            source = :projected_q_shell_sidecar_fixture,
            descriptor = descriptor,
            installed_in_fixed_block = false,
            public_api = false,
        ),
        diagnostics,
    )
end

function _cartesian_projected_q_shell_installed_sidecar_fixture(
    sidecar::_CartesianProjectedQShellSidecarFixture3D,
)
    diagnostics = merge(
        sidecar.diagnostics,
        (
            source = :projected_q_shell_fixed_block_sidecar_fixture,
            fixture_only = true,
            fixed_block_sidecar_installed = true,
            by_center_consumes = false,
            default_builder_consumes = false,
            qw_consumes = false,
            hamiltonian_consumes = false,
            production_supported = false,
            metric_capability = :pqs_low_order_support_local_reference,
            pqs_product_optimized_path_ready = false,
        ),
    )
    provenance = merge(
        sidecar.provenance,
        (
            source = :projected_q_shell_fixed_block_sidecar_fixture,
            installed_in_fixed_block = true,
            fixed_block_sidecar_slot = :staged_by_center_sidecar_fixture_only,
            by_center_consumes = false,
            public_api = false,
        ),
    )
    return _CartesianProjectedQShellSidecarFixture3D(
        sidecar.dims,
        sidecar.payloads,
        provenance,
        diagnostics,
    )
end

function _nested_projected_q_shell_sidecar_fixture(fixed_block::_NestedFixedBlock3D)
    sidecar = fixed_block.staged_by_center_sidecar[]
    isnothing(sidecar) && return nothing
    sidecar isa _CartesianProjectedQShellSidecarFixture3D || throw(
        ArgumentError("nested fixed block does not carry a projected q-shell sidecar fixture"),
    )
    dims = _nested_parent_axis_counts(fixed_block.parent_basis)
    sidecar.dims == dims || throw(
        DimensionMismatch("projected q-shell fixed-block sidecar dimensions must match parent basis"),
    )
    ncolumns = size(fixed_block.coefficient_matrix, 2)
    isempty(sidecar.payloads) && throw(
        ArgumentError("projected q-shell fixed-block sidecar requires at least one payload"),
    )
    ranges = UnitRange{Int}[payload.column_range for payload in sidecar.payloads]
    for range in ranges
        first(range) >= 1 && last(range) <= ncolumns || throw(
            ArgumentError("projected q-shell sidecar column range $range exceeds fixed dimension $ncolumns"),
        )
    end
    covered_columns = sort!(reduce(vcat, (collect(range) for range in ranges)))
    covered_columns == collect(1:ncolumns) || throw(
        ArgumentError("projected q-shell fixed-block sidecar payloads must cover every fixed column exactly once"),
    )
    return sidecar
end

function _cartesian_projected_q_shell_fixed_block_sidecar_fixture(
    parent_basis,
    layer::_CartesianNestedProjectedQShellLayer3D;
    gausslet_backend::Symbol = :unknown,
)
    descriptor = _nested_projected_q_shell_staged_unit_descriptor(layer)
    dims = _nested_parent_axis_counts(parent_basis)
    parent_dimension = prod(dims)
    size(layer.coefficient_matrix, 1) == parent_dimension || throw(
        DimensionMismatch("projected q-shell fixed-block fixture parent dimension does not match parent basis"),
    )
    ncolumns = size(layer.coefficient_matrix, 2)
    ncolumns == descriptor.retained_count || throw(
        DimensionMismatch("projected q-shell fixed-block fixture column count must match descriptor retained count"),
    )
    sidecar = _cartesian_projected_q_shell_sidecar_fixture(
        descriptor;
        column_range = 1:ncolumns,
        dims,
    )
    installed_sidecar = _cartesian_projected_q_shell_installed_sidecar_fixture(sidecar)
    packet = layer.packet
    fixed_centers = Matrix{Float64}(undef, ncolumns, 3)
    @inbounds for column in 1:ncolumns
        fixed_centers[column, 1] = packet.position_x[column, column]
        fixed_centers[column, 2] = packet.position_y[column, column]
        fixed_centers[column, 3] = packet.position_z[column, column]
    end
    return _NestedFixedBlock3D(
        parent_basis,
        layer,
        gausslet_backend,
        layer.coefficient_matrix,
        copy(layer.support_indices),
        packet.overlap,
        packet.kinetic,
        packet.position_x,
        packet.position_y,
        packet.position_z,
        packet.x2_x,
        packet.x2_y,
        packet.x2_z,
        packet.weights,
        packet.gaussian_sum,
        packet.pair_sum,
        fixed_centers,
        Base.RefValue{Any}(nothing),
        Base.RefValue{Any}(installed_sidecar),
    )
end

function _cartesian_interval_product_dimension(intervals::NTuple{3,UnitRange{Int}})
    return prod(length(interval) for interval in intervals)
end

function _cartesian_support_box_from_states(states::AbstractVector{<:NTuple{3,Int}})
    isempty(states) && throw(ArgumentError("raw product source support states cannot be empty"))
    return ntuple(axis -> begin
        values = Int[state[axis] for state in states]
        minimum(values):maximum(values)
    end, 3)
end

function _cartesian_raw_product_source(
    payload::_CartesianExecutableProjectedQShellPayload3D;
    source_id::Symbol = :projected_q_shell_raw_product_source,
    parent_dims::Union{Nothing,NTuple{3,Int}} = nothing,
)
    dims = isnothing(parent_dims) ? ntuple(axis -> last(payload.axis_intervals[axis]), 3) : parent_dims
    source_dimension = _cartesian_interval_product_dimension(payload.axis_intervals)
    source_dimension == payload.q * payload.q * payload.L || throw(
        DimensionMismatch("projected q-shell raw source dimension must equal q*q*L"),
    )
    support_summary = (
        support_indices_available = false,
        support_scope = :full_local_product_box,
        raw_product_support_count = source_dimension,
        boundary_support_count = length(payload.support_indices),
        retained_count = length(payload.column_range),
    )
    diagnostics = (
        source = :projected_q_shell_raw_product_source_adapter,
        fixture_only = true,
        production_supported = false,
        raw_source_contract = :full_local_product_block,
        raw_source_weight_role = :raw_source_positive,
        retained_transform_expected = :boundary_projection_lowdin,
        operator_matrices_built = false,
        q = payload.q,
        L = payload.L,
    )
    return _CartesianRawProductSource3D(
        source_id,
        isnothing(payload.role) ? :projected_q_shell : payload.role,
        dims,
        payload.current_box,
        payload.axis_intervals,
        source_dimension,
        nothing,
        support_summary,
        :raw_source_positive,
        (
            source = :projected_q_shell_executable_payload_fixture,
            payload = payload,
            descriptor = payload.descriptor,
        ),
        diagnostics,
    )
end

function _cartesian_retained_transform(
    payload::_CartesianExecutableProjectedQShellPayload3D;
    source_id::Symbol = :projected_q_shell_raw_product_source,
)
    diagnostics = (
        source = :projected_q_shell_retained_transform_adapter,
        fixture_only = true,
        production_supported = false,
        transform_contract = :factored_raw_product_to_retained,
        transform_matrix_scope = :factored_raw_product_to_retained,
        full_raw_to_retained_matrix_materialized = false,
        full_raw_to_retained_matrix_available = false,
        cleanup_stage_matrix_scope = :boundary_mode_to_retained,
        boundary_projection_stage = :raw_product_modes_to_boundary_rows,
        cleanup_stage = :full_rank_symmetric_lowdin_boundary_cleanup,
        raw_source_dimension = payload.q * payload.q * payload.L,
        boundary_support_indices = copy(payload.support_indices),
        boundary_support_count = length(payload.support_indices),
        boundary_mode_indices = copy(payload.boundary_mode_indices),
        boundary_column_indices = copy(payload.boundary_column_indices),
        boundary_mode_count = payload.descriptor.mode_count,
        retained_count = length(payload.column_range),
        cleanup_method = payload.cleanup_diagnostics.method,
        cleanup_rank_count = payload.cleanup_diagnostics.rank_count,
        cleanup_rank_drop_count = payload.cleanup_diagnostics.rank_drop_count,
        retained_column_weight_role = :debug_reference_only,
        retained_weight_semantics = :debug_reference_only,
        retained_weights_used_for_ida_division = false,
        retained_weight_positive_checked = false,
        ida_weight_division_allowed = false,
        quadrature_weight_semantics_claimed = false,
    )
    return _CartesianRetainedTransform3D(
        source_id,
        length(payload.column_range),
        :boundary_projection_lowdin,
        nothing,
        (
            :raw_product_modes,
            :raw_boundary_projection,
            :full_rank_symmetric_lowdin_cleanup,
            :retained_columns,
        ),
        payload.cleanup_diagnostics,
        :debug_reference_only,
        (
            source = :projected_q_shell_executable_payload_fixture,
            payload = payload,
            descriptor = payload.descriptor,
            cleanup_stage_matrix = payload.cleanup_transform,
        ),
        diagnostics,
    )
end

function _cartesian_raw_product_source_retained_transform(
    payload::_CartesianExecutableProjectedQShellPayload3D;
    source_id::Symbol = :projected_q_shell_raw_product_source,
    parent_dims::Union{Nothing,NTuple{3,Int}} = nothing,
)
    return (
        raw_source = _cartesian_raw_product_source(
            payload;
            source_id,
            parent_dims,
        ),
        retained_transform = _cartesian_retained_transform(payload; source_id),
    )
end

function _cartesian_raw_product_source_retained_transform(
    sidecar::_CartesianProjectedQShellSidecarFixture3D,
)
    return Tuple(
        _cartesian_raw_product_source_retained_transform(
            payload;
            source_id = Symbol(:projected_q_shell_raw_product_source_, index),
            parent_dims = sidecar.dims,
        ) for (index, payload) in pairs(sidecar.payloads)
    )
end

function _cartesian_raw_product_source_retained_transform(
    fixed_block::_NestedFixedBlock3D,
)
    sidecar = _nested_projected_q_shell_sidecar_fixture(fixed_block)
    isnothing(sidecar) && throw(
        ArgumentError("nested fixed block does not carry a projected q-shell sidecar fixture"),
    )
    return _cartesian_raw_product_source_retained_transform(sidecar)
end

function _cartesian_raw_product_source(
    unit::_CartesianNestedProductStagedByCenterUnit3D;
    source_id::Symbol = Symbol(unit.role, :_raw_product_source),
    parent_dims::NTuple{3,Int},
)
    axis_intervals = ntuple(axis -> begin
        staged_axis = unit.axes[axis]
        staged_axis.kind == :fixed && return staged_axis.fixed_index:staged_axis.fixed_index
        staged_axis.kind == :active && return staged_axis.interval
        throw(ArgumentError("unsupported product-staged axis kind $(staged_axis.kind)"))
    end, 3)
    local_box = _cartesian_support_box_from_states(unit.support_states)
    source_dimension = _cartesian_interval_product_dimension(axis_intervals)
    source_dimension == length(unit.support_indices) || throw(
        DimensionMismatch("product-staged raw source dimension must match unit support count"),
    )
    support_summary = _contraction_rule_support_summary(
        unit.support_indices;
        parent_dimension = prod(parent_dims),
    )
    raw_source_weight_role =
        unit.kind == :product_doside ? :raw_source_positive : :debug_reference_only
    return _CartesianRawProductSource3D(
        source_id,
        unit.role,
        parent_dims,
        local_box,
        axis_intervals,
        source_dimension,
        copy(unit.support_indices),
        support_summary,
        raw_source_weight_role,
        (
            source = :nested_product_staged_by_center_unit,
            unit = unit,
        ),
        (
            source = :product_staged_raw_product_source_adapter,
            fixture_only = true,
            production_supported = false,
            staged_by_center_kind = unit.kind,
            raw_source_weight_role = raw_source_weight_role,
            operator_matrices_built = false,
        ),
    )
end

function _cartesian_retained_transform(
    unit::_CartesianNestedProductStagedByCenterUnit3D;
    source_id::Symbol = Symbol(unit.role, :_raw_product_source),
)
    transform_kind = unit.kind == :product_doside ?
                     :product_axis_transform :
                     :support_dense_coefficients
    return _CartesianRetainedTransform3D(
        source_id,
        length(unit.column_range),
        transform_kind,
        unit.coefficient_matrix,
        (
            :raw_product_source,
            unit.kind == :product_doside ?
                :separable_axis_transform_metadata :
                :support_local_dense_transform,
            :retained_columns,
        ),
        (
            method = :none,
            rank_count = length(unit.column_range),
            rank_drop_count = 0,
        ),
        :debug_reference_only,
        (
            source = :nested_product_staged_by_center_unit,
            unit = unit,
        ),
        (
            source = :product_staged_retained_transform_adapter,
            fixture_only = true,
            production_supported = false,
            retained_column_weight_role = :debug_reference_only,
            retained_weight_positive_checked = false,
            ida_weight_division_allowed = false,
            transform_matrix_scope = :support_local_to_retained,
            full_raw_to_retained_matrix_materialized = true,
            fast_product_path_requires_separable_axis_transforms = true,
            separable_axis_transforms_available = unit.kind == :product_doside,
        ),
    )
end

function _cartesian_raw_product_source_retained_transform(
    unit::_CartesianNestedProductStagedByCenterUnit3D;
    source_id::Symbol = Symbol(unit.role, :_raw_product_source),
    parent_dims::NTuple{3,Int},
)
    return (
        raw_source = _cartesian_raw_product_source(
            unit;
            source_id,
            parent_dims,
        ),
        retained_transform = _cartesian_retained_transform(unit; source_id),
    )
end

function _cartesian_raw_product_source_pair_operator_packet(
    left::_CartesianRawProductSource3D,
    right::_CartesianRawProductSource3D;
    operator_kind::Symbol,
    supported_terms::Tuple{Vararg{Symbol}},
    symmetry_status::Symbol = :symmetric_upper_triangle_placeholder,
    backend::Symbol = :metadata_only,
    provenance = (; source = :raw_product_source_pair_operator_packet_placeholder),
)
    return _CartesianRawProductSourcePairOperatorPacket3D(
        left.source_id,
        right.source_id,
        operator_kind,
        supported_terms,
        symmetry_status,
        backend,
        nothing,
        provenance,
        (
            source = :raw_product_source_pair_operator_packet_placeholder,
            placeholder_only = true,
            production_supported = false,
            operator_matrices_built = false,
            ids_terms_provenance_only = true,
            all_pairs_inventory_built = false,
            left_right_sources_embedded = false,
            left_right_retained_transforms_embedded = false,
            future_inventory_must_resolve_sources_and_transforms = true,
            left_source_dimension = left.source_dimension,
            right_source_dimension = right.source_dimension,
            raw_operator_block_ready = false,
            retained_operator_block_built = false,
        ),
    )
end

function _cartesian_raw_product_source_pair_plan(
    records;
    operator_kind::Symbol,
    supported_terms::Tuple{Vararg{Symbol}},
    symmetry_status::Symbol = :symmetric_upper_triangle_placeholder,
    backend::Symbol = :metadata_only,
    source = :raw_product_source_pair_plan,
)
    raw_sources = Dict{Symbol,_CartesianRawProductSource3D}()
    retained_transforms = Dict{Symbol,_CartesianRetainedTransform3D}()
    for record in records
        raw_source = record.raw_source
        retained_transform = record.retained_transform
        raw_source.source_id == retained_transform.source_id || throw(
            ArgumentError("raw source and retained transform source ids must match"),
        )
        haskey(raw_sources, raw_source.source_id) && throw(
            ArgumentError("duplicate raw product source id $(raw_source.source_id)"),
        )
        raw_sources[raw_source.source_id] = raw_source
        retained_transforms[retained_transform.source_id] = retained_transform
    end

    source_ids = sort!(collect(keys(raw_sources)); by = string)
    pair_keys = NTuple{2,Symbol}[]
    pair_packets = _CartesianRawProductSourcePairOperatorPacket3D[]
    for right_index in eachindex(source_ids)
        right_id = source_ids[right_index]
        for left_index in 1:right_index
            left_id = source_ids[left_index]
            push!(pair_keys, (left_id, right_id))
            push!(
                pair_packets,
                _cartesian_raw_product_source_pair_operator_packet(
                    raw_sources[left_id],
                    raw_sources[right_id];
                    operator_kind,
                    supported_terms,
                    symmetry_status,
                    backend,
                    provenance = (
                        source = source,
                        left_source_id = left_id,
                        right_source_id = right_id,
                    ),
                ),
            )
        end
    end

    diagnostics = (
        source = source,
        fixture_only = true,
        production_supported = false,
        source_count = length(source_ids),
        retained_transform_count = length(retained_transforms),
        pair_count = length(pair_packets),
        expected_upper_triangle_pair_count =
            length(source_ids) * (length(source_ids) + 1) ÷ 2,
        upper_triangle_only = true,
        all_pairs_resolve_sources = all(
            packet -> haskey(raw_sources, packet.left_source_id) &&
                      haskey(raw_sources, packet.right_source_id),
            pair_packets,
        ),
        all_pairs_resolve_retained_transforms = all(
            packet -> haskey(retained_transforms, packet.left_source_id) &&
                      haskey(retained_transforms, packet.right_source_id),
            pair_packets,
        ),
        pair_packets_placeholder_only = all(
            packet -> isnothing(packet.operator_matrices),
            pair_packets,
        ),
        raw_operator_matrices_built = false,
        retained_operator_blocks_built = false,
        metric_execution_changed = false,
        qwhamiltonian_consumes = false,
        public_default_consumes = false,
        backend_policy_changed = false,
        quadrature_policy_changed = false,
        cr2_science_status_changed = false,
    )

    return _CartesianRawProductSourcePairPlan3D(
        operator_kind,
        supported_terms,
        symmetry_status,
        source_ids,
        raw_sources,
        retained_transforms,
        Tuple(pair_packets),
        pair_keys,
        diagnostics,
    )
end

function _cartesian_raw_product_source_pair_plan(
    sidecar::_CartesianProjectedQShellSidecarFixture3D;
    operator_kind::Symbol,
    supported_terms::Tuple{Vararg{Symbol}},
    symmetry_status::Symbol = :symmetric_upper_triangle_placeholder,
    backend::Symbol = :metadata_only,
)
    return _cartesian_raw_product_source_pair_plan(
        _cartesian_raw_product_source_retained_transform(sidecar);
        operator_kind,
        supported_terms,
        symmetry_status,
        backend,
        source = :projected_q_shell_sidecar_pair_plan,
    )
end

function _cartesian_raw_product_source_pair_plan(
    fixed_block::_NestedFixedBlock3D;
    operator_kind::Symbol,
    supported_terms::Tuple{Vararg{Symbol}},
    symmetry_status::Symbol = :symmetric_upper_triangle_placeholder,
    backend::Symbol = :metadata_only,
)
    return _cartesian_raw_product_source_pair_plan(
        _cartesian_raw_product_source_retained_transform(fixed_block);
        operator_kind,
        supported_terms,
        symmetry_status,
        backend,
        source = :projected_q_shell_fixed_block_pair_plan,
    )
end

function _cartesian_explicit_quadrature_weight_role(role::Symbol)
    return role in (
        :raw_source_positive,
        :positive_required,
        :not_quadrature_weight,
        :debug_reference_only,
    )
end

function _cartesian_retained_ida_weight_division_allowed(
    transform::_CartesianRetainedTransform3D,
)
    return hasproperty(transform.diagnostics, :ida_weight_division_allowed) ?
           Bool(transform.diagnostics.ida_weight_division_allowed) : false
end

function _cartesian_pair_plan_source_index(
    plan::_CartesianRawProductSourcePairPlan3D,
    source_id::Symbol,
)
    index = findfirst(==(source_id), plan.source_ids)
    isnothing(index) && throw(
        ArgumentError("source id $source_id is not present in raw-source pair plan"),
    )
    return Int(index)
end

function _cartesian_resolve_raw_product_source_pair(
    plan::_CartesianRawProductSourcePairPlan3D,
    pair_index::Integer,
)
    packet = plan.pair_packets[Int(pair_index)]
    return _cartesian_resolve_raw_product_source_pair(plan, packet, Int(pair_index))
end

function _cartesian_resolve_raw_product_source_pair(
    plan::_CartesianRawProductSourcePairPlan3D,
    packet::_CartesianRawProductSourcePairOperatorPacket3D,
    pair_index::Integer = something(
        findfirst(==(packet), collect(plan.pair_packets)),
        0,
    ),
)
    haskey(plan.raw_sources, packet.left_source_id) || throw(
        ArgumentError("left raw source $(packet.left_source_id) is missing from pair plan"),
    )
    haskey(plan.raw_sources, packet.right_source_id) || throw(
        ArgumentError("right raw source $(packet.right_source_id) is missing from pair plan"),
    )
    haskey(plan.retained_transforms, packet.left_source_id) || throw(
        ArgumentError("left retained transform $(packet.left_source_id) is missing from pair plan"),
    )
    haskey(plan.retained_transforms, packet.right_source_id) || throw(
        ArgumentError("right retained transform $(packet.right_source_id) is missing from pair plan"),
    )

    left_raw_source = plan.raw_sources[packet.left_source_id]
    right_raw_source = plan.raw_sources[packet.right_source_id]
    left_retained_transform = plan.retained_transforms[packet.left_source_id]
    right_retained_transform = plan.retained_transforms[packet.right_source_id]
    left_index = _cartesian_pair_plan_source_index(plan, packet.left_source_id)
    right_index = _cartesian_pair_plan_source_index(plan, packet.right_source_id)
    upper_triangular = left_index <= right_index
    placeholder_only = isnothing(packet.operator_matrices)
    raw_roles_explicit =
        _cartesian_explicit_quadrature_weight_role(left_raw_source.raw_source_weight_role) &&
        _cartesian_explicit_quadrature_weight_role(right_raw_source.raw_source_weight_role)
    retained_roles_explicit =
        _cartesian_explicit_quadrature_weight_role(
            left_retained_transform.retained_column_weight_role,
        ) &&
        _cartesian_explicit_quadrature_weight_role(
            right_retained_transform.retained_column_weight_role,
        )
    left_retained_ida_allowed =
        _cartesian_retained_ida_weight_division_allowed(left_retained_transform)
    right_retained_ida_allowed =
        _cartesian_retained_ida_weight_division_allowed(right_retained_transform)

    diagnostics = (
        source = :raw_product_source_pair_plan_resolution,
        pair_index = Int(pair_index),
        pair_key = (packet.left_source_id, packet.right_source_id),
        upper_triangular = upper_triangular,
        placeholder_only = placeholder_only,
        operator_matrices_built = false,
        raw_operator_block_ready = false,
        retained_operator_block_built = false,
        plan_operator_kind_matches = packet.operator_kind == plan.operator_kind,
        plan_supported_terms_match = packet.supported_terms == plan.supported_terms,
        plan_symmetry_status_matches = packet.symmetry_status == plan.symmetry_status,
        left_source_dimension = left_raw_source.source_dimension,
        right_source_dimension = right_raw_source.source_dimension,
        left_retained_dimension = left_retained_transform.retained_dimension,
        right_retained_dimension = right_retained_transform.retained_dimension,
        left_transform_kind = left_retained_transform.transform_kind,
        right_transform_kind = right_retained_transform.transform_kind,
        left_transform_stages = left_retained_transform.transform_stages,
        right_transform_stages = right_retained_transform.transform_stages,
        left_raw_source_weight_role = left_raw_source.raw_source_weight_role,
        right_raw_source_weight_role = right_raw_source.raw_source_weight_role,
        left_retained_column_weight_role =
            left_retained_transform.retained_column_weight_role,
        right_retained_column_weight_role =
            right_retained_transform.retained_column_weight_role,
        raw_weight_roles_explicit = raw_roles_explicit,
        retained_weight_roles_explicit = retained_roles_explicit,
        left_retained_ida_weight_division_allowed = left_retained_ida_allowed,
        right_retained_ida_weight_division_allowed = right_retained_ida_allowed,
        retained_ida_weight_division_allowed =
            left_retained_ida_allowed || right_retained_ida_allowed,
        pqs_factored_transform_present = any(
            stages -> stages == (
                :raw_product_modes,
                :raw_boundary_projection,
                :full_rank_symmetric_lowdin_cleanup,
                :retained_columns,
            ),
            (left_retained_transform.transform_stages, right_retained_transform.transform_stages),
        ),
        fixture_only = true,
        production_supported = false,
    )

    return _CartesianResolvedRawProductSourcePair3D(
        (packet.left_source_id, packet.right_source_id),
        Int(pair_index),
        left_raw_source,
        right_raw_source,
        left_retained_transform,
        right_retained_transform,
        packet,
        diagnostics,
    )
end

function _cartesian_resolved_raw_product_source_pairs(
    plan::_CartesianRawProductSourcePairPlan3D,
)
    return Tuple(
        _cartesian_resolve_raw_product_source_pair(plan, index)
        for index in eachindex(plan.pair_packets)
    )
end

function _cartesian_raw_product_source_pair_plan_audit(
    plan::_CartesianRawProductSourcePairPlan3D,
)
    resolved_pairs = _cartesian_resolved_raw_product_source_pairs(plan)
    diagnostics = (
        source = :raw_product_source_pair_plan_audit,
        fixture_only = true,
        production_supported = false,
        source_count = length(plan.source_ids),
        pair_count = length(plan.pair_packets),
        expected_upper_triangle_pair_count =
            length(plan.source_ids) * (length(plan.source_ids) + 1) ÷ 2,
        every_pair_resolves_raw_sources = all(
            pair -> pair.left_raw_source.source_id == pair.pair_key[1] &&
                    pair.right_raw_source.source_id == pair.pair_key[2],
            resolved_pairs,
        ),
        every_pair_resolves_retained_transforms = all(
            pair -> pair.left_retained_transform.source_id == pair.pair_key[1] &&
                    pair.right_retained_transform.source_id == pair.pair_key[2],
            resolved_pairs,
        ),
        every_pair_upper_triangular = all(
            pair -> pair.diagnostics.upper_triangular,
            resolved_pairs,
        ),
        every_pair_placeholder_only = all(
            pair -> pair.diagnostics.placeholder_only,
            resolved_pairs,
        ),
        raw_weight_roles_explicit = all(
            pair -> pair.diagnostics.raw_weight_roles_explicit,
            resolved_pairs,
        ),
        retained_weight_roles_explicit = all(
            pair -> pair.diagnostics.retained_weight_roles_explicit,
            resolved_pairs,
        ),
        retained_ida_weight_division_allowed = any(
            pair -> pair.diagnostics.retained_ida_weight_division_allowed,
            resolved_pairs,
        ),
        raw_operator_matrices_built = false,
        retained_operator_blocks_built = false,
        metric_execution_changed = false,
        qwhamiltonian_consumes = false,
        public_default_consumes = false,
        backend_policy_changed = false,
        quadrature_policy_changed = false,
        cr2_science_status_changed = false,
    )
    return _CartesianRawProductSourcePairPlanAudit3D(resolved_pairs, diagnostics)
end

function _cartesian_require_raw_pair_supported_term(
    resolved_pair::_CartesianResolvedRawProductSourcePair3D,
    term::Symbol;
    helper::Symbol,
)
    term in resolved_pair.pair_packet.supported_terms && return nothing
    throw(
        ArgumentError(
            "private raw packet helper $(helper) cannot build term $(repr(term)); " *
            "resolved pair $(resolved_pair.pair_key) advertises supported_terms = " *
            repr(resolved_pair.pair_packet.supported_terms),
        ),
    )
end

function _cartesian_raw_low_order_operator_packet(
    resolved_pair::_CartesianResolvedRawProductSourcePair3D;
    term::Symbol,
    backend::Symbol = :private_raw_product_reference,
)
    term in (:overlap, :axis_index_x) || throw(
        ArgumentError("private raw low-order packet currently supports only :overlap and :axis_index_x"),
    )
    _cartesian_require_raw_pair_supported_term(
        resolved_pair,
        term;
        helper = :private_raw_low_order_operator_packet,
    )
    resolved_pair.pair_key[1] == resolved_pair.pair_key[2] || throw(
        ArgumentError("private raw low-order packet currently supports only self-pairs"),
    )
    left_dimension = resolved_pair.left_raw_source.source_dimension
    right_dimension = resolved_pair.right_raw_source.source_dimension
    left_dimension == right_dimension || throw(
        DimensionMismatch("raw self-pair dimensions must match"),
    )
    matrix, raw_reference, reference_error = if term == :overlap
        identity_matrix = Matrix{Float64}(I, left_dimension, right_dimension)
        (
            identity_matrix,
            :orthonormal_raw_product_mode_overlap_identity,
            norm(identity_matrix - Matrix{Float64}(I, left_dimension, right_dimension), Inf),
        )
    else
        _cartesian_identity_product_slab_axis_index_x_packet_matrix(
            resolved_pair.left_raw_source,
        )
    end
    return _CartesianRawProductSourceLowOrderOperatorPacket3D(
        resolved_pair.pair_key[1],
        resolved_pair.pair_key[2],
        :low_order_metric,
        term,
        (left_dimension, right_dimension),
        resolved_pair.pair_packet.symmetry_status,
        backend,
        matrix,
        (
            source = :private_raw_low_order_operator_packet,
            resolved_pair = resolved_pair,
        ),
        (
            source = :private_raw_low_order_operator_packet,
            fixture_only = true,
            production_supported = false,
            term = term,
            raw_reference = raw_reference,
            raw_reference_error = reference_error,
            raw_operator_matrix_built = true,
            retained_operator_block_built = false,
            retained_transform_applied = false,
            separable_axis_metadata_used = term == :axis_index_x,
            axis_index_diagnostic = term == :axis_index_x,
            physical_position_operator = false,
            dense_parent_matrix_used = false,
            all_pair_matrices_built = false,
            metric_execution_changed = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _cartesian_identity_product_slab_axis_index_x_packet_matrix(
    raw_source::_CartesianRawProductSource3D,
)
    raw_source.source_id == :identity_product_slab_source || throw(
        ArgumentError("private :axis_index_x raw packet is restricted to the identity product/slab fixture"),
    )
    hasproperty(raw_source.provenance, :unit) || throw(
        ArgumentError("identity product/slab :axis_index_x packet requires staged-unit provenance"),
    )
    unit = raw_source.provenance.unit
    unit.kind == :product_doside || throw(
        ArgumentError("identity product/slab :axis_index_x packet requires a product_doside unit"),
    )
    length(unit.axis_function_indices) == raw_source.source_dimension || throw(
        DimensionMismatch("axis-function metadata must match the raw source dimension"),
    )
    matrix = zeros(Float64, raw_source.source_dimension, raw_source.source_dimension)
    x_values = _cartesian_raw_source_axis_index_values(unit, 1)
    for (left_column, left_indices) in pairs(unit.axis_function_indices)
        for (right_column, right_indices) in pairs(unit.axis_function_indices)
            left_indices[2] == right_indices[2] || continue
            left_indices[3] == right_indices[3] || continue
            left_indices[1] == right_indices[1] || continue
            matrix[left_column, right_column] = x_values[left_indices[1]]
        end
    end
    support_reference = zeros(Float64, raw_source.source_dimension, raw_source.source_dimension)
    for (index, state) in pairs(unit.support_states)
        support_reference[index, index] = Float64(state[1])
    end
    return (
        matrix,
        :identity_product_slab_axis_index_x,
        norm(matrix - support_reference, Inf),
    )
end

function _cartesian_raw_source_axis_index_values(
    unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis::Int,
)
    staged_axis = unit.axes[axis]
    if staged_axis.kind == :fixed
        return Float64[Int(staged_axis.fixed_index)]
    elseif staged_axis.kind == :active
        return Float64[index for index in staged_axis.interval]
    end
    throw(ArgumentError("unsupported staged axis kind $(staged_axis.kind)"))
end

function _cartesian_physical_position_axis(term::Symbol)
    term == :position_x && return 1
    term == :position_y && return 2
    term == :position_z && return 3
    throw(
        ArgumentError(
            "private physical raw low-order packet currently supports only :position_x, :position_y, and :position_z",
        ),
    )
end

function _cartesian_axis_metric_source(axis_metrics, axis_name::Symbol)
    data = getproperty(axis_metrics, axis_name)
    return hasproperty(data, :source) ? Symbol(getproperty(data, :source)) :
           :explicit_axis_metrics
end

function _cartesian_axis_metric_matrix(axis_metrics, axis::Int, kind::Symbol)
    axis_name = (:x, :y, :z)[axis]
    data = getproperty(axis_metrics, axis_name)
    return Matrix{Float64}(getproperty(data, kind))
end

function _cartesian_validate_overlap_axis_metrics(axis_metrics, dims::NTuple{3,Int})
    for (axis, axis_name) in enumerate((:x, :y, :z))
        overlap = _cartesian_axis_metric_matrix(axis_metrics, axis, :overlap)
        size(overlap) == (dims[axis], dims[axis]) || throw(
            DimensionMismatch("$(axis_name)-axis overlap metric size must match raw source parent dimensions"),
        )
    end
    return nothing
end

function _cartesian_validate_position_axis_metrics(axis_metrics, dims::NTuple{3,Int})
    for (axis, axis_name) in enumerate((:x, :y, :z))
        overlap = _cartesian_axis_metric_matrix(axis_metrics, axis, :overlap)
        position = _cartesian_axis_metric_matrix(axis_metrics, axis, :position)
        size(overlap) == (dims[axis], dims[axis]) || throw(
            DimensionMismatch("$(axis_name)-axis overlap metric size must match raw source parent dimensions"),
        )
        size(position) == (dims[axis], dims[axis]) || throw(
            DimensionMismatch("$(axis_name)-axis position metric size must match raw source parent dimensions"),
        )
    end
    return nothing
end

function _cartesian_validate_factorized_low_order_axis_metrics(
    axis_metrics,
    dims::NTuple{3,Int},
    term::Symbol,
)
    _cartesian_validate_overlap_axis_metrics(axis_metrics, dims)
    term == :overlap && return nothing
    position_axis = _cartesian_physical_position_axis(term)
    axis_name = (:x, :y, :z)[position_axis]
    position = _cartesian_axis_metric_matrix(axis_metrics, position_axis, :position)
    size(position) == (dims[position_axis], dims[position_axis]) || throw(
        DimensionMismatch("$(axis_name)-axis position metric size must match raw source parent dimensions"),
    )
    return nothing
end

function _cartesian_product_doside_axis_interval(axis)
    axis.kind == :fixed && return axis.fixed_index:axis.fixed_index
    axis.kind == :active && return axis.interval
    throw(ArgumentError("unsupported product/doside staged axis kind $(axis.kind)"))
end

function _cartesian_product_doside_axis_factor_matrix(left_axis, right_axis, axis_matrix)
    parent_count = size(axis_matrix, 1)
    size(axis_matrix, 2) == parent_count || throw(
        ArgumentError("product/doside axis factor matrix must be square"),
    )
    left_interval = _cartesian_product_doside_axis_interval(left_axis)
    right_interval = _cartesian_product_doside_axis_interval(right_axis)
    first(left_interval) >= 1 && last(left_interval) <= parent_count || throw(
        ArgumentError("left product/doside staged axis interval exceeds axis metric dimensions"),
    )
    first(right_interval) >= 1 && last(right_interval) <= parent_count || throw(
        ArgumentError("right product/doside staged axis interval exceeds axis metric dimensions"),
    )
    return Matrix{Float64}(axis_matrix[left_interval, right_interval])
end

function _cartesian_product_doside_axis_local_index(axis, index::Int)
    interval = _cartesian_product_doside_axis_interval(axis)
    index in interval || throw(
        ArgumentError("raw product source support state index $(index) lies outside staged axis interval $(interval)"),
    )
    return index - first(interval) + 1
end

function _cartesian_factorized_low_order_axis_kind(term::Symbol, axis::Int)
    term == :overlap && return :overlap
    position_axis = _cartesian_physical_position_axis(term)
    return axis == position_axis ? :position : :overlap
end

function _cartesian_factorized_product_doside_raw_reference(term::Symbol)
    term == :overlap && return :product_doside_factorized_overlap
    _cartesian_physical_position_axis(term)
    return Symbol(:product_doside_factorized_, term)
end

function _cartesian_raw_source_support_states(raw_source::_CartesianRawProductSource3D)
    isnothing(raw_source.support_indices) && throw(
        ArgumentError("private physical raw packet requires explicit raw source support indices"),
    )
    length(raw_source.support_indices) == raw_source.source_dimension || throw(
        DimensionMismatch("raw source support count must match raw source dimension"),
    )
    nx, ny, nz = raw_source.parent_dims
    parent_dimension = nx * ny * nz
    states = NTuple{3,Int}[]
    sizehint!(states, length(raw_source.support_indices))
    for index in raw_source.support_indices
        1 <= index <= parent_dimension || throw(
            ArgumentError("raw source support index lies outside parent dimensions"),
        )
        push!(states, _cartesian_unflat_index(index, raw_source.parent_dims))
    end
    return states
end

function _cartesian_require_product_doside_raw_packet_source(
    raw_source::_CartesianRawProductSource3D;
    side::Symbol,
)
    hasproperty(raw_source.provenance, :unit) || throw(
        ArgumentError("product/doside raw packet requires staged-unit provenance on $(side) source"),
    )
    unit = raw_source.provenance.unit
    unit.kind == :product_doside || throw(
        ArgumentError("product/doside raw packet requires a product_doside unit on $(side) source"),
    )
    return unit
end

function _cartesian_product_doside_axis_factor_low_order_matrix(
    left_raw_source::_CartesianRawProductSource3D,
    right_raw_source::_CartesianRawProductSource3D;
    term::Symbol,
    axis_metrics,
)
    term in (:overlap, :position_x, :position_y, :position_z) || throw(
        ArgumentError("private factorized product/doside raw packet currently supports only :overlap and :position_x/y/z"),
    )
    left_unit = _cartesian_require_product_doside_raw_packet_source(
        left_raw_source;
        side = :left,
    )
    right_unit = _cartesian_require_product_doside_raw_packet_source(
        right_raw_source;
        side = :right,
    )
    left_raw_source.parent_dims == right_raw_source.parent_dims || throw(
        DimensionMismatch("factorized product/doside low-order packet requires matching parent dimensions"),
    )
    _cartesian_validate_factorized_low_order_axis_metrics(
        axis_metrics,
        left_raw_source.parent_dims,
        term,
    )
    axis_factors = ntuple(
        axis -> _cartesian_product_doside_axis_factor_matrix(
            left_unit.axes[axis],
            right_unit.axes[axis],
            _cartesian_axis_metric_matrix(
                axis_metrics,
                axis,
                _cartesian_factorized_low_order_axis_kind(term, axis),
            ),
        ),
        3,
    )
    left_support_states = _cartesian_raw_source_support_states(left_raw_source)
    right_support_states = _cartesian_raw_source_support_states(right_raw_source)
    matrix = zeros(Float64, length(left_support_states), length(right_support_states))
    for col in eachindex(right_support_states)
        right_state = right_support_states[col]
        jx = _cartesian_product_doside_axis_local_index(right_unit.axes[1], right_state[1])
        jy = _cartesian_product_doside_axis_local_index(right_unit.axes[2], right_state[2])
        jz = _cartesian_product_doside_axis_local_index(right_unit.axes[3], right_state[3])
        for row in eachindex(left_support_states)
            left_state = left_support_states[row]
            ix = _cartesian_product_doside_axis_local_index(left_unit.axes[1], left_state[1])
            iy = _cartesian_product_doside_axis_local_index(left_unit.axes[2], left_state[2])
            iz = _cartesian_product_doside_axis_local_index(left_unit.axes[3], left_state[3])
            matrix[row, col] =
                axis_factors[1][ix, jx] *
                axis_factors[2][iy, jy] *
                axis_factors[3][iz, jz]
        end
    end
    return (
        matrix,
        _cartesian_factorized_product_doside_raw_reference(term),
        0.0,
    )
end

function _cartesian_factorized_product_doside_raw_low_order_operator_packet(
    resolved_pair::_CartesianResolvedRawProductSourcePair3D;
    term::Symbol,
    axis_metrics,
    backend::Symbol = :private_raw_product_factor_reference,
)
    term in (:overlap, :position_x, :position_y, :position_z) || throw(
        ArgumentError("private factorized product/doside raw packet currently supports only :overlap and :position_x/y/z"),
    )
    _cartesian_require_raw_pair_supported_term(
        resolved_pair,
        term;
        helper = :private_factorized_product_doside_raw_low_order_operator_packet,
    )
    left_dimension = resolved_pair.left_raw_source.source_dimension
    right_dimension = resolved_pair.right_raw_source.source_dimension
    matrix, raw_reference, reference_error =
        _cartesian_product_doside_axis_factor_low_order_matrix(
            resolved_pair.left_raw_source,
            resolved_pair.right_raw_source;
            term,
            axis_metrics,
        )
    size(matrix) == (left_dimension, right_dimension) || throw(
        DimensionMismatch("factorized raw overlap matrix dimensions must match left/right raw source dimensions"),
    )
    axis_metric_sources = (
        x = _cartesian_axis_metric_source(axis_metrics, :x),
        y = _cartesian_axis_metric_source(axis_metrics, :y),
        z = _cartesian_axis_metric_source(axis_metrics, :z),
    )
    position_operator = term != :overlap
    position_axis = position_operator ? _cartesian_physical_position_axis(term) : nothing
    return _CartesianRawProductSourceLowOrderOperatorPacket3D(
        resolved_pair.pair_key[1],
        resolved_pair.pair_key[2],
        :low_order_metric,
        term,
        (left_dimension, right_dimension),
        resolved_pair.pair_packet.symmetry_status,
        backend,
        matrix,
        (
            source = :private_factorized_product_doside_raw_low_order_operator_packet,
            resolved_pair = resolved_pair,
        ),
        (
            source = :private_factorized_product_doside_raw_low_order_operator_packet,
            fixture_only = true,
            production_supported = false,
            term = term,
            raw_reference = raw_reference,
            raw_reference_error = reference_error,
            cross_pair = resolved_pair.pair_key[1] != resolved_pair.pair_key[2],
            left_source_dimension = left_dimension,
            right_source_dimension = right_dimension,
            factorized_axis_path_used = true,
            axis_factor_scope = :staged_local_axis_intervals,
            support_row_reference_used = false,
            raw_basis_scope = :raw_product_source_rows,
            raw_operator_matrix_built = true,
            retained_operator_block_built = false,
            retained_transform_applied = false,
            axis_metric_sources = axis_metric_sources,
            physical_position_operator = position_operator,
            position_axis = position_operator ? (:x, :y, :z)[position_axis] : nothing,
            dense_parent_matrix_used = false,
            all_pair_matrices_built = false,
            metric_execution_changed = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            cr2_science_status_changed = false,
            retained_positive_weight_claim = false,
            ida_weight_division_allowed = false,
        ),
    )
end

function _cartesian_product_doside_physical_position_packet_matrix(
    left_raw_source::_CartesianRawProductSource3D,
    right_raw_source::_CartesianRawProductSource3D;
    term::Symbol,
    axis_metrics,
)
    _cartesian_require_product_doside_raw_packet_source(
        left_raw_source;
        side = :left,
    )
    _cartesian_require_product_doside_raw_packet_source(
        right_raw_source;
        side = :right,
    )
    left_raw_source.parent_dims == right_raw_source.parent_dims || throw(
        DimensionMismatch("product/doside physical position raw packet requires matching parent dimensions"),
    )
    _cartesian_validate_position_axis_metrics(axis_metrics, left_raw_source.parent_dims)
    position_axis = _cartesian_physical_position_axis(term)
    axis_matrices = ntuple(
        axis -> _cartesian_axis_metric_matrix(
            axis_metrics,
            axis,
            axis == position_axis ? :position : :overlap,
        ),
        3,
    )
    left_support_states = _cartesian_raw_source_support_states(left_raw_source)
    right_support_states = _cartesian_raw_source_support_states(right_raw_source)
    matrix = zeros(Float64, length(left_support_states), length(right_support_states))
    for col in eachindex(right_support_states)
        jx, jy, jz = right_support_states[col]
        for row in eachindex(left_support_states)
            ix, iy, iz = left_support_states[row]
            matrix[row, col] =
                axis_matrices[1][ix, jx] *
                axis_matrices[2][iy, jy] *
                axis_matrices[3][iz, jz]
        end
    end
    return (
        matrix,
        Symbol(:product_doside_physical_, term),
        0.0,
        position_axis,
    )
end

function _cartesian_product_doside_physical_position_packet_matrix(
    raw_source::_CartesianRawProductSource3D;
    term::Symbol,
    axis_metrics,
)
    return _cartesian_product_doside_physical_position_packet_matrix(
        raw_source,
        raw_source;
        term,
        axis_metrics,
    )
end

function _cartesian_physical_raw_low_order_operator_packet(
    resolved_pair::_CartesianResolvedRawProductSourcePair3D;
    term::Symbol,
    axis_metrics,
    backend::Symbol = :private_raw_product_reference,
)
    position_axis = _cartesian_physical_position_axis(term)
    _cartesian_require_raw_pair_supported_term(
        resolved_pair,
        term;
        helper = :private_physical_raw_low_order_operator_packet,
    )
    left_dimension = resolved_pair.left_raw_source.source_dimension
    right_dimension = resolved_pair.right_raw_source.source_dimension
    matrix, raw_reference, reference_error, _ =
        _cartesian_product_doside_physical_position_packet_matrix(
            resolved_pair.left_raw_source,
            resolved_pair.right_raw_source;
            term,
            axis_metrics,
        )
    size(matrix) == (left_dimension, right_dimension) || throw(
        DimensionMismatch("raw physical operator matrix dimensions must match left/right raw source dimensions"),
    )
    axis_metric_sources = (
        x = _cartesian_axis_metric_source(axis_metrics, :x),
        y = _cartesian_axis_metric_source(axis_metrics, :y),
        z = _cartesian_axis_metric_source(axis_metrics, :z),
    )
    return _CartesianRawProductSourceLowOrderOperatorPacket3D(
        resolved_pair.pair_key[1],
        resolved_pair.pair_key[2],
        :low_order_metric,
        term,
        (left_dimension, right_dimension),
        resolved_pair.pair_packet.symmetry_status,
        backend,
        matrix,
        (
            source = :private_physical_raw_low_order_operator_packet,
            resolved_pair = resolved_pair,
        ),
        (
            source = :private_physical_raw_low_order_operator_packet,
            fixture_only = true,
            production_supported = false,
            term = term,
            raw_reference = raw_reference,
            raw_reference_error = reference_error,
            cross_pair = resolved_pair.pair_key[1] != resolved_pair.pair_key[2],
            left_source_dimension = left_dimension,
            right_source_dimension = right_dimension,
            raw_operator_matrix_built = true,
            retained_operator_block_built = false,
            retained_transform_applied = false,
            separable_axis_metadata_used = false,
            axis_index_diagnostic = false,
            physical_position_operator = true,
            position_axis = (:x, :y, :z)[position_axis],
            axis_metric_sources = axis_metric_sources,
            raw_basis_scope = :raw_product_source_rows,
            dense_parent_matrix_used = false,
            all_pair_matrices_built = false,
            metric_execution_changed = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _cartesian_materialized_retained_transform_matrix(
    transform::_CartesianRetainedTransform3D,
)
    isnothing(transform.transform_matrix) && throw(
        ArgumentError("private retained low-order block requires a materialized retained transform"),
    )
    return Matrix{Float64}(transform.transform_matrix)
end

function _cartesian_retained_low_order_operator_block(
    raw_packet::_CartesianRawProductSourceLowOrderOperatorPacket3D,
    left_transform::_CartesianRetainedTransform3D,
    right_transform::_CartesianRetainedTransform3D = left_transform,
)
    raw_packet.term in (:overlap, :axis_index_x, :position_x, :position_y, :position_z) || throw(
        ArgumentError(
            "private retained low-order block currently supports only :overlap, :axis_index_x, and :position_x/y/z",
        ),
    )
    raw_packet.left_source_id == left_transform.source_id || throw(
        ArgumentError("left retained transform source id does not match raw packet"),
    )
    raw_packet.right_source_id == right_transform.source_id || throw(
        ArgumentError("right retained transform source id does not match raw packet"),
    )
    left_matrix = _cartesian_materialized_retained_transform_matrix(left_transform)
    right_matrix = _cartesian_materialized_retained_transform_matrix(right_transform)
    size(raw_packet.raw_operator_matrix) == (size(left_matrix, 1), size(right_matrix, 1)) ||
        throw(
            DimensionMismatch("raw operator matrix dimensions must match retained transform source dimensions"),
        )
    retained_matrix =
        transpose(left_matrix) * raw_packet.raw_operator_matrix * right_matrix
    all(isfinite, retained_matrix) || throw(
        ArgumentError("private retained low-order block entries must be finite"),
    )
    return _CartesianRetainedLowOrderOperatorBlock3D(
        raw_packet.left_source_id,
        raw_packet.right_source_id,
        raw_packet.operator_kind,
        raw_packet.term,
        (size(left_matrix, 2), size(right_matrix, 2)),
        raw_packet.symmetry_status,
        retained_matrix,
        (
            source = :private_retained_low_order_operator_block,
            raw_packet = raw_packet,
            left_retained_transform = left_transform,
            right_retained_transform = right_transform,
        ),
        (
            source = :private_retained_low_order_operator_block,
            fixture_only = true,
            production_supported = false,
            term = raw_packet.term,
            left_transform_kind = left_transform.transform_kind,
            right_transform_kind = right_transform.transform_kind,
            left_transform_materialized = true,
            right_transform_materialized = true,
            retained_transform_applied = true,
            retained_operator_block_built = true,
            retained_reference = :explicit_materialized_transform,
            raw_operator_matrix_source = raw_packet.provenance.source,
            all_pair_matrices_built = false,
            metric_execution_changed = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            cr2_science_status_changed = false,
            ida_weight_division_allowed =
                _cartesian_retained_ida_weight_division_allowed(left_transform) ||
                _cartesian_retained_ida_weight_division_allowed(right_transform),
        ),
    )
end

function _cartesian_resolved_contraction_payloads(
    sidecar::_CartesianProjectedQShellSidecarFixture3D,
)
    parent_dimension = prod(sidecar.dims)
    fixed_block_sidecar_installed =
        hasproperty(sidecar.diagnostics, :fixed_block_sidecar_installed) ?
        Bool(sidecar.diagnostics.fixed_block_sidecar_installed) : false
    return [
        _cartesian_resolved_contraction_payload(
            payload;
            parent_dimension,
            fixed_block_sidecar_installed,
        )
        for payload in sidecar.payloads
    ]
end

function _packet_build_column_coverage(
    column_ranges::AbstractVector{<:UnitRange{Int}},
    contracted_dimension::Integer,
)
    columns = Int[]
    for range in column_ranges
        append!(columns, collect(range))
    end
    return _entry_counts(columns, 1:Int(contracted_dimension))
end

function _packet_build_source_unit_descriptor(payload::_CartesianResolvedContractionPayload3D)
    return (
        payload_kind = payload.payload_kind,
        metric_path = payload.metric_path,
        ready_for_metric_execution = payload.ready_for_metric_execution,
        column_range = payload.column_range,
        support_count = isnothing(payload.support_indices) ? nothing : length(payload.support_indices),
        block_role = hasproperty(payload.diagnostics, :block_role) ?
                     payload.diagnostics.block_role : nothing,
        metric_capability = hasproperty(payload.diagnostics, :metric_capability) ?
                            payload.diagnostics.metric_capability : nothing,
        numerical_packet_matrices_built = false,
    )
end

function _cartesian_packet_build_source_from_resolved_payloads(
    resolved_payloads;
    parent_dimension::Integer,
    contracted_dimension::Union{Nothing,Integer} = nothing,
    provenance = (; source = :resolved_contraction_payloads),
)
    payloads = collect(resolved_payloads)
    isempty(payloads) && throw(
        ArgumentError("Cartesian packet build source requires at least one resolved payload"),
    )
    all(payload -> payload.ready_for_metric_execution, payloads) || throw(
        ArgumentError("Cartesian packet build source requires executable resolved payloads"),
    )
    any(payload -> isnothing(payload.column_range), payloads) && throw(
        ArgumentError("Cartesian packet build source requires every payload to have a column range"),
    )
    column_ranges = UnitRange{Int}[payload.column_range for payload in payloads]
    resolved_contracted_dimension = isnothing(contracted_dimension) ?
        maximum(last(range) for range in column_ranges) :
        Int(contracted_dimension)
    support_indices = Int[]
    for payload in payloads
        isnothing(payload.support_indices) && throw(
            ArgumentError("Cartesian packet build source requires every payload to have support indices"),
        )
        append!(support_indices, payload.support_indices)
    end
    column_coverage = _packet_build_column_coverage(
        column_ranges,
        resolved_contracted_dimension,
    )
    support_union_summary = _contraction_rule_support_summary(
        support_indices;
        parent_dimension = Int(parent_dimension),
    )
    payload_kind_counts = _symbol_count_pairs(payload.payload_kind for payload in payloads)
    candidate_packet_fields = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :weights,
        :first_moments,
        :kinetic,
    )
    missing_packet_fields = (
        :nuclear_one_body,
        :local_coulomb_one_body,
        :local_ecp_one_body,
        :gaussian_local_terms,
        :gaussian_sum,
        :pair_sum,
        :mwg_interaction,
        :interaction,
    )
    axis_operator_requirements = (
        overlap = (:overlap,),
        position = (:position,),
        x2 = (:x2,),
        weights = (:weights,),
        centers = (:centers,),
        kinetic = (:overlap, :kinetic),
    )
    unit_descriptors = [
        _packet_build_source_unit_descriptor(payload) for payload in payloads
    ]
    return _CartesianPacketBuildSource3D(
        Int(parent_dimension),
        resolved_contracted_dimension,
        payloads,
        column_ranges,
        column_coverage,
        support_union_summary,
        payload_kind_counts,
        candidate_packet_fields,
        missing_packet_fields,
        axis_operator_requirements,
        unit_descriptors,
        (
            source = :cartesian_packet_build_source,
            metadata_only = true,
            descriptor_does_not_drive_builder = true,
            current_builder_authority = :nested_shell_packet,
            numerical_packet_matrices_built = false,
            operator_data_available = false,
            packet_operator_data_checked = false,
            overlap_matrix_built = false,
            position_matrices_built = false,
            weights_built = false,
            kinetic_matrix_built = false,
            gaussian_terms_built = false,
            interaction_terms_built = false,
            column_ranges_cover_contract =
                column_coverage.missing_count == 0 &&
                column_coverage.outside_count == 0 &&
                column_coverage.duplicate_count == 0,
            column_layout_ready_for_candidate_fields =
                column_coverage.missing_count == 0 &&
                column_coverage.outside_count == 0 &&
                column_coverage.duplicate_count == 0,
            support_union_summary_informational = true,
            parent_support_complete_required = false,
            overlapping_payload_support_allowed = true,
            packet_construction_consumes_source = false,
            source_object_builds_packet_matrices = false,
            nested_shell_packet_remains_authoritative = true,
            full_packet_builder_ready = false,
        ),
        provenance,
    )
end

function _cartesian_endcap_panel_pre_packet_build_source(
    owned_units::_CartesianNestedEndcapPanelOwnedUnits3D,
    coefficient_matrix::AbstractMatrix{<:Real},
    unit_column_ranges::AbstractVector{<:UnitRange{Int}},
    support_indices::AbstractVector{<:Integer},
    dims::NTuple{3,Int};
    provenance = (;),
)
    owned_units.coefficient_contract == :product_doside || throw(
        ArgumentError("endcap/panel pre-packet build source requires product_doside owned units"),
    )
    owned_units.audit.coverage_ok || throw(
        ArgumentError("endcap/panel pre-packet build source requires exact owned-unit support coverage"),
    )
    length(owned_units.units) == length(unit_column_ranges) || throw(
        ArgumentError("endcap/panel pre-packet build source unit column ranges must match owned units"),
    )
    parent_dimension = prod(dims)
    size(coefficient_matrix, 1) == parent_dimension || throw(
        DimensionMismatch("endcap/panel pre-packet coefficient rows must match parent dimension"),
    )
    expected_support = sort(Int[index for index in owned_units.expected_support_indices])
    actual_support = sort(Int[Int(index) for index in support_indices])
    actual_support == expected_support || throw(
        ArgumentError("endcap/panel pre-packet support indices must match owned-unit expected support"),
    )

    column_ranges = UnitRange{Int}[
        first(range):last(range) for range in unit_column_ranges
    ]
    column_coverage = _packet_build_column_coverage(
        column_ranges,
        size(coefficient_matrix, 2),
    )
    column_coverage.missing_count == 0 &&
        column_coverage.outside_count == 0 &&
        column_coverage.duplicate_count == 0 || throw(
            ArgumentError("endcap/panel pre-packet unit column ranges must exactly cover coefficient columns"),
        )

    staged_units = [
        _nested_product_staged_unit_from_owned_unit(
            owned_unit;
            column_range,
            dims,
        ) for (owned_unit, column_range) in zip(owned_units.units, column_ranges)
    ]
    resolved_payloads = [
        _cartesian_resolved_contraction_payload(unit; parent_dimension)
        for unit in staged_units
    ]
    support_counts = Int[unit.diagnostics.support_count for unit in staged_units]
    return _cartesian_packet_build_source_from_resolved_payloads(
        resolved_payloads;
        parent_dimension,
        contracted_dimension = size(coefficient_matrix, 2),
        provenance = (
            source = :endcap_panel_pre_packet_build_source,
            support_contract = owned_units.support_contract,
            coefficient_contract = owned_units.coefficient_contract,
            current_box = owned_units.current_box,
            inner_box = owned_units.inner_box,
            bond_axis = owned_units.bond_axis,
            q = owned_units.q,
            L = owned_units.L,
            parent_dimension = parent_dimension,
            contracted_dimension = size(coefficient_matrix, 2),
            unit_count = length(staged_units),
            product_unit_count = length(staged_units),
            generic_unit_count = 0,
            support_counts = support_counts,
            max_support_count = isempty(support_counts) ? 0 : maximum(support_counts),
            packet_construction_consumes_source = false,
            source_object_builds_packet_matrices = false,
            nested_shell_packet_remains_authoritative = true,
            user_provenance = provenance,
        ),
    )
end

function _cartesian_nested_sequence_pre_packet_build_source(
    dims::NTuple{3,Int},
    core_indices::AbstractVector{<:Integer},
    core_coefficients::AbstractMatrix{<:Real},
    core_column_range::UnitRange{Int},
    shell_layers::AbstractVector,
    layer_column_ranges::AbstractVector{<:UnitRange{Int}},
    coefficient_matrix::AbstractMatrix{<:Real},
    support_indices::AbstractVector{<:Integer};
    atol::Real = 1.0e-14,
    provenance = (;),
)
    length(shell_layers) == length(layer_column_ranges) || throw(
        ArgumentError("nested sequence pre-packet source layer ranges must match shell layers"),
    )
    parent_dimension = prod(dims)
    size(coefficient_matrix, 1) == parent_dimension || throw(
        DimensionMismatch("nested sequence pre-packet coefficient rows must match parent dimension"),
    )
    size(core_coefficients, 1) == parent_dimension || throw(
        DimensionMismatch("nested sequence pre-packet core coefficient rows must match parent dimension"),
    )
    length(core_column_range) == size(core_coefficients, 2) || throw(
        ArgumentError("nested sequence pre-packet core column range must match core coefficient columns"),
    )

    sequence_column_ranges = UnitRange{Int}[core_column_range]
    append!(
        sequence_column_ranges,
        UnitRange{Int}[first(range):last(range) for range in layer_column_ranges],
    )
    column_coverage = _packet_build_column_coverage(
        sequence_column_ranges,
        size(coefficient_matrix, 2),
    )
    column_coverage.missing_count == 0 &&
        column_coverage.outside_count == 0 &&
        column_coverage.duplicate_count == 0 || throw(
            ArgumentError("nested sequence pre-packet source column ranges must exactly cover coefficient columns"),
        )

    staged_units = _CartesianNestedProductStagedByCenterUnit3D[]
    product_count = 0
    generic_count = 0

    core_unit = _nested_product_staged_generic_unit(
        :core,
        core_coefficients,
        core_column_range,
        dims;
        atol,
        provenance = (; source = :sequence_core_pre_packet),
    )
    sort(core_unit.support_indices) == sort(Int[index for index in core_indices]) || throw(
        ArgumentError("nested sequence pre-packet core support must match core indices"),
    )
    push!(staged_units, core_unit)
    generic_count += 1

    for (layer, layer_range) in zip(shell_layers, layer_column_ranges)
        if hasproperty(layer, :owned_units) &&
           hasproperty(layer, :unit_column_ranges) &&
           getproperty(layer, :owned_units).coefficient_contract == :product_doside
            owned_units = getproperty(layer, :owned_units)
            unit_column_ranges = getproperty(layer, :unit_column_ranges)
            length(owned_units.units) == length(unit_column_ranges) || throw(
                ArgumentError("nested sequence pre-packet product layer unit ranges must match owned units"),
            )
            layer_first = first(layer_range)
            for (owned_unit, unit_range) in zip(owned_units.units, unit_column_ranges)
                global_range =
                    (layer_first + first(unit_range) - 1):(layer_first + last(unit_range) - 1)
                push!(
                    staged_units,
                    _nested_product_staged_unit_from_owned_unit(
                        owned_unit;
                        column_range = global_range,
                        dims,
                    ),
                )
                product_count += 1
            end
        else
            push!(
                staged_units,
                _nested_product_staged_generic_unit(
                    Symbol(:layer_, length(staged_units)),
                    coefficient_matrix[:, layer_range],
                    layer_range,
                    dims;
                    atol,
                    provenance = (
                        source = :sequence_layer_pre_packet,
                        layer_type = nameof(typeof(layer)),
                    ),
                ),
            )
            generic_count += 1
        end
    end

    product_count == 0 && throw(
        ArgumentError("nested sequence pre-packet source currently requires at least one product_doside layer"),
    )
    staged_support = sort(unique(reduce(vcat, (unit.support_indices for unit in staged_units))))
    sequence_support = sort(unique(Int[Int(index) for index in support_indices]))
    staged_support == sequence_support || throw(
        ArgumentError("nested sequence pre-packet staged support must match sequence support indices"),
    )

    resolved_payloads = [
        _cartesian_resolved_contraction_payload(unit; parent_dimension)
        for unit in staged_units
    ]
    support_counts = Int[unit.diagnostics.support_count for unit in staged_units]
    return _cartesian_packet_build_source_from_resolved_payloads(
        resolved_payloads;
        parent_dimension,
        contracted_dimension = size(coefficient_matrix, 2),
        provenance = (
            source = :nested_sequence_pre_packet_build_source,
            parent_dimension = parent_dimension,
            contracted_dimension = size(coefficient_matrix, 2),
            unit_count = length(staged_units),
            product_unit_count = product_count,
            generic_unit_count = generic_count,
            support_counts = support_counts,
            max_support_count = isempty(support_counts) ? 0 : maximum(support_counts),
            packet_construction_consumes_source = false,
            source_object_builds_packet_matrices = false,
            nested_shell_packet_remains_authoritative = true,
            fixed_block_construction_changed = false,
            qwhamiltonian_changed = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            cr2_science_status_changed = false,
            user_provenance = provenance,
        ),
    )
end

function _cartesian_packet_build_source(
    sidecar::_CartesianNestedProductStagedByCenterSidecar3D,
)
    parent_dimension = prod(sidecar.dims)
    resolved_payloads = [
        _cartesian_resolved_contraction_payload(unit; parent_dimension)
        for unit in sidecar.units
    ]
    contracted_dimension = hasproperty(sidecar.diagnostics, :final_dimension) ?
        Int(sidecar.diagnostics.final_dimension) :
        maximum(last(unit.column_range) for unit in sidecar.units)
    return _cartesian_packet_build_source_from_resolved_payloads(
        resolved_payloads;
        parent_dimension,
        contracted_dimension,
        provenance = (
            source = :nested_product_staged_by_center_sidecar,
            sidecar,
            original_provenance = sidecar.provenance,
        ),
    )
end

function _cartesian_packet_build_plan(
    sidecar::_CartesianNestedProductStagedByCenterSidecar3D,
)
    source = _cartesian_packet_build_source(sidecar)
    return _CartesianPacketBuildPlan3D(
        source,
        :nested_shell_packet,
        false,
        false,
        (
            source = :cartesian_packet_build_plan,
            metadata_only = true,
            descriptor_does_not_drive_builder = true,
            current_nested_shell_packet_authoritative = true,
            fixed_block_construction_changed = false,
            metric_packet_execution_changed = false,
            qwhamiltonian_changed = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            ida_positive_weight_semantics_changed = false,
            cr2_science_status_changed = false,
            numerical_packet_matrices_built = false,
            operator_data_available = false,
            packet_operator_data_checked = false,
            candidate_packet_fields = source.candidate_packet_fields,
            missing_packet_fields = source.missing_packet_fields,
        ),
        (
            source = :nested_product_staged_by_center_sidecar,
            packet_build_source = source,
        ),
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

function _cartesian_resolved_contraction_payload(
    unit::CartesianContractionUnit3D;
    parent_dimension::Union{Nothing,Int} = nothing,
)
    if hasproperty(unit.metadata, :staged_by_center_unit)
        payload = unit.metadata.staged_by_center_unit
        rule = hasproperty(unit.metadata, :contraction_rule) ?
               unit.metadata.contraction_rule :
               cartesian_contraction_rule(payload; parent_dimension)
        return _cartesian_resolved_contraction_payload(
            payload;
            parent_dimension,
            rule,
        )
    end
    rule = cartesian_contraction_rule(unit; parent_dimension)
    return _CartesianResolvedContractionPayload3D(
        :unsupported_missing_payload,
        false,
        rule.kind,
        rule.column_range,
        copy(rule.support_indices),
        nothing,
        nothing,
        (:staged_by_center_unit,),
        (
            source = :cartesian_contraction_unit,
            rule_family = rule.rule_family,
            rule_kind = rule.kind,
            metric_capability = rule.metric_capability,
            linear_vector_path = :unsupported,
            block_role = :unsupported,
            unsupported = true,
            prototype = _rule_prototype_only(rule),
            payload_ready_for_current_metric_execution = false,
        ),
        (
            source = :cartesian_contraction_unit,
            contraction_rule = rule,
        ),
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
