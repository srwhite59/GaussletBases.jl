module CartesianContractedParentMetrics
import LinearAlgebra
import SparseArrays
import ..GaussletBases: CoulombGaussianExpansion,
                         _NestedFixedBlock3D,
                         _BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
                         _CartesianNestedAxisBundles3D,
                         _CartesianNestedProjectedQShellStagedUnitDescriptor3D,
                         _CartesianNestedProductStagedByCenterUnit3D,
                         _MappedOrdinaryGausslet1DBundle,
                         _MappedOrdinaryPGDGIntermediate1D,
                         _cartesian_flat_index,
                         _cartesian_raw_product_box_plan,
                         _cartesian_raw_product_box_source_mode_indices,
                         _cartesian_unflat_index,
                         _nested_axis_bundle,
                         _nested_axis_lengths,
                         _nested_product_axis_function_indices,
                         _nested_product_staged_active_axis,
                         _nested_product_staged_fixed_axis,
                         _nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block,
                         _nested_projected_q_shell_staged_unit_descriptor,
                         _nested_projected_q_shell_descriptor_seed_coefficients,
                         _qwrg_bond_aligned_axis_bundles,
                         _require_analytic_primitive_backend,
                         bond_aligned_homonuclear_qw_basis,
                         centers,
                         contract_primitive_matrix,
                         coulomb_gaussian_expansion,
                         gaussian_factor_matrices,
                         integral_weights,
                         overlap_matrix,
                         position_matrix,
                         primitive_set
import ..GaussletBases.CartesianContractedParents:
    CartesianContractedParent3D,
    _CartesianExecutableProjectedQShellPayload3D,
    _CartesianProjectedQShellSidecarFixture3D,
    _cartesian_resolved_contraction_payload,
    _cartesian_resolved_contraction_payloads,
    cartesian_contracted_parent,
    contracted_parent_contraction_rules,
    contracted_parent_basis,
    contracted_parent_coefficients,
    contracted_parent_dimension,
    contracted_parent_parent_dimension,
    contracted_parent_units
import ..GaussletBases.CartesianParentGaussletBases:
    CartesianParentGaussletBasis3D,
    axis_basis,
    parent_axis_counts,
    parent_dimension

export CartesianContractedParentMetricPacket3D,
       cartesian_contracted_parent_metric_packet,
       cartesian_contracted_parent_metric_packet_dense_reference,
       contracted_parent_metric_packet_parent,
       contracted_parent_metric_packet_overlap,
       contracted_parent_metric_packet_weights,
       contracted_parent_metric_packet_centers,
       contracted_parent_metric_packet_diagnostics
include("cartesian_contracted_parent_metrics/core.jl")
include("cartesian_contracted_parent_metrics/product_staged_metric_fallbacks.jl")
include("cartesian_contracted_parent_metrics/source_box_pair_shadow.jl")
include("cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl")
include("cartesian_contracted_parent_metrics/current_route_metadata_export.jl")
include("cartesian_contracted_parent_metrics/component_smoke_sidecars.jl")
end
