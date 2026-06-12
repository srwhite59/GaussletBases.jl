"""
    CartesianPairBlockMaterialization

Internal metadata/preflight layer for Cartesian pair-block materialization.

This module owns the transition:

    pair-operator plan metadata -> numerical pair blocks

It currently provides metadata/preflight records for pair-block
materialization, plus local direct/direct one-body pair-block pilots for:

    :overlap
    :position_x, :position_y, :position_z
    :x2_x, :x2_y, :x2_z
    :kinetic

The direct/direct one-body selector is local pair-block materialization only.
Signs and prefactors for supplied 1D factors are owned by the caller-provided
parent/axis blocks.

PQS/PQS raw-source pair support currently has metadata preflight plus raw
source-space safe one-body helpers for overlap, position, x2, kinetic, and a
selector over those terms. It can also contract a raw source-space result to
retained PQS source modes by source-mode boundary selector columns. Caller-
supplied 1D factors own signs and prefactors. It also provides metadata-only
shell-realization bridge summaries for those source-space blocks and batches,
plus metadata-only readiness summaries for future final PQS pair blocks. This
path does not build shell projection, Lowdin realization, or final shell-
realized PQS pair blocks.

The private route-shaped global one-body adapter currently composes existing
local block collections, placement plans, and global safe one-body matrix pilots
for the term-separated safe one-body terms overlap, kinetic, position_x/y/z,
x2_x/y/z, plus uncharged by-center electron-nuclear matrices over decomposed
White--Lindsey unit-pair inventories. The module also exposes the narrow
shellification -> lowering -> retained-unit -> unit-pair source constructor
needed by those inventories. It still does not assemble Hamiltonians, Ham
bundles, exports, artifacts, final retained PQS pair blocks, or complete
White--Lindsey route drivers.
"""
module CartesianPairBlockMaterialization

using ..CartesianCPB
using ..CartesianShellification
using ..CartesianTerminalLowering
using ..CartesianUnitPairs
using ..CartesianRetainedUnits
using ..CartesianRouteCore
using ..CartesianPairOperatorPlans
using ..CartesianRetainedUnitTransformContracts
using ..CartesianRawProductSources
using ..TimeG: @timeg
using LinearAlgebra
using SparseArrays

const CPB = CartesianCPB
const CSH = CartesianShellification
const CTL = CartesianTerminalLowering
const CUP = CartesianUnitPairs
const CRU = CartesianRetainedUnits
const CRC = CartesianRouteCore
const CPOP = CartesianPairOperatorPlans
const CRTC = CartesianRetainedUnitTransformContracts
const CRPS = CartesianRawProductSources
const ParentGaussletBases = Base.parentmodule(@__MODULE__)

export PairBlockMaterializationPolicy,
       MetadataOnlyPairBlockMaterialization,
       PairBlockMaterializationRecord,
       PairBlockMaterializationResult,
       PairBlockMaterializationBatchResult,
       PairBlockMaterializationPlan,
       pair_block_materialization_plan,
       pair_block_materialization_records,
       pqs_source_pair_overlap_block,
       pqs_source_pair_position_block,
       pqs_source_pair_x2_block,
       pqs_source_pair_kinetic_block,
       pqs_source_pair_overlap_blocks,
       pqs_source_pair_position_blocks,
       pqs_source_pair_x2_blocks,
       pqs_source_pair_kinetic_blocks,
       pqs_source_pair_retained_one_body_block,
       pqs_source_pair_retained_one_body_blocks,
       pqs_source_pair_retained_overlap_block,
       pqs_source_pair_retained_kinetic_block,
       pqs_retained_source_one_body_matrix,
       pqs_source_pair_one_body_block,
       pqs_source_pair_one_body_blocks,
       pqs_source_pair_shell_realization_bridge_summary,
       pqs_source_pair_final_block_readiness_summary,
       white_lindsey_boundary_stratum_adapter_summary,
       white_lindsey_boundary_stratum_unit_adapter_descriptor,
       white_lindsey_boundary_stratum_pair_adapter_descriptor,
       white_lindsey_materialized_seed_oracle_summary,
       white_lindsey_boundary_stratum_unit_coefficients,
       white_lindsey_boundary_stratum_unit_coefficient_context,
       white_lindsey_boundary_stratum_pair_unit_coefficients,
       white_lindsey_shellification_decomposed_unit_pair_inventory,
       white_lindsey_decomposed_unit_pair_inventory,
       white_lindsey_boundary_stratum_overlap_block,
       white_lindsey_boundary_stratum_position_block,
       white_lindsey_boundary_stratum_x2_block,
       white_lindsey_boundary_stratum_kinetic_block,
       white_lindsey_boundary_stratum_electron_nuclear_by_center_block,
       route_global_decomposed_wl_density_density_matrix,
       white_lindsey_boundary_stratum_one_body_block,
       white_lindsey_boundary_stratum_one_body_blocks,
       white_lindsey_boundary_stratum_one_body_adapter_summary,
       one_body_electron_nuclear_by_center_placement_plan,
       one_body_global_electron_nuclear_by_center_matrix,
       route_global_combined_gto_basis_layout,
       route_global_mixed_gto_blocks_from_decomposed_units,
       route_global_combined_gto_one_electron_matrices,
       route_global_combined_gto_residual_moment_matrices,
       route_global_combined_gto_final_basis_projection,
       route_global_residual_gto_mwg_representation,
       route_global_combined_gto_final_basis_density_density_matrix,
       direct_direct_overlap_block,
       direct_direct_overlap_blocks,
       direct_direct_position_block,
       direct_direct_position_blocks,
       direct_direct_x2_block,
       direct_direct_x2_blocks,
       direct_direct_kinetic_block,
       direct_direct_kinetic_blocks,
       direct_direct_one_body_block,
       direct_direct_one_body_blocks,
       unavailable_summary,
       summary

# File organization:
#
# records.jl
#     Pair-block materialization policy, record, and plan types.
#
# preflight.jl
#     Metadata-only pair-operator to pair-block readiness rules.
#
# summaries.jl
#     Compact summaries for tests and reports.
#
# one_body_terms.jl
#     Shared internal one-body term descriptors and selector-surface audit for
#     the future mixed one-body consumer.
#
# one_body_factor_inputs.jl
#     Caller-supplied one-body factor input normalization for future mixed
#     dispatch.
#
# one_body_dispatch.jl
#     Mixed one-body dispatch preflight and local selector orchestration.
#
# one_body_block_collection.jl
#     Private local one-body block collection entry vocabulary.
#
# route_one_body_adapter.jl
#     Private route-shaped adapter from pair-block plans to local one-body
#     block collections.
#
# route_global_one_body_adapter.jl
#     Private route-shaped adapter from pair-block plans to one global retained
#     one-body term matrix for the safe one-body terms only.
#
# route_global_combined_gto_basis_layout.jl
#     Metadata-only layout for appending a small GTO supplement sector to the
#     decomposed WL retained gausslet sector.
#
# route_global_mixed_gto_blocks.jl
#     Narrow route-global mixed gausslet/GTO block placement over decomposed
#     WL retained units.
#
# route_global_combined_gto_matrix_assembly.jl
#     Narrow combined gausslet+GTO overlap and one-electron Hamiltonian matrix
#     assembly from already materialized gausslet route-global and provider
#     GTO blocks.
#
# route_global_combined_gto_final_basis.jl
#     Final-basis orthogonalization/projection of the GTO supplement residual
#     against the decomposed WL gausslet sector.
#
# route_global_combined_gto_density_density.jl
#     Residual GTO -> MWG effective Gaussian representation and first
#     final-basis combined gausslet+residual-supplement density-density
#     readiness surface. It blocks rather than treating raw GTO density-density
#     as final MWG data.
#
# route_global_atom_gto_final_basis_route.jl
#     Private one-center atom+GTO construction seam that wires mapped parent
#     axes, shellification-backed WL inventory, combined-GTO one-electron
#     assembly, final-basis projection, and residual-MWG density-density
#     materialization for driver-facing probes.
#
# one_body_placement_plan.jl
#     Metadata-only local one-body placement records for future global retained
#     operator assembly.
#
# one_body_global_matrix_helpers.jl
#     Shared private validation and symmetric insertion helpers for dense
#     global one-body matrix pilots.
#
# one_body_global_overlap.jl
#     First dense global retained-overlap assembly pilot from placeable local
#     placement records only.
#
# one_body_global_kinetic.jl
#     Kinetic-only dense global retained-matrix assembly pilot from placeable
#     local placement records only.
#
# one_body_global_position.jl
#     Position-axis dense global retained-matrix assembly pilots from placeable
#     local placement records only.
#
# one_body_global_x2.jl
#     X2-axis dense global retained-matrix assembly pilots from placeable local
#     placement records only.
#
# one_body_global_electron_nuclear.jl
#     By-center electron-nuclear dense global retained-matrix pilot from
#     placeable local placement records only.
#
# pqs_source_safe_terms.jl
#     PQS/PQS raw source-space safe-term helpers.
#
# pqs_source_shell_bridge.jl
#     Metadata-only bridge summaries for future PQS shell realization.
#
# pqs_source_final_readiness.jl
#     Metadata-only readiness summaries for future final PQS pair blocks.
#
# white_lindsey_adapter_summary.jl
#     Metadata-only old-kernel reuse guidance and unit descriptors for future
#     LW adapters.
#
# white_lindsey_seed_oracle_summary.jl
#     Compact old-seed oracle summary for LW adapter validation.
#
# white_lindsey_unit_coefficients.jl
#     Narrow LW unit coefficient adapters.
#
# white_lindsey_pair_unit_coefficients.jl
#     Pair-level gathering of LW unit coefficient maps, without operator
#     blocks.
#
# white_lindsey_decomposed_unit_pair_inventory.jl
#     Route-owned metadata validator for decomposed LW unit pairs with retained
#     column ranges.
#
# white_lindsey_overlap.jl
#     First overlap-only LW pair-block pilot.
#
# white_lindsey_position.jl
#     Position-only LW pair-block pilot.
#
# white_lindsey_x2.jl
#     X2-only LW pair-block pilot.
#
# white_lindsey_kinetic.jl
#     Kinetic-only LW pair-block pilot.
#
# white_lindsey_electron_nuclear.jl
#     By-center electron-nuclear LW pair-block pilot.
#
# white_lindsey_density_density.jl
#     Route-global decomposed WL electron-electron density-density interaction
#     matrix assembly from pair-factor local blocks.
#
# white_lindsey_one_body.jl
#     Selector and compact summaries over LW one-body safe-term pilots.
#
# direct_overlap.jl
#     First tiny direct/direct overlap pair-block pilot.
#
# direct_position.jl
#     Direct/direct position pair-block pilot.
#
# direct_x2.jl
#     Direct/direct x2 pair-block pilot.
#
# direct_kinetic.jl
#     Direct/direct kinetic pair-block pilot.
#
# direct_one_body.jl
#     Direct/direct one-body term selector.
#
# pqs_source_one_body.jl
#     PQS/PQS raw/retained source-mode one-body term selectors and one-unit
#     retained source-mode matrix helper.
include("records.jl")
include("preflight.jl")
include("summaries.jl")
include("one_body_terms.jl")
include("one_body_factor_inputs.jl")
include("one_body_dispatch.jl")
include("one_body_block_collection.jl")
include("route_one_body_adapter.jl")
include("route_global_one_body_adapter.jl")
include("route_global_combined_gto_basis_layout.jl")
include("route_global_mixed_gto_blocks.jl")
include("route_global_combined_gto_matrix_assembly.jl")
include("route_global_combined_gto_final_basis.jl")
include("route_global_combined_gto_density_density.jl")
include("route_global_atom_gto_final_basis_route.jl")
include("one_body_placement_plan.jl")
include("one_body_global_matrix_helpers.jl")
include("one_body_global_overlap.jl")
include("one_body_global_kinetic.jl")
include("one_body_global_position.jl")
include("one_body_global_x2.jl")
include("one_body_global_electron_nuclear.jl")
include("direct_overlap.jl")
include("pqs_source_safe_terms.jl")
include("pqs_source_shell_bridge.jl")
include("pqs_source_final_readiness.jl")
include("white_lindsey_adapter_summary.jl")
include("white_lindsey_seed_oracle_summary.jl")
include("white_lindsey_unit_coefficients.jl")
include("white_lindsey_pair_unit_coefficients.jl")
include("white_lindsey_decomposed_unit_pair_inventory.jl")
include("white_lindsey_overlap.jl")
include("white_lindsey_position.jl")
include("white_lindsey_x2.jl")
include("white_lindsey_kinetic.jl")
include("white_lindsey_electron_nuclear.jl")
include("white_lindsey_density_density.jl")
include("white_lindsey_one_body.jl")
include("direct_position.jl")
include("direct_x2.jl")
include("direct_kinetic.jl")
include("direct_one_body.jl")
include("pqs_source_one_body.jl")

end # module CartesianPairBlockMaterialization
