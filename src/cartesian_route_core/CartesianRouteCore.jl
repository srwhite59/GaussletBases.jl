"""
    CartesianRouteCore

Internal typed contract layer for Cartesian shellification, CPB lowering, and
retained-unit planning.

This module is intentionally not a public user API. It exists to keep the
Cartesian route code from becoming a flat collection of named-tuple fields. The
module defines small typed objects and constructors for the route concepts that
must remain distinct:

    parent geometry
        -> ShellificationRegion / OwnedSupport
        -> LoweringSource / CoordinateProductBox
        -> IntermediateRetainedSpace
        -> ShellRealization
        -> FinalRetainedUnit

Responsibilities:
- re-export Coordinate Product Box (CPB) primitives from `CartesianCPB` for
  compatibility with existing route code;
- represent shellification-owned support without confusing shells with CPBs;
- represent lowering sources for White--Lindsey and PQS routes;
- represent intermediate retained spaces and shell-realization plans;
- represent final retained units consumed by terminal realization.

Non-responsibilities:
- do not build numerical coefficient matrices;
- do not build Hamiltonian matrices;
- do not write artifacts;
- do not perform dense parent-space operator construction;
- do not encode route-driver report plumbing.

Important invariants:
- `B_outer \\ B_inner` is shell/owned support, not a CPB.
- White--Lindsey shell lowering may use facet/edge/corner CPBs.
- PQS shell lowering uses a filled source CPB plus boundary-mode selection and
  later shell realization.

Metadata policy:
- `metadata` fields are for diagnostics, provenance, labels, and transitional
  adapter notes.
- Required route semantics should be typed fields or constructor arguments, not
  hidden in metadata.
"""
module CartesianRouteCore

using ..CartesianCPB

import ..CartesianCPB: role, support_count

const CPB = CartesianCPB
const CoordinateProductBox = CartesianCPB.CoordinateProductBox
const cpb = CartesianCPB.cpb
const filled_cpb = CartesianCPB.filled_cpb
const slab_cpb = CartesianCPB.slab_cpb
const complete_shell_boundary_strata = CartesianCPB.complete_shell_boundary_strata
const intervals = CartesianCPB.intervals
const shape = CartesianCPB.shape
const codimension = CartesianCPB.codimension

# Export policy:
# These exports are internal-to-GaussletBases conveniences, not public package
# API. External package users should not rely on this module. Route-driver code
# may either import these names or call them as `CartesianRouteCore.name`.
#
export CoordinateProductBox,
       OwnedSupport,
       ShellificationRegion,
       LoweringSource,
       IntermediateRetainedSpace,
       ShellRealization,
       FinalRetainedUnit,
       cpb,
       filled_cpb,
       slab_cpb,
       complete_shell_support,
       complete_shell_boundary_strata,
       shellification_region,
       lowering_source,
       white_lindsey_boundary_strata_lowering,
       pqs_filled_source_lowering,
       intermediate_retained_space,
       boundary_product_mode_count,
       shell_realization,
       trivial_shell_realization,
       pqs_shell_realization,
       final_retained_unit,
       intervals,
       shape,
       codimension,
       role,
       support_count,
       owned_support,
       source_cpbs,
       lowering_recipe,
       unit_keys,
       final_units

# File organization:
#
# shellification_regions.jl
#     OwnedSupport and ShellificationRegion. Shells live here as owned support,
#     not as CPBs.
#
# lowering_sources.jl
#     Recipe-specific source CPBs for LW, PQS, direct, or future lowerings.
#
# retained_spaces.jl
#     IntermediateRetainedSpace, ShellRealization, FinalRetainedUnit.
include("shellification_regions.jl")
include("lowering_sources.jl")
include("retained_spaces.jl")

# Example: White--Lindsey complete-shell lowering
#
# outer = filled_cpb(1:5, 1:5, 1:5; role = :outer_box)
# inner = filled_cpb(2:4, 2:4, 2:4; role = :inner_box)
# support = complete_shell_support(outer, inner)
# region = shellification_region(:atom_local_shell, support)
# strata = complete_shell_boundary_strata(outer, inner)
# lowering = white_lindsey_boundary_strata_lowering(region, strata.all_strata)
#
# Example: PQS complete-shell lowering
#
# source = filled_cpb(1:5, 1:5, 1:5; role = :pqs_source_box)
# lowering = pqs_filled_source_lowering(region, source)
# space = intermediate_retained_space(
#     lowering;
#     retained_rule = :pqs_boundary_comx_product_modes,
#     source_mode_dims = (5, 5, 5),
# )
# realization = pqs_shell_realization(space, region)

end # module CartesianRouteCore
