# Integration/slow nested runner. Do not include in default nested runner.
# Run manually for slow/materializer/audit/probe coverage.

using Test
using LinearAlgebra
using SparseArrays
using GaussletBases
using JLD2

include("pqs_source_metadata_real_artifact_acceptance_runtests.jl")
include("pqs_explicit_core_spacing_parent_axis_probe_runtests.jl")
include("pqs_raw_product_box_plan_probe_runtests.jl")
include("pqs_source_box_route_driver_report_runtests.jl")
include("cartesian_terminal_shellification_geometry_runtests.jl")
include("cartesian_shellification_plan_runtests.jl")
include("cartesian_ham_builder_one_center_config_smoke_runtests.jl")
include("cartesian_ham_builder_diatomic_config_smoke_runtests.jl")
include("cartesian_route_diatomic_materializer_probe_runtests.jl")
include("white_lindsey_materialized_seed_runtests.jl")

include("pqs_projected_q_shell_local_layer_integration_runtests.jl")

include("bond_aligned_diatomic_high_order_recipe_realization_audit_runtests.jl")
include("bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl")
include("bond_aligned_diatomic_endcap_panel_shared_shell_source_policy_runtests.jl")

include("one_center_atomic_fixed_block_timing_surface_runtests.jl")
