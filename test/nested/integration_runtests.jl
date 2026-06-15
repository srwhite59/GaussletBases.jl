# Integration/slow nested runner. Do not include in default nested runner.
# Run through the top-level harness with:
#   GAUSSLETBASES_TEST_GROUPS=nested GAUSSLETBASES_SLOW_TESTS=1 julia --project=. test/runtests.jl
#
# Some fixture-bound tests below depend on helpers defined by test/runtests.jl.
# Do not use this runner, docs builds, route-driver tests, or oracle/gate tests
# as routine per-pass validation for small mixed one-body consumer edits.

using Test
using LinearAlgebra
using SparseArrays
using GaussletBases
using JLD2

isdefined(Main, :_cached_fixture) ||
    error("nested integration tests require top-level test/runtests.jl fixture helpers")

# Fixture-bound bundle/projector and nested fixed-block tests.
include("cartesian_nested_face_first_primitive_runtests.jl")
include("cartesian_endcap_panel_owned_shell_producer_runtests.jl")
include("cartesian_nested_shell_first_packet_runtests.jl")
include("cartesian_nested_support_immediate_contraction_helpers_runtests.jl")
include("cartesian_nested_shell_interface_runtests.jl")
include("cartesian_contracted_parent_metric_packet_runtests.jl")

include("one_center_atomic_full_parent_nested_contract_runtests.jl")
include("one_center_atomic_legacy_profile_nested_contract_runtests.jl")

include("cartesian_basis_representation_direct_product_qw_bases_runtests.jl")
include("cartesian_basis_representation_nested_fixed_blocks_runtests.jl")

include("nested_coefficient_maps_sparse_storage_runtests.jl")

include("cartesian_basis_representation_atomic_qw_residual_bases_runtests.jl")
include("cartesian_basis_representation_cross_overlap_runtests.jl")
include("cartesian_basis_projector_orbital_transfer_runtests.jl")
include("cartesian_basis_bundle_export_runtests.jl")
include("cartesian_basis_bundle_overlap_projector_runtests.jl")

include("atomic_direct_product_he_extent_change_runtests.jl")
include("atomic_hybrid_he_orbital_transfer_runtests.jl")
include("mapped_ordinary_cartesian_1d_working_representation_runtests.jl")
include("one_center_atomic_factorized_direct_packet_kernel_runtests.jl")
include("one_center_atomic_legacy_profile_residual_completion_runtests.jl")
include("one_center_atomic_ns9_legacy_profile_residual_stabilization_runtests.jl")

include("cartesian_nested_shell_sequence_fixed_block_runtests.jl")
include("cartesian_nested_complete_shell_layer_runtests.jl")
