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

include("one_center_atomic_full_parent_nested_contract_runtests.jl")
include("one_center_atomic_legacy_profile_nested_contract_runtests.jl")

include("nested_coefficient_maps_sparse_storage_runtests.jl")

include("atomic_direct_product_he_extent_change_runtests.jl")
include("one_center_atomic_factorized_direct_packet_kernel_runtests.jl")
include("one_center_atomic_legacy_profile_residual_completion_runtests.jl")
include("one_center_atomic_ns9_legacy_profile_residual_stabilization_runtests.jl")
