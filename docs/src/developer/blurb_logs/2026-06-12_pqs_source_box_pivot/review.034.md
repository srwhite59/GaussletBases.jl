Review 034: accepted.

The extraction matches the pass boundary. `CartesianFinalBasisRealization` now
owns the shell-realization final-basis object, shell-support oracle projection,
and retained-boundary overlap/kinetic transfer. CPBM keeps direct aliases so
existing CPBM-qualified callers still work, and the CPBM file now contains only
the helpers that still depend on CPBM result types.

Independent manager validation:

```text
julia --project=. -e 'using GaussletBases; ... alias assertions ...'
julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl
git diff --check
```

The CPBM contract test passed in about 13 seconds.

Remaining cleanup before adding new PQS kernels:

- The final-basis behavior checks still live inside the large CPBM contract
  file even though the implementation moved.
- Move/shrink those checks into compact `CartesianFinalBasisRealization`
  module-contract coverage.
- Leave CPBM with only small alias compatibility and CPBM-owned
  nuclear-by-center/final-Hamiltonian checks.

Do not proceed to direct retained-boundary kernels until this test ownership
cleanup is done.

-- repo-manager@macmini
