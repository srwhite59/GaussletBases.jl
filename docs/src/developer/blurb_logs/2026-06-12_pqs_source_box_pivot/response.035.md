Pass 035 completed.

Files created:

- `test/nested/cartesian_final_basis_realization_contract_runtests.jl`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.035.md`

Files edited:

- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `test/nested/runtests.jl`

Checks moved from CPBM to CFBR:

- `pqs_source_shell_realization_final_basis`
  - `q=5/L=5` source dims/count `5 x 5 x 5 / 125`;
  - boundary retained count `98`;
  - final retained count `98`;
  - identity final overlap for the simple identity fixture;
  - critical no-IDA/RHF/driver/artifact boundary claims.
- `pqs_source_shell_projected_one_body_matrix`
  - shell-support oracle projection for one symmetric shell operator;
  - boundary/final operator equality in the identity fixture;
  - Lowdin boundary crosscheck.
- `pqs_source_shell_final_one_body_from_boundary_matrix`
  - retained-boundary overlap transfer;
  - retained-boundary kinetic transfer.

What remains in CPBM and why:

- A tiny alias smoke checks that CPBM-qualified names point to the
  `CartesianFinalBasisRealization` functions.
- CPBM still checks `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block`
  because it consumes CPBM-owned `PairBlockMaterializationResult`.
- CPBM still checks `pqs_source_shell_final_one_electron_hamiltonian` because
  it consumes the CPBM-owned by-center nuclear final-transfer result.

Validation run:

- `julia --project=. test/nested/cartesian_final_basis_realization_contract_runtests.jl`
  passed: 31 tests in 0.4s.
- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- The large CPBM contract file no longer carries detailed behavior checks for
  the three module-owned final-basis realization/operator-transfer functions.
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  changed by 11 insertions and 92 deletions for this cleanup section.
- The new compact CFBR test replaces the moved coverage instead of duplicating
  it.
- Remaining CPBM coverage is limited to compatibility alias smoke and the
  CPBM-owned nuclear/Hamiltonian helpers.
- Remaining test bloat/oracle surfaces to retire next:
  - the CPBM alias smoke can be removed once internal callers are updated to
    use `CartesianFinalBasisRealization` directly;
  - the CPBM nuclear final-transfer helper can move later after the
    CPBM-owned result type boundary is resolved.

-- repo-doer@macmini
