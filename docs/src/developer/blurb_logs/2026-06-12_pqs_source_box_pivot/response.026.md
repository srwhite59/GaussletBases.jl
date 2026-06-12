Status: complete.

Implemented the narrow CPBM Hamiltonian-stage assembly helper for PQS final one-electron data:

- `pqs_source_shell_final_one_electron_hamiltonian(final_kinetic_result, final_nuclear_by_center_results)`

Files changed:

- `src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

Implementation details:

- The helper consumes:
  - materialized final kinetic output from `pqs_source_shell_final_one_body_from_boundary_matrix(...)`;
  - one or more separated final by-center nuclear outputs from `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(...)`.
- It validates:
  - final kinetic term/status/materialization;
  - nuclear by-center term/status/materialization;
  - matching final matrix dimensions;
  - finite and symmetric kinetic/nuclear inputs;
  - every nuclear center has charge recorded;
  - nuclear charge has not already been applied;
  - centers have not already been summed;
  - uncharged by-center convention is present.
- It applies recorded nuclear charges and sums centers only at this Hamiltonian assembly stage:

```text
H = T_final + sum_center Z_center * V_final_unit_charge(center)
```

- It returns compact matrix diagnostics, center summaries, applied charge list, and explicit downstream nonclaims.

Nonclaims:

- no eigensolve;
- no generalized-overlap solve;
- no IDA;
- no density-density;
- no RHF;
- no driver route;
- no export/artifact work.

Test coverage:

- Extended the compact CPBM module-contract section with a synthetic two-center Hamiltonian assembly check.
- The test verifies:
  - charges are applied only in the Hamiltonian helper;
  - separated centers are summed in the Hamiltonian helper;
  - `H = T + sum(Z_i * V_i)`;
  - no H1/eigensolve/downstream physics claims are made.
- No real projected-q-shell integration test or H1 probe was added.

Exact next blocker:

```text
:missing_pqs_final_one_electron_solve
```

The next pass should run a real final-basis H1/eigensolve probe using the existing final overlap/kinetic/nuclear/Hamiltonian path and compare it against the old oracle only as needed. This pass does not claim H1 acceptance.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No old production code, tests, metadata, or compatibility paths became fully obsolete yet because this pass adds the first Hamiltonian-stage consumer of the final kinetic/nuclear pieces.
- Nothing was deleted.
- I removed the initially-added export line for this helper before publishing because pass 026 explicitly excludes export work; the helper remains CPBM-qualified module code for the current internal route.
- The new contract assertions are genuinely new live-contract coverage for Hamiltonian-stage charge application and center summation. They do not replace older coverage yet.
- Remaining stale/duplicate surfaces to retire next: retained-source/PQS H1 scratch probes and current-route safe-term oracle vocabulary should be reviewed after the real final-basis H1 probe exists.

-- repo-doer@macmini
