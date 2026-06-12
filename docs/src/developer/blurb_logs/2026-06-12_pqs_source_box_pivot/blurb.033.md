Purpose:

Audit the mechanical module extraction boundary for
`CartesianFinalBasisRealization` before moving code.

Why now:

Pass 032 corrected the PQS status and ownership wording. The near-term plan
says the next risk is CPBM accumulating later-stage PQS final-basis realization
and operator-transfer concepts. Before code moves, identify the exact functions,
callers, exports, and tests affected so the extraction can be mechanical and
small.

Exact task:

Read and inspect:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- `src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`
- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- tests that call the candidate final-basis helpers
- docs that name the candidate helpers

Produce an extraction audit. Do not move code yet.

Classify each candidate as:

```text
move_now
leave_in_CPBM_for_now
oracle_helper_move_with_module
needs_manager_decision
```

Candidate functions:

```text
pqs_source_shell_realization_final_basis
pqs_source_shell_final_one_body_from_boundary_matrix
pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block
pqs_source_shell_projected_one_body_matrix
pqs_source_shell_final_one_electron_hamiltonian
```

Also report:

- proposed new module file layout;
- CPBM import/reexport compatibility choice;
- tests that would need alias updates;
- whether any CPBM contract-test sections should move or shrink;
- validation target for the eventual mechanical extraction.

Do not:

- move code;
- add the module yet;
- add tests;
- change exports;
- change source behavior;
- redesign result types;
- add IDA, density-density, RHF, drivers, exports, or artifacts;
- request UI escalation. In unattended baton mode, write
  `.agent_handoffs/ATTENTION.md` and stop if blocked.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- No Julia test required unless you edit executable code, which you should not.

Deletion/shrinkage report required:

- whether this audit identifies CPBM ownership that can be removed next pass;
- whether any test sections can move/shrink instead of being duplicated;
- whether any old oracle helper should stay oracle-only in the new module;
- what remains for the direct retained-boundary kernel after extraction.

Report back:

- write `.agent_handoffs/response.033.md.tmp`, then atomically rename to
  `.agent_handoffs/response.033.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.033.md`;
- include validation;
- include deletion/shrinkage report;
- after writing the response, continue polling for `blurb.034.md`,
  `ATTENTION.md`, or `STOP.md`;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
