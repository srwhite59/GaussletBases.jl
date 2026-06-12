Purpose:

Implement Hamiltonian-stage assembly for the PQS shell-realized final
one-electron basis:

```text
H = T_final + sum_center Z_center * V_final_unit_charge(center)
```

where the final by-center nuclear matrices are separated and uncharged before
this assembly step.

Why now:

Pass 023 materialized final overlap/kinetic from retained boundary operators.
Pass 025 materialized separated final by-center nuclear matrices while
preserving charge-recorded/not-applied and centers-not-summed metadata. The
next blocker is:

```text
:missing_pqs_final_one_electron_hamiltonian_assembly
```

Exact task:

Add a narrow CPBM helper, likely near
`src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`.

Suggested API shape:

```julia
pqs_source_shell_final_one_electron_hamiltonian(
    final_kinetic_result,
    final_nuclear_by_center_results,
)
```

or include `final_basis` if the implementation needs it for dimension checks.

The helper should:

- require a materialized final kinetic one-body result from
  `pqs_source_shell_final_one_body_from_boundary_matrix(...)`;
- require one or more materialized final nuclear by-center results from
  `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(...)`;
- validate all final matrix shapes agree;
- validate finite and symmetric inputs;
- validate every nuclear input has:
  - nuclear charge recorded;
  - nuclear charge not applied;
  - centers not summed;
  - uncharged by-center convention;
- apply the recorded nuclear charge at this Hamiltonian assembly stage;
- sum the charged nuclear center contributions;
- build `hamiltonian_matrix = kinetic + charged_nuclear_sum`;
- report center count, applied charges, final dimension, finite/symmetry
  diagnostics, and explicit nonclaims:
  no eigensolve, no generalized overlap solve, no IDA, no density-density, no
  RHF, no driver route, no export, no artifact.

Trust boundary:

Do not run H1/eigensolve in this pass. Do not build a physics acceptance probe.
Do not implement IDA, density-density, RHF, drivers, exports, or artifacts. Do
not call `_pqs_current_route_safe_term_matrices(...)`.

Test policy:

Add at most one compact module-contract test. It should use synthetic final
kinetic and nuclear-by-center result objects and verify:

- charge application happens only in this helper;
- separated centers are summed here;
- `H = T + sum(Z_i * V_i)`;
- nonclaims remain false for solves and downstream physics.

Do not add a real projected-q-shell integration test in this pass. A real H1
probe should be the next pass after assembly exists.

Expected blocker after this pass:

The next blocker should advance to:

```text
:missing_pqs_final_one_electron_solve
```

or a more precise final-H1 probe blocker. Do not claim H1 is accepted until the
ordinary final-basis eigensolve is run and compared.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- focused module-contract test if one is added or edited;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Report back:

- write `.agent_handoffs/response.026.md.tmp`, then atomically rename to
  `.agent_handoffs/response.026.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.026.md`;
- include files changed;
- include exact next blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
