Purpose:

Run the first real projected-q-shell final-basis H1 probe through the new PQS
source-box/final-basis seams.

Why now:

The route now has:

- final-basis object from shell projection and Lowdin cleanup;
- final overlap/kinetic from retained boundary operators;
- final by-center nuclear matrices with charge still unapplied;
- Hamiltonian-stage charge application and center summation.

The next blocker is:

```text
:missing_pqs_final_one_electron_solve
```

Exact task:

Create a `tmp/work` probe. Do not change production code or tests unless the
probe exposes a tiny missing seam.

Use the real projected-q-shell fixture from passes 020, 022, and 024:

- `current_box = (1:5, 1:5, 1:5)`
- `inner_box = (2:4, 2:4, 2:4)`
- `q = 5`, `L = 5`
- CPBM/CRPS raw source plan and PQS boundary retained rule
- `pqs_source_shell_realization_final_basis(...)`

Build the final one-electron path using the new CPBM seams:

1. retained-source overlap and kinetic from existing CPBM source-box helpers;
2. final overlap/kinetic from `pqs_source_shell_final_one_body_from_boundary_matrix(...)`;
3. one origin center, charge `1.0`, using the retained centered
   electron-nuclear by-center helper;
4. final by-center nuclear from
   `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(...)`;
5. final Hamiltonian from
   `pqs_source_shell_final_one_electron_hamiltonian(...)`;
6. ordinary symmetric eigensolve of the final Hamiltonian matrix.

Do not use a generalized overlap solve. The final overlap should be identity to
roundoff before the eigensolve.

Oracle comparison:

Use shell-support/oracle data only as a comparison, not as route authority.
Compare at least:

- final overlap identity error;
- final Hamiltonian matrix versus a shell-support oracle final Hamiltonian, if
  straightforward;
- lowest H1 eigenvalue from the new path versus the oracle eigensolve, if
  straightforward.

If a direct oracle comparison is awkward, report the exact missing comparison
seam and still report final overlap identity, Hamiltonian symmetry/finite
diagnostics, and lowest eigenvalue.

Trust boundary:

Probe only. Do not add a permanent acceptance test. Do not implement IDA,
density-density, RHF, drivers, exports, or artifacts. Do not call
`_pqs_current_route_safe_term_matrices(...)` unless it is strictly an oracle
comparison and clearly labeled as such; prefer the direct shell-support oracle
construction used in earlier probes.

Questions to answer:

- Does the final overlap remain identity to roundoff?
- Is the Hamiltonian finite and symmetric?
- What is the lowest ordinary final-basis H1 eigenvalue?
- Does it match the shell-support oracle, if available?
- What exact implementation or cleanup should pass 028 do?

Test policy:

No permanent test in this pass. Use `tmp/work` artifacts only.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new probe artifact was added and why it earned its temporary
  carrying cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. <tmp/work probe>`;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Report back:

- write `.agent_handoffs/response.027.md.tmp`, then atomically rename to
  `.agent_handoffs/response.027.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.027.md`;
- include probe artifact path;
- include energy/matrix comparison diagnostics;
- include recommended pass-028 target;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
