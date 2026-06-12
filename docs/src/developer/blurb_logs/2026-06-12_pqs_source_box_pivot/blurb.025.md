Purpose:

Implement the narrow CPBM final-basis by-center nuclear helper for PQS:

```text
V_final(center) = L' * V_boundary(center) * L
```

where `V_boundary(center)` is the existing retained-source
electron-nuclear-by-center boundary block.

Why now:

Pass 024 proved on the real projected-q-shell fixture that the existing
retained-source centered nuclear block equals the shell-projected boundary
nuclear operator `P' V_shell P`, for both origin and off-origin centers, while
preserving charge-recorded/not-applied and centers-not-summed metadata.

Exact task:

Add a small by-center-specific helper in CPBM, likely near
`src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`.

Preferred API shape:

```julia
pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(
    final_basis,
    retained_boundary_result,
)
```

where `retained_boundary_result` is the result returned by:

```text
pqs_source_pair_retained_centered_electron_nuclear_by_center_block(...)
```

or the equivalent retained-source electron-nuclear by-center result.

The helper should:

- require `final_basis.status == :available_pqs_shell_realization_final_basis`;
- require a retained-source electron-nuclear by-center input:
  - term `:retained_source_electron_nuclear_by_center`;
  - block space `:retained_pqs_source_modes`;
  - materialized retained-source operator block;
  - by-center metadata present;
- validate boundary block shape against `final_basis.boundary_source_mode_count`;
- validate finite and symmetric boundary entries;
- compute `final_operator = L' * V_boundary * L`;
- preserve center metadata:
  - center key;
  - center index;
  - center location;
  - nuclear charge recorded;
  - nuclear charge value;
  - nuclear charge applied is false;
  - centers summed is false;
  - uncharged by-center convention is true;
- return a compact result with final matrix diagnostics and explicit nonclaims:
  no charge summing, no H1 solve, no Hamiltonian assembly, no IDA,
  no density-density, no RHF, no driver route, no export, no artifact.

Trust boundary:

Do not assemble H1 in this pass. Do not apply nuclear charge. Do not sum
centers. Do not run a one-electron solve. Do not call
`_pqs_current_route_safe_term_matrices(...)`. Do not implement IDA,
density-density, RHF, drivers, exports, or artifacts.

Test policy:

Add at most one compact module-contract test for this helper. It should use a
small synthetic retained-source nuclear-by-center result or reuse the existing
contract fixture; it should focus on:

- accepted term/space;
- final matrix transformation;
- center metadata preserved;
- charge unapplied and centers unsummed.

Do not add a real projected-q-shell integration test. If useful, run a
`tmp/work` real-fixture probe to confirm the final nuclear transform works for
the origin/off-origin blocks from pass 024, but keep it temporary.

Expected blocker after this pass:

The next blocker should advance to Hamiltonian-stage assembly:

```text
:missing_pqs_final_one_electron_hamiltonian_assembly
```

or a more precise charge-application/summing blocker if the implementation
exposes one. H1 is still not ready until charges are applied and centers are
summed at the Hamiltonian assembly stage.

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

- write `.agent_handoffs/response.025.md.tmp`, then atomically rename to
  `.agent_handoffs/response.025.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.025.md`;
- include files changed;
- include any probe artifact path;
- include exact next blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
