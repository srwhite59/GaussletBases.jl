Purpose:

Audit the PQS shell-support one-body operator source boundary before adding
more implementation. Determine whether existing CPBM retained-source
overlap/kinetic blocks already represent the boundary operators:

```text
O_boundary = P' * O_shell_support * P
```

for the real projected-q-shell fixture, or whether a new shell-support operator
source is truly required.

Why now:

Pass 021 added the final-basis projection seam for a caller-supplied
shell-support one-body operator. The next blocker is:

```text
:missing_pqs_shell_support_one_body_operator_source
```

Do not guess this seam. Earlier drift came from treating raw retained-source
operators as if Lowdin cleanup alone made them final. This pass should establish
whether the existing retained-source blocks are actually boundary-projected
shell operators on the real fixture.

Exact task:

Create a `tmp/work` probe. Do not change production code or tests unless the
audit exposes a tiny missing introspection seam.

Use the real projected-q-shell setup from pass 020:

- `current_box = (1:5, 1:5, 1:5)`
- `inner_box = (2:4, 2:4, 2:4)`
- `q = 5`, `L = 5`
- old `_pqs_shell_realization_plan(...)` as oracle/input only
- `CPBM.pqs_source_shell_realization_final_basis(...)`

Audit overlap first:

1. Build the CPBM raw product source plan and PQS boundary retained rule.
2. Build the CPBM retained-source overlap matrix using the existing
   source-box-first retained-source helper.
3. Build shell support overlap using
   `_pqs_product_box_support_overlap_matrix(...)` as oracle/input.
4. Compute:

```text
S_boundary_oracle = P' * S_shell_support * P
S_final_from_boundary = L' * S_boundary_oracle * L
```

5. Compare:

```text
retained_source_overlap vs S_boundary_oracle
S_final_from_boundary vs final_basis.final_overlap
```

Then audit kinetic if the fixture setup has a clean existing CPBM retained-source
kinetic path and old shell/support oracle for kinetic. If kinetic requires
calling `_pqs_current_route_safe_term_matrices(...)`, that is allowed only as an
oracle comparison in the probe and must be labeled as debug/oracle, not route
authority.

Do not call old fixed-block/current-route helpers as production inputs.

Questions to answer in the response:

- Does retained-source overlap equal the shell-projected boundary overlap?
- Does Lowdin-transforming that boundary overlap reproduce the final identity
  overlap?
- Is the same statement currently demonstrated for kinetic?
- If yes, what is the smallest implementation target for pass 023?
- If no, what exact source object is missing?

Trust boundary:

No H1 solve, no Hamiltonian assembly, no charge summing, no electron-nuclear
implementation, no IDA, no density-density, no RHF, no driver route, no exports,
and no artifacts.

Test policy:

No new permanent test in this pass. Use `tmp/work` artifacts only.

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

- write `.agent_handoffs/response.022.md.tmp`, then atomically rename to
  `.agent_handoffs/response.022.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.022.md`;
- include probe artifact path;
- include overlap and kinetic comparison numbers or exact blockers;
- include the recommended pass-023 implementation target;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
