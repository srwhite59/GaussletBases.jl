Purpose:

Implement the narrow CPBM helper that converts a proven PQS retained boundary
one-body operator into the shell-realized final basis:

```text
O_final = L' * O_boundary * L
```

for overlap and kinetic only.

Why now:

Pass 022 showed on the real projected-q-shell fixture that existing CPBM
retained-source overlap and kinetic blocks equal the shell-projected boundary
operators `P' O_shell_support P` to roundoff. That clears the algebraic and
provenance concern for overlap/kinetic. It does not clear electron-nuclear or
H1.

Exact task:

Add a small CPBM helper, likely near
`src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`.

Preferred API shape:

```julia
pqs_source_shell_final_one_body_from_boundary_matrix(
    final_basis,
    retained_boundary_result;
    term = nothing,
)
```

where `retained_boundary_result` is the result returned by
`pqs_retained_source_one_body_matrix(...)`.

If a matrix-only method is useful internally, make it explicit that the matrix
space must be `:retained_pqs_source_modes` / boundary source modes. Do not make
the helper look like it accepts arbitrary raw source-space operators.

The helper should:

- require `final_basis.status == :available_pqs_shell_realization_final_basis`;
- require `retained_boundary_result.object_kind ==
  :pqs_retained_source_one_body_matrix`;
- require `retained_boundary_result.matrix_materialized == true`;
- require `retained_boundary_result.matrix_space == :retained_pqs_source_modes`;
- accept only retained-source overlap/kinetic terms in this first pass;
- validate shape against `final_basis.boundary_source_mode_count`;
- validate finite and symmetric boundary operator entries;
- compute `final_operator = L' * O_boundary * L`;
- return a compact result with:
  - status/blocker;
  - term;
  - boundary and final operator shapes;
  - final operator matrix;
  - finite/symmetry diagnostics;
  - provenance that the input was a retained boundary operator, not a raw source
    operator;
  - explicit nonclaims: no electron-nuclear, no charge summing, no H1 solve, no
    Hamiltonian assembly, no IDA, no density-density, no RHF, no driver route,
    no export, no artifact.

Trust boundary:

Do not generate shell-support operators in this pass. Do not call
`_pqs_current_route_safe_term_matrices(...)`. Do not run H1. Do not assemble
Hamiltonians. Do not implement electron-nuclear, IDA, density-density, RHF,
drivers, exports, or artifacts.

Test policy:

Add at most one compact module-contract test for the helper. It can use the
existing synthetic final-basis fixture and a synthetic retained-source
overlap/kinetic result if the current test setup already has a clean way to
build one.

Use a `tmp/work` real-fixture probe if useful to confirm:

- final overlap from retained-source boundary overlap is identity to roundoff;
- final kinetic from retained-source boundary kinetic matches the pass-022
  shell-support projection result to roundoff.

Do not add the real-fixture probe as a permanent test.

Expected blocker after this pass:

The exact next blocker should become:

```text
:missing_pqs_shell_boundary_electron_nuclear_operator_source
```

or a more precise name if the implementation exposes one. Do not claim H1 is
available.

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

- write `.agent_handoffs/response.023.md.tmp`, then atomically rename to
  `.agent_handoffs/response.023.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.023.md`;
- include files changed;
- include any probe artifact path;
- include exact next blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
