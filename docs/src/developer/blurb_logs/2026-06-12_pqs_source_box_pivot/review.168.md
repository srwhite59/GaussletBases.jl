Pass 168 review - accepted

Reviewed implementation commit:

```text
f679ba86 Add Be2 CR2 artifact generator
```

Verdict: accepted.

This pass adds the repo-owned private generator path that CR2 requested:

```text
tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

The script is tracked with `git add -f` because `tmp/` is ignored. It writes
ignored runtime outputs under the same directory:

```text
tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
```

The generated outputs remain ignored and are not committed.

The generator builds the same probe-enabled Be2/PQS assembly as the focused
test, calls the private pass-167 payload/writer, adds cheap producer and fixture
provenance, writes JLD2/TSV, and prints artifact paths plus a short readiness
summary.

Validation:

```text
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

passed and printed the JLD2/TSV paths with:

```text
pqs_status=available_diatomic_complete_core_shell_hamiltonian_handoff_payload
cr2_read_only_inspector_ready=true
cr2_solver_ready=false
white_lindsey_status=unavailable
```

JLD2/TSV readback confirmed schema name/version, both route groups, fixture
`q=5`, fixture `n_s=5`, and the TSV header.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --cached --check
```

passed.

No `src` or `test` files were edited, so the source/test line-budget rule was
not triggered.

Important note:

The pre-commit generated artifact reported `dirty=true` because the generator
script was staged/uncommitted at generation time. After publishing this pass,
the manager should rerun the generator once on clean HEAD so the ignored local
artifact records the committed revision with a clean dirty marker.

Remaining blockers:

- WL route remains placeholder-only.
- CR2 solver/export handoff remains blocked by
  `:missing_cr2_solver_handoff_format`.
- Downstream solver/HamV6/HFDMRG readiness remains blocked by
  `:missing_hfdmrg_density_density_contract` and related downstream contracts.

-- repo-manager@macmini
