Pass 168 response - private Be2 PQS CR2 artifact generator path

Status: complete
Role: repo-doer@macmini
Commit: uncommitted; repo-manager owns commits.

Generator path:
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- The script was force-added with `git add -f` because `tmp/` is ignored.

Generated output paths:
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2`
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv`

Generator behavior:
- Builds the same probe-enabled Be2/PQS assembly used by the focused Be2 fingerprint test.
- Calls the private pass-167 payload/writer helpers.
- Writes the ignored JLD2 and TSV outputs under `tmp/work/be2_wl_pqs_cr2_inspection_artifact/`.
- Prints both paths and a short readiness summary:
  - `pqs_status=available_diatomic_complete_core_shell_hamiltonian_handoff_payload`
  - `cr2_read_only_inspector_ready=true`
  - `cr2_solver_ready=false`
  - `white_lindsey_status=unavailable`

Top-level JLD2/fingerprint confirmation:
- `schema_name=be2_wl_pqs_handoff_inspection_bundle`
- `schema_version=1`
- `has_pqs=true`
- `has_wl=true`
- fingerprint header includes `route_label`, `status`, `blocker`, dimensions, shapes, H1 scalar, finite flags, CR2 readiness flags, and system scalar columns.

Producer/fixture provenance included:
- `producer/package`
- `producer/repo_commit`
- `producer/dirty`
- `producer/generated_at`
- `producer/generator_entrypoint`
- `fixture/atom_symbols`
- `fixture/q`
- `fixture/n_s`
- `fixture/reference_spacing`
- `fixture/tail_spacing`
- `fixture/parent_axis_counts`
- `fixture/map_backend`
- `fixture/route_kind`

Validation results:
- `julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl` passed and wrote both artifacts.
- `ls -l tmp/work/be2_wl_pqs_cr2_inspection_artifact` showed:
  - `be2_wl_pqs_handoff_fingerprint.tsv`
  - `be2_wl_pqs_handoff_inspection_bundle.jld2`
  - `generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- JLD2/fingerprint readback command passed and confirmed schema, routes, producer commit/dirty, fixture `q=5`, and fixture `n_s=5`.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed: `load ok`.
- `git diff --check` passed.
- Additional staged check `git diff --cached --check` passed for the force-added ignored script.
- `git diff --numstat -- src test` produced no output; no `src` or `test` files were edited, so the source/test line-budget rule was not triggered.

Git status before writing this response:
`## main...origin/main`
`A  tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`

Generated-output tracking check:
- The generated `.jld2` and `.tsv` outputs did not appear in `git status --short --branch`.
- Only the force-added generator script appeared.

Deletion/shrinkage report:
- deleted: none.
- simplified: CR2 now has a repo-owned private generator path instead of needing to call underscored writer helpers directly.
- quarantined: generated JLD2/TSV artifacts remain ignored under `tmp/work`; WL remains unavailable/not-applicable placeholder data in the existing schema.
- not deleted because: this pass only added the private generator path and cheap provenance; no source/test cleanup was needed.
- exact remaining caller/blocker: WL arrays are still not populated; CR2 solver/export handoff remains blocked by `:missing_cr2_solver_handoff_format`; downstream solver/HamV6/HFDMRG readiness remains blocked by `:missing_hfdmrg_density_density_contract` and downstream contract objects.

-- repo-doer@macmini
