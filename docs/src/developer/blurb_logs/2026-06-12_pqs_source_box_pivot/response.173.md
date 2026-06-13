Pass 173 generator shrink complete.

Files changed:
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`

Tracked implementation line budget:
```text
2	3	tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```
- Added: 2
- Deleted: 3
- Net: -1
- Requirement satisfied: deleted > added.

Artifact size:
- Before: `242M`
- After: `202M`
- Shrink: about `40M`

Generator runtime:
- Measured with Julia-level `@elapsed` wrapper:
  - `generator_elapsed_s=240.012615125`

WL overlap convention/readback:
- `routes/white_lindsey/final_basis/overlap_convention = orthonormal_identity_diagnostic_checked_not_stored`
- `routes/white_lindsey/final_basis/overlap_matrix_stored = false`
- `routes/white_lindsey/final_basis/overlap_identity_defect = 3.9895900571894835e-14`
- `routes/white_lindsey/one_body/overlap` present: `false`

PQS overlap readback:
- `routes/pqs_source_box/final_basis/overlap_convention = orthonormal_identity_by_contract`
- `routes/pqs_source_box/final_basis/overlap_matrix_stored = false`
- `routes/pqs_source_box/one_body/overlap` present: `false`

H1/Vee readback:
- PQS H1 shape: `(221, 221)`
- PQS Vee shape: `(221, 221)`
- WL H1 shape: `(2287, 2287)`
- WL Vee shape: `(2287, 2287)`
- PQS solver/export: `false` / `false`
- WL solver/export: `false` / `false`

Generated artifact provenance/readback:
- `producer_commit=e94dff01ec9a60479190a0e85d2f8f271f48b438`
- `producer_dirty=true` before commit, expected

Validation commands/results:
- `julia --project=. -e 't = @elapsed include("tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl"); println("generator_elapsed_s=", t)'`
  - Passed.
  - Printed:
    - `pqs_status=available_diatomic_complete_core_shell_hamiltonian_handoff_payload`
    - `cr2_read_only_inspector_ready=true`
    - `cr2_solver_ready=false`
    - `white_lindsey_status=available_route_configured_diatomic_atom_growth_ham_payload`
    - `generator_elapsed_s=240.012615125`
- JLD2 readback command passed; values summarized above.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`
- `git diff --check`
  - Passed.
- `git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - Passed line-budget rule; details above.

Generated outputs ignored status:
```text
 M tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
```

Final `git status --short --branch`:
```text
## main...origin/main [ahead 1]
 M tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Deletion/shrinkage report:
deleted: dense downstream WL overlap matrix from the generated artifact payload
simplified: WL now matches PQS as orthonormal final-basis downstream `H1 + Vee`, with overlap retained only as an identity-defect diagnostic
quarantined: WL overlap matrix remains a GaussletBases-side diagnostic/check, not CR2 working data
not deleted because: WL still computes overlap to validate the identity defect; PQS/WL H1 and final-basis Vee matrices remain required CR2 inspection data
exact remaining caller/blocker: HF/solver/export promotion remains blocked by `:missing_reviewed_density_density_hf_fock_energy_convention`

-- repo-doer@macmini

Manager post-commit regeneration note:

After committing the generator shrink as `adb55213`, the generator was rerun.
The ignored JLD2 artifact now reports:

```text
producer_commit=adb55213b51264aa1956fd46a83a4e74e34a2d6c
producer_dirty=false
bundle_size=202M
```

The WL overlap matrix remains absent from `routes/white_lindsey/one_body`, and
PQS/WL H1 and final-basis Vee shapes remain `(221, 221)` and `(2287, 2287)`.
