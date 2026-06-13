Pass 164 response - draft CR2-facing read-only Hamiltonian inspection handoff

Role: repo-doer@macmini

Task type: no-edit handoff/drafting pass. I edited no source, tests, docs, or
runtime behavior, and ran no CR2/HFDMRG/Juila solver commands.

Draft CR2-doer blurb:

```text
You are `cr2-doer` for the CR2 lane on macmini.

Purpose:
Inspect the current GaussletBases Be2/PQS private Hamiltonian handoff and tell
GaussletBases what concrete downstream handoff format CR2 wants for a first
Be2 WL-vs-PQS comparison.

Context:
GaussletBases has a private Be2/PQS diagnostic/handoff path. This is not a
public API, not an export, and not a solver-ready CR2/HFDMRG object.

Current GaussletBases facts:
- The read-only inspection view lives at
  `consumer_contract_payload.readiness`.
- `cr2_read_only_inspector_ready = true`.
- `cr2_solver_ready = false`.
- `cr2_export_ready = false`.
- `cr2_handoff_blocker = :missing_cr2_solver_handoff_format`.
- `two_body_representation_kind = :pre_final_density_interaction`.
- `density_gauge = :pre_final_localized_positive_weight`.
- `raw_pair_factor_convention = :raw_numerator`.
- Overall GaussletBases readiness still blocks on
  `:missing_hfdmrg_density_density_contract`.

Available GaussletBases handoff pieces:
- final one-body Hamiltonian reference in the PQS final basis;
- pre-final density-interaction representation;
- final-to-pre-final coefficients;
- support/pre-final weights and raw pair numerator/provenance;
- Be2 nuclear charges, coordinates, and nuclear repulsion;
- electron count 8 and closed-shell singlet spin convention.

Boundaries:
- Do not ask GaussletBases to run CR2 or HFDMRG.
- Do not treat this as HamV6, dense `Vee`, `V6`, `Vblocks`, public API,
  JLD2/export, artifact, RHF/SCF, H1/J, or production readiness.
- CR2 normally owns/runs HF. GaussletBases should not become the HF runner.
- Qiu-White Hamiltonian corrections require atom-local HF; identify what CR2
  needs for that path instead of pushing HF into GaussletBases.

Questions for CR2:
1. For the first Be2 WL-vs-PQS comparison, what exact Hamiltonian/basis data
   should GaussletBases expose?
2. Does CR2 want a HamV6-like file, an in-memory Julia object, a JLD2/JSON
   inspection artifact, or just a printed fingerprint first?
3. Is the current representation
   `one_body_hamiltonian + pre_final_density_interaction + final_to_pre_final`
   useful to CR2, or should GaussletBases convert it before handoff?
4. If conversion is needed, what representation and ordering should
   GaussletBases target: density-density `H,V`, sliced `V6`/`Vblocks`,
   HamV6-like blocks, or a CR2-specific inspection bundle?
5. What atom-local HF inputs are required for Qiu-White corrections?
6. What is the smallest Be2 comparison fixture CR2 would use first?
7. What exact fields should GaussletBases include in a read-only artifact so
   CR2 can compare WL and PQS dimensions, H1 spectrum, nuclear metadata,
   electron/spin convention, and missing downstream blockers?

Expected CR2 output:
- Recommend one first handoff format.
- State whether the current pre-final density-interaction representation is
  directly useful or must be converted.
- List exact required fields and ordering conventions.
- State the smallest Be2 fixture and comparison metrics CR2 would inspect
  before any solver run.
- State what remains blocked for solver/export readiness.
```

GaussletBases source/test fields used:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - consumer readiness fields:
    `cr2_read_only_inspector_ready`, `cr2_solver_ready`, `cr2_export_ready`,
    `cr2_handoff_blocker`, `two_body_representation_kind`, `density_gauge`,
    `raw_pair_factor_convention`;
  - overall readiness blocker:
    `:missing_hfdmrg_density_density_contract`;
  - handoff-owned objects: final H1 reference, pre-final density interaction,
    final-to-pre-final coefficients, support/pre-final weights, raw numerator,
    and Be2 nuclear/electron/spin metadata.
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - verifies the CR2 read-only view, false solver/export readiness, the
    pre-final representation labels, and the unchanged overall blocker.

Orientation files inspected:

- `/Users/srw/Dropbox/codexhome/work/cr2/AGENTS.md`
- `/Users/srw/Dropbox/codexhome/work/cr2/answers.md`
- `/Users/srw/Dropbox/codexhome/work/hfdmrg/README.md`

I used them only for lightweight orientation. I did not run CR2, HFDMRG, HF,
DMRG, or any downstream solver.

Git status:

```text
git status --short --branch
## main...origin/main
```

-- repo-doer@macmini
