Pass 164 - draft CR2-facing read-only Hamiltonian inspection handoff

Role: repo-doer@macmini

This is a no-edit handoff/drafting pass. Do not change source, tests, docs, or
runtime behavior.

Purpose:

Pass 163 made the Be2 PQS Hamiltonian consumer readiness tuple CR2-inspectable
without making it solver/export ready. The next useful move is not another
GaussletBases payload. Prepare a concise, user-pastable blurb for a future
`cr2-doer` agent so CR2 can inspect the current GaussletBases Be2/PQS line and
tell us what concrete handoff format it wants.

Read:

- `AGENTS.md`
- `.agent_handoffs/response.163.md`
- `.agent_handoffs/review.163.md`
- `src/pqs_source_box_diatomic_complete_core_shell.jl` around the diatomic
  Hamiltonian handoff and consumer readiness helpers
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

If readable without escalation, you may also inspect only lightweight top-level
orientation files in:

- `~/Dropbox/codexhome/work/cr2`
- `~/Dropbox/codexhome/work/hfdmrg`

Use that only to make the CR2-facing prompt concrete. Do not run CR2, HFDMRG,
HF, DMRG, Julia package tests outside this repo, or any solver.

Current facts to preserve in the CR2-facing blurb:

- This is GaussletBases Be2/PQS private diagnostic/handoff work, not a public
  API or export.
- Current read-only inspection view lives at
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
- CR2 normally runs HF itself. GaussletBases should not try to become the HF
  runner.
- Qiu-White Hamiltonian corrections require atom-local HF; CR2 should tell us
  what it needs for that path.

Draft the CR2-doer blurb around these questions:

1. For a first Be2 WL-vs-PQS comparison, what exact Hamiltonian/basis data
   should GaussletBases expose?
2. Does CR2 want a HamV6-like file, an in-memory Julia object, a JLD2/JSON
   inspection artifact, or just a printed fingerprint first?
3. Is the current representation
   `one_body_hamiltonian + pre_final_density_interaction + final_to_pre_final`
   useful to CR2, or should GaussletBases convert it before handoff?
4. What atom-local HF inputs are required for Qiu-White corrections?
5. What is the smallest Be2 comparison fixture CR2 would use first?

Boundaries:

- Do not edit files.
- Do not add a new GaussletBases payload.
- Do not add source/test/docs.
- Do not run CR2/HFDMRG.
- Do not implement HamV6, dense `Vee`, sliced integrals, HF, RHF, DMRG,
  Qiu-White corrections, exports, artifacts, or public APIs.
- Do not ask the user for permission from inside the unattended loop. If you hit
  a permission or external-state blocker, write `.agent_handoffs/ATTENTION.md`
  and stop.

Report back with:

- the drafted CR2-doer blurb, ready for the user to paste;
- exact GaussletBases source/test fields you used;
- whether you inspected any CR2/HFDMRG orientation files, and which ones;
- git status.

-- repo-manager@macmini
