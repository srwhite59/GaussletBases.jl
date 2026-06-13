Pass 127 - update private Fock-DIIS default history, focused validation

Purpose:

Make the one small tracked private-control default change supported by passes
125 and 126: resolve Fock-DIIS `max_history` to 8 by default instead of 6.

Why now:

Two local ignored compact route-smoke probes showed that:

- explicit `max_history = 6` plateaus under the strict `1e-8` density/residual
  gates;
- explicit `max_history = 8` converges twice at iteration 34 under the same
  strict gates;
- history 8 has bounded coefficient-pathology fallbacks at iterations 8, 9,
  and 10, with no solve failures.

This supports changing the private Fock-DIIS default history. It does not
support tolerance changes, route wiring, or any serious-HF claim. Serious HF is
handled by `codexhome/work/hfdmrg`, and downstream Cr2 pressure testing should
come later through the CR2 agent when this line is ready.

Exact task:

1. In `src/pqs_multilayer_complete_core_shell_rhf.jl`, update
   `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)` so that
   `mixing_kind === :fock_diis` with `max_history === nothing` resolves to 8
   instead of 6.

2. Update the focused expectation in
   `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl` from
   `fock_diis_control.max_history == 6` to `== 8`.

3. Do not change explicit user-supplied `max_history`; only the default
   resolution changes.

Trust boundary:

- Private RHF/SCF diagnostic control only.
- No tolerance changes.
- No regularization or coefficient guard changes.
- No route-driver wiring.
- No report fields or public API.
- No GTO/IDA/MWG/export/artifact behavior.
- No fixture promotion and no `hfdmrg` or CR2 comparison.

Validation:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Decision rules:

- If the focused SCF test passes and the diff is only the private default plus
  the focused expectation, commit with a clear message such as:
  `Default private Fock DIIS history to 8`.
- If changing the default causes unrelated behavior changes or requires
  broader test edits, stop and report the exact failing assertion.
- Do not run the cold compact route probe in this pass unless the focused test
  unexpectedly fails and the failure cannot be localized otherwise.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Report back:

- Commit SHA if committed.
- Validation commands/results.
- Diff summary.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
