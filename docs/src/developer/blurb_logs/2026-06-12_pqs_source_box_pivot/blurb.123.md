Pass 123 - tighten RHF SCF final convergence gating

Baseline:

- Pass 122 changes are present in the working tree but are not accepted or
  committed yet.
- Current pushed HEAD should still include `a807fa24 Publish PQS blurb 122`.
- Manager review for pass 122 is `.agent_handoffs/review.122.md`.

Task:

Fix the convergence gating issue in the pass-122 private Fock-DIIS SCF
implementation.

Problem:

- `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` currently decides
  convergence from iteration-input density/Fock checks.
- It then recomputes the final one-step payload and residual diagnostics for
  the returned `final_density`.
- But it can still return `rhf_converged = true` even if the recomputed final
  residual diagnostics are not converged.

Required correction:

- In the converged branch, after recomputing the final one-step payload and
  residual diagnostics, require the final diagnostics to pass:
  - residual `<= residual_atol`;
  - trace error `<= trace_atol`;
  - idempotency error `<= idempotency_atol`;
  - final density match/shape checks already present should remain.
- Only then return
  `:materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload` with
  `rhf_converged = true`.
- If iteration-input convergence passes but recomputed final diagnostics fail,
  return the blocked `:scf_not_converged` payload, carrying:
  - returned/final density;
  - recomputed final-one-step payload;
  - residual diagnostics;
  - compact counters and iteration history.
- Keep iteration records unchanged as history.
- Keep route/report/public/export/artifact nonclaims.

Tests:

- Update `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`.
- Existing tests should still pass.
- Add a small assertion that the converged synthetic case has final
  `summary.convergence_diagnostics.residual_converged`,
  `trace_converged`, and `idempotency_converged` all true for the recomputed
  final residual diagnostics.
- Add a blocked-path test only if a tiny synthetic fixture can exercise final
  residual failure without overbuilding.

Exclusions:

- Do not add route wiring.
- Do not run compact route probes.
- Do not add public API/report fields/exports/artifacts/GTO/IDA/MWG.
- Do not change the Fock-DIIS formula unless needed for this convergence-gating
  fix.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- Exact final convergence gating behavior.
- Any additional status/blocker labels.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
