Pass 111 - no-edit review of private RHF SCF payload contract

Baseline:

- Current pushed HEAD should include `672886c1 Record PQS pass 110 audit`.
- Private SCF helper added in pass 109:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`
  in `src/pqs_multilayer_complete_core_shell_rhf.jl`.
- Do not implement route wiring from pass 110 yet.

Task:

No code. Audit the private SCF payload against the pass-102 RHF contract and the
pass-109 implementation intent.

Focus questions:

1. Energy and density consistency:
   - Does `final_one_step_payload` correspond to the returned `final_density`,
     or to the previous density used to build the Fock matrix?
   - Is that acceptable under the current convergence rule, or should the
     payload labels be tightened before any real probe?
   - Is `final_total_energy` clearly the energy of the density used in the
     final one-step payload?

2. Convergence semantics:
   - Is the first-iteration rule (`energy_change = nothing`, energy convergence
     treated as satisfied) acceptable when density convergence is still
     required?
   - Does `:scf_not_converged` return the last attempted one-step payload and
     final density consistently?

3. Density/occupation convention:
   - Does the loop preserve spin-summed final density with trace equal to
     `electron_count`?
   - Is the ordinary final-basis metric assumption explicit enough?
   - Are open-shell/fractional occupations still blocked by the input contract?

4. Payload carrying cost:
   - Are large matrices kept only in top-level payload fields, not duplicated
     in summaries/metadata?
   - Are iteration records compact enough?

5. Test pressure:
   - Are the four current RHF test files still justified as separate fast
     contract tests for now?
   - If not, recommend consolidation, but do not edit files in this pass.

Deliverables:

- State whether the SCF payload is ready for a local real compact PQS probe.
- If not ready, give the exact blocker and the smallest corrective pass.
- If ready, recommend the next pass as a `tmp/work` local real compact PQS SCF
  probe, not a tracked test and not route wiring.

Decision rules:

- Do not edit files.
- Do not run Julia tests unless a tiny grep/read-only inspection is
  insufficient. This is primarily a code review/audit.
- Do not run the heavy source-box dry-run.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- `git status --short --branch` is enough.

Report back:

- Audit result.
- Ready/not ready for real compact local probe.
- Any exact corrective task.
- Test consolidation recommendation.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
