Pass 120 - design private RHF SCF-control/acceleration contract, no code

Baseline:

- Current pushed HEAD should include `b9f7e36a Record PQS RHF residual payload probe`.
- Private compact RHF SCF reports residual diagnostics but remains blocked by
  `:scf_not_converged` with commutator residual around `1.34e-5`.

Task:

No code. Design the smallest private SCF-control/acceleration contract for the
diagnostic RHF path.

Context:

- Simple scalar damping did not converge the compact route-smoke fixture within
  100 iterations.
- Fock/energy convention audit found no formula bug.
- Residual-aware private payload now reports:
  - fixed-point density change;
  - energy change;
  - trace/idempotency errors;
  - ordinary final-basis commutator residual.

Questions to answer:

1. Control object:
   - What compact private control object or NamedTuple should configure SCF?
   - Fields might include `mixing_kind`, `max_history`, `mixing_alpha`,
     `residual_metric`, `density_atol`, `energy_atol`, `residual_atol`,
     `max_iterations`.
   - Which fields are required versus optional defaults?

2. First acceleration strategy:
   - Should the first implementation be:
     - residual-aware scalar damping with better convergence labels;
     - density DIIS;
     - Fock DIIS;
     - or a local-only DIIS experiment before production private code?
   - Explain the choice.

3. Residual/history definition:
   - What vector should be used for DIIS/residual history?
     Examples: flattened commutator `[F, P]`, flattened density delta,
     flattened Fock delta.
   - What history size is appropriate for the compact route-smoke fixture?
   - How should trace/idempotency be preserved or checked after mixing?

4. Convergence labels:
   - What does `rhf_converged = true` require?
   - Should convergence require density delta, energy delta, and commutator
     residual, or only a subset?
   - What blocker/status should be used for max-iteration failure?

5. Implementation boundary:
   - Which file/helper should be touched first if implementation is approved?
   - Should route-driver integration remain blocked? I expect yes.
   - Should default tests remain synthetic only? I expect yes.

Decision rules:

- Do not edit files.
- Do not run the compact route probe.
- Do not add damping/DIIS/mixing code.
- Do not route-wire RHF.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- `git status --short --branch` is enough.

Report back:

- Recommended first implementation pass.
- Exact control fields and labels.
- Exact tests to add/update.
- Explicit exclusions.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
