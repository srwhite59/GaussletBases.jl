Pass 133 - fingerprint Be2 PQS Ham payload blocker

Purpose:

Add a focused Be2/diatomic PQS Ham readiness fingerprint. The goal is to prove
exactly where the current Be2 PQS route blocks relative to the new private Ham
payload seam.

Why now:

Pass 132 showed that the medium-term CR2 Be2 WL/PQS comparison is blocked first
by the PQS Be2/diatomic side: one-center PQS has a private Ham payload, but Be2
PQS still lacks a route-owned final-basis/H1/Ham payload seam. We need a small
tracked fingerprint before designing or coding the missing seam.

Target:

```text
Be2 PQS staged route
-> cartesian_assembly(...)
-> assembly.complete_core_shell_diagnostic_route_payload
-> complete_core_shell_ham_payload status/blocker/missing_inputs
```

Exact task:

- Add one focused test file, for example:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Do not add it to the default nested runner.
- Build the same Be2 PQS route family used by existing route-driver report
  tests, but stop at `cartesian_assembly(...)` rather than using
  `cartesian_report(...)` aliases.
- Reuse the existing Be2 PQS route parameters from
  `test/nested/pqs_source_box_route_driver_report_runtests.jl`:
  - `route_family = :pqs_source_box`
  - route kind like `:be2_cartesian_nesting_route_driver_spine`
  - `atom_symbols = ("Be", "Be")`
  - `nuclear_charges = (4, 4)`
  - `atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))`
  - `parent_axis_counts = (x = 9, y = 7, z = 9)`
  - `q = 5`, `n_s = 5`
  - `route_shape = (:pqs_left, :product, :pqs_right)`
  - standard source-box retained rules and terms already used there.

What to assert:

- `assembly.object_kind == :cartesian_assembly`
- `assembly.route_skeleton.route_family === :pqs_source_box`
- `assembly.complete_core_shell_diagnostic_route_payload` exists.
- Inspect `payload = assembly.complete_core_shell_diagnostic_route_payload`.
- Inspect `ham = payload.complete_core_shell_ham_payload`.
- Assert the observed current status/blocker/missing inputs for Be2 PQS.

Decision rule for expected status:

- If `ham.status` is blocked, commit the fingerprint with the exact observed
  blocker and `missing_inputs`.
- If `ham.status` unexpectedly materializes, assert only compact facts:
  final dimension, one-body Hamiltonian status, density interaction status,
  electron-electron representation, and nonclaim flags. Do not add route/export
  behavior.
- If constructing the Be2 `cartesian_assembly(...)` path takes too long or
  requires broad report/materialization plumbing, stop and report the exact
  bottleneck instead of adding scaffolding.

Trust boundary:

- Test/fingerprint only.
- No production code unless the fingerprint exposes a tiny typo.
- No WL payload implementation.
- No export/artifact writing.
- No public API.
- No RHF/SCF work.
- No hfdmrg or CR2 execution.
- No fixture promotion.
- No IDA/MWG semantic change.
- Do not parse broad report aliases as authority.

Validation:

- Run the focused new test directly.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- If the focused test takes more than 60 seconds cold, report elapsed time and
  whether precompilation dominated.

Report back:

- Commit SHA if committed.
- Observed Be2 PQS Ham payload status/blocker/missing inputs.
- Files changed.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
