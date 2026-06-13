Pass 114 response - compact PQS RHF iteration trace

Existing artifact decision:
- Existing pass-113 artifacts were not enough because they only saved the last SCF iteration, not `scf.iteration_records`.
- Reran a local ignored trace probe.

Local ignored artifacts:
- Updated `tmp/work/pqs_compact_route_owned_rhf_scf_probe.jl` so it can be included without immediately running.
- Added `tmp/work/pqs_compact_route_owned_rhf_scf_trace_probe.jl`.
- Wrote:
  - `tmp/work/pqs_compact_route_owned_rhf_scf_trace.tsv`
  - `tmp/work/pqs_compact_route_owned_rhf_scf_trace_summary.txt`

Command:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_scf_trace_probe.jl`

Elapsed time:
- `96.844649625` seconds, measured with Julia-level `@elapsed` inside the probe.
- This exceeded 60 seconds as expected and was reported before launch.

SCF status:
- status: `blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- blocker: `scf_not_converged`
- iteration count: `50`

Trace summary:
- last total energy: `-10.032175867411672`
- density change min: `4.4092162672509927e-05`
- density change max: `0.047547563330606474`
- density change last: `4.5438494082283842e-05`
- energy change min: `1.5979158796142201e-07`
- energy change max: `0.5849762487415795`
- energy change last: `1.750599221139737e-07`
- total energy monotone decreasing: `true`
- density change generally decreasing: `false`
- energy two-cycle signal: `false`
- density two-cycle signal: `false`

Qualitative diagnosis:
- Not unstable divergence: total energy decreases monotonically across 50 iterations.
- No obvious two-cycle in the last several density or energy values.
- Not simple slow monotone density convergence: density change drops rapidly through iteration 10-11, bottoms near `4.409e-5`, then slowly increases to `4.544e-5` by iteration 50.
- This looks like a stable-but-not-converging fixed-point update or residual plateau/drift, not a clear bookkeeping bug from the trace alone.

Recommended next pass:
- Run a local ignored controlled damping/mixing probe first, not production code.
- Suggested shape: reuse the same compact route-owned fixture, apply simple density mixing in the local probe only, compare a few damping factors such as `0.25`, `0.5`, and maybe `0.75`, and record whether density residual falls below tolerance.
- If local damping works, then consider a bounded private SCF-control seam. If it does not, audit the Fock/density residual definition before changing production code.

Validation/status:
- `git status --short --branch`

```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none
- simplified: none; this was a local trace probe only.
- quarantined: RHF SCF remains private route-smoke diagnostic output.
- not deleted because: ignored trace artifacts are useful local evidence for the next decision.
- exact remaining caller/blocker: compact route-owned private RHF SCF remains blocked by nonconvergence without damping/mixing; no route-driver caller should be added yet.

-- repo-doer@macmini
