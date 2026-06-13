Purpose:
Measure/design a compact driver H1/J smoke before adding any permanent test.

Why now:
The ignored stale probes are out of active circulation, and the driver-owned
private H1/J diagnostic path exists. However, the focused one-center source-box
route dry-run used in manager validation took about 124 seconds in a fresh Julia
process. We should not add that as a permanent smoke without understanding
runtime and whether a cheaper fixture or warm-process validation shape exists.

Governing framework:
Use `AGENTS.md` test runtime policy and
`docs/src/developer/pqs_source_box_operator_framework.md`.

Keep these boundaries sharp:

- source-box-first PQS is the algorithmic framing;
- shell/support-row contraction is oracle/debug;
- retained diagnostic weights are not IDA/quadrature weights;
- H1/J remains diagnostic/private until explicitly promoted.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
approved commands and focused local probes. If approval would be required, write
`ATTENTION.md` with the exact command, reason, and blocker, then stop.

Exact task:
Do a local-only timing/design probe for the future compact driver H1/J smoke.

Use `tmp/work/` for any probe script. Do not edit source, tests, docs, or tracked
fixtures.

Measure or determine:

1. Fresh-process elapsed time for the accepted one-center source-box dry-run
   path that reaches materialized private H1/J.
2. If practical, warm-process second-run elapsed time for the same dry-run in
   one Julia process.
3. Whether the 124-second cold time is mainly first-call compilation or repeated
   route construction.
4. Whether a smaller source-box route fixture can still exercise the same
   accepted path without changing semantics into a toy-only nonrepresentative
   check.
5. Whether a permanent tracked smoke should be:
   - skipped for now;
   - added but marked slow/manual;
   - added as a generated/local-only probe;
   - or added after route construction is made cheaper.

Do not add the smoke test in this pass.
Do not add caching/checkpointing.
Do not change lattice size in source/tests.
Do not add RHF/SCF/Fock.
Do not add GTO, IDA/MWG, exports, artifacts, fixture promotion, or production
route behavior.

Validation:
Use Julia-level timing with `@elapsed`, not `/usr/bin/time`.

Because the accepted dry-run is expected to exceed 60 seconds cold, it is
acceptable to run it once for this timing/design pass. Do not run broad suites.

Run:

`git status --short --branch`

Report back:

- timing commands and results;
- final dimension, H1 energy, self-Coulomb, and density gauge from the measured
  accepted dry-run;
- warm-vs-cold conclusion if measured;
- recommendation for permanent smoke shape;
- whether smaller lattice is advisable or would weaken the route contract too
  much;
- git status;
- deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.101.md`, continue polling for
`blurb.102.md`, `STOP.md`, or `ATTENTION.md`.
