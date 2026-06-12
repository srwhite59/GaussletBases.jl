Review 075: accepted.

The response gives the needed boundary decision. `pqs_multilayer_shell_source_plan(...)`
should remain a tactical bridge, while the next route-owned boundary should be
a shellification/lowering-backed PQS region plan. That plan should consume
`CartesianShellification.ShellificationPlan` and
`CartesianTerminalLowering.lower_terminal_regions(...)` output and carry core
owned region, ordered shell layer records, source CPBs, support coverage, and
disjointness/duplicate fingerprints. It should not carry PQS descriptors,
Lowdin matrices, support operators, H1, IDA, RHF, or final-basis transfer data.

The fact split is specific enough to guide the next implementation pass:
geometry ownership, box-depth/layer construction, coverage, duplicate checks,
and provenance should move out of `pqs_multilayer_shell_source_plan(...)`;
PQS source realization should keep projected-q descriptors, shell realization
plans, collapsed shell-sector coefficients, retained counts, and support
ordering.

Additional manager guardrail:

- support-space dense operator helpers are accepted as H1 seam/oracle
  machinery for the current complete core/shell gate;
- they should not become the general PQS operator algorithm or the path for
  future density/RHF scaling without an explicit scale contract;
- future implementation should keep private probe/H1 support glue from
  outrunning the driver spine.

Validation:

- doer reported `git diff --check` passed;
- manager inspected the docs diff and found it concise and aligned with the
  boundary decision.

Deletion/shrinkage:

- no deletion expected in this docs-only audit;
- next implementation should make explicit-box source planning less central by
  introducing the shellification/lowering-backed entry point;
- after that lands, the explicit-box entry point should be documented or
  quarantined as a bridge, and probe/test code that treats `core_box` /
  `outer_box` as PQS route authority should shrink.

-- repo-manager@macmini
