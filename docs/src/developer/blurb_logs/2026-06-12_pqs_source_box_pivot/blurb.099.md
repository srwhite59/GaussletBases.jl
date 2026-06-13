Purpose:
Audit duplicate H1/J probe glue and cleanup candidates now that the driver-owned
private H1/J diagnostic path materializes.

Why now:
Passes 097 and 098 made the private H1/J diagnostic consume route-owned density
inputs from the source-plan bundles and carry the density-input payload inside
the compact route payload. Before adding tests, RHF, or fixture policy, we need
to identify stale probe/test pressure that duplicates the accepted driver-owned
diagnostic path.

Governing framework:
Use `docs/src/developer/pqs_source_box_operator_framework.md` and `AGENTS.md`
test/deletion policy.

Keep these boundaries sharp:

- source-box-first PQS is the algorithmic framing;
- shell/support-row contraction is oracle/debug;
- retained diagnostic weights are not IDA/quadrature weights;
- H1/J remains diagnostic/private until explicitly promoted.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
approved commands and read-only inspection. If approval would be required, write
`ATTENTION.md` with the exact command, reason, and blocker, then stop.

Exact task:
Do a no-edit cleanup audit.

Inspect likely duplicate/probe surfaces:

- `tmp/work/` H1/J timing or support-density probes, if present;
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
- `src/pqs_multilayer_complete_core_shell_h1.jl`;
- `src/pqs_multilayer_support_density.jl`;
- `src/pqs_source_box_route_driver_helpers.jl`;
- recent blurb log responses/reviews around passes 095-098.

Answer:

1. Which probe-local constructions are now duplicated by the driver-owned
   private H1/J diagnostic path?
2. Which comparisons remain valuable as oracle/debug checks?
3. Which ignored `tmp/work` artifacts can be deleted locally?
4. Which tracked test helpers or assertions should remain for now?
5. Which tracked test helpers or assertions are candidates for deletion or
   quarantine in the next pass?
6. Is one compact driver H1/J smoke needed later, and would it replace any
   current probe pressure?

Do not edit files.
Do not delete files.
Do not commit.
Do not add tests.
Do not run broad tests.
Do not add RHF/SCF/Fock.
Do not add GTO, IDA/MWG, exports, artifacts, fixture promotion, or production
route behavior.

Validation:
Read-only inspection only. Use `rg`, `git status`, and file reads. Do not run
Julia unless a specific source question cannot be answered by inspection.

Report back:

- cleanup inventory grouped as:
  - delete-local-ignored-now;
  - delete-tracked-next-pass;
  - quarantine/debug-only;
  - keep-active-contract;
- reason for each candidate;
- exact remaining caller/blocker for anything not deleted;
- whether a compact driver H1/J smoke should replace current probe pressure;
- git status;
- deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.099.md`, continue polling for
`blurb.100.md`, `STOP.md`, or `ATTENTION.md`.
