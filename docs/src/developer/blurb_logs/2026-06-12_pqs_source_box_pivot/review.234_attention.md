# Pass 234 attention review - exception granted

The `ATTENTION.md` line-budget blocker is accepted as a legitimate stop.

Manager decision:

- Grant a one-pass line-budget exception for the compact support-region
  materializer implementation described in `ATTENTION.md`.
- The allowed positive budget is approximately the reported compact diff:
  `119` added / `10` deleted, net `+109` in `src + test + bin`.
- This exception applies only to pass 234.

Reason:

- The doer reports that the implementation generated the support counts
  `(275, 578, 362)` from route geometry/shellification grouping.
- This is real route-authority progress: it replaces target constants with a
  generated independent H2 PQS support-region plan.
- The safe local deletion candidates are exhausted, and forcing unrelated
  deletion into this pass would be higher risk than accepting the small debt.

Future handling:

- A separate old-flat-Cartesian-path audit is being started to identify deletion
  candidates for future passes.
- Future passes should use that audit to pay down line debt; pass 234 should not
  bundle that broader cleanup.

-- repo-manager@macmini
