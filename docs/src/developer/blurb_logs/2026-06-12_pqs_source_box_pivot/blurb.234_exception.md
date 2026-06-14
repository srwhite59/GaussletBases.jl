# Pass 234 exception addendum - proceed despite local line-budget blocker

Role: repo-doer.

Manager decision:

Proceed with the compact pass-234 implementation described in `ATTENTION.md`.
The line-budget rule is waived for this pass only, up to approximately the
reported compact implementation size:

```text
about 119 added / 10 deleted, net +109 in src + test + bin
```

Reason:

- The blocker is line budget, not technical route correctness.
- The reported implementation generates the independent H2 PQS support counts
  `(275, 578, 362)` from route geometry/shellification grouping.
- This replaces target constants with real support-region authority, which
  advances the independent H2 PQS route.
- Forcing unrelated deletion into this pass risks unsafe cleanup.

Constraints still apply:

- Implement only the support-region materializer/fingerprint.
- Keep retained transforms, source-plan authority beyond support regions, final
  basis, H1, H1-J, RHF, supplements, CR2, export, and public API blocked.
- Do not use fake-PQS/WL coefficient matrices or fixed-source retained
  transforms.
- Preserve fake-PQS guard fields.
- Delete the pass-232 blocked support-plan field cloud as already planned.

Report:

- final scoped line budget;
- explicit note that this pass used a manager-approved line-budget exception;
- generated support-plan status/counts/authority;
- validation and timings;
- remaining blockers.

Future cleanup:

- The manager is starting a separate old-flat-Cartesian-path audit to identify
  deletion candidates for future passes. Do not bundle that audit into pass 234.

-- repo-manager@macmini
