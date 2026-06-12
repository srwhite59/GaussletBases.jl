Accepted.

Pass 028 was doc-only and did the intended cleanup: it scoped
`:shell_realization_not_materialized` to bridge-level source-space readiness,
not the current status of the explicit PQS final-basis H1 route. It also
recorded in the retirement ledger that:

```text
source-box/final-basis H1 probe succeeds
_pqs_current_route_safe_term_matrices(...) remains oracle/debug vocabulary
```

Manager validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

No tests were run because only developer docs changed.

Next target:

Do a caller-driven audit of `_pqs_current_route_safe_term_matrices(...)` and
related current-route authority comparison helpers. Delete or shrink only the
callers that are clearly obsolete now that final-basis H1 exists. If deletion
is not safe, produce a precise keep/delete map.

-- repo-manager@macmini
