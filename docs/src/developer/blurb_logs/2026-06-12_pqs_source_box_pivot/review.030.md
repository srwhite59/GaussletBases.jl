Review 030

Accepted.

The one-line `@nospecialize` boundary is the right minimal fix for the
validation blocker. It keeps the existing route-driver test semantics, avoids
adding new tests or metadata, and directly addresses the stack overflow caused
by specializing on the giant route-report `NamedTuple`.

Validation:

- The previously failing route-driver file now passes directly: 1189 checks in
  about 13m43s.
- Load check passed.
- `git diff --check` passed.

Deletion/shrinkage:

- Nothing became obsolete from this fix.
- No new tests were added.
- Remaining report-shape pressure is real, but this pass intentionally stopped
  at the validation blocker.

Next:

- Rerun the slow nested harness with `TMPDIR=/private/tmp` to validate that the
  pass 029 edited integration file is no longer blocked by the earlier
  route-driver stack overflow.

-- repo-manager@macmini
