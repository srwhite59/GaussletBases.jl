# Pass 241 manager review - accepted with line-budget exception

Accepted under the explicitly recorded pass-241 line-budget exception.

What landed:

- Added the independent H2 PQS shared-shell realization payload.
- Materialized only the two shared-shell realization coefficient blocks.
- Kept `:atom_contact_core` as descriptor/identity-like only.
- Kept the complete source plan blocked at:

```text
:missing_independent_pqs_complete_core_shell_source_plan_assembly
```

Acceptance checks:

- The route remains fake-free:
  `fake_pqs/enabled = false`,
  `route/source_backed_fixed_source_oracle_used = false`.
- `physics/endpoint_ready` remains false.
- No final basis, H1, H1-J, RHF, supplements, CR2, export, public API, or
  fake-PQS/WL fixed-source coefficients were added.
- Shared-shell retained counts are `(98, 98)`.
- Reported realized overlap identity errors are small:
  `(1.2667660635192445e-14, 5.2966890053049083e-14)`.

Validation:

- Doer: package load passed; focused independent artifact/readiness check
  passed in `78.1625405s`; `git diff --check` passed.
- Manager: reviewed the diff, ran `git diff --check`, and ran package load;
  load passed in `0.645702292s`.

Line budget:

```text
src + test + bin: 280 added / 2 deleted, net +278
```

Exception reason:

- This was a narrow numerical seam. Forcing net-negative here would have
  required broader cleanup not assigned in the pass.
- Pass 239 had just removed the large projected-shell integration file, so the
  current exception is recorded rather than hidden.

Next:

- Return immediately to cleanup/shrink pressure. A subagent is auditing
  `test/nested/pqs_source_metadata_real_artifact_acceptance_support.jl` as the
  likely next carefully bounded cleanup target.

-- repo-manager@macmini
