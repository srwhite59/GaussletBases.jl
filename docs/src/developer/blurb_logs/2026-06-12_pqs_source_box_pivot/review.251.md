# Pass 251 manager review

Decision: accepted.

Commit reviewed:

- pending commit: clarify independent H2 PQS input taxonomy

Scope reviewed:

- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl`
- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl`
- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

Findings:

- No blocking findings.
- The existing independent H2 PQS input remains the no-physics readiness input.
- The three new input files are minimal include/override variants for
  final-basis, H1, and H1-J diagnostics. None enables private RHF.
- The endpoint manifest now distinguishes readiness, final-basis, H1, and H1-J
  diagnostic roles and removes stale rows for absent H2 diagnostic input/test
  files.
- Fake-PQS remains a separate fake/source-backed WL/QW reproduction row.
- The deletion offset removes stale route-core pair/materialization mirror
  assertions while preserving compact terminal route/scaffold counts.

Validation accepted:

- Doer ran include/flag smoke for all three new input variants; it passed.
- Doer ran `git diff --check`; it passed.
- Doer ran direct trailing-whitespace search over touched input/docs/test files;
  it found no matches.
- Manager reviewed the input/manifest/test diffs and accepted the lightweight
  validation. No slow H2 route run was needed because passes 248-250 already
  validated the route seams.

Line budget:

- Scoped `src + test + bin`, counting new driver inputs: `24` added / `27`
  deleted, net `-3`.
- Manifest/input taxonomy plus deletion offset: `28` added / `29` deleted,
  net `-1`.

Remaining blocker / next:

- Independent H2 PQS stage taxonomy is clearer. Next substantive physics choice
  is private RHF diagnostic contract versus additional cleanup; do not promote
  public/export readiness yet.

-- repo-manager@macmini
