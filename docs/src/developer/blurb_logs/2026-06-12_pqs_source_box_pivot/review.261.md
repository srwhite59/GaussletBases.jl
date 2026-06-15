# Pass 261 manager review

Decision: accepted with a line-count exception.

Commit reviewed:

- pending commit: materialize independent H2 PQS supplement support partition

Scope reviewed:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`

Findings:

- No blocking findings.
- The new private support-partition payload is the right boundary before
  provider blocks. It carries atom-contact per-piece support tiles and
  shared-shell outer-minus-inner support tiles/row maps, rather than letting the
  next pass infer support from filled shared-shell source CPBs.
- The reported partition facts are internally consistent: support counts
  `(275, 578, 362)`, retained counts `(275, 98, 98)`, total support count
  `1215`, total tile count `55`, and zero duplicate/missing/outside parent
  rows.
- Provider blocks, mixed/GTO matrices, route-global matrices, residual MWG,
  combined density readiness, supplemented values, RHF, CR2/export, HamV6,
  public API, and fake/WL comparison paths remain blocked.
- Caveat: the current local-row mapping assumes tile rows belong to the unit
  support. That is fine for the validated H2 route, but before generalizing this
  helper or making it a provider-block error boundary, prefer blocking with an
  outside-row count over throwing on an unexpected row.

Validation accepted:

- Doer ran package load, focused support-partition smoke, and `git diff --check`.
  The focused smoke passed in about 164 seconds and did not run H1, H1-J, RHF,
  provider blocks, or supplemented values.
- Manager reran package load and `git diff --check`; both passed.

Line budget:

- Scoped `src + test + bin`: `435` added / `0` deleted, net `+435`.
- Exception accepted because this creates a route-authority payload required to
  prevent provider-block misuse. The next cleanup-capable pass should pay this
  down with mature deletion candidates rather than widening this implementation
  pass.

Remaining blocker / next:

- The next provider-block pass may consume
  `independent_h2_pqs_supplement_support_partition_payload`, but it should
  still keep matrices local/provider-level and avoid route-global supplemented
  values. A small hardening/paydown pass is also acceptable first.

-- repo-manager@macmini
