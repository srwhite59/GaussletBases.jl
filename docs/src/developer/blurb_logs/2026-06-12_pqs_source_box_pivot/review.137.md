Pass 137 manager review

Accepted.

Plain-language state:

- The existing complete-core/shell final-basis and H1 helpers consume an
  available `:pqs_multilayer_shell_source_plan` with `bundles`, `metrics`,
  core/shell support states and indices, `shell_final_coefficients`, and
  compatible support ordering.
- The Be2/PQS route owns source-box route structure, retained units/ranges,
  pair entries, center metadata, and, in the probe-enabled path, the parent
  axis bundle.
- It does not yet own a diatomic source realization object that honestly
  satisfies the existing complete-core/shell consumer shape.

Decision:

- Add a private blocked diatomic complete-core/shell source-plan payload next.
- The first pass should not fake or materialize a
  `:pqs_multilayer_shell_source_plan`.
- The payload should narrow the Be2 blocker from "missing producer" to the more
  precise missing source-realization contract/materializer while preserving all
  nonclaim flags.
- Reuse existing readiness summaries/helpers where possible to avoid duplicating
  field groups.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the next pass can move the source-plan blocker into its own
  route-owned payload instead of carrying it only as a Ham-readiness missing
  object.
- quarantined: one-center shellification/support-row semantics, WL adapters,
  final-basis/H1/H1-J materialization, RHF/SCF, public APIs, exports, artifacts,
  hfdmrg, and CR2 execution remain out of scope.
- not deleted because: Be2 fingerprint tests are still the compact validation
  for parent-axis-bundle availability and source-plan readiness.
- exact remaining caller/blocker: no diatomic source realization contract
  currently supplies the existing final-basis/H1 source-plan consumer shape.

-- repo-manager@macmini
