# GaussletBases Cartesian/PQS Manager Running Log

This is a manager-level decision ledger. It does not replace doer responses,
baton reviews, or `state.md`. It records the strategic interpretation of
accepted work:

- what changed;
- why it was accepted;
- how it advances the real goals;
- what remains blocked;
- what future agents should not misunderstand.

Keep entries compact but substantive. Use the detailed pass response for
implementation detail.

## Initiation Snapshot

GaussletBases has several development lines, and they should not all be treated
as equally unsettled. The radial and angular basis work is comparatively stable
and well developed: those lines contain working basis constructions, atomic
operator machinery, radial/angular tests, and export paths that are not the main
source of current architectural uncertainty. The active management focus is the
Cartesian/high-order line. That line is where the repo is being pushed toward
consumer-ready Hamiltonian generation, CR2/HFDMRG/DMRG handoff, public driver
documentation, efficient operator construction, and a unified treatment of WL,
PQS, high-order slab/endcap/panel variants, MWG residual-Gaussian supplements,
and EGOI-type corrections.

The Cartesian code currently contains multiple generations of machinery. Older
ordinary/QW/WL, nested fixed-source, high-order doside, and residual-Gaussian
paths may mostly work and remain valuable as references, oracles, diagnostics,
and regression tests. But they are not the desired long-term architecture. The
desired direction is the newer layered Cartesian route spine: shellification
owns disjoint parent support; lowering chooses recipe-specific coordinate
product boxes; construction produces intermediate and final retained spaces;
pair planning starts from final retained units; pair-block materialization builds
operator blocks; and the driver/artifact layer exposes stable consumer
contracts. Future work should therefore avoid propagating old flat route/report
surfaces as the final design. Supplant them gradually with typed, modular,
route-authority-explicit construction while preserving old paths only where they
serve clear reference or migration purposes.

Recent cleanup confirms that this is not only architectural rhetoric: the
current Cartesian/PQS lane has already removed thousands of source/test/bin
lines while adding driver-owned H2 route checks, fake-PQS guard fields, and
supplement request/representation boundaries. Future work should preserve this
direction. New scaffolding should usually pay for itself by deleting stale
probes, duplicate report-key assertions, or obsolete route-shadow surfaces;
line count is not the only metric, but persistent line growth is a warning
sign.

For PQS in particular, the important distinction is between the intended
source-box algorithm and compatibility/oracle paths. PQS should be built from
filled source boxes, one-dimensional source transforms, boundary product-mode
retained rules, source-space operator blocks, and final shell realization where
needed. Shell-row projection, dense parent construction, and old fixed-source
WL/QW retained transforms may be useful for checks, but they must not be
mistaken for PQS route authority. The H2 463 fake-PQS episode is the current
warning example: a useful source-backed WL/QW reproduction was almost allowed to
propagate as a PQS physics endpoint. The manager log begins after that route has
been quarantined as fake-PQS. The next Cartesian/PQS work should recover an
independent H2 PQS route, keep fake/source-backed paths temporary, stage
supplements only after retained-transform authority is clear, and continue
reducing bloat while moving the Cartesian driver line toward public, efficient,
documented use.

## Long-Term Goals

LT1. Reliable Hamiltonian construction and driver workflow.
Build robust Hamiltonian creator functions and a visible, script-like driver
workflow for WL/QW, PQS, slab/endcap/panel and related high-order routes,
MWG residual-Gaussian supplements, EGOI-type corrections, and future CR2,
HFDMRG, or DMRG consumers.

LT2. Reduce repo bloat and complexity.
Maintain a long-term trend toward fewer source/test/bin lines, fewer stale
probes, fewer duplicated report aliases, less route-shadow vocabulary, and
cleaner JuliaStyle-compliant modular code.

LT3. Public Cartesian driver line.
Make the Cartesian driver line public-facing, documented, and scientifically
understandable, with algorithm pages, examples, supported/unsupported status,
and stable artifact meanings.

LT4. High efficiency.
Make Hamiltonian construction fast and memory-efficient, using source-first,
factorized, and one-dimensional contraction strategies where possible, while
separating diagnostic dense paths from production paths.

LT5. Scientific correctness and provenance.
Every route must expose its construction authority: WL/QW, independent PQS,
source-backed oracle, diagnostic fixture, or production method. Prevent
route-label drift.

LT6. Stable consumer contracts.
Driver artifacts should provide stable basis/operator/supplement/provenance
contracts for downstream consumers, especially CR2 and future public APIs.

LT7. Stratified testing.
Maintain a clear distinction among kernel tests, module-contract tests, driver
artifact tests, physics endpoint tests, slow/manual tests, and performance
benchmarks.

LT8. Extensible architecture.
Keep common concepts common: physical support plan, retained transform,
operator construction path, supplement request/representation, and final
artifact contract. Allow implementation divergence only where needed for speed
or memory.

## Anti-Goals

AG1. Do not let fake-PQS or source-backed oracle routes propagate as real PQS.

AG2. Do not add broad report-key clouds when a compact summary object would do.

AG3. Do not make public APIs out of private route scaffolds.

AG4. Do not compare gausslet-only endpoints to supplemented WL/QW references.

AG5. Do not use dense parent/operator construction as a production claim unless
it is explicitly labeled diagnostic.

AG6. Do not keep stale tests merely because they once caught something.

AG7. Do not let old flat Cartesian paths become the public architecture merely
because they still work. Preserve them only as references, oracles,
diagnostics, or migration scaffolds until typed route-owned replacements exist.

## Current Medium-Term Goals: PQS H2 Recovery And Supplement Staging

MT1. Quarantine fake-PQS H2 463.
Keep the current H2 463 source-backed WL/QW reproduction labeled as fake-PQS.
It is a temporary driver golden regression only.

MT2. Recover independent H2 PQS.
Create a separate H2 route whose retained transform is generated by independent
PQS source-box construction, not imported from the WL/QW fixed-source oracle.

MT3. Preserve common physical support vocabulary.
Use the common H2 physical support plan where appropriate: atom-contact core
plus shared shells, with retained-transform authority explicit.

MT4. Stage MWG/GTO supplements only after authority is clear.
The supplement request/representation boundary may remain, but provider-block
and supplemented-value work should not be promoted as real PQS until the base
retained transform is real PQS.

MT5. Keep cleanup pressure.
Any added scaffolding should be paid for by deleting stale probes, route-shadow
assertions, or duplicate metadata unless a safety-critical exception is
recorded.

MT6. Audit and classify old Cartesian flat paths.
Classify old ordinary/QW/WL, high-order doside, fixed-source, route-global, and
supplement paths as one of: reference/oracle, temporary migration scaffold,
active route-owned implementation, or deletion candidate. Do not promote a flat
path to public driver architecture without an explicit review.

## Medium-Term Goal Maintenance

Medium-term goals are current-lane goals, not permanent project goals. The
manager should review them after major milestones, strategic corrections, or
lane changes, and every 5 accepted Cartesian/PQS passes in a long loop.

When a pass completes, invalidates, or materially changes a medium-term goal,
the manager should update this section in the same commit as the running-log
entry or in a small follow-up commit. The entry should say that the
medium-term goals were updated and why.

Every 5 accepted Cartesian/PQS passes, add a medium-term goal checkpoint entry.
Classify each current MT goal as active, completed, blocked, stale, or needing
refinement. Update MT wording only when evidence warrants it.

Do not churn medium-term goals for every small pass. Prefer stable goals over
pass-local task lists. Ask the user before changing the broad direction of the
lane, but do not ask for permission merely to mark a completed goal as done or
to clarify a blocker already established by review.

## Entry Policy

Every accepted substantive Cartesian/PQS pass should get a compact but real
running-log entry, usually 100-250 words or a short bullet list using the
template below. The entry should not duplicate the doer response, but it must
preserve strategic interpretation: what changed, which LT/MT goals advanced,
what future agents should not misread, and the next blocker.

For purely mechanical passes with no strategic change, a 1-3 sentence tick is
allowed. The tick should explicitly say "no strategic change" and name the
current goal or guardrail it leaves unchanged.

Every 10-20 accepted passes, or after a major correction, add a strategic
compression entry. Summarize durable decisions, stale stories, false starts,
updated MT goals, and next lane direction. Do not prune prior running-log
entries by default; prefer append-only compression unless the user asks for
archival reorganization.

## Entry Template

```markdown
## Pass NNN - Short Title

Commit(s):
- `sha` - message

Summary:
- ...

Validation:
- Doer: ...
- Manager: ...

Goal advancement:
- LT#: ...
- MT#: ...

Medium-goal update:
- none / updated MT# because ...

Risk / guardrail:
- ...

Remaining blocker / next:
- ...

Line-count / complexity note:
- ...
```

## Pass 229 - Fake-PQS Quarantine

Commit(s):
- `208df0ad` - Relabel H2 source-backed route as fake PQS
- `26ecc737` - Quarantine fake PQS endpoint readiness
- `5f0bc97b` - Record pass 229 fake PQS audit

Summary:
- The H2 463 source-backed route was demoted from a "PQS physical endpoint" to
  a fake-PQS source-backed WL/QW reproduction.
- The route and test/input names now say fake-PQS.
- The artifact records fake-PQS provenance:
  `fake_pqs/enabled = true`,
  `fake_pqs/source = :source_backed_fixed_source_oracle`, and
  `fake_pqs/independent_pqs_transform = false`.
- The corrective pass set `physics/endpoint_ready = false` and added same-group
  fake/independent markers under `route/*` and `physics/*`, so downstream
  consumers should not classify this artifact as an independent PQS endpoint.

Validation:
- Manager ran the renamed fake-PQS endpoint test for `208df0ad`; after the
  final test shrink it passed `82/82` in about `141.8s`.
- Doer reported the corrective fake-PQS endpoint test passed `88/88` in about
  `3m22s`; load and `git diff --check` also passed.
- Manager reviewed the corrective diff and pushed without rerunning the long
  endpoint test.

Goal advancement:
- LT5: corrected construction provenance by marking the H2 463 route as
  fake-PQS instead of independent PQS.
- LT6: strengthened the artifact contract by placing fake/independent markers
  in stable groups that downstream consumers may read.
- MT1: quarantined the source-backed WL/QW reproduction as a temporary driver
  golden regression.
- MT2: protects independent PQS recovery by requiring future real PQS work to
  use a separate route.

Risk / guardrail:
- Do not treat the H2 463 fake-PQS route as an independent PQS physics
  endpoint.
- Do not build MWG/GTO provider blocks against this fake route unless the pass
  explicitly labels the work as fake-PQS supplement plumbing.
- A future independent H2 PQS target must have `fake_pqs/enabled = false` and a
  retained transform generated by PQS source-box construction.

Remaining blocker / next:
- Independent H2 PQS source-box construction is still missing.
- Some internal type/status names still contain `physical_gausslet` because the
  fake reproduction reuses that implementation path. This remains cleanup
  pressure after an independent PQS route exists.

Line-count / complexity note:
- The initial fake relabel stayed net-negative in `src + test` by shrinking
  stale test assertions.
- The corrective relabel was intentionally allowed to be small net-positive.
  Stable fake-PQS guard fields were safety-critical, so the normal line-reduction
  rule was temporarily ignored for that correction. The rule resumes after this
  documented exception.

## Pass 230 - Independent H2 PQS Recovery Audit

Commit(s):
- `bc61a46b` - Publish independent H2 PQS audit blurb
- this commit - Accept independent H2 PQS audit

Summary:
- Doer performed a no-edit audit of the independent H2 PQS source-box surfaces.
  The audit found useful lower-level PQS machinery, including raw product
  source-box facts and q=5 boundary retained rules, but no independent route
  authority for the H2 463 retained transform.
- The fake-PQS H2 463 path still imports its retained transform from the WL/QW
  fixed-source oracle. The retained counts `atom_contact_core = 251` and
  `shared_shell_2 = 114` are not currently explained by independent PQS
  retained rules.

Validation:
- Doer: read-only inspection only; no Julia commands or tests.
- Manager: confirmed response files match, checked the named live source
  surfaces exist, and verified the current H2 PQS driver input/test surface is
  the fake-PQS route rather than an independent route.

Goal advancement:
- LT5: clarified retained-transform authority and prevented WL/QW fixed-source
  data from being normalized as PQS route authority.
- MT2: refined the independent H2 PQS blocker to missing atom-contact-core and
  shared-shell retained rules.
- MT3: preserved the common H2 support vocabulary while separating support
  metadata from retained-transform authority.

Medium-goal update:
- none. MT2 remains active, and MT1 remains a maintained guardrail.

Risk / guardrail:
- Do not treat `251`, `98`, `114` as independent PQS retained counts merely
  because they are useful fake-PQS/WL reference values. Only `shared_shell_1`
  has a plausible current q=5 PQS boundary-count explanation.

Remaining blocker / next:
- Add a separate `fake_pqs=false` H2 PQS target/readiness surface that blocks on
  `:missing_independent_pqs_atom_contact_core_retained_rule` before source-plan,
  final-basis, H1, H1-J, RHF, supplement, or CR2 work.

Line-count / complexity note:
- No source/test/bin changes in this audit. The next implementation pass should
  pay for new target/readiness scaffolding by deleting stale scaffold pressure
  or stop with `ATTENTION.md`.

## Pass 231 - Independent H2 PQS Target Readiness

Commit(s):
- this commit - Accept independent H2 PQS target/readiness surface

Summary:
- Added a separate driver input and route kind for the independent H2 PQS target:
  `h2_pqs_q5_independent_source_box_r4.jl` and
  `:bond_aligned_diatomic_independent_pqs_source_box_core_shell`.
- The artifact/readiness path is explicitly fake-free:
  `fake_pqs/enabled = false`,
  `source_backed_fixed_source_oracle_used = false`, and
  `retained_transform_authority = :pqs_source_box_construction`.
- The target records H2 support metadata only, with support counts
  `(275, 578, 362)`, retained counts empty, and no expected final dimension.
  It blocks before source-plan/final-basis/H1 materialization.

Validation:
- Doer: package load passed; focused readiness artifact check passed in about
  `128.7s` after two scoping failures caught during development; `git diff
  --check` passed.
- Manager: reviewed the diff and input file, confirmed response files match,
  ran `git diff --check`, and ran package load in about `0.66s`. Manager did
  not rerun the long focused driver check.

Goal advancement:
- MT2: created a separate independent-PQS target so future work no longer has
  to mutate the fake-PQS reproduction.
- MT3: preserved the common H2 physical support vocabulary while keeping
  retained-transform authority blocked.
- LT5: kept route authority explicit and prevented fake-PQS retained counts
  from becoming PQS claims.

Medium-goal update:
- none.

Risk / guardrail:
- The support counts are still target/readiness metadata, not independently
  generated support-plan authority. Retained counts `(251, 98, 114)` remain
  fake-PQS/WL reference values unless independently generated by PQS.

Remaining blocker / next:
- Build or block a route-owned support/region plan for the independent route.
  Do not implement retained rules, final basis, H1, H1-J, RHF, supplements, or
  CR2 until support-plan authority is clear.

Line-count / complexity note:
- Pass 231 stayed barely net-negative in `src + test + bin`: tracked `src`
  added/deleted `154/166`, plus the 11-line new driver input, for net `-1`.
  This is acceptable but leaves little margin; pass 232 should continue deleting
  redundant readiness/report scaffold in the same surface.

## Pass 232 - Independent H2 PQS Support-Plan Blocker

Commit(s):
- this commit - Accept independent H2 PQS support-plan blocker

Summary:
- Added a compact `support_plan` fingerprint to the independent H2 PQS target.
  The pass intentionally did not generate support regions; it records
  `support_plan_status = :blocked_independent_pqs_support_region_plan` and
  `support_plan_blocker = :missing_independent_pqs_support_region_materializer`.
- The support counts `(275, 578, 362)` remain target constants with
  `support_counts_generated = false` and
  `support_counts_source = :target_constants_pending_support_region_materializer`.
- The pass removed stale source-plan-candidate aliases from the visible blocked
  target report rather than adding another field cloud.

Validation:
- Doer: focused independent-input artifact/readiness check passed in about
  `67.6s`; package load passed; `git diff --check` passed.
- Manager: reviewed the diff, ran `git diff --check`, and ran package load in
  about `0.66s`. Manager did not rerun the focused driver check.

Goal advancement:
- MT2: refined the next blocker from retained-rule work to the missing
  support-region materializer.
- MT3: preserved the common H2 support vocabulary while refusing to call target
  constants generated support.
- LT2: kept line pressure by deleting stale source-plan-candidate report aliases.

Medium-goal update:
- none.

Risk / guardrail:
- Do not treat `(275, 578, 362)` as route-generated until a materializer derives
  actual support regions from geometry/shellification/lowering.

Remaining blocker / next:
- Audit the support-region materializer seam. The next pass should identify
  whether existing shellification/lowering/PQS multilayer machinery can generate
  atom-contact core and shared-shell support rows without WL/QW coefficient
  matrices.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `23` added / `24` deleted, net `-1`.

## Pass 233 - Support-Region Materializer Audit

Commit(s):
- this commit - Record independent H2 PQS support-region audit

Summary:
- No-edit audit identified a geometry-owned path for the independent H2 PQS
  support plan. Existing terminal shellification can emit atom-local cores,
  midpoint slabs, and shared molecular shells without WL/QW coefficient
  matrices.
- The audit derived the target support counts geometrically:
  `atom_contact_core = 125 + 125 + 25 = 275`, inner shared shell
  `7*7*13 - 5*5*11 = 362`, and outer shared shell
  `9*9*15 - 7*7*13 = 578`.
- It also identified the ordering issue: raw terminal geometry may emit shared
  shells inside-out, but the H2 target vocabulary wants outside-in order
  `578, 362`.

Validation:
- Doer: read-only inspection only; no Julia commands or tests.
- Manager: reviewed the audit and confirmed the worktree stayed clean.

Goal advancement:
- MT2: located the next implementation seam for independent H2 PQS.
- MT3: made the common support vocabulary concrete without importing retained
  transforms.
- LT5: reinforced that fixed-source coefficients and fake-PQS adapters are not
  support authority.

Medium-goal update:
- none.

Risk / guardrail:
- Do not use `bond_aligned_diatomic_nested_fixed_source(...)`,
  fake-PQS `source.sequence.coefficient_matrix`, or WL/QW retained-transform
  data as independent support authority.

Remaining blocker / next:
- Implement a compact private support-region plan materializer that groups
  shellification-owned primitive support into `:atom_contact_core`,
  `:shared_shell_1`, and `:shared_shell_2`, while keeping retained transforms
  and physics blocked.

Line-count / complexity note:
- No source/test/bin changes in this audit. The implementation pass should
  delete the pass-232 hard-coded blocked support-plan field cloud once generated
  support-plan authority exists.

## Pass 234 - Independent H2 PQS Support-Region Plan

Commit(s):
- `433d684e` - Authorize H2 PQS support-plan line exception
- this commit - Accept independent H2 PQS support-region plan

Summary:
- Implemented the independent H2 PQS support-region materializer/fingerprint
  under the pass-234 line-budget exception. The support plan now derives its
  counts from `CartesianShellification.raw_terminal_geometry(...)`, not from
  fake-PQS/WL coefficient matrices.
- Generated support units match the intended route order and counts:
  `:atom_contact_core = 275`, `:shared_shell_1 = 578`, and
  `:shared_shell_2 = 362`. Shared shells are ordered outside-in.
- Artifact fields now report generated support authority:
  `support_plan_status = :available_independent_pqs_support_region_plan`,
  `support_plan_authority = :cartesian_shellification_route_geometry`, and
  `support_counts_generated = true`.

Validation:
- Doer: package load passed; focused independent-input artifact check passed in
  about `66.6s`; `git diff --check` passed.
- Manager: reviewed the diff, ran `git diff --check`, and ran package load in
  about `0.65s`. Manager did not rerun the long focused driver check.

Goal advancement:
- MT2: advanced independent H2 PQS from target constants to generated
  support-region authority.
- MT3: preserved common physical support vocabulary while keeping retained
  transforms blocked.
- LT5: kept construction authority explicit and fake-free.

Medium-goal update:
- none.

Risk / guardrail:
- This is support-region authority only. Source plan, retained rules, final
  basis, H1, H1-J, RHF, supplements, CR2, export, and public API remain blocked.

Remaining blocker / next:
- Pay down the pass-234 line-budget exception using the old-flat-path deletion
  audit. Then return to independent PQS retained-rule/source-plan work.

Line-count / complexity note:
- Pass 234 used the approved exception: `118` added / `10` deleted, net `+108`
  in `src + test + bin`.

## Medium-Term Goal Checkpoint - Passes 230-234

- MT1 Fake-PQS quarantine: active/maintained. Fake-PQS remains a golden
  regression only and was not mutated in passes 230-234.
- MT2 Independent H2 PQS recovery: active. The route now has a fake-free
  target and generated support-region authority; retained rules/source-plan
  remain blocked.
- MT3 Common physical support vocabulary: active. The H2 support vocabulary is
  now generated from shellification geometry, not just copied as constants.
- MT4 Supplement staging after authority: active/blocked. Supplement work
  remains intentionally deferred until retained-transform authority exists.
- MT5 Cleanup pressure: active with one recorded exception. Pass 234 incurred
  net `+108`; pass 235 should pay this down.
- MT6 Audit/classify old Cartesian flat paths: active. A read-only old-flat-path
  retirement audit identified route-shadow density-density fixture pressure as
  the best near-term deletion target.

## Pass 235 - Route-Shadow Density-Density Cleanup

Commit(s):
- this commit - Delete route-shadow density-density fixture pressure

Summary:
- Deleted the old PQS/PQS/product route-shaped density-density producer/consumer
  fixture path from `current_route_metadata_export.jl` and removed its slow
  integration-test call sites.
- Preserved compact lower-level math checks for boundary retained counts,
  density-normalized versus raw-weighted pair conventions, and nuclear
  charge/sign conventions.
- Independent H2 PQS support-plan code, fake-PQS guard fields, and the fake-PQS
  golden regression were untouched.

Validation:
- Doer: package load passed after precompile; focused integration check passed
  `3820/3820`, with test-reported time `6m04.6s` and Julia elapsed
  `845.4s`; `git diff --check` passed.
- Manager: reviewed the diff, reran deleted-name caller search with no matches,
  ran `git diff --check`, and ran package load in about `0.65s`. Manager did
  not rerun the long integration test.

Goal advancement:
- LT2: removed a large route-shadow fixture surface and old metadata pressure.
- MT5: paid down the pass-234 line-budget exception.
- MT6: used the old-flat-path audit to retire a classified deletion candidate.

Medium-goal update:
- none.

Risk / guardrail:
- Do not reintroduce the route-shaped density-density producer/consumer just to
  satisfy helper-name or private metadata tests. If a convention is needed,
  preserve it in compact lower-level tests.

Remaining blocker / next:
- Return to independent H2 PQS retained-rule/source-plan authority. The next
  step should audit `atom_contact_core`, `shared_shell_1`, and
  `shared_shell_2` retained-rule ownership before implementation.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `19` added / `1587` deleted, net `-1568`,
  more than paying down the pass-234 `+108` exception.

## Pass 236 - Independent H2 PQS Retained-Rule Audit

Commit(s):
- this commit - Record independent H2 PQS retained-rule audit

Summary:
- No-edit audit found that the fake/WL retained counts `(251, 98, 114)` are not
  independent PQS authority. Existing route-owned authority supports
  `:atom_contact_core => 275` via direct source modes and q=5 boundary product
  mode counts `98` for each shared shell.
- The resulting independent-PQS readiness target is retained counts
  `(275, 98, 98)` and final target dimension `471`, not the fake-PQS/WL `463`.
- The audit recommends a compact private retained-rule/readiness plan before
  any coefficient/source-plan materialization.

Validation:
- Doer: read-only inspection only; no Julia commands or tests.
- Manager: reviewed the audit and accepted `(275, 98, 98)` as the next
  readiness target.

Goal advancement:
- MT2: clarified the next independent-PQS seam and rejected fake/WL retained
  counts as route authority.
- MT3: kept common support vocabulary but allowed PQS retained counts to differ
  from WL/fake counts where independent rules differ.
- LT5: strengthened construction provenance by separating support authority,
  retained-rule authority, and fake/WL reference values.

Medium-goal update:
- none.

Risk / guardrail:
- Do not reintroduce `251` or `114` as independent-PQS expectations unless a
  real PQS retained rule is introduced and reviewed.

Remaining blocker / next:
- Add retained-rule/readiness metadata for `(275, 98, 98)` and final target
  dimension `471`. Do not materialize source coefficients, final basis, H1,
  H1-J, RHF, supplements, CR2, export, or public API.

Line-count / complexity note:
- No source/test/bin changes in this audit.

## Pass 237 - Independent H2 PQS Retained-Rule Readiness

Commit(s):
- this commit - Add independent H2 PQS retained-rule readiness

Summary:
- Added a compact private retained-rule readiness plan for the independent H2
  PQS route. The route now reports generated support counts `(275, 578, 362)`,
  retained counts `(275, 98, 98)`, and expected readiness dimension `471`.
- Per-unit authority is explicit: `:atom_contact_core` uses direct source modes,
  while both shared shells use q=5 PQS boundary COMX product-mode retained
  rules.
- The remaining endpoint/source-plan blocker is now the single blocker
  `:missing_independent_pqs_physical_source_plan_materializer`.

Validation:
- Doer: package load passed; focused independent-input artifact/readiness check
  passed in about `82.2s`; `git diff --check` passed.
- Manager: reviewed the diff, ran `git diff --check`, and ran package load in
  about `0.65s`. Manager did not rerun the focused artifact check.

Goal advancement:
- MT2: moved independent H2 PQS from retained-rule audit to retained-rule
  readiness metadata.
- MT3: preserved common support vocabulary while using independent PQS retained
  counts where WL/fake values are not route authority.
- LT2/MT5: deleted stale route-skeleton and axis-count scaffold tests while
  adding the new readiness plan.
- LT5: kept fake-PQS/WL fixed-source coefficients out of the independent route.

Medium-goal update:
- none.

Risk / guardrail:
- This is still metadata/readiness only. Do not treat `471` as a final basis or
  H1-ready endpoint; the source-plan materializer remains missing.

Remaining blocker / next:
- Audit the source-plan materializer seam before implementation, including the
  atom-contact-core direct source-mode representation and the shared-shell q=5
  source-box retained-rule descriptors.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `173` added / `248` deleted, net `-75`.

## Pass 238 - Source-Plan Materializer Audit

Commit(s):
- this commit - Record independent H2 PQS source-plan audit

Summary:
- No-edit audit identified the next source-plan seam: add a descriptor-only
  independent H2 PQS source-plan payload, not coefficient or final-basis
  materialization.
- The audit found the key hazard in the current assembly path: the generic
  source-plan candidate calls `bond_aligned_diatomic_nested_fixed_source(...)`,
  which is source-backed WL/QW data and must be gated away from the independent
  route.
- Recommended descriptor contents are compact per-unit source/retained-rule
  facts for `:atom_contact_core`, `:shared_shell_1`, and `:shared_shell_2`,
  while keeping coefficient materialization blocked.

Validation:
- Doer: read-only inspection only; no Julia commands or tests.
- Manager: reviewed the audit and confirmed the worktree was clean except for
  the tracked response file.

Goal advancement:
- MT2: defined the next implementation seam for independent H2 PQS.
- LT5: separated independent source-plan authority from the fake/source-backed
  WL/QW candidate path.
- MT5/MT6: identified same-surface cleanup candidates for the descriptor pass.

Medium-goal update:
- none.

Risk / guardrail:
- Do not let the independent route call `bond_aligned_diatomic_nested_fixed_source(...)`.
  That path belongs to fake-PQS/WL reproduction only.

Remaining blocker / next:
- Implement the descriptor-only source-plan payload, gate out the source-backed
  candidate, and leave the next numerical blocker explicit.

Line-count / complexity note:
- No source/test/bin changes in this audit.

## Pass 239 - Independent H2 PQS Source-Plan Descriptor

Commit(s):
- this commit - Add independent H2 PQS source-plan descriptor

Summary:
- Added a descriptor-only independent H2 PQS physical source-plan payload. The
  descriptor records support counts `(275, 578, 362)`, retained counts
  `(275, 98, 98)`, expected dimension `471`, and compact per-unit source/retained
  rule facts.
- Gated the source-backed WL/QW candidate away from the independent route, so
  `bond_aligned_diatomic_nested_fixed_source(...)` is no longer invoked for
  `:bond_aligned_diatomic_independent_pqs_source_box_core_shell`.
- Advanced the blocker to
  `:missing_independent_pqs_source_plan_numerical_materialization` while keeping
  `source_coefficients_materialized = false`.

Validation:
- Doer: package load passed; focused independent artifact/readiness check
  passed in about `128.9s`; `git diff --check` passed.
- Manager: reviewed the diff, ran `git diff --check`, and ran package load.
  Manager package load passed after precompile in about `58.7s`.

Goal advancement:
- MT2: created the first independent source-plan descriptor without fake/WL
  retained transforms.
- LT5: strengthened route authority by gating out the source-backed candidate.
- MT5/MT6/LT2: deleted the obsolete projected-q-shell integration file and its
  include, removing a large route-shadow test surface.

Medium-goal update:
- none.

Risk / guardrail:
- Descriptor availability is not numerical source-plan materialization. Do not
  treat the route as final-basis, H1, H1-J, RHF, supplement, or CR2 ready.

Remaining blocker / next:
- Audit the numerical materialization seam before implementation. The next
  decision is whether to materialize atom-contact direct identity facts,
  shared-shell source transforms, shared-shell realization coefficients, or a
  narrower blocked numerical payload.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `185` added / `5671` deleted, net `-5486`.

## Medium-Term Goal Checkpoint - Passes 235-239

- MT1 Fake-PQS quarantine: active/maintained. Fake-PQS remains separate and was
  not used as independent route authority.
- MT2 Independent H2 PQS recovery: active. The route now has support authority,
  retained-rule readiness, and a descriptor-only source-plan payload. Numerical
  materialization remains blocked.
- MT3 Common physical support vocabulary: active. Support vocabulary remains
  common, while independent PQS retained counts differ from fake/WL where route
  authority differs.
- MT4 Supplement staging after authority: active/blocked. Supplement work
  remains intentionally deferred.
- MT5 Cleanup pressure: active. Passes 235 and 239 removed large stale
  route-shadow/test surfaces and more than paid for new readiness descriptors.
- MT6 Audit/classify old Cartesian flat paths: active. Classified flat/source-
  backed paths are being gated, quarantined, or deleted rather than promoted.

## Pass 240 - Numerical Source-Plan Audit

Commit(s):
- this commit - Record independent H2 PQS numerical source-plan audit

Summary:
- No-edit audit recommended the next numerical object: shared-shell realization
  coefficients for `:shared_shell_1` and `:shared_shell_2`, not a full combined
  source plan or final basis.
- The audit keeps `:atom_contact_core` as a direct identity-like descriptor for
  now, avoiding a dense identity matrix before final-basis assembly.
- Recommended the next blocker
  `:missing_independent_pqs_shared_shell_realization_coefficients`, followed by
  complete-core/shell source-plan assembly after the shared-shell payload exists.

Validation:
- Doer: read-only inspection only; no Julia commands or tests.
- Manager: reviewed the audit and confirmed the worktree was clean except for
  the tracked response file.

Goal advancement:
- MT2: identified the next numerical seam for independent H2 PQS.
- LT5: kept old projected-shell machinery demoted to possible internal math
  adapter, not route authority.
- LT7: preserved staged validation rather than jumping to final basis/H1.

Medium-goal update:
- none.

Risk / guardrail:
- If `_nested_projected_q_shell_layer(...)` is reused, it must be fed by
  route-owned independent support/source boxes and must not import fake/WL
  fixed-source data.

Remaining blocker / next:
- Implement the shared-shell realization payload only. Keep complete source-plan
  assembly, final basis, H1, H1-J, RHF, supplements, CR2, export, and public API
  blocked.

Line-count / complexity note:
- No source/test/bin changes in this audit.

## Pass 241 - Shared-Shell Realization Payload

Commit(s):
- `4f3f2a66` - Authorize pass 241 line-budget exception
- this commit - Add independent H2 PQS shared-shell realization

Summary:
- Added the narrow shared-shell realization payload for independent H2 PQS.
  The route now materializes realization coefficients for `:shared_shell_1` and
  `:shared_shell_2`, each with retained count `98`.
- The atom-contact core remains descriptor/identity-like; no dense core identity
  matrix was materialized.
- The complete source plan remains blocked at
  `:missing_independent_pqs_complete_core_shell_source_plan_assembly`. Final
  basis, H1, H1-J, RHF, supplements, CR2, export, and public API remain blocked.

Validation:
- Doer: package load passed; focused independent artifact/readiness check
  passed in `78.1625405s`; `git diff --check` passed.
- Manager: reviewed the diff, ran `git diff --check`, and ran package load.
  Manager package load passed in `0.645702292s`.

Goal advancement:
- MT2: advanced independent H2 PQS from descriptor-only source-plan readiness
  to shared-shell numerical realization.
- LT5: kept fake/WL fixed-source data out of the independent route.
- LT7: preserved the staged validation line by keeping final basis and physics
  blocked.

Medium-goal update:
- none.

Risk / guardrail:
- The pass used old projected-shell machinery only as an internal mathematical
  adapter fed by independent route-owned support/source boxes. Do not promote
  that old machinery as route authority.

Remaining blocker / next:
- Assemble or audit the complete core/shell source plan only after the cleanup
  exception is paid down. The immediate next pass should return to targeted
  cleanup/shrink pressure.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `280` added / `2` deleted, net `+278`.
  This exception was explicitly granted and recorded in pass-241 exception
  files. A subagent was dispatched to prepare the next cleanup target.

## Pass 242 - Source Metadata Acceptance Shrink

Commit(s):
- this commit - Shrink PQS source metadata acceptance scaffold

Summary:
- Deleted the stale explicit missing-artifact export wrapper, path-writing
  branch, path/privacy assertions, fixed-column source-relation inventory
  checks, and no-op operator/Hamiltonian/postprocess rows from the Be2 PQS
  source metadata acceptance scaffold.
- Preserved the live private source-shell/source-mode metadata acceptance
  checks: TSV presence/headers, source-local axis labels, explicit parent-axis
  fields, unavailable ray/shell/radial labels, and no route/Hamiltonian/public
  consumption changes.

Validation:
- Doer: ran the focused opt-in source-metadata test without
  `BE2_PQS_Q5_ARTIFACT_DIR`; it skipped cleanly with one broken/skipped test in
  `0.0s`; ran `git diff --check`.
- Manager: reviewed the test-only diff and reran `git diff --check`.

Goal advancement:
- MT5/LT2: paid down part of the pass-241 line-budget exception by deleting
  obsolete acceptance scaffolding instead of broadening the independent H2 PQS
  implementation pass.
- LT7: reinforced the rule that stale scaffold tests should not be kept alive
  merely to prove obsolete branches still run.
- MT6: classified this metadata export wrapper/source-relation cloud as
  migration-era test pressure, while leaving the live metadata contract intact.

Medium-goal update:
- none.

Risk / guardrail:
- Do not delete the remaining source-shell/source-mode acceptance contract until
  there is a separate decision on which source count/category assertions are
  migration detail versus live acceptance.

Remaining blocker / next:
- Continue using the deletion-candidate audit queue for near-term cleanup, or
  return to independent H2 PQS complete source-plan assembly after enough
  cleanup paydown.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `11` added / `143` deleted, net `-132`.

## Pass 243 - Retire Legacy Source-Box Comparison Oracles

Commit(s):
- this commit - Delete uncalled PQS route-shadow comparison oracles

Summary:
- Deleted three uncalled safe-term/operator comparison wrapper families from
  `legacy_source_box_fixtures.jl`: contact-cap, outer-mismatch, and atom-box.
- Preserved the live retained-unit fixture builders used by
  `current_route_metadata_export.jl`:
  `_pqs_contact_cap_product_doside_unit`,
  `_pqs_outer_mismatch_product_doside_units`, and
  `_pqs_atom_box_support_dense_units`.
- Updated developer docs to say the old comparison wrappers are retired and
  that current authority is the retained-unit builders plus active route
  metadata, CPB, and source-box contracts.

Validation:
- Doer: package load passed and printed `load ok`; `git diff --check` passed;
  no non-slow metadata/report smoke was available without a Be2 artifact, so
  doer used the active-caller audit allowed by the blurb.
- Manager: reviewed source/doc diff, reran deleted-symbol search, and reran
  `git diff --check`.

Goal advancement:
- MT5/LT2: removed a large block of route-shadow oracle code instead of adding
  new scaffolding around independent H2 PQS.
- MT6/AG7: classified these old flat/source-box comparison paths as retired,
  while preserving the live migration builders that current-route metadata
  still uses.
- LT5: kept route authority clearer by removing wrappers that compared old
  product/support oracles without serving active route ownership.

Medium-goal update:
- none.

Risk / guardrail:
- Do not delete the retained-unit fixture builders until
  `current_route_metadata_export.jl` no longer calls them or a replacement
  inventory authority exists.

Remaining blocker / next:
- Additional cleanup candidates remain queued from sidecar audits. Good next
  candidates are report-stage low-order alias/test shrink or legacy-default
  low-order policy vocabulary shrink.

Line-count / complexity note:
- Total diff was `39` added / `850` deleted, net `-811`.
- Scoped source deletion was `0` added / `749` deleted in
  `legacy_source_box_fixtures.jl`.

## Pass 244 - Low-Order Report CRC Alias Shrink

Commit(s):
- this commit - Shrink low-order report CRC alias assertions

Summary:
- Shrank `cartesian_report_stage_low_order_policy_runtests.jl` by deleting
  duplicate CRC typed pair-operator count/source-path/materialization assertion
  clouds and the helper used only by those assertions.
- Kept compact coverage for atom-growth RouteCore final-unit/pair counts, the
  aggregate-subtree typed pair-operator blocker, and short CRC print
  compatibility substrings.

Validation:
- Doer: package load passed; the compact CRC print-line test passed; `git diff
  --check` passed.
- Doer attempted the broad report-stage test, stopped after it exceeded the 60s
  threshold and exposed stale/unrelated terminal-shellification exact-field
  failures, and used the blurb fallback.
- Manager: reviewed the diff, confirmed the broad report-stage file is not in
  default/integration runners, and reran `git diff --check`.

Goal advancement:
- MT5/LT2: deleted report-alias assertion bloat while preserving compact route
  summary/print coverage.
- LT7: reinforced that stale broad manual gates should not be forced as
  per-pass validation when the live contract has a smaller smoke.
- MT6: identified the same report-stage file's terminal-shellification
  exact-field assertions as a future cleanup candidate.

Medium-goal update:
- none.

Risk / guardrail:
- Do not treat `cartesian_report_stage_low_order_policy_runtests.jl` as a
  focused per-pass validation gate until its unrelated stale terminal/report
  assertions are retired or split.

Remaining blocker / next:
- Add the required medium-term checkpoint for passes 240-244. Then either use
  another audited cleanup target or resume independent H2 PQS source-plan
  assembly with a bundled deletion target if new source is needed.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `4` added / `121` deleted, net `-117`.

## Medium-Term Goal Checkpoint - Passes 240-244

- MT1 Fake-PQS quarantine: active/maintained. None of passes 240-244 touched
  fake-PQS endpoint authority or promoted fake data as independent PQS.
- MT2 Independent H2 PQS recovery: active. Pass 240 identified the next
  numerical seam; pass 241 added shared-shell realization coefficients. The
  route still lacks complete source-plan assembly, final basis, H1, H1-J, RHF,
  supplements, CR2, export, and public API.
- MT3 Common physical support vocabulary: active. The independent route still
  uses the common support vocabulary while keeping retained-transform authority
  separate from fake/WL source-backed data.
- MT4 Supplement staging after authority: active/blocked. Supplement work
  remains intentionally deferred until the independent retained-transform path
  is authoritative.
- MT5 Cleanup pressure: active. Pass 241 required an explicit `+278` exception;
  passes 242-244 then deleted or simplified stale source/test surfaces with
  scoped net changes of `-132`, `-811`, and `-117`.
- MT6 Audit/classify old Cartesian flat paths: active. The sidecar audits
  produced concrete deletion candidates and the loop has begun retiring old
  route-shadow comparison wrappers and report-alias test pressure.

## Pass 245 - Synthetic RouteCore Sidecar Blocker Test Shrink

Commit(s):
- this commit - Shrink synthetic RouteCore sidecar blocker tests

Summary:
- Deleted the synthetic missing-RouteCore-sidecar unit-pair helper/testset and
  collapsed pair-operator assertions that only preserved private
  metadata-only/blocker-count vocabulary.
- Preserved active unit-pair shape/order/family checks and pair-operator
  shape/path/transform matching checks.

Validation:
- Doer: `cartesian_unit_pairs_contract_runtests.jl` passed;
  `cartesian_pair_operator_plans_contract_runtests.jl` passed; `git diff
  --check` passed.
- Manager: reviewed the two focused diffs and accepted the doer validation.

Goal advancement:
- MT5/LT2: removed another 100 lines of stale contract-test pressure.
- LT7: kept focused module-contract tests on active behavior instead of
  synthetic obsolete blocker vocabularies.
- MT6: further retired old RouteCore sidecar transition scaffolding while
  preserving active module contracts.

Medium-goal update:
- none.

Risk / guardrail:
- Do not delete `_pair_ops_count(...)` while active path/final-block summary
  assertions still use it.

Remaining blocker / next:
- Cleanup candidates remain available, especially legacy-default low-order
  policy test vocabulary and broader flat report/status alias surfaces. The
  independent H2 PQS source-plan line can resume once manager chooses to spend
  the accumulated deletion budget.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `0` added / `100` deleted, net `-100`.

## Pass 246 - Legacy-Default Low-Order Vocabulary Shrink

Commit(s):
- this commit - Shrink legacy-default low-order policy vocabulary tests

Summary:
- Collapsed exact default `:legacy_diatomic_source*`, deferred legacy, and
  not-selected legacy-source-pair vocabulary assertions across shell, unit,
  transform, assembly, and report staged low-order policy tests.
- Kept compact legacy/reference smoke in the edited default blocks and left
  active atom-growth plus terminal-shellification sections untouched.

Validation:
- Doer: shell-stage low-order policy test passed; package load passed; `git
  diff --check` passed.
- Doer attempted unit/transform/assembly focused gates, but they exceeded 60s
  and failed on unrelated stale terminal-shellification exact-field assertions
  outside the edited blocks; fallback validation was used per blurb.
- Manager: reviewed five diffs, reran stale-symbol search for removed legacy
  vocabulary in the edited files, and reran `git diff --check`.

Goal advancement:
- MT5/LT2: deleted another legacy-vocabulary assertion layer from staged tests.
- LT7: continued moving tests from exact private status clouds toward compact
  live-contract checks.
- MT6/AG7: classified more old flat/legacy low-order route status vocabulary as
  reference/migration-only rather than public route architecture.

Medium-goal update:
- none.

Risk / guardrail:
- The broad unit/transform/assembly/report staged tests still include unrelated
  stale terminal-shellification exact-field assertions. Do not use them as
  per-pass validation gates until those blocks are retired or split.

Remaining blocker / next:
- The deletion budget is now substantial. The loop can either continue
  retiring stale terminal/report assertions or resume independent H2 PQS
  complete source-plan assembly with any new source offset by a targeted
  deletion.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `0` added / `176` deleted, net `-176`.

## Pass 247 - Independent H2 PQS Source-Plan Payload

Commit(s):
- this commit - Assemble independent H2 PQS source-plan payload

Summary:
- Assembled an available route-owned source-plan payload for the independent H2
  PQS target from atom-contact support rows plus the pass-241 shared-shell
  realization payloads.
- The artifact now reports source-plan availability with support counts
  `(275, 578, 362)`, retained counts `(275, 98, 98)`, and final dimension `471`.
- Final basis, H1, H1-J, RHF, supplements, CR2, export, and public API remain
  blocked.

Validation:
- Doer: package load passed; the focused independent H2 PQS driver/artifact
  check passed in about 78 seconds; `git diff --check` passed.
- Manager: reviewed the source/report/test diffs, reran `git diff --check`, and
  accepted the focused driver validation instead of duplicating the slow route
  run.

Goal advancement:
- MT2/LT5: advanced independent H2 PQS recovery from descriptor plus
  shared-shell payload to a route-owned source-plan payload with PQS authority.
- MT3/LT8: preserved the common physical support vocabulary while keeping the
  retained-transform authority explicitly independent of fake/WL source-backed
  data.
- MT5/LT2: offset new source-plan assembly with stale terminal-shellification
  exact-field assertion deletions.

Medium-goal update:
- none.

Risk / guardrail:
- Do not enable final basis or physics yet. The downstream final-basis helper
  still carries old physical-gausslet assumptions, including stale retained
  counts and source-backed metadata labels.

Remaining blocker / next:
- Review and correct the independent H2 PQS final-basis seam so it consumes the
  new `(275, 98, 98)` source plan without reintroducing fake/source-backed
  authority.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `267` added / `287` deleted, net `-20`.

## Pass 248 - Independent H2 PQS Final Basis

Commit(s):
- this commit - Materialize independent H2 PQS final basis

Summary:
- Updated the independent H2 PQS final-basis seam so retained counts come from
  the route-owned source-plan summary instead of the stale fake/WL-era
  `(251, 98, 114)` tuple.
- The focused route now materializes a 471-dimensional final basis with retained
  counts `(275, 98, 98)`, full rank, and final overlap identity error about
  `1.3e-13`.
- Final-basis overlap/rank fields are exposed as diagnostics only. H1, H1-J,
  RHF, supplements, CR2, export, and public API remain off.

Validation:
- Doer: package load passed; focused independent H2 PQS driver/artifact check
  with `run_final_basis=true` passed in about 78 seconds; `git diff --check`
  passed.
- Manager: reviewed the source/report/test diffs, reran `git diff --check`, and
  accepted the focused driver validation without duplicating the slow run.

Goal advancement:
- MT2/LT5: advanced the independent H2 PQS route from source-plan availability
  to a route-owned final basis without reusing fake/WL authority.
- LT6/LT7: exposed compact diagnostic final-basis contract fields while keeping
  physics and export surfaces blocked.
- MT5/LT2: offset the source/report additions by deleting stale
  terminal-shellification exact-field test mirrors.

Medium-goal update:
- none.

Risk / guardrail:
- Do not treat dense support-overlap/final cleanup as production performance
  evidence. It is first final-basis diagnostic materialization.
- Do not attach H1-J/RHF/supplements before H1 and density conventions are
  reviewed for this independent basis.

Remaining blocker / next:
- Begin the independent H2 PQS H1 one-body seam: kinetic and by-center
  electron-nuclear construction for the 471-dimensional final basis. Keep
  H1-J, RHF, supplements, CR2, export, and public API blocked.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `87` added / `100` deleted, net `-13`.

## Pass 249 - Independent H2 PQS H1

Commit(s):
- this commit - Expose independent H2 PQS H1 artifact facts

Summary:
- The existing physical-gausslet H1 path worked for the independent H2 PQS final
  basis when requested; this pass exposed its status facts in the artifact.
- The focused route reports a finite, symmetric 471-dimensional H1 with lowest
  one-electron energy `-0.7946037173365885`.
- H1-J, RHF, supplements, CR2, export, and public API remain off. The endpoint
  blocker is now the expected H1-J/density seam.

Validation:
- Doer: package load passed; focused independent H2 PQS H1 driver/artifact
  check passed in about 80 seconds; `git diff --check` passed.
- Manager: reviewed the focused diff and accepted doer validation without
  duplicating the slow route run.

Goal advancement:
- MT2/LT5: advanced independent H2 PQS from final-basis readiness to one-body
  H1 readiness without fake/source-backed authority.
- LT6: made H1 status, energy, finiteness, symmetry, and support/final one-body
  statuses visible to the artifact consumer.
- MT5/LT2: kept the pass net-negative by deleting stale selected-terminal
  sidecar count/status mirror assertions.

Medium-goal update:
- See the passes 245-249 checkpoint below.

Risk / guardrail:
- Do not jump directly from H1 to RHF. The next seam must review H1-J/density
  convention, because earlier PQS density attempts exposed bad weight/gauge
  behavior.

Remaining blocker / next:
- H1-J/density convention review for the independent H2 PQS basis. Keep RHF,
  supplements, CR2, export, and public API blocked until that seam is coherent.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `17` added / `26` deleted, net `-9`.

## Medium-Term Goal Checkpoint - Passes 245-249

- MT1 Fake-PQS quarantine: active/maintained. None of these passes promoted the
  fake-PQS source-backed reproduction; pass outputs continue to assert
  `fake_pqs_enabled = false` and
  `source_backed_fixed_source_oracle_used = false` for independent H2 PQS.
- MT2 Independent H2 PQS recovery: active and substantially advanced. Passes
  247-249 moved the route from source-plan blocker to available source plan,
  available 471-dimensional final basis, and materialized one-body H1.
  H1-J/density, RHF, supplements, CR2, export, and public API remain blocked.
- MT3 Common physical support vocabulary: active. The independent route keeps
  the shared physical support vocabulary while using independent PQS authority
  and retained counts `(275, 98, 98)`.
- MT4 Supplement staging after authority: active/blocked. Base gausslet
  authority is clearer now, but supplement work remains deferred until H1-J and
  density conventions are reviewed.
- MT5 Cleanup pressure: active. Passes 245-249 were all scoped net-negative in
  `src + test + bin`, even while source-plan/final-basis/H1 capability was
  added.
- MT6 Audit/classify old Cartesian flat paths: active. These passes continued
  retiring stale staged-test mirror assertions and old flat vocabulary while
  preserving compact active route smoke.

## Pass 250 - Independent H2 PQS H1-J Diagnostic

Commit(s):
- this commit - Validate independent H2 PQS H1-J diagnostic

Summary:
- The existing H1-J/density diagnostic path materialized for the independent H2
  PQS basis without source edits.
- The focused route reports density gauge
  `:pre_final_localized_positive_weight`, raw pair convention `:raw_numerator`,
  positive support weights, finite/symmetric `(471, 471)` pre-final pair matrix,
  and H1-J self-Coulomb `0.4569117646737236`.
- This remains diagnostic only. RHF/private RHF, supplements, CR2, export, and
  public API remain off.

Validation:
- Doer: package load passed; focused independent H2 PQS H1-J driver/artifact
  check passed in about 79 seconds; `git diff --check` passed.
- Manager: reviewed the deletion-only diff and accepted doer validation without
  duplicating the slow route run.

Goal advancement:
- MT2/LT5: advanced independent H2 PQS from H1 readiness to a coherent H1-J
  density diagnostic while preserving fake/source-backed guardrails.
- LT6: confirmed the artifact already exposes the key density/H1-J convention
  facts needed for consumer review.
- MT5/LT2: deleted another stale selected-terminal-sidecar mirror block.

Medium-goal update:
- none.

Risk / guardrail:
- H1-J is not solver readiness. Do not promote RHF/export/public readiness
  until a private RHF/solver-contract pass explicitly reviews the density/Fock
  convention.

Remaining blocker / next:
- Choose the next lane deliberately: either private RHF diagnostic contract for
  this independent H2 PQS basis, or cleanup of driver input naming so
  readiness/final-basis/H1/H1-J variants are not confused.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `0` added / `12` deleted, net `-12`.

## Pass 251 - Independent H2 PQS Input Taxonomy

Commit(s):
- this commit - Clarify independent H2 PQS input taxonomy

Summary:
- Added three tiny include/override driver inputs for independent H2 PQS
  final-basis, H1, and H1-J diagnostic stages while leaving the existing input
  as the no-physics readiness input.
- Updated the endpoint manifest to distinguish readiness, final-basis, H1, and
  H1-J diagnostic roles, and removed stale rows for absent H2 diagnostic input
  and test files.
- Fake-PQS remains a separate fake/source-backed WL/QW reproduction row.

Validation:
- Doer: include/flag smoke passed for the three new input variants; `git diff
  --check` passed; direct trailing-whitespace search found no matches.
- Manager: reviewed the input, manifest, and cleanup diffs. No slow H2 route run
  was needed because passes 248-250 already validated the route seams.

Goal advancement:
- LT6/LT7: made the driver input taxonomy match the demonstrated artifact
  stages and reduced future confusion between readiness, final basis, H1, and
  H1-J.
- LT5/MT1: preserved the fake-PQS quarantine as a separate manifest row.
- MT5/LT2: offset the small input additions by deleting stale route-core
  pair/materialization mirror assertions.

Medium-goal update:
- none.

Risk / guardrail:
- These variants are diagnostic inputs, not public endpoints. Do not use them
  to claim solver/export/public readiness.

Remaining blocker / next:
- Choose deliberately between a private RHF diagnostic contract pass and another
  cleanup/classification pass. Supplements, CR2, export, and public API remain
  blocked.

Line-count / complexity note:
- Scoped `src + test + bin`, counting new driver inputs, was `24` added / `27`
  deleted, net `-3`.
- Manifest/input taxonomy plus deletion offset was `28` added / `29` deleted,
  net `-1`.

## Pass 252 - Independent H2 PQS Private RHF Diagnostic

Commit(s):
- this commit - Materialize independent H2 PQS private RHF diagnostic

Summary:
- Repaired the existing private RHF diagnostic seam so it accepts the
  independent H2 PQS route/role family, not only the old physical/fake-PQS role.
- The focused route now materializes private RHF on the independent 471-D H2 PQS
  H1-J basis with `fake_pqs_enabled=false`,
  `source_backed_fixed_source_oracle_used=false`, and
  `retained_transform_authority=:pqs_source_box_construction`.
- RHF converged in 8 iterations with total energy `-1.1589735957658853`,
  one-body energy `-1.5609752182694798`, two-body energy
  `0.4020016225035945`, density trace `1.9999999999999993`, idempotency
  residual `3.1e-17`, and commutator residual `4.0e-9`.

Validation:
- Doer: package load passed; focused independent H2 PQS private RHF
  driver/artifact probe passed in about 84 seconds; `git diff --check` passed.
- Manager: reviewed the role-gate/RHF-helper diff, reran package load, and
  accepted doer validation without duplicating the slow RHF route run.

Goal advancement:
- MT2/LT5: advanced independent H2 PQS from H1-J diagnostic to private RHF
  diagnostic while preserving route authority and fake-PQS quarantine.
- LT6/LT7: confirmed the driver artifact exposes the RHF input contract,
  execution status, convergence, and residual diagnostics needed for review.
- MT5/LT2: offset the RHF source changes by deleting stale terminal RouteCore
  mirror assertions.

Medium-goal update:
- none.

Risk / guardrail:
- Private RHF is still diagnostic-only. Do not promote export/public solver
  readiness, CR2 readiness, or supplement work from this result.
- Do not compare this gausslet-only independent PQS value to supplemented WL/QW
  references.

Remaining blocker / next:
- Decide whether to add a tiny explicit private-RHF input/manifest row, perform
  a strategic checkpoint after pass 254, or resume old-flat cleanup before any
  supplement/provider-block work.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `29` added / `31` deleted, net `-2`.

## Pass 253 - Independent H2 PQS Private RHF Input Taxonomy

Commit(s):
- this commit - Add independent H2 PQS private RHF input

Summary:
- Added an explicit include/override input for the independent H2 PQS private
  RHF diagnostic:
  `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl`.
- The input enables final basis, H1, H1-J, and private RHF with
  `private_rhf_electron_count = 2`, while keeping
  `physics_endpoint_ready = false` and blocker
  `:private_rhf_diagnostic_not_public_solver_contract`.
- Updated the endpoint manifest and compact independent-H2 artifact-role
  classifier so the new role keeps fake-PQS/source-backed authority disabled.

Validation:
- Doer: include/flag smoke passed; package load passed; `git diff --check` and
  trailing-whitespace search passed.
- Manager: reran the include/flag smoke and package load. Both passed. No slow
  H2 RHF route run was needed because pass 252 already validated the route.

Goal advancement:
- LT6/LT7: made the private RHF diagnostic stage explicit in driver inputs and
  the endpoint manifest.
- MT2/LT5: preserved independent-PQS authority and fake-PQS quarantine for the
  new RHF role.
- MT5/LT2: offset the taxonomy addition by deleting stale terminal RouteCore
  status mirror assertions.

Medium-goal update:
- none.

Risk / guardrail:
- The new input is a private diagnostic input, not a public solver endpoint.

Remaining blocker / next:
- Pass 254 should either do a small cleanup/classification pass or a checkpoint
  style pass. After accepting pass 254, add the required medium-term goal
  checkpoint for passes 250-254.

Line-count / complexity note:
- Scoped `src + test + bin`, counting the new input, was `11` added / `17`
  deleted, net `-6`.

## Pass 254 - Terminal Report-Stage Alias Cleanup

Commit(s):
- this commit - Shrink terminal report-stage low-order aliases

Summary:
- Deleted stale terminal-shellification report-stage alias assertions from
  `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`.
- Removed exact object-identity mirrors, flat
  `low_order_terminal_shellification_*` alias mirrors, and deferred terminal
  `low_order_route_core_*` status mirrors.
- Preserved compact coverage for terminal selection, summary availability,
  region/unit counts, pair/assembly deferral statuses, materialization blocker,
  no operator blocks materialized, central counts, and print-line content.

Validation:
- Doer: Julia parse smoke passed; `git diff --check` passed.
- Manager: reviewed the deletion-only diff and reran the Julia parse smoke. It
  passed. The broad report-stage low-order test was intentionally not used as a
  per-pass gate for stale-alias cleanup.

Goal advancement:
- MT5/LT2: removed another 36 lines of stale flat report alias pressure.
- MT6/AG7: classified more old terminal-shellification flat report vocabulary
  as compatibility/test scaffolding rather than public architecture.

Medium-goal update:
- Added the required checkpoint below for passes 250-254.

Risk / guardrail:
- Do not replace these deleted aliases with new flat report-key assertions.

Remaining blocker / next:
- Continue with either another classified cleanup candidate or a deliberate
  supplement-staging decision. Public/export/CR2 readiness remains blocked.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `0` added / `36` deleted, net `-36`.

## Medium-Term Goal Checkpoint - Passes 250-254

- MT1 Fake-PQS quarantine: active/maintained. None of passes 250-254 promoted
  fake-PQS/source-backed WL/QW data as independent PQS evidence.
- MT2 Independent H2 PQS recovery: active and substantially advanced. These
  passes moved the independent route from H1-J/density diagnostic to
  materialized private RHF diagnostic, then added explicit stage inputs through
  private RHF. The route remains private/diagnostic-only.
- MT3 Common physical support vocabulary: active. The independent route keeps
  the H2 support vocabulary while retaining independent PQS counts and route
  authority.
- MT4 Supplement staging after authority: active/blocked. The base gausslet
  independent route now has private RHF diagnostics, but supplement/provider
  work is still blocked until manager deliberately opens that lane.
- MT5 Cleanup pressure: active. Passes 250-254 were all scoped net-negative:
  `-12`, `-3`, `-2`, `-6`, and `-36`.
- MT6 Audit/classify old Cartesian flat paths: active. These passes continued
  retiring stale selected-terminal, RouteCore, and report-stage flat alias
  assertions while preserving compact active smoke checks.

## Pass 255 - Atom-Growth Report-Stage RouteCore Mirror Cleanup

Commit(s):
- this commit - Shrink atom-growth report-stage RouteCore mirrors

Summary:
- Deleted stale atom-growth report-stage RouteCore mirror assertions from
  `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`.
- Removed nested preflight/plan object field mirrors and one flat
  report-to-summary family-count mirror.
- Preserved compact checks for atom-growth route selection, materialization
  deferral, pair inventory/count, RouteCore inventory/status, pair-operator
  readiness, preflight/plan status, and LW complete-shell enumeration.

Validation:
- Doer: Julia parse smoke passed; `git diff --check` passed.
- Manager: reviewed the deletion-only diff and reran the Julia parse smoke. It
  passed. The broad report-stage test was intentionally not used as a per-pass
  gate.

Goal advancement:
- MT5/LT2: removed another 13 lines of stale report-stage mirror pressure.
- MT6/AG7: classified nested RouteCore preflight/plan object-shape assertions
  as old flat scaffolding, not public route authority.

Medium-goal update:
- none.

Risk / guardrail:
- Do not reintroduce nested preflight/plan object-shape assertions as a
  substitute for compact readiness facts.

Remaining blocker / next:
- Continue old-flat cleanup from mature candidates or deliberately decide when
  to open independent-PQS supplement staging.

Line-count / complexity note:
- Scoped `src + test + bin` diff was `0` added / `13` deleted, net `-13`.

## Pass 256 - Independent H2 PQS Supplement Staging Audit

Commit(s):
- this commit - Record independent H2 PQS supplement staging audit

Summary:
- No source/test/bin implementation edits. The pass audited whether existing
  supplement request, representation, and preflight helpers can attach to the
  independent H2 PQS target.
- Finding: the helpers are target-driven enough to attach mechanically to the
  independent route, but their function/type/status vocabulary is still old
  physical-gausslet shaped and not yet a clean independent-PQS supplement seam.
- The first safe next step is an independent supplement-preflight input/artifact
  role only, with provider blocks and supplemented values still blocked.

Validation:
- Doer: `git diff --check` passed; no Julia run because this was no-edit
  inspection.
- Manager: reviewed the audit response and accepted scoped line impact `0`.

Goal advancement:
- MT4/LT5: opened supplement staging deliberately after independent retained
  authority was established, without starting provider-block work.
- LT6: clarified the artifact facts expected for independent supplement
  preflight: fake disabled, source-backed oracle disabled, retained authority
  PQS, counts `(275, 578, 362)` and `(275, 98, 98)`, final dimension `471`,
  H/cc-pVTZ lmax-1 supplement representation, and provider-block blockers.

Medium-goal update:
- none.

Risk / guardrail:
- Do not copy fake-PQS retained counts `(251, 98, 114)`, final dimension `463`,
  source-backed labels, or WL/QW scalar comparisons into independent supplement
  work.

Remaining blocker / next:
- Add only the independent supplement-preflight input/artifact role. Provider
  blocks, mixed matrices, residual MWG representation, combined density
  readiness, supplemented values, CR2/export, and public API remain blocked.

Line-count / complexity note:
- Scoped `src + test + bin` impact was `0`.

## Pass 257 - Independent H2 PQS Supplement Preflight Input

Commit(s):
- this commit - Add independent H2 PQS supplement preflight input

Summary:
- Added
  `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl`
  as a tiny include/override input based on the H1-J diagnostic input.
- The input sets `supplement_policy = :mwg_residual_gto`, marks comparison as
  preflight-only, keeps private RHF off, and leaves endpoint readiness blocked
  by `:missing_provider_gto_supplement_blocks`.
- Added the new artifact role to the independent-H2 classifier and the endpoint
  manifest. No provider blocks or supplemented values were implemented.

Validation:
- Doer: include/flag smoke, package load, report-stage parse smoke, classifier
  smoke, and `git diff --check` passed.
- Manager: reran include/flag smoke and classifier smoke. Both passed.

Goal advancement:
- MT4/LT5: created the explicit independent supplement-preflight entry point
  while preserving independent-PQS route authority.
- LT6/LT7: made the preflight role visible in the driver input taxonomy and
  manifest.
- MT5/LT2: offset the new input/classifier lines by deleting stale report-level
  RouteCore mirror assertions.

Medium-goal update:
- none.

Risk / guardrail:
- This is preflight only. Do not interpret it as provider-block availability,
  residual-MWG readiness, supplemented values, or export/public readiness.

Remaining blocker / next:
- The next supplement pass should either run/verify the independent preflight
  artifact facts or start a narrow provider-block design audit. Actual
  materialization remains blocked by `:missing_provider_gto_supplement_blocks`.

Line-count / complexity note:
- Scoped `src + test + bin`, counting the new input, was `13` added / `24`
  deleted, net `-11`.

## Pass 258 - Independent H2 PQS Supplement Preflight Verification

Commit(s):
- this commit - Record independent H2 PQS supplement preflight verification

Summary:
- Ran the new independent H2 PQS supplement-preflight input through the driver
  and inspected the JLD2 artifact.
- The artifact reported fake-PQS disabled, source-backed oracle disabled,
  retained authority PQS, support counts `(275, 578, 362)`, retained counts
  `(275, 98, 98)`, final dimension `471`, and H/cc-pVTZ lmax-1 supplement
  representation with 18 orbitals.
- Preflight correctly remains blocked on
  `:missing_provider_gto_supplement_blocks` with mixed/GTO/MWG/density missing
  facts. No source fix was needed.

Validation:
- Doer: focused route/artifact probe passed on rerun in about 78 seconds;
  `git diff --check` passed.
- Manager: accepted doer validation without duplicating the slow focused route
  run.

Goal advancement:
- MT4/LT5: verified the independent supplement-preflight artifact without
  provider-block materialization or fake-PQS evidence.
- LT6: confirmed the preflight artifact exposes the expected consumer-facing
  blockers and counts.

Medium-goal update:
- none.

Risk / guardrail:
- Preflight is not supplement readiness. Provider blocks, residual MWG
  representation, combined density readiness, supplemented values, CR2/export,
  and public API remain blocked.

Remaining blocker / next:
- Open only a narrow provider-block seam or design audit next. Do not jump to
  supplemented energies.

Line-count / complexity note:
- Scoped `src + test + bin` tracked impact was `0`.

## Pass 259 - Independent H2 PQS Provider-Block Seam Audit

Commit(s):
- this commit - Record independent H2 PQS provider-block seam audit

Summary:
- No implementation edits. Audited the first route-owned provider-block seam for
  independent H2 PQS MWG/GTO supplement staging.
- The first object should be private and route-owned, carrying compact
  provider-block status/count/fingerprint summaries, not route-global
  supplemented operators.
- Key hazard: CPB providers operate on rectangular CPBs, but independent shared
  shells are outer-minus-inner support. The implementation must carry explicit
  support CPB tiling/row ownership or equivalent row maps before provider calls.

Validation:
- Doer: `git diff --check` passed; no Julia run because this was no-edit audit.
- Manager: reviewed the audit and accepted scoped line impact `0`.

Goal advancement:
- MT4/LT4/LT5: identified the performance/authority seam for supplement
  provider blocks without building route-global matrices or supplemented
  values.
- LT8: clarified that support partition/row ownership is the common object the
  next implementation needs before local provider kernels can be reused.

Medium-goal update:
- Added the required checkpoint below for passes 255-259.

Risk / guardrail:
- Do not call provider blocks blindly on filled shared-shell source CPBs.
- Keep matrices local/provider-level in the first implementation; route-global
  matrices, residual MWG, combined density readiness, supplemented values,
  CR2/export, and public API remain blocked.

Remaining blocker / next:
- Implement only the private independent-H2 provider-block payload and compact
  support-partition summary, or do a narrower support-tiling audit if that seam
  is not yet ready.

Line-count / complexity note:
- Scoped `src + test + bin` tracked impact was `0`.

## Medium-Term Goal Checkpoint - Passes 255-259

- MT1 Fake-PQS quarantine: active/maintained. None of passes 255-259 used
  fake-PQS or source-backed WL/QW data as independent-PQS evidence.
- MT2 Independent H2 PQS recovery: active. The base independent gausslet route
  remains available through private RHF diagnostics; these passes did not alter
  the base route.
- MT3 Common physical support vocabulary: active. Supplement preflight now uses
  the independent support/retained counts `(275, 578, 362)` and `(275, 98, 98)`.
- MT4 Supplement staging after authority: active and advanced. The lane moved
  from audit to explicit preflight input/artifact verification and then to a
  provider-block seam audit. Provider-block implementation remains next, but
  supplemented values remain blocked.
- MT5 Cleanup pressure: active. Passes 255 and 257 were scoped net-negative;
  passes 256, 258, and 259 were no-edit/audit or verification records with
  scoped impact `0`.
- MT6 Audit/classify old Cartesian flat paths: active. Cleanup continued for
  old report-stage RouteCore mirrors, and the supplement audit classified
  fake-PQS supplement preflight as schema/history only.

## Pass 260 - Independent H2 PQS Supplement Support-Tiling Audit

Commit(s):
- this commit - Record independent H2 PQS supplement support-tiling audit

Summary:
- No implementation edits. The pass audited where support rows, source boxes,
  parent-row indices, retained ranges, and support-unit facts are available for
  the future independent-H2 supplement provider-block payload.
- Finding: a support-partition summary is ready to implement, but provider
  calls must wait until atom-contact per-piece descriptors and shared-shell
  outer-minus-inner support tiles are carried as route-owned row maps.
- The next implementation should expose a private support-partition payload and
  compact row/tile/coverage fingerprints only. No provider blocks or
  route-global matrices should be materialized yet.

Validation:
- Doer: `git diff --check` passed; no Julia run because this was no-edit audit.
- Manager: reviewed the audit and accepted scoped line impact `0`.

Goal advancement:
- MT4/LT5: tightened the supplement-provider boundary without advancing into
  provider blocks or supplemented values.
- LT8: identified the support-partition object needed to keep common physical
  support vocabulary separate from provider implementation details.
- LT4: preserved the intended local/provider-level path instead of building
  diagnostic route-global matrices.

Medium-goal update:
- none.

Risk / guardrail:
- Do not call CPB provider helpers on filled shared-shell source CPBs as if
  they were support. Shared-shell support is `outer_box \ inner_box` and must
  be tiled or row-mapped explicitly.

Remaining blocker / next:
- Implement only the private support-partition payload/summary. If the existing
  complete-shell boundary-strata helper does not fit a shared shell, add only a
  narrow rectangular-difference tiler and validate it against stored support
  rows.

Line-count / complexity note:
- Scoped `src + test + bin` tracked impact was `0`.

## Pass 261 - Independent H2 PQS Supplement Support Partition

Commit(s):
- this commit - Materialize independent H2 PQS support partition

Summary:
- Implemented the private support-partition payload needed before independent
  H2 PQS supplement provider blocks.
- The payload exposes atom-contact per-piece tiles and shared-shell
  outer-minus-inner support tiles/row maps. It reports support counts
  `(275, 578, 362)`, retained counts `(275, 98, 98)`, final support count
  `1215`, total tile count `55`, and zero duplicate/missing/outside rows.
- No provider blocks, mixed/GTO matrices, residual MWG, route-global matrices,
  supplemented values, CR2/export, HamV6, public API, or fake/WL comparison
  paths were materialized.

Validation:
- Doer: package load, focused support-partition smoke, and `git diff --check`
  passed. The focused smoke took about 164 seconds and did not run H1, H1-J,
  RHF, provider blocks, or supplemented values.
- Manager: reran package load and `git diff --check`; both passed.

Goal advancement:
- MT4/LT5: created the route-authority boundary required before provider-block
  work can proceed without abusing filled shared-shell CPBs as support.
- LT8: made the common support vocabulary concrete as row-owned support tiles
  and unit partitions.
- LT4: preserved the local/provider-first path by keeping route-global matrices
  unmaterialized.

Medium-goal update:
- none.

Risk / guardrail:
- This is still pre-provider infrastructure. Provider blocks must consume this
  support partition and keep matrices local/provider-level; supplemented values
  remain blocked.

Remaining blocker / next:
- Either harden/pay down this support-partition seam or implement only the first
  private provider-block payload against these tiles. Do not jump to residual
  MWG or supplemented energies.

Line-count / complexity note:
- Scoped `src + test + bin` was `435` added / `0` deleted, net `+435`.
- Exception accepted because the payload prevents a worse provider-block
  authority error. The next cleanup-capable pass should pay this down with a
  mature deletion candidate.

## Pass 262 - Terminal Assembly Flat-Mirror Paydown

Commit(s):
- this commit - Shrink terminal assembly low-order flat mirrors

Summary:
- Tests-only cleanup in
  `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`.
- Removed duplicated terminal-shellification flat assembly/summary mirror
  assertions, deferred-status aliases, duplicate false matrix/materialization
  checks, and empty helper-family mirrors.
- Preserved compact live checks for terminal policy selection, summary/deferred
  status, region/unit counts, midpoint/distorted-box facts, non-CPB support
  records, and major stage-object preservation.

Validation:
- Doer: package-load/parse validation and `git diff --check` passed.
- Manager: reran package-load/parse validation and `git diff --check`; both
  passed. The stale low-order integration file was intentionally not used as a
  per-pass gate.

Goal advancement:
- MT5/LT2: paid down part of the pass-261 implementation exception by removing
  stale test pressure.
- MT6/AG7: continued classifying old flat terminal-shellification report keys
  as migration scaffolding rather than public route authority.

Medium-goal update:
- none.

Risk / guardrail:
- The test still protects compact terminal-policy behavior. Do not re-expand it
  with exact flat report-key mirrors.

Remaining blocker / next:
- The pass-261 line-count debt is only partially paid down. Continue with
  mature cleanup candidates, or proceed to provider-block payload work only if
  it consumes the support-partition payload and keeps matrices local/provider
  level.

Line-count / complexity note:
- Scoped `src + test + bin` was `6` added / `51` deleted, net `-45`.

## Pass 263 - Current-Route Metadata Export Retirement

Commit(s):
- this commit - Retire current-route metadata export stack

Summary:
- Deleted the old current-route metadata export stack:
  `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`.
- Removed its include from `CartesianContractedParentMetrics`, deleted the old
  Be2 source-metadata acceptance support/test files, and removed that test from
  the slow nested integration runner.
- Caller audit found no remaining `src/test/bin` users of `_pqs_current_route_*`,
  `_be2_pqs_q5_source_metadata_*`,
  `pqs_source_metadata_real_artifact_acceptance_*`, or
  `current_route_metadata_export`.

Validation:
- Doer: package load, integration-runner parse, caller-proof `rg` checks, and
  `git diff --check` passed.
- Manager: reran the same checks; all passed. The slow integration runner was
  intentionally not executed.

Goal advancement:
- LT2/MT5: removed 5,981 scoped lines, more than paying down the pass-261
  support-partition implementation exception.
- MT6/AG7: retired old flat Cartesian/PQS migration scaffolding instead of
  preserving it as public architecture.
- LT5: removed a stale source-metadata route-authority surface that was not the
  independent H2 PQS route.

Medium-goal update:
- none.

Risk / guardrail:
- Historical docs and old pass logs still mention the retired helpers as
  history. Do not revive them through compatibility wrappers; if a needed fact
  reappears, route it through current typed route-owned objects.

Remaining blocker / next:
- Re-audit `legacy_source_box_fixtures.jl` and `source_box_pair_shadow.jl`
  pockets that were previously kept alive by current-route metadata callers, or
  resume supplement provider-block staging against the new support-partition
  payload.

Line-count / complexity note:
- Scoped `src + test + bin` was `0` added / `5,981` deleted, net `-5,981`.

## Pass 264 - Legacy Source-Box Fixture Builder Retirement

Commit(s):
- this commit - Retire legacy source-box fixture builders

Summary:
- Deleted the now-unblocked legacy retained-unit fixture builder tail from
  `src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl`.
- Removed the old contact-cap, outer-mismatch, atom-box support-dense fixture
  builders and their private route-fact audit/helper scaffolding.
- The deletion was enabled by pass 263, which removed the only old metadata
  export path that had kept these helpers alive.

Validation:
- Doer: package load, caller-proof `rg`, and `git diff --check` passed.
- Manager: reran package load, caller-proof `rg`, and `git diff --check`; all
  passed.

Goal advancement:
- LT2/MT5: deleted another 999 scoped lines from old migration scaffolding.
- MT6/AG7: classified these fixture builders as old flat Cartesian route-shadow
  scaffolding, not active route-owned implementation.
- LT5: kept independent H2 PQS route authority untouched while removing stale
  alternate fixture authority.

Medium-goal update:
- Added the required checkpoint below for passes 260-264.

Risk / guardrail:
- Active product/doside source-box pair-plan, density-density, nuclear/local
  Gaussian, and raw-plan diagnostics remain in place. Do not delete them without
  a fresh caller audit and replacement-authority review.

Remaining blocker / next:
- The pass-261 implementation exception is now fully paid down. The next cleanup
  candidate can be either remaining legacy source-box pockets or tracked
  docs-history compression for old PQS blurb logs.

Line-count / complexity note:
- Scoped `src + test + bin` was `0` added / `999` deleted, net `-999`.

## Medium-Term Goal Checkpoint - Passes 260-264

- MT1 Fake-PQS quarantine: active/maintained. None of passes 260-264 promoted
  fake-PQS or source-backed WL/QW data as independent PQS evidence.
- MT2 Independent H2 PQS recovery: active. The base independent route remains
  available through final-basis, H1, H1-J, and private RHF diagnostics; these
  passes did not alter those endpoint seams.
- MT3 Common physical support vocabulary: active. Passes 260-261 clarified and
  materialized support tiling/row ownership for the independent H2 supplement
  lane.
- MT4 Supplement staging after authority: active. The lane advanced only to
  support-partition pre-provider infrastructure; provider blocks and
  supplemented values remain blocked.
- MT5 Cleanup pressure: active and strong. Pass 261 took a `+435` implementation
  exception; passes 262-264 deleted 7,025 scoped lines combined, more than
  paying it down.
- MT6 Audit/classify old Cartesian flat paths: active. Passes 263-264 retired
  old current-route metadata export and legacy retained-unit fixture-builder
  scaffolding after caller audits.

## Docs-History Compression - PQS Passes 101-229

Commit(s):
- this commit - Compress PQS blurb logs 101-229

Summary:
- Consolidated raw PQS blurb/response/review logs for passes 101-229 into
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/summary.md`.
- Kept durable decisions, authority changes, validation landmarks, cleanup
  milestones, and fake-PQS guardrails for the three windows 101-150, 151-200,
  and 201-229.
- Deleted the corresponding tracked raw `blurb.NNN.md`, `response.NNN.md`, and
  `review.NNN.md` files. The active 230-current tail remains intact.

Validation:
- Manager: verified no tracked raw blurb/response/review files remain in the
  101-229 range; `git diff --cached --check` passed.

Goal advancement:
- LT2/MT5: removed historical log bulk while preserving strategic content.
- MT6/AG7: kept the old flat-path/fake-PQS interpretation in curated summaries
  rather than requiring agents to reread obsolete raw handoffs.

Medium-goal update:
- none.

Risk / guardrail:
- Do not delete the active 230-current raw logs until their decisions have aged
  enough to compress safely.

Remaining blocker / next:
- Continue normal independent-H2-PQS/supplement-staging work. Future docs-log
  compression should start from 230-current only after a newer active tail
  exists.

Line-count / complexity note:
- Added 197 summary lines and deleted 38,867 raw markdown lines, net
  `-38,670` tracked docs lines.

## Pass 265 - H2 PQS H1 Common Operator Path

Commit(s):
- this commit - Route H2 PQS H1 through complete core shell path

Summary:
- Added a narrow adapter from the current H2 independent physical source plan
  into the existing `pqs_multilayer_shell_source_plan` / complete core-shell
  operator contract.
- Routed the H2 PQS physical H1 payload through
  `pqs_multilayer_complete_core_shell_h1_payload` while preserving the current
  driver artifact behavior.
- Extended the common final H1 solve with the lowest orbital coefficients so
  the existing H1-J and RHF consumers continue to receive real computational
  data.
- Deleted the now-dead H2-local duplicate support kinetic, support
  electron-nuclear, final one-body, final H1 Hamiltonian, and H1 solve helper
  family.

Validation:
- Doer: `git diff --check`, package load, and
  `tools/run_cartesian_line_ladder.jl --line=pqs_diatomic` passed. The ladder
  checked the H2 PQS base route, residual-GTO preflight, and materialized
  sidecar/provider artifact route.

Goal advancement:
- LT2/LT8: removed duplicate diatomic-local one-body authority and moved H2
  PQS H1 onto the shared complete core-shell operator path.
- LT5: preserved independent H2 PQS route authority while making the operator
  authority explicit in the H1 payload metadata.
- MT4: keeps residual-GTO/provider work staged above a common one-body path
  rather than adding another route-specific operator stack.

Medium-goal update:
- none.

Risk / guardrail:
- This is an H2 adapter into the existing common path, not the final atomic /
  diatomic unification. Do not generalize it into a provider registry or revive
  status/readiness payloads.

Remaining blocker / next:
- Atomic PQS/WL still has its own direct materializer. The next unification
  pass should design or implement the atomic source/support/final-basis adapter
  below geometry construction and above one-body operators.

Line-count / complexity note:
- Scoped source diff before the log was `188` added / `266` deleted, net
  `-78`; the net direction remains simplification while adding the adapter seam.

## Pass 266 - Atomic WL H1 Common Operator Adapter

Commit(s):
- this commit - Route atomic WL H1 through complete core shell path

Summary:
- Added a narrow one-center atomic adapter from the existing full-parent shell
  sequence into the shared `pqs_multilayer_shell_source_plan` / complete
  core-shell operator contract.
- Changed the WL atomic pure-gausslet materializer so final basis, kinetic,
  charged nuclear, H1, and H1 lowest energy come from
  `pqs_multilayer_complete_core_shell_h1_payload`.
- Kept the old fixed-block path only as an optional artifact sidecar source for
  legacy pair-sum/fixed-center fields when artifact saving is requested; it is
  no longer the H1 authority for the active driver materialization.
- Added explicit adapter checks for shell duplicate support, core/shell
  disjointness, sequence support agreement, retained count agreement, and shell
  coefficient row agreement.

Validation:
- Doer: `git diff --check`, package load,
  `tools/run_cartesian_line_ladder.jl --line=wl_atomic`, and
  `tools/run_cartesian_line_ladder.jl --line=pqs_atomic` passed.

Goal advancement:
- LT8: atomic and H2 now both have a path through the same complete core-shell
  one-body authority above geometry/source construction.
- LT5: preserved route provenance by keeping the atomic shell sequence as the
  source/support producer while moving operator authority to the common path.
- LT2/MT5: this is a positive-line adapter pass, but it avoids creating a new
  framework and sets up later deletion of atomic fixed-block H1 use.

Medium-goal update:
- none.

Risk / guardrail:
- Atomic GTO/supplement materialization remains intentionally blocked/not
  materialized. Do not interpret this pass as full atomic supplement support or
  as a generalized provider registry.

Remaining blocker / next:
- The next common-path cleanup should either retire old atomic fixed-block H1
  artifact assumptions after a replacement pair/density sidecar exists, or
  route atomic PQS materialization through the same adapter once that endpoint
  becomes active.

Line-count / complexity note:
- Scoped source diff before the log was `214` added / `22` deleted, net `+192`.
  The addition is accepted as a narrow adapter seam; the carrying cost should be
  paid down by deleting fixed-block H1 authority once pair/density sidecars are
  replaced.

## Pass 267 - Atomic WL Artifact Sidecar Without Fixed Block

Commit(s):
- this commit - Remove atomic WL fixed-block artifact dependency

Summary:
- Removed the remaining `_nested_fixed_block(common.sequence, common.bundle)`
  call from the WL atomic pure-gausslet materializer.
- Added a compact artifact-sidecar projection that consumes the already-built
  atomic shell sequence packet and the common complete core-shell final-basis
  cleanup transform.
- Basis/Ham artifact sidecars now write final-gauge weights, fixed centers, and
  the density-density pair matrix directly from common-path data, without
  constructing a `_NestedFixedBlock3D`.

Validation:
- Doer: `git diff --check`, package load,
  `tools/run_cartesian_line_ladder.jl --line=wl_atomic`, and
  `tools/run_cartesian_line_ladder.jl --line=pqs_atomic` passed.
- Doer also ran a direct WL atomic artifact-save smoke with temporary JLD2
  files and reloaded the sidecar fields: basis coefficients `(1331, 419)`,
  fixed centers `(419, 3)`, one-body Hamiltonian `(419, 419)`, and density
  pair matrix `(419, 419)`.

Goal advancement:
- LT8: keeps atomic WL on the common complete core-shell one-body authority
  while removing the last fixed-block conversion from the active artifact path.
- LT2/MT5: narrows old fixed-block authority to historical/reference code
  rather than current driver materialization.
- LT5: preserves the atomic shell sequence as the route-owned source/support
  producer while final-gauge artifact facts come from common final-basis data.

Medium-goal update:
- none.

Risk / guardrail:
- This is not a new density/pair provider implementation. The pair sidecar is
  the existing shell packet pair matrix projected through the common final-basis
  cleanup transform. Do not treat it as supplemented atomic pair-provider
  support.

Remaining blocker / next:
- Atomic PQS materialization is still not implemented as an endpoint. The next
  unification step should either add the atomic PQS endpoint through the same
  common-path adapter or continue retiring old atomic fixed-block helper
  pressure after caller audit.

Line-count / complexity note:
- Scoped source diff before the log was `46` added / `8` deleted, net `+38`.
  The extra helper replaces a conceptually heavier legacy block conversion and
  should make later fixed-block retirement simpler.

## Pass 268 - H2 PQS Residual-GTO Density Descriptor

Commit(s):
- this commit - Add H2 PQS residual GTO density descriptor

Summary:
- Added a narrow density-provider descriptor for the H2 independent PQS
  residual-GTO materialized route.
- The descriptor explicitly fixes the next density-provider gauge as
  `(:pre_final_pqs, :residual_gto)` with
  `density_gauge = :pre_final_localized_positive_weight`, not `[F, R]`.
- It carries the concrete residual-orbital carrier data needed by the next
  provider pass: `p_projection_of_g = A*S_FG` and
  `[-A*S_FG; I_G] * L`, where `A` is the common final-basis cleanup and `L` is
  the residual-GTO transform.
- The H2 sidecar artifact roundtrip now checks the descriptor shape and
  provenance, while still leaving residual density moments, P-R/R-R blocks,
  supplemented values, and RHF untouched.

Validation:
- Doer: `git diff --check`, package load, and
  `tools/run_cartesian_line_ladder.jl --line=pqs_diatomic` passed.
- Doer inspected the materialized artifacts and confirmed
  `density_provider_kind = :descriptor_only`,
  `augmented_density_space = (:pre_final_pqs, :residual_gto)`,
  `p_projection_of_g` shape `(471, 18)`, residual carrier shape `(489, 18)`,
  P-P pair matrix size `(471, 471)`, residual mode source
  `:requires_mwg_residual_moments`, and
  `fixed_block_pair_data_authority_used = false`.

Goal advancement:
- MT4: starts P1e at the correct density-gauge boundary before building
  density/pair provider values.
- LT8: keeps the common complete core-shell positive-weight pre-final density
  convention as authority for the PQS sector.
- LT5: prevents a scientifically wrong shortcut through `[F, R]` or
  `L' * V_GG * L` by making the projected residual orbital carrier explicit.

Medium-goal update:
- none.

Risk / guardrail:
- This is descriptor-only. It does not materialize P-R/R-R provider blocks,
  supplemented H1-J, density/pair values, RHF, or CR2 artifacts. Do not
  interpret the new descriptor as provider-block completion.

Remaining blocker / next:
- Implement MWG residual moments for projected residual orbitals, then build
  the first real P-R and R-R density/provider blocks in the
  `[pre_final PQS, residual GTO]` gauge.

Line-count / complexity note:
- Scoped source diff before the log was `318` added / `0` deleted. Most of
  the increase is local artifact readback validation for the new contract. This
  is a deliberate guardrail pass; the next value-producing pass should avoid
  expanding it into a general provider registry.

## Pass 269 - H2 PQS Residual-GTO Density Provider Blocks

Commit(s):
- this commit - Add H2 PQS residual GTO density provider blocks

Summary:
- Advanced the P1e lane from descriptor-only to the first concrete
  residual-GTO density/provider values for the H2 independent PQS materialized
  route.
- The materialized case now requests
  `residual_gto_provider_blocks = :one_body_and_density_provider`.
- Built residual MWG moment descriptors for projected residual orbitals in the
  explicit `[pre_final_pqs, residual_gto]` density gauge, then used the MWG
  interaction kernel to form P-R, R-R, and assembled augmented density/provider
  pair matrices.
- The artifact contract now writes and roundtrips concrete P-R/R-R fields:
  `v_pr_pair_matrix`, `v_rr_pair_matrix`, `augmented_pair_matrix`,
  residual centers/widths, source labels, and symmetry/shape facts.
- The pass does not compute supplemented H1-J scalar values, RHF, CR2, or a
  full two-body Hamiltonian.

Validation:
- Manager/doer: `git diff --check` passed.
- Manager/doer: package load passed with
  `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- Manager/doer: `tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases.
- Manager/doer inspected the written artifacts and confirmed
  `provider_blocks_included = :one_body_and_density_provider`,
  `density_provider_kind = :augmented_residual_gto_pair_matrix`,
  density space `(:pre_final_pqs, :residual_gto)`, residual center/width shapes
  `(18, 3)`, P-R size `(471, 18)`, R-R size `(18, 18)`, augmented pair size
  `(489, 489)`, and augmented pair symmetry error
  `2.55351295663786e-15`.

Goal advancement:
- MT4/LT5: converts the residual-GTO density boundary from "known missing
  object" to concrete provider blocks without entering supplemented values.
- LT8: preserves the common positive-weight pre-final PQS density convention:
  the augmented density gauge is `[P, R]`, not `[F, R]`.
- MT1: keeps fake-PQS/source-backed WL/QW authority quarantined; the new values
  use direct source/support moment kernels and MWG interaction components, not
  receipt wrappers.

Medium-goal update:
- none.

Risk / guardrail:
- The residual density sector uses MWG moment descriptors for projected
  residual orbitals. That is the intended first provider value, but it still
  needs a numerical/science audit before being treated as a production
  supplemented-H1-J authority.
- The artifact field is honestly named
  `:one_body_and_density_provider`; it must not be interpreted as full
  supplemented density/pair/H1-J/RHF completion.

Remaining blocker / next:
- Next P1 work should decide how to consume the augmented density/provider
  matrix: first as a supplemented H1-J diagnostic or as the minimal object
  needed by the existing RHF path. Do not skip directly to a broad provider
  registry.

Line-count / complexity note:
- Scoped source/input diff before the log was `609` added / `28` deleted. This
  is a positive-line value-producing pass, with most new bulk in local kernel
  plumbing and artifact roundtrip checks. No new tests, status/readiness
  payloads, QW receipt wrappers, or provider registry were added.

## Medium-Term Goal Checkpoint - Passes 265-269

- MT1 Fake-PQS quarantine: active/maintained. Passes 265-269 moved H2 H1 and
  residual-GTO density/provider work through common complete-core-shell and
  source/support kernels without promoting fake-PQS or old source-backed WL/QW
  routes as independent PQS authority.
- MT2 Independent H2 PQS recovery: active and now extended. The base route
  remains available through final-basis, H1, H1-J, private RHF diagnostics, and
  the materialized residual-GTO Ham/Basis artifact.
- MT3 Common physical support vocabulary: active and strengthened. H2 H1 now
  uses the common complete-core-shell H1 path, and atomic WL has a common
  packet/sidecar cleanup path, but atomic PQS is still not an endpoint.
- MT4 Supplement staging after authority: active and substantially advanced.
  The lane progressed from sidecar/roundtrip and one-body provider fields to a
  concrete `[pre_final_pqs, residual_gto]` density/provider matrix. Remaining
  work is supplemented H1-J/RHF consumption and a science/performance audit.
- MT5 Cleanup pressure: active but temporarily line-positive. These passes were
  architecture-paydown and P1 feature-slice work rather than demolition; pass
  267 deleted old fixed-block pressure, while passes 268-269 added the guarded
  density-provider contract needed to retire donor lines later.
- MT6 Audit/classify old Cartesian flat paths: active. No new helper/schema
  tests were added; the driver/line ladders remain the validation authority for
  this route lane.

## Pass 270 - Residual-GTO P-R Density Gauge Correction

Commit(s):
- this commit - Correct H2 PQS residual GTO P-R density normalization

Summary:
- Corrected the P-R density-provider block for the H2 independent PQS
  residual-GTO materialized route.
- The prior implementation contracted the pre-final weighted PQS coefficients
  against density-normalized support-residual MWG rows directly. That missed
  the support-row weight factor required by the same raw-numerator convention
  used to build the trusted P-P pre-final pair matrix.
- The corrected path recomputes support-row weights from the route source-plan
  axis PGDG weights in `packet.support_states` order, forms
  `diag(support_weights) * gausslet_residual`, and only then applies
  `pre_final_coefficients ./ pre_final_weights`.
- Added two narrow ham-artifact provenance fields:
  `support_residual_input_convention = :density_normalized` and
  `p_r_weight_application = :support_row_weights_before_pre_final_projection`.

Validation:
- Manager/doer: `git diff --check` passed.
- Manager/doer: package load passed with
  `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- Manager/doer: `tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases.
- Manager/doer reran the artifact PSD/Schur audit after the ladder rewrote the
  artifacts. The augmented matrix eigenvalue minimum changed from the previous
  failing `-6.7885` to `3.1943e-9`, and the regularized Schur eigenvalue range
  is now `(3.1943e-9, 0.01952)`.

Goal advancement:
- MT4/LT8: repairs the density gauge so the `[pre_final_pqs, residual_gto]`
  provider matrix is internally consistent enough for the next private H1-J
  diagnostic attempt.
- LT5: preserves source/support authority and the common pre-final
  positive-weight convention rather than introducing a final-gauge shortcut.

Medium-goal update:
- none.

Risk / guardrail:
- This is still provider-block readiness, not a supplemented value. No H1-J
  scalar, RHF, CR2, public API, or provider registry was added.

Remaining blocker / next:
- The next pass can attempt a private supplemented H1-J diagnostic using the
  corrected augmented one-body and density-provider objects, with coarse
  semantic checks only.

Line-count / complexity note:
- Scoped source diff before the log was `35` added / `2` deleted. The patch is
  intentionally narrow and fixes an existing gauge error rather than adding a
  new feature surface.

## Pass 271 - Private Residual-GTO Augmented H1-J Diagnostic

Commit(s):
- this commit - Add H2 PQS residual GTO private H1-J diagnostic

Summary:
- Added the first private consumer of the corrected H2 residual-GTO
  density-provider blocks.
- The diagnostic solves the augmented one-body Hamiltonian in `[F, R]`, maps the
  `F` component through the existing `final_to_pre_final_coefficients`, keeps
  the residual component in the `R` sector, and contracts the resulting
  `[P, R]` density coefficient vector with the augmented pair matrix.
- The ham artifact now records
  `augmented_h1_j_diagnostic_kind = :private_augmented_h1_j_self_coulomb`,
  `augmented_h1_j_self_coulomb`, density coefficient length, and orbital
  source. The artifact still records `supplemented_values_kind = :not_computed`.
- No RHF, CR2, public API, provider registry, or generalized two-body framework
  was added.

Validation:
- Manager/doer: `git diff --check` passed.
- Manager/doer: package load passed with
  `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- Manager/doer: `tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases after the final roundtrip-shape patch.
- Manager/doer inspected the written ham artifact and confirmed
  `augmented_h1_j_self_coulomb = 0.457435475059184`,
  base `h1_j_self_coulomb = 0.4569117646737236`,
  `augmented_h1_lowest = -0.79590283450777`, and augmented density coefficient
  length `489`.

Goal advancement:
- MT4/LT8: moves the residual-GTO provider lane from block readiness to a
  private H1-J diagnostic consumer while preserving the `[P, R]` density gauge.
- LT5: continues to consume route-owned source/support artifacts and the common
  complete-core-shell density convention.

Medium-goal update:
- none.

Risk / guardrail:
- This is a private scalar diagnostic, not a production supplemented H1-J value
  and not an RHF input contract. The next step should audit the scalar and
  decide whether it is adequate as the minimal bridge to private augmented RHF.

Remaining blocker / next:
- The immediate next pass can be a narrow science/performance audit of the
  private augmented H1-J diagnostic, or, if accepted, a private augmented RHF
  smoke that consumes the same `[F, R]` one-body and `[P, R]` density objects.

Line-count / complexity note:
- Scoped source diff before the log was `117` added / `0` deleted. This is
  positive-line feature work with no new tests or registries; a later cleanup
  pass should consolidate local scalar-contraction helpers if more consumers
  appear.

## Pass 272 - Private Augmented H1-J Diagnostic Audit

Commit(s):
- this commit - Audit H2 PQS residual GTO private H1-J diagnostic

Summary:
- Performed a read-only numerical audit of the private H2 PQS residual-GTO
  augmented H1-J diagnostic introduced in pass 271.
- Recomputed the augmented one-body lowest orbital, mapped its F sector through
  the pre-final PQS density transform, appended the residual R coefficients,
  and decomposed the `[P, R]` self-Coulomb diagnostic into P-P, P-R, and R-R
  contributions.
- The audit did not change source behavior, did not add RHF, CR2, public API,
  tests, or provider registries, and left `supplemented_values_kind =
  :not_computed`.

Validation:
- Manager/doer inspected the current artifacts at
  `/Users/srw/dmrgtmp/h2_pqs_q5_independent_source_box_r4_gto_basis.jld2` and
  `/Users/srw/dmrgtmp/h2_pqs_q5_independent_source_box_r4_gto_ham.jld2`.
- Recomputed `augmented_h1_lowest = -0.7959028345077876`, matching the artifact
  value `-0.79590283450777`.
- Recomputed `augmented_h1_j_self_coulomb = 0.457435475059184`, matching the
  artifact. Base `h1_j_self_coulomb` remains `0.4569117646737236`.
- Augmented orbital weights are F `0.9998627116279778` and R
  `0.0001372883720226732`; the density coefficient vector has length `489`.
- Self-Coulomb decomposition:
  P-P `0.4573636390600626`, P-R `7.183224722295418e-5`, R-R
  `3.75189826063245e-9`, summing to `0.4574354750591838`.
- The augmented pair matrix remains positive to audit tolerance:
  eigenvalue range `(3.194313211426547e-9, 141.59598639666223)`, regularized
  Schur range `(3.194312980933833e-9, 0.01952029093239803)`.

Goal advancement:
- MT4/LT8: classifies the private augmented H1-J diagnostic as numerically sane
  enough for a first private augmented RHF smoke.
- LT5: confirms that the scalar is dominated by the trusted P-P sector with a
  small P-R correction and negligible R-R contribution for this H2 fixture.

Medium-goal update:
- none.

Risk / guardrail:
- This is still not a production supplemented H1-J value. The diagnostic is a
  private consumer sanity check over the current provider blocks. Any RHF pass
  must keep the same private/diagnostic label until a separate science and
  performance review accepts the route.

Remaining blocker / next:
- The next pass may attempt a private augmented RHF smoke using the augmented
  one-body Hamiltonian and corrected `[P, R]` density provider, with only
  convergence/trace/residual diagnostics. Do not add CR2, public API, or broad
  provider abstractions in that pass.

Line-count / complexity note:
- No source-code behavior changed. The only tracked change is this audit log
  entry.

## Pass 273 - Private Residual-GTO Augmented RHF Smoke

Commit(s):
- this commit - Add H2 PQS residual GTO private RHF smoke

Summary:
- Added a private closed-shell RHF smoke for the H2 PQS residual-GTO materialized
  route. The smoke consumes the augmented one-body Hamiltonian in `[F, R]` and
  the corrected augmented density provider in `[P, R]`.
- The implementation is intentionally local and bounded: a 25-iteration maximum
  SCF, first-three-iteration density damping, convergence on energy, density
  change, and Fock-density commutator, and no DIIS/general solver framework.
- The ham artifact now carries compact private RHF facts and the roundtrip
  check validates convergence, finite energies, trace, idempotency, and
  commutator residual.
- The route remains private/diagnostic only: `private_augmented_rhf_public_api =
  false`, no CR2/export consumer contract, no provider registry, and no
  production supplemented RHF claim.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Artifact readback from
  `/Users/srw/dmrgtmp/h2_pqs_q5_independent_source_box_r4_gto_ham.jld2`:
  converged `true`, iterations `15`, electronic energy
  `-1.1611254039289651`, nuclear repulsion `0.25`, total with nuclear
  repulsion `-0.9111254039289651`, density trace `2.0000000000000018`, trace
  error `1.7763568394002505e-15`, idempotency error
  `1.249000902703301e-16`, commutator residual `9.90647328058536e-10`, F/R
  orbital weights `0.9998235792815352` / `0.0001764207184651116`.

Goal advancement:
- MT4/LT8: advances the residual-GTO P1 lane from a private H1-J scalar
  diagnostic to the first private augmented RHF smoke over the same provider
  objects.
- LT5: keeps the authority on route-owned augmented one-body and density
  provider data, with no fake WL/QW receipt wrappers or readiness/status
  payloads.

Medium-goal update:
- none.

Risk / guardrail:
- This is not a public solver, not a production supplemented RHF value, and not
  a broad RHF framework. It is a narrow private smoke for one H2 PQS residual-GTO
  fixture. Science and performance review are still required before broadening
  it or using it as a downstream consumer contract.

Remaining blocker / next:
- Audit the private RHF smoke numerically and operationally: energy
  decomposition, stability under tolerances, timing/allocation cost, and whether
  the local Fock construction is the right long-lived consumer. After that,
  choose between cleanup/paydown or a broader supplemented-consumer step.

Line-count / complexity note:
- The source diff before this log entry was `369` added / `0` deleted. This is
  positive-line feature work on an active physics target; it should be followed
  by audit and later paydown once the consumer shape is accepted.

## Pass 274 - Private RHF Smoke Payload Paydown

Commit(s):
- this commit - Trim H2 PQS residual GTO RHF smoke payload

Summary:
- Reduced bloat introduced by the private H2 PQS residual-GTO RHF smoke without
  changing the RHF computation.
- Simplified the nuclear-repulsion displacement expression and removed an
  unused final-orbital sector decomposition/eigensolve from the RHF helper.
- Stripped duplicate artifact-roundtrip assertions for the private RHF smoke.
  The SCF still enforces convergence before writing; reload now stays focused
  on the existing artifact shape rather than re-testing every convergence
  scalar.
- Removed nonessential saved/returned fields such as private-public status,
  electron-count mirror, delta bookkeeping, and F/R orbital weights. The ham
  artifact still carries the useful private smoke facts: kind, convergence,
  iteration count, density trace, idempotency, commutator, one-/two-body energy,
  electronic energy, nuclear repulsion, and total energy.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Artifact readback confirms the retained private RHF facts still exist:
  converged `true`, iterations `15`, density trace `2.0000000000000018`,
  idempotency error `1.249000902703301e-16`, commutator residual
  `9.90647328058536e-10`, electronic energy `-1.1611254039289651`, total with
  nuclear repulsion `-0.9111254039289651`.

Goal advancement:
- MT4/LT8: keeps the private RHF smoke working while reducing the carrying cost
  of the P1 residual-GTO lane.
- LT5: reinforces the driver-ladder validation policy and avoids rebuilding
  artifact/schema test pressure inside the runtime materializer.

Medium-goal update:
- none.

Risk / guardrail:
- This remains a private diagnostic, not a production supplemented RHF contract.
  The cleanup intentionally did not remove the core convergence checks inside
  the SCF helper.

Remaining blocker / next:
- The same science/performance audit remains the next decision point before
  broadening the private RHF smoke or using it as a downstream consumer
  contract.

Line-count / complexity note:
- Net cleanup before this log entry was `1` added / `114` deleted in
  `src/pqs_source_box_low_order_materialization.jl`; the full commit is net
  negative despite this log entry.

## Pass 275 - Residual-GTO Provider Payload Paydown

Commit(s):
- this commit - Trim H2 PQS residual GTO provider payload labels

Summary:
- Audited recent residual-GTO provider additions for payload/schema bloat and
  removed label-only mirrors from the descriptor, density-provider, H1-J
  diagnostic, JLD2 writer, and artifact readback path.
- Deleted fields and reload checks such as density-provider kind mirrors,
  duplicated sector labels, residual-mode source labels, final-gauge/fixed-block
  negative flags, component-source labels, P-R weight-application labels, and
  H1-J diagnostic/source labels.
- Kept concrete data and physics facts: residual transform, residual carrier
  coefficients, augmented density gauge/space, dimensions, moment centers and
  widths, one-body/provider matrices, augmented pair matrix, H1-J self-Coulomb,
  and private RHF convergence/energy facts.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Artifact readback still shows the retained private RHF facts: converged
  `true`, iterations `15`, density trace `2.0000000000000018`, commutator
  residual `9.90647328058536e-10`, total with nuclear repulsion
  `-0.9111254039289651`.

Goal advancement:
- MT4/LT8: keeps the H2 PQS residual-GTO provider lane functional while
  reducing runtime artifact/schema carrying cost.
- LT5: reinforces that the thin route should carry real provider data and
  compact physics facts, not status/provenance mirrors or helper-schema tests.

Medium-goal update:
- none.

Risk / guardrail:
- This pass deliberately did not change the numerical provider construction or
  private RHF smoke. The remaining P1 risks are still science/performance
  validation and deciding which consumer contract should survive.

Remaining blocker / next:
- Continue with the private RHF science/performance audit or a further paydown
  pass if review finds more non-physics payload duplication.

Line-count / complexity note:
- Source cleanup before this log entry was `210` deleted / `0` added in
  `src/pqs_source_box_low_order_materialization.jl`.

## Pass 276 - Private Residual-GTO RHF Science/Performance Audit

Commit(s):
- this commit - Audit H2 PQS residual GTO private RHF smoke

Summary:
- Performed a read-only audit of the private H2 PQS residual-GTO augmented RHF
  smoke after the payload paydown passes.
- Re-read the H2 PQS residual-GTO Ham/Basis artifacts and checked the augmented
  one-body, density-provider, H1-J, and RHF scalar facts. No source behavior,
  artifact schema, tests, public API, or provider registry changed.
- Compared against the base H2 PQS private-RHF driver path where the current
  driver exposes data: the base run completes and reports commutator residual
  `3.965093933744335e-9`; the compact summary does not currently expose a base
  total energy, so that remains a visibility gap rather than a reason to
  rebuild reporting.

Validation / audit results:
- Base driver probe:
  `julia --project=. -e 'empty!(ARGS); append!(ARGS,
  ["test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl",
  "save_artifact=false", "save_tsv=false"]); include("bin/cartesian_ham_builder.jl")'`
  completed through `cartesian_print/save`. Assembly took `23.466620` seconds
  with `53.17 M` allocations / `6.825 GiB`, mostly compilation.
- Artifact readback/eigen audit elapsed `1.486418084` seconds.
- Artifact facts: provider blocks `:one_body_and_density_provider`, final
  dimension `471`, residual rank `18`, augmented dimension `489`, augmented
  density dimension `489`, final overlap identity error
  `1.295907825493714e-13`.
- One-body/H1-J facts: base H1 lowest `-0.7946037173365885`, augmented H1
  lowest `-0.79590283450777`, augmented H1 symmetry error
  `1.0658141036401503e-14`, base H1-J self-Coulomb
  `0.4569117646737236`, augmented H1-J self-Coulomb
  `0.457435475059184`, delta `5.237103854604519e-4`.
- Residual density facts: residual overlap identity error
  `3.753268124085306e-15`, residual moment overlap error
  `1.976174779372286e-11`, residual width range
  `(0.7873020609606813, 10.51582331678006)`.
- Density-provider PSD/symmetry facts: V_PP eigen range
  `(0.2006880365831959, 138.44830404378297)`, V_RR eigen range
  `(2.1365833223596135e-8, 3.5225668776085977)`, augmented pair eigen range
  `(3.1943132473945437e-9, 141.59598639666217)`, regularized Schur eigen range
  `(3.194312938499561e-9, 0.01952029093239794)`, augmented pair symmetry error
  `2.55351295663786e-15`.
- Private augmented RHF facts: converged `true`, iterations `15`, density trace
  `2.0000000000000018`, idempotency error `1.249000902703301e-16`,
  commutator residual `9.90647328058536e-10`, one-body energy
  `-1.5634981066427645`, two-body energy `0.4023727027137994`, electronic
  energy `-1.1611254039289651`, nuclear repulsion `0.25`, total with nuclear
  repulsion `-0.9111254039289651`.

Goal advancement:
- MT4/LT8: classifies the private augmented RHF smoke as coherent enough for a
  next private supplemented-consumer step. The provider matrix is positive to
  audit tolerance, the self-Coulomb correction is small and positive, and the
  RHF fixed point is numerically tight.
- LT5: confirms that the thin route remains anchored in route-owned augmented
  one-body and `[P, R]` density-provider data, not QW/WL receipt wrappers.

Medium-goal update:
- none.

Risk / guardrail:
- This is still private/prototype-grade. Performance is compilation-heavy and
  allocation-heavy at the driver level, and the Fock construction is local to
  the slice rather than a reviewed production solver. Base private RHF energy
  is not exposed in the current compact print/report path.

Remaining blocker / next:
- It is reasonable to move to the next private supplemented-consumer step, but
  not to public RHF or CR2. Before any production claim, do a performance pass
  on provider construction/Fock assembly and decide whether the private RHF
  consumer should be refactored into a reusable module-owned contract.

Line-count / complexity note:
- No source behavior changed. The audit used a temporary `tmp/work` script that
  was removed before commit; this commit is log-only.

## Pass 277 - Print Private Residual-GTO Consumer Facts

Commit(s):
- this commit - Print H2 PQS residual GTO private consumer facts

Summary:
- Promoted the existing private H2 PQS residual-GTO consumer facts to driver
  glass-box output in the `[route_materialization]` section.
- Added no new runtime payloads, tests, artifact fields, solver framework,
  provider registry, public API, or CR2/export path. This pass only makes the
  already-computed materialization facts visible to the driver ladder and human
  inspection.
- Printed facts now include provider mode, residual rank, augmented dimension,
  augmented H1 lowest/symmetry, augmented density dimension/symmetry,
  augmented H1-J self-Coulomb, and private augmented RHF convergence, iteration
  count, energy, total with nuclear repulsion, and commutator residual.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- The materialized H2 PQS residual-GTO case now prints:
  `provider_blocks_included = :one_body_and_density_provider`,
  residual rank `18`, augmented dimension `489`, augmented H1 lowest
  `-0.79590283450777`, augmented H1-J self-Coulomb `0.457435475059184`,
  private augmented RHF converged `true`, iterations `15`, electronic energy
  `-1.1611254039289651`, total with nuclear repulsion
  `-0.9111254039289651`, and commutator residual `9.90647328058536e-10`.

Goal advancement:
- MT4/LT8: turns the private supplemented-consumer lane into visible driver
  diagnostics, which is the intended validation mechanism for this Cartesian
  route.
- LT5: keeps validation at the driver/materialization boundary instead of
  adding helper-schema tests or status/report payloads.

Medium-goal update:
- none.

Risk / guardrail:
- This is still private diagnostic output, not a production supplemented RHF or
  public artifact-consumer contract. Exact print strings must not become tests.

Remaining blocker / next:
- Next useful work is either a small performance/paydown pass on the provider
  construction/Fock assembly, or a carefully scoped private consumer step that
  uses these visible facts without broadening to public RHF/CR2.

Line-count / complexity note:
- Source change before this log entry was `13` added / `0` deleted in the
  simple print helper.

## Pass 278 - Private Residual-GTO Ham Handoff Payload

Commit(s):
- this commit - Define H2 PQS residual GTO Ham handoff payload

Summary:
- Added a private, solver-neutral in-memory Ham handoff payload for the H2 PQS
  residual-GTO materialized route.
- The handoff explicitly separates gauges: one-body Hamiltonian in orbital
  basis `[F, R] = (:final_pqs, :residual_gto)`, density interaction in
  provider basis `[P, R] = (:pre_final_pqs, :residual_gto)`, and
  `orbital_to_density` as the required map from orbital coefficients to the
  density-provider gauge.
- The payload also carries electron count `2`, spin sectors `(nup = 1, ndn =
  1)`, nuclear repulsion, and compact diagnostics. No public API, external
  solver call, CR2 path, provider registry, or new physics was added.
- To avoid accidental report/TSV bloat, the generic durable materialization
  helper now elides the heavy `ham_handoff` object while preserving scalar
  handoff facts in the driver output.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- The materialized case prints `ham_handoff_kind =
  :pqs_h2_residual_gto_ham_handoff`, orbital basis
  `(:final_pqs, :residual_gto)`, density basis
  `(:pre_final_pqs, :residual_gto)`, orbital dimension `489`, and density
  dimension `489`.
- Direct durable-materialization probe confirmed `ham_handoff` is elided from
  generic saved materialization output:
  `d.ham_handoff = nothing` and `heavy_materialization_objects_elided = true`.

Goal advancement:
- MT4/LT8: creates the first explicit solver-neutral Hamiltonian handoff
  boundary for the residual-GTO slice, using already-validated one-body,
  density-provider, H1-J, and private RHF data.
- LT5: preserves the corrected `[F,R]`/`[P,R]` gauge distinction as a concrete
  consumer contract, preventing downstream code from misusing the pair matrix.

Medium-goal update:
- none.

Risk / guardrail:
- This is an in-memory private experimental payload only. The next artifact
  pass must decide how to write/read a compact handoff without dumping route
  payload trees. This is not an HFDMRG integration and not a public Ham API.

Remaining blocker / next:
- Add a narrow JLD2 handoff artifact/group and readback check for this payload:
  H in `[F,R]`, V in `[P,R]`, `T: [F,R] -> [P,R]`, electron/spin counts,
  nuclear repulsion, dimensions, finiteness, and symmetry. Do not add solver
  execution in that pass.

Line-count / complexity note:
- Source diff before this log entry was `130` added / `2` deleted. The added
  lines are a live private contract plus save-path elision; a later extraction
  pass should move the provider/handoff family out of
  `pqs_source_box_low_order_materialization.jl` after the artifact handoff shape
  is accepted.

## Pass 279 - Handoff Summary and Durable Save Paydown

Commit(s):
- this commit - Simplify H2 PQS Ham handoff summary fields

Summary:
- Factored the repeated nullable `ham_handoff` flat-field extraction into
  `_pqs_source_box_route_driver_ham_handoff_summary`.
- The materialization summary and returned materialization object now both
  splice that compact summary while still exposing the same flat
  `ham_handoff_*` fields for the current driver print path.
- Fixed two durable-save details introduced with the handoff: the
  `heavy_materialization_objects_elided` marker now reports true only when a
  non-`nothing` handoff object was actually elided, and TSV output iterates the
  durable materialization keys so that the marker is written.
- Added the missing `augmented_h1_j` contract check to the private Ham handoff
  constructor, so a missing H1-J diagnostic fails at the handoff boundary with
  a local error.
- No physics, artifact fields, handoff object contents, or driver behavior were
  intentionally changed.
- A local scan of the recent materialization code found no other duplicated
  handoff nullable cloud. The larger flat-summary smell for optional
  one-body/density/RHF facts remains; this pass intentionally treated that as a
  deferred carrying-cost cleanup rather than broadening into a summary rewrite.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Direct durable-materialization probe confirmed `(ham_handoff = nothing,
  heavy_materialization_objects_elided = false)` when no object is present and
  `(ham_handoff = nothing, heavy_materialization_objects_elided = true)` when a
  real handoff object is elided.
- The materialized residual-GTO case still printed
  `ham_handoff_kind = pqs_h2_residual_gto_ham_handoff`, orbital/density bases,
  and matching orbital/density dimensions `489`.

Goal advancement:
- LT5/LT8: keeps the new private handoff boundary while reducing flat-field
  duplication and lowering carrying cost in the materialization layer.

Medium-goal update:
- none.

Risk / guardrail:
- The flat `ham_handoff_*` mirrors are still present for current report/print
  compatibility. Longer-term they should probably shrink behind the compact
  `ham_handoff` object or a single `ham_handoff_summary` field when the print
  path can consume that directly. The same principle likely applies to the
  wider optional provider/RHF summary cloud, but that should be a focused
  cleanup pass, not an incidental edit.

Remaining blocker / next:
- Continue with the planned handoff artifact/readback pass, or first extract
  the provider/handoff helpers out of
  `pqs_source_box_low_order_materialization.jl` if carrying cost becomes the
  dominant blocker.

Line-count / complexity note:
- Source/reporting diff before this log entry was `36` added / `45` deleted.
  The log entry itself adds documentation lines, but the code path is net
  smaller.

## Pass 280 - Residual-GTO Artifact Writer Paydown

Commit(s):
- this commit - Trim H2 PQS residual GTO artifact writing

Summary:
- Accepted the doer cleanup that made the residual-GTO provider mode explicit
  as `provider_block_mode` while continuing to write the legacy
  `provider_blocks_included` artifact key for compatibility.
- Replaced repeated nullable summary extraction in the residual-GTO
  materialization summary with small optional property/size helpers.
- Collapsed the new artifact-writer helper family into one local
  `_pqs_source_box_route_driver_write_jld2_fields!` plus two named field lists
  for repeated residual-transform and density-moment artifact fields.
- The long inline basis/Ham JLD2 write blocks remain decomposed, but without a
  dozen one-off setter helpers. No physics or artifact key meaning was changed.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- The materialized case still reports final dimension `471`, residual rank
  `18`, augmented dimension `489`, augmented H1-J self-Coulomb
  `0.457435475059184`, private RHF convergence in `15` iterations, and the
  private Ham handoff dimensions `489`/`489`.

Goal advancement:
- LT5/LT8: preserves the residual-GTO Ham/Basis/sidecar/handoff lane while
  reducing carrying cost in the route-specific materialization writer.
- MT4: keeps the H2 PQS residual-GTO artifact path as the active downstream
  consumer boundary without promoting it to a public solver API.

Medium-goal checkpoint:
- Active: residual-GTO Ham handoff artifact/readback remains the next
  consumer-boundary goal.
- Active: provider/handoff helper extraction from
  `pqs_source_box_low_order_materialization.jl` remains desirable once the
  handoff artifact layout is stable.
- Active: performance review for the residual-GTO provider construction/Fock
  path is still required before broadening beyond the H2 private lane.
- Deferred: atomic PQS unification and density/pair generalization should not
  be mixed into the current handoff-artifact paydown lane.
- Completed for now: current private H2 PQS residual-GTO materialized ladder
  still exercises source plan, final basis, H1, H1-J, one-body provider,
  density provider, private RHF smoke, and private Ham handoff facts.

Risk / guardrail:
- `provider_blocks_included` is now explicitly legacy compatibility beside
  `provider_block_mode`. It should not become a second authority. Delete the
  legacy key once artifact consumers and roundtrip readers no longer need it.
- The flat materialization summary still carries optional provider facts; this
  pass reduced repeated extraction but did not fully replace the flat summary
  with compact nested summaries.

Remaining blocker / next:
- Continue with the narrow H2 residual-GTO Ham handoff artifact/readback pass:
  write/read `H` in `[F,R]`, `V` in `[P,R]`, `T: [F,R] -> [P,R]`,
  electron/spin counts, nuclear repulsion, dimensions, finiteness, and
  symmetry. No solver execution.

Line-count / complexity note:
- Source diff before this log entry was `266` added / `242` deleted in
  `pqs_source_box_low_order_materialization.jl`, net `+24`.

## Pass 281 - Residual-GTO Ham Handoff Artifact Readback

Commit(s):
- this commit - Add H2 PQS residual GTO Ham handoff readback

Summary:
- Extended the H2 PQS residual-GTO Ham artifact with the compact private Ham
  handoff contract needed by downstream consumers: handoff kind/visibility,
  model, orbital basis `[F,R]`, density basis `[P,R]`, `T: [F,R] -> [P,R]`,
  electron count, spin sectors, and nuclear repulsion.
- Reused the existing artifact matrices for `H` in `[F,R]` and `V` in `[P,R]`
  rather than duplicating those large matrices under a second handoff tree.
- Extended the narrow sidecar roundtrip to require and validate the handoff
  fields when `provider_block_mode === :one_body_and_density_provider`,
  including shape/finiteness for the orbital-to-density map and exact labels
  for the private handoff convention.
- No solver execution, CR2 integration, public API, provider registry, or route
  payload tree was added.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- The materialized case still reports final dimension `471`, residual rank
  `18`, augmented dimension `489`, augmented H1-J self-Coulomb
  `0.457435475059184`, private RHF convergence in `15` iterations, and private
  Ham handoff orbital/density dimensions `489`/`489`.

Goal advancement:
- MT4/LT8: upgrades the residual-GTO Ham artifact from route-local fields to a
  small solver-neutral handoff boundary with explicit gauge transform.
- LT5: keeps the critical distinction that the one-body Hamiltonian lives in
  `[F,R]`, while the density interaction lives in `[P,R]` and requires the
  stored orbital-to-density transform.

Medium-goal update:
- The narrow H2 residual-GTO Ham handoff artifact/readback goal is complete for
  the private H2 lane. Broader consumer coverage and provider-block artifact
  coverage remain separate work.

Risk / guardrail:
- This remains private experimental data. It is not an HFDMRG integration and
  not a public Hamiltonian API. The stored handoff fields must not grow into a
  route payload dump.

Remaining blocker / next:
- A minimal consumer smoke can now read the artifact and reconstruct invariants
  from `H`, `V`, and `T` without solving. Performance/paydown and helper
  extraction remain open before broadening this lane.

Line-count / complexity note:
- Source diff before this log entry was `57` added / `0` deleted. The added
  lines are explicit artifact readback checks and compact handoff fields; a
  later extraction pass should move this route-specific handoff logic out of
  `pqs_source_box_low_order_materialization.jl` if the file continues to grow.

## Pass 282 - Residual-GTO Handoff Label Paydown

Commit(s):
- this commit - Drop H2 PQS residual GTO handoff labels

Summary:
- Removed label-only fields from the private H2 residual-GTO handoff surface:
  `ham_handoff_kind`, `ham_handoff_visibility`, `ham_handoff_model`, and
  `private_augmented_rhf_kind`.
- The Ham artifact and readback now rely on the concrete consumer facts:
  orbital basis, density basis, orbital-to-density map, electron count, spin
  sectors, nuclear repulsion, and the existing artifact/result discriminators.
- The materialization print surface no longer shows `ham_handoff_kind`; it
  still shows orbital/density basis and dimensions as the human-facing handoff
  facts.
- No numerical path, matrix, driver input, or construction mode changed.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- The materialized case still reports final dimension `471`, residual rank
  `18`, augmented dimension `489`, augmented H1-J self-Coulomb
  `0.457435475059184`, private RHF convergence in `15` iterations, and handoff
  orbital/density dimensions `489`/`489`.

Goal advancement:
- LT5/LT8: reduces route-specific label carrying cost while preserving the
  actual solver-neutral handoff contract.

Medium-goal update:
- none.

Risk / guardrail:
- Keep artifact/result discriminators and mathematical convention labels, but
  do not add new constant warning/schema labels unless they drive behavior or
  protect a real consumer contract.

Remaining blocker / next:
- The next useful step remains a minimal consumer invariant smoke or a
  provider/handoff helper extraction pass; avoid broadening into solver
  execution.

Line-count / complexity note:
- Source/reporting diff before this log entry was `0` added / `22` deleted.

## Pass 283 - Residual-GTO Ham Handoff Consumer Invariant

Commit(s):
- this commit - Add H2 PQS residual GTO handoff consumer invariant

Summary:
- Added a minimal private consumer invariant to the H2 PQS residual-GTO
  sidecar artifact roundtrip.
- The readback now consumes only artifact-level handoff facts: augmented
  one-body `H` in `[F,R]`, augmented density interaction `V` in `[P,R]`, and
  the stored `T: [F,R] -> [P,R]` orbital-to-density map.
- It diagonalizes artifact `H`, maps the lowest orbital through artifact `T`,
  contracts artifact `V`, and checks that the reconstructed one-orbital H1-J
  self-Coulomb matches the stored augmented H1-J diagnostic.
- No RHF solve, DMRG, CR2, public API, provider registry, or route payload tree
  was added.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Direct artifact readback returned consumer H1-J self-Coulomb
  `0.4574354750591831`, delta `9.43689570931383e-16` against the stored
  diagnostic, and orbital-to-density size `(489, 489)`.

Goal advancement:
- MT4/LT8: proves the residual-GTO Ham artifact is usable by a downstream
  consumer without route-internal payload knowledge.
- LT5: reinforces the explicit gauge contract: orbital coefficients in
  `[F,R]` must be mapped by `T` before using the density interaction in
  `[P,R]`.

Medium-goal update:
- The narrow H2 residual-GTO Ham handoff artifact now has both readback checks
  and a minimal consumer invariant. Further consumer work should either be a
  similarly small invariant smoke or an external integration, not more internal
  solver development.

Risk / guardrail:
- This remains a private route-specific consumer smoke. It should not become a
  hidden solver test or a general artifact registry.

Remaining blocker / next:
- Reasonable next work is either route-specific helper extraction/paydown from
  `pqs_source_box_low_order_materialization.jl`, or a performance review of the
  residual-GTO provider construction before broadening.

Line-count / complexity note:
- Source diff before this log entry was `26` added / `0` deleted. The added
  lines are direct artifact-consumer invariant checks inside the existing
  roundtrip helper.

## Pass 284 - Residual-GTO Ham Handoff Decoupling

Commit(s):
- this commit - Decouple H2 PQS residual GTO Ham handoff

Summary:
- Changed the private H2 residual-GTO Ham handoff constructor so it no longer
  depends on the augmented H1-J diagnostic or private RHF smoke.
- The handoff is now built directly from producer-side data: one-body provider
  blocks, density provider blocks, the density interaction transform, and
  nuclear repulsion from route metadata.
- Removed H1-J and private RHF facts from `ham_handoff.diagnostics`; those
  remain separate consumer/audit results in the materialization path.
- Reordered the materialization path so the Ham handoff is produced before the
  private H1-J/RHF diagnostics. No artifact keys, matrix values, or driver
  output facts were intentionally changed.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.

Goal advancement:
- MT4/LT8: makes the Ham handoff a producer contract rather than a wrapper
  around private solver/audit results.
- LT5: keeps the handoff boundary focused on `H`, `V`, `T`, electron/spin
  counts, nuclear repulsion, and minimal dimensions/symmetry diagnostics.

Medium-goal update:
- The next slimming pass can target artifact-field pruning: remove default Ham
  artifact keys that are provider-development decomposition data rather than
  consumer contract data.

Risk / guardrail:
- Private H1-J and RHF facts are still present in the materialization artifact
  and print surface as audits. They should not be allowed to define the Ham
  handoff contract.

Remaining blocker / next:
- Prune intermediate provider block fields from the default Ham artifact and
  readback path while preserving the minimal consumer invariant.

Line-count / complexity note:
- Source diff before this log entry was `14` added / `24` deleted, net `-10`.
