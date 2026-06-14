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
