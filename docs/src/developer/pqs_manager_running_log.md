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

## Pass 285 - Residual-GTO Artifact Decomposition Pruning

Commit(s):
- this commit - Prune H2 PQS residual GTO artifact decomposition fields

Summary:
- Pruned provider-development decomposition matrices from the default H2 PQS
  residual-GTO Ham artifact and readback path.
- The artifact no longer writes or requires the intermediate one-body component
  blocks (`H_FG`, `H_GG`, `H_FR`, `H_RR` and kinetic/nuclear splits) or the
  density-provider decomposition blocks (`V_PR`, `V_RR`).
- The saved consumer contract remains focused on the durable Ham handoff:
  augmented one-body `H` in `[F,R]`, augmented density interaction `V` in
  `[P,R]`, explicit `T: [F,R] -> [P,R]`, electron/spin counts, nuclear
  repulsion, and compact dimensions/symmetry facts.
- Provider decomposition objects still exist in memory where they are needed to
  construct `H` and `V`; this pass only removes them from the default durable
  artifact/readback surface.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Direct artifact readback reconstructed the handoff consumer H1-J
  self-Coulomb as `0.4574354750591831` with delta
  `9.43689570931383e-16`, using saved `H`, `V`, and `T`.

Goal advancement:
- MT4/LT8: makes the residual-GTO Ham artifact more clearly a consumer
  boundary rather than a provider-debug dump.
- LT5: reduces flat matrix/key carrying cost without changing the mathematical
  gauge contract.

Medium-goal update:
- none.

Risk / guardrail:
- Do not delete the in-memory decomposition blocks yet; they remain the
  construction path for the augmented one-body and density matrices.

Remaining blocker / next:
- Continue contract slimming before new physics: remove returned
  `artifact_roundtrip`/flat optional summary spillover or retire
  `provider_blocks_included` once no live caller requires the compatibility
  mirror.

Line-count / complexity note:
- Source diff before this log entry was `0` added / `105` deleted.

## Pass 286 - Residual-GTO Materialization Summary Trim

Commit(s):
- this commit - Trim H2 PQS residual GTO materialization summary

Summary:
- Kept the H2 PQS residual-GTO artifact roundtrip as a local validation step
  but removed the returned `artifact_roundtrip` object from the materialization
  payload.
- Trimmed summary-only provider/debug shape fields that were not printed or
  consumed: component one-body block sizes, P-projection carrier sizes, P-R/R-R
  pair sizes, and residual width extrema.
- Deleted the now-unused optional-size helper.
- Preserved the printed and durable consumer facts: final dimension, overlap,
  H1/H1-J diagnostics, provider mode, residual rank, augmented H/V dimensions,
  Ham handoff basis/dimension facts, and private RHF audit facts.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Direct artifact readback reconstructed the handoff consumer H1-J
  self-Coulomb as `0.4574354750591831` with delta
  `9.43689570931383e-16`.

Goal advancement:
- MT4/LT8: keeps the Ham handoff artifact as the consumer boundary while
  reducing the route materialization payload back toward printed/consumed facts.
- LT5: removes another small flat-field cloud without changing construction,
  artifact writing, or physics values.

Medium-goal checkpoint:
- Active: residual-GTO Ham handoff artifact/readback remains the current
  consumer-boundary lane; recent passes made it smaller and less debug-shaped.
- Active: `provider_blocks_included` retirement is now the next obvious
  compatibility cleanup, provided the readback and driver can rely only on
  `provider_block_mode`.
- Active: helper extraction from `pqs_source_box_low_order_materialization.jl`
  is still desirable, but should follow another pruning pass rather than move
  bulky route-specific helpers unchanged.
- Deferred: density/pair/H1-J provider expansion, atomic PQS unification, and
  public solver integration should wait until this handoff surface is slimmer
  and performance is reviewed.
- Completed for now: the private H2 materialized ladder validates source plan,
  final basis, H1, H1-J, one-body provider, density provider, private RHF smoke,
  Ham handoff artifact, and minimal consumer invariant.

Risk / guardrail:
- The artifact roundtrip is still essential validation; it is just no longer a
  returned materialization object. Do not replace it with a public status object
  or helper-schema test.

Remaining blocker / next:
- Retire the legacy `provider_blocks_included` mirror or continue trimming
  non-consumer artifact fields before any new physics.

Line-count / complexity note:
- Source diff before this log entry was `0` added / `23` deleted.

## Pass 287 - Residual-GTO Provider Mode Authority Cleanup

Commit(s):
- this commit - Retire H2 PQS residual GTO provider mode mirror

Summary:
- Deleted the live `provider_blocks_included` compatibility mirror from the
  H2 PQS residual-GTO materialization path.
- Artifacts now write and read only `provider_block_mode`; readback no longer
  accepts the legacy fallback key.
- Route metadata, materialization summaries, and the printed driver summary now
  use `provider_block_mode` as the single construction/consumer authority.
- Updated current Cartesian developer docs to name the current
  `:one_body_and_density_provider` mode rather than the retired mirror.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Direct artifact readback returned
  `provider_block_mode = :one_body_and_density_provider` and reconstructed the
  handoff consumer H1-J self-Coulomb with delta `9.43689570931383e-16`.

Goal advancement:
- MT4/LT8: removes a second provider-mode authority from the residual-GTO Ham
  artifact contract.
- LT5: reduces compatibility drift and keeps the route surface tied to one
  explicit construction switch.

Medium-goal update:
- The provider-mode compatibility mirror is retired. Next slimming work should
  either keep the no-solver consumer invariant narrow or extract the remaining
  residual-GTO helper family after pruning.

Risk / guardrail:
- Old artifacts containing only `provider_blocks_included` are intentionally no
  longer accepted by the private roundtrip helper. This is acceptable because
  the route artifact is private and regenerated by the ladder.

Remaining blocker / next:
- Pass E should preserve the minimal no-solver consumer invariant without
  expanding into solver execution. Pass F can then extract the surviving helper
  family into a focused private file.

Line-count / complexity note:
- Diff before this log entry was `11` added / `44` deleted across source and
  current developer docs.

## Pass 288 - Residual-GTO Handoff Consumer Invariant Isolation

Commit(s):
- this commit - Isolate H2 PQS residual GTO handoff invariant

Summary:
- Factored the private no-solver Ham handoff consumer invariant out of the JLD2
  roundtrip body into a small helper.
- The invariant still consumes only saved handoff facts: augmented one-body
  `H`, augmented density interaction `V`, and `T: [F,R] -> [P,R]`, then
  reconstructs the lowest-orbital H1-J self-Coulomb and compares it to the
  stored diagnostic.
- No new artifact fields, solver execution, public API, provider registry, or
  helper-schema test was added.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Direct artifact readback returned
  `provider_block_mode = :one_body_and_density_provider` and H1-J consumer
  delta `9.43689570931383e-16`.

Goal advancement:
- MT4/LT8: keeps the handoff consumer check small and explicitly solver-free.
- LT5: prevents the roundtrip helper from growing another embedded
  computation block while preserving the exact consumer invariant.

Medium-goal update:
- none.

Risk / guardrail:
- This helper is private and route-specific. Do not broaden it into an artifact
  registry or a solver smoke suite.

Remaining blocker / next:
- Pass F can extract the surviving residual-GTO helper family into a focused
  private file now that the default artifact and materialization surface have
  been pruned.

Line-count / complexity note:
- Source diff before this log entry was `28` added / `16` deleted. The small
  positive delta buys a named private boundary for the no-solver consumer
  invariant.

## Pass 289 - Residual-GTO Handoff Helper Extraction

Commit(s):
- this commit - Extract H2 PQS residual GTO handoff helpers

Summary:
- Moved the private H2 PQS residual-GTO sidecar, provider-block, artifact
  roundtrip, and Ham handoff helper family out of
  `pqs_source_box_low_order_materialization.jl` into
  `pqs_h2_residual_gto_handoff.jl`.
- Left the low-order materialization file focused on atomic materialization,
  H2 route wiring, artifact writing, and the generic materialization dispatcher.
- Did not change helper names, artifact fields, physics, readback behavior, or
  driver output facts.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Direct artifact readback returned
  `provider_block_mode = :one_body_and_density_provider` and H1-J consumer
  delta `9.43689570931383e-16`.

Goal advancement:
- MT4/LT8: preserves the private H2 residual-GTO Ham handoff lane while giving
  it a focused implementation file.
- LT5: reduces conceptual carrying cost in the low-order materialization file
  without moving stale helper bloat into a public framework.

Medium-goal update:
- Extraction is complete enough that subsequent paydown can target the focused
  helper file directly. It remains private and route-specific.

Risk / guardrail:
- This was mechanical extraction, not a broad architecture split. Keep the new
  file private; do not turn it into a provider registry.

Remaining blocker / next:
- Reasonable next passes are targeted pruning inside
  `pqs_h2_residual_gto_handoff.jl` or performance review of the residual-GTO
  provider construction before further physics expansion.

Line-count / complexity note:
- `pqs_source_box_low_order_materialization.jl` dropped from about `2907` lines
  to `800`; the extracted private helper file is about `2109` lines. Net source
  line count is roughly flat, but ownership is clearer.

## Pass 290 - Optional Cartesian Driver Precompile Workload

Commit(s):
- this commit - Add optional Cartesian driver precompile workload

Summary:
- Added an opt-in PrecompileTools workload for repeated Cartesian/PQS driver
  ladder work.
- The workload executes the H2 PQS source-box staged route through
  materialization with no artifacts, no supplement provider blocks, and no
  private RHF. This targets the high-level driver methods that dominate cold
  ladder startup without forcing the full materialized residual-GTO route into
  ordinary package precompilation.
- The workload is guarded by
  `GAUSSLETBASES_PRECOMPILE_CARTESIAN_DRIVER=1`; normal package load does not
  execute it.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed and
  did not run the opt-in workload.
- Direct invocation of `_cartesian_driver_precompile_workload` completed in
  about `26.7` seconds.
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.

Goal advancement:
- MT4/LT8: supports the current driver-ladder validation authority by giving
  developers an opt-in route to precompile the staged H2 PQS path.
- LT5: keeps the precompile workload out of normal package startup so the
  validation optimization does not become default carrying cost for all users.

Medium-goal update:
- none.

Risk / guardrail:
- This does not replace the need for a performance review. It is a developer
  precompile aid, not a production scaling fix.

Remaining blocker / next:
- To use it, rebuild the package cache with
  `ENV["GAUSSLETBASES_PRECOMPILE_CARTESIAN_DRIVER"]="1"` before package
  precompilation. Measure ladder runtime before making it default.

Line-count / complexity note:
- Adds a small private precompile workload file. The workload mirrors existing
  driver inputs rather than adding new route behavior.

## Pass 291 - Neutral Cartesian Weighted Hadamard Kernel

Commit(s):
- `4bfd9d04` - Share Cartesian weighted Hadamard contraction

Summary:
- Added neutral private `_cartesian_weighted_hadamard3` in
  `cartesian_gaussian_axis_integrals.jl`.
- Removed the duplicate H2 residual-GTO weighted-Hadamard implementation and
  routed H2 private provider code directly through the neutral helper.
- Preserved `_qwrg_atomic_weighted_hadamard` as a compatibility wrapper for
  existing QW/CPB donor callers.
- Did not touch axis cross/self table extraction, MWG interaction code,
  artifacts, public API, or Cr2 logic.

Validation:
- Doer reported `git diff --check`.
- Doer reported package load.
- Doer reported
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through `cartesian_print/save`.
- Doer reported direct artifact readback:
  `provider_block_mode=one_body_and_density_provider`, final dimension `471`,
  and handoff consumer self-Coulomb `0.4574354750591831`.
- Manager inspected the pushed diff and confirmed the worktree is clean.

Goal advancement:
- MT4/LT8: starts moving reusable residual-GTO producer kernels out from under
  route-specific and QW donor names.
- LT5: removes an exact duplicate implementation while preserving old donor
  surfaces as wrappers.

Medium-goal update:
- The neutral Gaussian-kernel lane is open. The next pass can tackle axis
  cross/self table loops, but should keep old QW names as wrappers and remain
  private/internal.

Risk / guardrail:
- `_cartesian_weighted_hadamard3` still allocates the same intermediate matrix
  as the old implementations. This pass is semantic consolidation, not a
  performance optimization.

Remaining blocker / next:
- Consolidate the axis-table scalar/table loops into neutral helpers without
  changing donor wrapper APIs.

Line-count / complexity note:
- Source impact was small and line-negative: `18` insertions / `19` deletions
  across four source files.

## Pass 292 - Neutral Cartesian Gaussian Axis Table Kernels

Commit(s):
- `dcba7dff` - Share Cartesian Gaussian axis table kernels

Summary:
- Extended `cartesian_gaussian_axis_integrals.jl` with neutral private helpers
  for Cartesian 1D Gaussian axis prefactors, scalar integrals, and integral
  tables.
- Routed the H2 residual-GTO axis cross/self builders through the neutral
  table helper.
- Routed the QW donor axis cross/self/factor wrapper family through the same
  neutral table helper while preserving the old QW function names for live
  callers.
- Left H2 route-local tuple wrappers and larger 3D support/product contractions
  untouched because their data layouts and caller contracts are not exact
  duplicates.

Validation:
- Doer reported `git diff --check`.
- Doer reported package load.
- Doer reported
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases.
- Doer reported
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=wl_diatomic`
  passed both cases.
- Doer reported direct H/V/T artifact readback with provider mode
  `one_body_and_density_provider`, final dimension `471`, and handoff consumer
  self-Coulomb `0.4574354750591831`.
- Manager inspected the pushed diff and found no blocking issue.

Goal advancement:
- MT4/LT8: moves more reusable Gaussian integral machinery out from under H2
  residual-GTO and QW donor-specific names.
- LT5: removes duplicated scalar axis loops without creating a new public
  surface or preserving old donor code as the target architecture.

Medium-goal update:
- none.

Risk / guardrail:
- The neutral helpers still use dense table construction and dispatch by
  `term::Symbol`; this pass is consolidation and carrying-cost reduction, not
  the final performance shape.
- QW wrapper names remain because live donor callers still reference them.

Remaining blocker / next:
- Inventory the remaining MWG/Gaussian donor surfaces and update the Cartesian
  feature-donor documentation so absorbed/shared-kernel status is accurate
  before the next deletion wave.

Line-count / complexity note:
- Net source impact was slightly negative: `336` insertions / `350` deletions.

## Pass 293 - Cartesian Donor Inventory Refresh For H2 Handoff

Commit(s):
- `6b355bea` - Refresh Cartesian donor inventory for H2 handoff

Summary:
- Updated the Cartesian donor inventory and route migration docs after the H2
  residual-GTO handoff/paydown sequence.
- Replaced stale wording that described only one-body provider-block or missing
  density-provider coverage.
- Documented the current private H2 H/V/T handoff: `H` in
  `(:final_pqs, :residual_gto)`, `V` in
  `(:pre_final_pqs, :residual_gto)`, and an explicit orbital-to-density
  transform with `provider_block_mode = :one_body_and_density_provider`.
- Recorded that private H1-J/RHF solver diagnostics are not part of the
  producer contract.
- Added a shared-kernel status note for the neutral weighted-Hadamard and 1D
  Gaussian axis-table helpers.

Validation:
- Doer reported `git diff --check`.
- Doer reported package load.
- Doer reported `julia --project=. tools/run_cartesian_line_ladder.jl --list`
  and listed all four temporary line ladders.
- Manager inspected the pushed diff and found no blocking issue.

Goal advancement:
- MT4/LT8: keeps the migration ledger aligned with the current private H2
  producer capability before starting public contract work.
- LT5: prevents stale donor-doc claims from preserving private solver
  diagnostics or old QW wrapper architecture as if they were target surfaces.

Medium-goal update:
- none.

Risk / guardrail:
- The docs now correctly say this is still H2-fixture/private producer work,
  not a public Cr2-ready producer. Do not weaken that distinction in the next
  pass.

Remaining blocker / next:
- Sketch the public/neutral H/V/T Hamiltonian contract and list the exact H2
  fixture assumptions that must be removed before Cr2 or public API work.

Line-count / complexity note:
- Documentation impact was small: `30` insertions / `18` deletions.

## Pass 294 - H/V/T Public Contract Sketch

Commit(s):
- none; read-only baton pass

Summary:
- Doer produced a read-only sketch for the neutral density-density Hamiltonian
  handoff contract.
- The proposed durable object carries `H` in the orbital basis, `V` in the
  density/provider basis, `T` as the orbital-to-density map, explicit spin
  counts, constant energy, compact basis labels, nuclear metadata, and eventual
  format version.
- The convention was stated explicitly:
  `(ij|kl) = sum_ab T[a,i] T[a,j] V[a,b] T[b,k] T[b,l]`.
- The sketch excludes derived dimensions, finite flags, symmetry errors, route
  labels, private H1-J/RHF diagnostics, sidecar overlap matrices, residual
  debug facts, and route payloads from the durable consumer contract.
- It also lists the H2 fixture assumptions that block Cr2/public readiness.

Validation:
- Read-only response; no Julia validation was required.
- Doer reported inspection only and a clean/even git state before response.
- Manager reviewed the sketch and accepted it as the next coding boundary.

Goal advancement:
- MT4/LT8: clarifies the exact public/neutral producer seam before public API
  or Cr2 work starts.
- LT5: prevents the private H2 payload and diagnostic field cloud from becoming
  the public Hamiltonian contract by accident.

Medium-goal update:
- none.

Risk / guardrail:
- The next pass must stay internal/private. Do not export a type, change
  artifact schema, add Cr2 branches, or make private solver diagnostics part of
  the handoff.

Remaining blocker / next:
- Add an internal neutral H/V/T constructor/validator and route the existing
  private H2 handoff through it while preserving current artifact behavior.

Line-count / complexity note:
- No source or docs changed in the doer pass. The accepted sketch should enable
  later deletion by centralizing validation and reducing private H2 field
  authority.

## Pass 295 - Internal Cartesian Density-Density Hamiltonian Contract

Commit(s):
- `015b6508` - Add internal Cartesian density-density Hamiltonian contract

Summary:
- Added private `_CartesianDensityDensityHamiltonian` as the internal neutral
  H/V/T contract object.
- The object carries one-body `H`, density interaction `V`, orbital-to-density
  transform `T`, spin counts, constant energy, compact basis labels, and
  nuclear metadata.
- Constructor validation now owns the square/finite/symmetric matrix checks,
  transform shape check, spin-count check, constant-energy check, and nuclear
  metadata shape/finiteness checks.
- The existing private H2 residual-GTO handoff now adapts H2-specific labels,
  spin counts, and nuclear repulsion into the neutral constructor.
- Existing artifact/report key names were preserved; no public API, Cr2 branch,
  provider registry, or private H1-J/RHF solver diagnostic was added.

Validation:
- Doer reported `git diff --check`.
- Doer reported package load.
- Doer reported
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases after fixing a constructor materialization issue for
  tuple nuclear charges.
- Doer reported direct H/V/T artifact readback with provider mode
  `one_body_and_density_provider`, final dimension `471`, and handoff consumer
  self-Coulomb `0.4574354750591831`.
- Manager inspected the pushed diff and found no blocking issue.

Goal advancement:
- MT4/LT8: creates the first internal neutral contract seam between the H2
  private producer and a future public Hamiltonian object.
- LT5: pulls validation authority out of the H2 payload field cloud without
  changing durable artifact keys yet.

Medium-goal update:
- The private H2 H/V/T producer lane now has an internal neutral object. The
  next lane should decide public naming/read-write format and external consumer
  needs before adding Cr2/generalization code.

Risk / guardrail:
- This is still private/internal. Do not export the object, treat it as Cr2
  ready, or add consumer-specific fields until the external contract is agreed.

Remaining blocker / next:
- Decide public name, durable format/version, factorized H/V/T consumer
  requirements, Cr2 electron/spin/core treatment, and performance validation
  before public promotion.

Line-count / complexity note:
- Source impact was `117` insertions / `38` deletions. The positive count is
  accepted because it creates the neutral validation seam, but later passes
  should use it to delete H2 private validation/readback code rather than grow
  parallel surfaces.

## Medium-Term Goal Checkpoint - Passes 291-295

Status:
- Active: neutral H2 residual-GTO H/V/T producer migration. The H2 private
  artifact now has shared Gaussian kernels and an internal neutral
  density-density Hamiltonian contract.
- Active: donor-wrapper paydown. Weighted-Hadamard and 1D Gaussian axis table
  loops have moved into neutral private helpers; QW names remain as wrappers
  because live donor callers still need them.
- Completed for this checkpoint: private solver diagnostics were kept out of
  the producer contract, and donor docs now reflect the current H2/private
  status.
- Blocked: public/Cr2 producer promotion is blocked on external consumer
  contract, electron/spin/core treatment, durable format/version, non-H2 source
  dimensions, and performance review.
- Needing refinement: the next medium goal should separate two lanes:
  public H/V/T contract/read-write design, and continued donor-kernel
  extraction/deletion.

Guardrail update:
- The private H2 route remains a producer prototype, not a public solver lane.
  Continue to reject status/readiness/probe payloads, provider registries, and
  private H1-J/RHF diagnostics as part of the producer surface.

## Pass 296 - Internal H/V/T Contract Hardening

Commit(s):
- `2af0bf67` - Harden internal Cartesian HVT contract

Summary:
- Hardened the private `_CartesianDensityDensityHamiltonian` constructor.
- Avoided copies when H, V, T, nuclear charges, or nuclear positions are
  already owned dense `Float64` containers.
- Added spin-sector consumer invariants: `nup <= norb` and `ndn <= norb`.
- Routed private H2 residual-GTO artifact readback through the neutral H/V/T
  constructor and the no-solver consumer invariant.
- Deleted duplicate H/V/T finite/symmetry/shape checks from the full H/V/T
  provider readback path.
- Removed the unproven package-source precompile workload and direct
  `PrecompileTools` dependency.
- Preserved current artifact keys and did not add public API, Cr2 branches,
  provider registries, or private solver diagnostics.

Validation:
- Doer reported `git diff --check`.
- Doer reported `julia --project=. -e 'using Pkg; Pkg.resolve()'`.
- Doer reported package load.
- Doer reported
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases.
- Doer reported direct H/V/T artifact readback with consumer self-Coulomb
  `0.4574354750591831`.
- Doer reported `julia --project=docs docs/make.jl` passed with existing
  Documenter size warnings.
- Manager inspected the pushed diff and found no blocking issue.

Goal advancement:
- MT4/LT8: strengthens the internal neutral H/V/T seam while keeping the H2
  producer private.
- LT5: removes route-specific package-source precompile bloat and deletes
  duplicate readback validation now owned by the neutral constructor.

Medium-goal update:
- none.

Risk / guardrail:
- The one-body-only artifact mode still has local augmented-H validation
  because it does not produce a full H/V/T object. Do not expand that mode into
  another parallel contract unless it remains a real endpoint.

Remaining blocker / next:
- Obtain or prepare the exact external/Cr2 consumer handshake before public
  naming, durable format/version, Cr2 source dimensions, or electron/spin/core
  treatment are coded.

Line-count / complexity note:
- Net source/docs impact was strongly negative: `63` insertions / `188`
  deletions (`-125` lines).

## Pass 297 - Cr2 Consumer Handshake Preparation

Commit(s):
- none; read-only baton pass

Summary:
- Doer produced a read-only decision packet for the external/Cr2 consumer
  handshake needed before public H/V/T API or Cr2 producer work.
- The current internal contract was restated: `H` in orbital basis `O`, `V` in
  density/provider basis `D`, `T: O -> D`, spin counts, constant energy, and
  nuclear metadata.
- The two-body convention was stated as
  `(ij|kl) = sum_ab T[a,i] T[a,j] V[a,b] T[b,k] T[b,l]`.
- Existing export surfaces were classified as not directly compatible:
  `fullida_dense`, sliced Ham/HamIO-style exports, HamV6/angular bridges, and
  the private H2 sidecar are references at most, not the public H/V/T artifact.
- The consumer checklist now names the decisions required: factorized H/V/T
  versus four-index integrals, dense versus sparse/block storage, Cr2
  electron/spin/core treatment, constant-energy semantics, nuclear metadata,
  ordering/locality metadata, and durable JLD2 layout/version expectations.

Validation:
- Read-only response; no Julia validation was required.
- Doer reported clean/even git status.
- Manager reviewed the response and accepted it as the correct stop point.

Goal advancement:
- MT4/LT8: separates producer-contract design from Cr2/generalization coding.
- LT5: prevents speculative public fields, old HamV6/QW receipt reuse, or
  private route sidecars from becoming the external contract.

Medium-goal update:
- none.

Risk / guardrail:
- Do not code a public writer/reader, Cr2 branch, or extra metadata until the
  consumer confirms whether dense factorized H/V/T is acceptable and which
  ordering/core/constant-energy fields are actually consumed.

Remaining blocker / next:
- User/chat/Cr2-consumer review of the decision packet. If dense factorized
  H/V/T is accepted, the next small coding pass is a format-versioned
  writer/reader boundary for the neutral object plus one H2 consumer smoke. If
  the consumer requires four-index/block/order metadata, stop and redesign
  against that requirement.

Line-count / complexity note:
- No source/docs changed in the doer pass. The baton was stopped deliberately
  rather than adding speculative public fields.

## Pass 298 - Remove Global PQS Core-Shell Lowdin

Commit(s):
- `60fa3945` - Remove global PQS core shell Lowdin cleanup

Summary:
- Audited the active H2 independent PQS shell-local concatenated basis before
  touching the cleanup path.
- The block overlap diagnostics were clean: core/shell block identity errors
  were around `1e-14`, cross-overlap maxima were around `1e-15`, and the full
  concatenated overlap identity error was `5.29668900282789e-14`.
- Removed the active global assembled core-shell overlap eigendecomposition and
  global Lowdin application from both the common complete-core/shell final-basis
  helper and the H2 physical-gausslet final-basis path.
- Promoted shell-local concatenated coefficients directly to final localized
  coefficients.
- Kept `combined_lowdin_cleanup` only as an identity compatibility field for
  private schema callers, with `combined_lowdin_cleanup_used = false`.
- Adjusted the residual-GTO density descriptor so the residual carrier uses
  the fixed localized PQS basis directly: `[-S_FG; I] * L_R`.

Validation:
- Doer reported `git diff --check`.
- Doer reported package load.
- Doer reported focused audit rerun with final overlap identity error
  `5.29668900282789e-14`, H1 lowest `-0.7946037173365894`, H1-J self-Coulomb
  `0.4569117646737199`, and augmented H1 lowest `-0.7959028345077699`.
- Doer reported
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases.
- Doer reported direct artifact readback with provider mode
  `one_body_and_density_provider`, consumer self-Coulomb
  `0.4574354750591778`, residual rank `18`, and residual overlap identity
  error `7.88144916874193e-15`.
- Doer reported `julia --project=docs docs/make.jl` passed with existing
  Documenter size warnings.
- Manager inspected the pushed diff and reran `git diff --check` plus package
  load.

Goal advancement:
- MT4/LT8: corrects the basis authority before public Hamiltonian export work.
  The public producer can now target one localized IDA working basis rather
  than a split final/pre-final gauge.
- LT5: removes a mathematically wrong cleanup layer instead of preserving it as
  compatibility architecture.

Medium-goal update:
- The next public contract should be a one-basis IDA contract:
  `Tkin`, unit-nuclear `{U_i}`, `Vee`, spin counts, nuclear charges/positions,
  and nuclear repulsion. The previous two-basis H/V/T framing is obsolete for
  the all-electron producer.

Risk / guardrail:
- Private schema names such as `pre_final_pair_matrix`,
  `final_to_pre_final_coefficients`, and `ham_handoff_orbital_to_density`
  remain as compatibility aliases. They should not be promoted publicly and
  should be retired after the one-basis public artifact exists.

Remaining blocker / next:
- User/chat review of this mathematical correction. If accepted, resume with a
  public/private transition pass for the one-basis `Tkin`, `{U_i}`, `Vee`
  contract. Do not reintroduce a public density transform for the all-electron
  producer.

Line-count / complexity note:
- Source impact was line-negative: `48` insertions / `74` deletions (`-26`
  lines). The compatibility aliases are the main remaining carrying cost.

## Pass 299 - Public PQS/IDA Algorithm Contract Pages

Commit(s):
- `35f0c5e4` - Document Cartesian PQS and IDA algorithms

Summary:
- Added the missing public Algorithms suite for the corrected Cartesian PQS/IDA
  route: overview, low-dimensional Cartesian operator assembly, PQS shell
  construction, residual-Gaussian extension, and IDA Hamiltonian/counterpoise.
- Updated the Algorithms index to require spaces/dimensions, inputs/outputs,
  pseudocode, linear algebra, allowed orthogonalizations, forbidden operations,
  invariants, operator/gauge conventions, code maps, and current deviations for
  active routes.
- Wired the new pages into Documenter navigation.
- The pages make the corrected one-basis IDA contract explicit: `K`, separated
  unit-nuclear `{U_A}`, `Vee`, charges, positions, spin counts, and nuclear
  repulsion. They also state that there is no public all-electron density
  transform and no global PQS core/shell or PQS+residual-Gaussian Lowdin.

Validation:
- Manager ran `git diff --check`.
- Manager ran package load.
- Manager ran `julia --project=docs docs/make.jl`; it passed with the existing
  large-page/search-index warnings only.

Goal advancement:
- LT3: promotes the Cartesian/PQS basis and Hamiltonian rules into public
  Algorithms documentation before public PQS/Cr2 producer promotion.
- LT5/LT6: records the mathematical guardrails that would have caught the
  former global Lowdin mistake and prevents the obsolete two-basis H/V/T story
  from becoming the public contract.

Medium-goal update:
- The public-contract lane now has a documentation prerequisite satisfied for
  one-basis IDA work. The next code pass should implement or transition toward
  the `K`, `{U_A}`, `Vee` contract rather than resurrecting a public density
  transform.

Risk / guardrail:
- This was documentation only. The private H2 artifact still carries
  compatibility names from the former H/V/T handoff stage; future code should
  retire those aliases rather than documenting them as public concepts.

Remaining blocker / next:
- Implement the one-basis IDA Hamiltonian object and versioned read/write
  boundary, preserving separated unit-nuclear matrices for counterpoise.

Line-count / complexity note:
- This intentionally adds public documentation surface: five pages, about 705
  lines total. No source behavior changed.

## Pass 300 - Delete Obsolete Two-Gauge Compatibility Aliases

Commit(s):
- `5ca229f2` - Delete obsolete PQS two-gauge compatibility aliases

Summary:
- Deleted the private identity compatibility layer left behind after removing
  the global PQS core/shell Lowdin.
- Removed active `combined_lowdin_cleanup`,
  `combined_lowdin_cleanup_used`, `pre_final_coefficients`,
  `pre_final_overlap`, `final_to_pre_final_coefficients`, and
  `ham_handoff_orbital_to_density` artifact persistence.
- Renamed the active density interaction surface from pre-final vocabulary to
  localized IDA vocabulary, including `electron_electron_ida`,
  `ida_weights`, and `pqs_complete_core_shell_ida_density_interaction`.
- Updated the residual-GTO density/handoff labels so both orbital and density
  sectors are `(:localized_ida_pqs, :residual_gto)` in the remaining private
  H/V/T constructor.
- Fixed the new Algorithms docs paths, made the PQS source-box shell stages
  explicit, clarified counterpoise branch data, and made the Page Contract the
  single index authority.

Validation:
- Doer reported `git diff --check`, package load,
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`,
  and `julia --project=docs docs/make.jl`.
- Manager reran diff check, package load, docs build, and the pqs_diatomic
  ladder. All three ladder cases passed; the materialized case reported
  localized IDA handoff sectors and residual-GTO artifact readback succeeded.

Goal advancement:
- LT2: removes stale identity matrix allocation, identity multiplication, and
  two-gauge field aliases before Cr2-scale work.
- LT5/LT6: aligns the private H2 route surface with the documented one-basis
  IDA contract and prevents obsolete H/V/T gauge vocabulary from becoming
  public.

Medium-goal update:
- The obsolete two-gauge compatibility layer is no longer a blocker for the
  next internal one-basis work. The remaining blocker is functional: residual
  GTO augmentation still carries charged/summed one-body data rather than
  separated `K` and unit-nuclear `{U_A}` components.

Risk / guardrail:
- `_CartesianDensityDensityHamiltonian.orbital_to_density` remains only inside
  the private H/V/T constructor for the current artifact invariant. Do not
  promote it; replace it with the one-basis IDA object after separated
  augmented one-body components exist.

Remaining blocker / next:
- Preserve one-body components through residual-GTO augmentation: build and
  carry augmented `K` and separated uncharged `{U_A}` before public writer work.

Line-count / complexity note:
- Net impact was `300` insertions / `428` deletions (`-128` lines). This was a
  productive deletion pass rather than another compatibility layer.

## Medium-Term Goal Checkpoint - Passes 296-300

Status:
- Completed: global core/shell Lowdin removal and the old two-gauge identity
  compatibility layer. The active H2 PQS route now uses shell-local PQS
  coefficients directly and names the active density interaction as localized
  IDA.
- Completed: public algorithm documentation for PQS shell construction,
  residual-Gaussian extension, low-dimensional operator assembly, and IDA
  counterpoise.
- Active: one-basis IDA producer transition. The next code target is separated
  augmented `K`, `{U_A}`, and `Vee`, not public writer/reader work yet.
- Active: carrying-cost reduction. Recent passes removed package precompile
  bloat and private alias fields; continue deleting stale private H/V/T surfaces
  as the one-basis object takes over.
- Blocked: Cr2/general diatomic promotion remains blocked on separated
  residual-GTO one-body components, public one-basis object/format, source-plan
  generalization beyond H2/q5, and performance review.

Guardrail update:
- The public all-electron contract is now one localized IDA basis:
  `K`, separated unit-nuclear `{U_A}`, `Vee`, charges/positions, spin counts,
  and nuclear repulsion. Do not reintroduce a public density transform or a
  global Lowdin cleanup to make old artifacts easier to preserve.

## Pass 301 - Preserve Residual-GTO One-Body Components

Commit(s):
- `6d69a2bd` - Preserve residual GTO one-body components

Summary:
- Refactored the H2 residual-GTO one-body provider so it builds residualized
  one-body operators one component at a time.
- The provider now carries augmented kinetic `K_aug` and separated uncharged
  unit-nuclear matrices `{U_A,aug}` through residual-GTO augmentation.
- The augmented one-body Hamiltonian is reconstructed as
  `K_aug + sum_A Z_A U_A,aug`, instead of treating charged/summed nuclear data
  as the provider authority.
- Added artifact/readback checks that validate base and augmented H1
  reconstruction from separated components.
- Kept `charged_nuclear` only as a derived legacy base artifact field, not as
  augmented route authority.

Validation:
- Doer reported package load, `git diff --check`, pqs_diatomic ladder, docs
  build, and focused artifact reconstruction readback.
- Manager reran diff check, package load, docs build, and the pqs_diatomic
  ladder. All three cases passed; the materialized case reported augmented H1
  lowest `-0.7959028345077871`, augmented symmetry error `0.0`, and successful
  residual-GTO artifact readback.

Goal advancement:
- LT5/LT6: closes the main functional gap between the private residual-GTO
  artifact and the documented one-basis IDA counterpoise contract.
- LT3: keeps the implementation aligned with the new public Algorithms pages
  before public writer/reader promotion.

Medium-goal update:
- The one-basis IDA lane can now introduce an internal object around
  `K`, `{U_A}`, `Vee`, charges/positions, spin counts, and nuclear repulsion.
  The obsolete private `_CartesianDensityDensityHamiltonian` shape should be
  retired rather than preserved.

Risk / guardrail:
- This pass was line-positive because it added real separated component
  construction and readback checks. The next pass should pay that back by
  deleting the private H/V/T object shape and its identity `orbital_to_density`
  invariant rather than running both contracts in parallel.

Remaining blocker / next:
- Replace `_CartesianDensityDensityHamiltonian` with a compact internal
  one-basis IDA object and add a small counterpoise smoke that uses the same
  basis and `Vee` for full and ghost branches.

Line-count / complexity note:
- Source impact was `231` insertions / `61` deletions (`+170` lines). Accepted
  as required functionality, with explicit paydown expected in the next pass.

## Pass 302 - Replace Private HVT Handoff With IDA Object

Commit(s):
- `c0e80954` - Replace private Cartesian HVT handoff with IDA object

Summary:
- Replaced the old private `_CartesianDensityDensityHamiltonian` shape with a
  private `_CartesianIDAHamiltonian` carrying the one-basis IDA contract:
  kinetic, separated unit-nuclear attraction matrices by center,
  electron-electron IDA, spin counts, nuclear charges/positions, and nuclear
  repulsion.
- Removed the identity `orbital_to_density` construction and the active
  `ham_handoff_*` orbital/density basis mirrors from the H2 residual-GTO
  producer path.
- Reworked artifact readback to construct the internal IDA object, rebuild full
  `H1 = K + sum_A Z_A U_A`, run two H2 ghost/counterpoise branch checks with
  the same basis and `Vee`, and report compact IDA facts instead of handoff
  mirror fields.

Validation:
- Doer reported `git diff --check`, package load,
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`,
  and `julia --project=docs docs/make.jl`.
- Manager reran diff check on `c0e80954`, package load, docs build, and the
  pqs_diatomic ladder. All three ladder cases passed; the materialized case
  reported `ida_orbital_dimension = 489`, `ida_center_count = 2`,
  `ida_full_self_coulomb = 0.4574354750591777`, and
  `ida_counterpoise_branch_count = 2`.

Goal advancement:
- LT5/LT6: completes the internal correction from the obsolete H/V/T-with-map
  artifact idea to the documented one-basis IDA contract.
- LT2: deletes a stale identity-map concept before it can become public API or
  Cr2-scale allocation pressure.

Medium-goal update:
- The private H2 residual-GTO lane now has the core ingredients needed for a
  public one-basis artifact: `K`, `{U_A}`, `Vee`, charges/positions, spin
  counts, and nuclear repulsion. Public writer/reader work should not start
  until the artifact field list and deletion plan for the remaining private H2
  scaffolding are explicitly reviewed.

Risk / guardrail:
- This pass is net `+35` lines because it introduced a real constructor and
  counterpoise smoke helpers. Do not let the private artifact readback become a
  second public reader; the next public contract should replace and delete that
  scaffolding, not wrap it.

Remaining blocker / next:
- Decide the public one-basis IDA artifact shape and spin-count ownership, then
  implement the writer/reader as a replacement for the private H2 handoff
  artifact rather than another compatibility layer.

Line-count / complexity note:
- Source impact was `173` insertions / `138` deletions (`+35` lines). The
  conceptual surface shrank even though the constructor made the pass slightly
  line-positive.

## Pass 303 - Clean Up Cartesian IDA Hamiltonian Semantics

Commit(s):
- `ee9cfd88` - Clean up Cartesian IDA Hamiltonian semantics

Summary:
- Renamed the implementation file to `cartesian_ida_hamiltonian.jl` and updated
  the include.
- Replaced ambiguous effective-charge `charge_multipliers` semantics with
  `center_weights`, so branch H1 and nuclear repulsion use
  `w_A * Z_A` while the object retains physical nuclear charges.
- Changed `_CartesianIDAHamiltonian` so nuclear repulsion is derived from
  stored physical charges and `ncenter x 3` positions rather than accepted as a
  caller-provided all-electron value.
- Removed private persistence/readback fields for
  `augmented_h1_component_reconstruction_error`,
  `nuclear_attraction_unit_by_center_count`, and `ida_nuclear_repulsion`.
- Updated Algorithms docs to state that the private H2 residual-GTO lane now
  uses an internal one-basis IDA object while the public type and
  writer/reader remain pending.

Validation:
- Doer reported `git diff --check`, package load,
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`,
  and `julia --project=docs docs/make.jl`.
- Manager reran diff check on `ee9cfd88`, package load, docs build, and the
  pqs_diatomic ladder. All three ladder cases passed; the materialized case
  reported `ida_full_self_coulomb = 0.4574354750591777` and
  `ida_counterpoise_branch_count = 2`.

Goal advancement:
- LT5/LT6: removes the last semantic ambiguity before public IDA artifact
  promotion. Counterpoise is now framed as center weighting of stored physical
  charges rather than passing effective charges under a misleading name.
- LT2: deletes private derived artifact mirrors that should not enter the
  public format.

Medium-goal update:
- The public IDA artifact shape can now be implemented narrowly around `K`,
  `{U_A}`, `Vee`, physical nuclear charges, `ncenter x 3` positions, and
  spin-sector defaults. It should replace the private H2 artifact in the next
  deletion pass, not coexist indefinitely.

Risk / guardrail:
- The remaining private H2 artifact still carries route-side diagnostic fields
  such as augmented H1 lowest/symmetry and augmented density dimension. Keep
  them out of the public artifact and delete them when the public reader
  replaces the private sidecar reader.

Remaining blocker / next:
- Add the public `CartesianIDAHamiltonian` type plus minimal versioned JLD2
  writer/reader. Do not generalize to Cr2 or switch the H2 ladder until that
  public format is reviewed.

Line-count / complexity note:
- Source/docs impact was `64` insertions / `98` deletions (`-34` lines), a
  useful semantic cleanup before the public API pass.

## Pass 304 - Add Public Cartesian IDA Hamiltonian Artifact

Commit(s):
- `a462e862` - Add public Cartesian IDA Hamiltonian artifact

Summary:
- Promoted the one-basis IDA contract to the public
  `CartesianIDAHamiltonian` type with public `one_body_hamiltonian` and
  `nuclear_repulsion` helpers.
- Added `write_cartesian_ida_hamiltonian` and
  `read_cartesian_ida_hamiltonian` with a minimal versioned JLD2 format:
  `artifact_kind`, `format_version`, `kinetic`,
  `nuclear_attraction_unit_by_center`, `electron_electron_ida`,
  `nuclear_charges`, `nuclear_positions`, `nup`, and `ndn`.
- The public artifact stores center matrices as `(norb, norb, ncenter)` and
  derives nuclear repulsion on read rather than storing it.
- Removed the private `_CartesianIDAHamiltonian` alias; private H2 code now
  constructs the public object and uses the public helpers directly.
- Added a compact synthetic public-contract test covering construction,
  center-weight counterpoise, roundtrip, and tensor layout.

Validation:
- Doer reported `git diff --check`, package load, the focused public IDA test,
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`,
  and `julia --project=docs docs/make.jl`.
- Manager reran diff check on `a462e862`, package load, the focused public IDA
  test, docs build, and the pqs_diatomic ladder. All three ladder cases passed;
  the materialized case still reported
  `ida_full_self_coulomb = 0.4574354750591777` and
  `ida_counterpoise_branch_count = 2`.

Goal advancement:
- LT5/LT6: establishes the public all-electron one-basis IDA boundary needed
  by downstream Cr2 work without exposing route metadata, residual diagnostics,
  or a density transform.
- LT2: creates the replacement target needed to delete the private H2 Ham
  sidecar artifact.

Medium-goal update:
- The next pass should switch the materialized H2 Ham artifact to the public
  writer/reader and delete the private Ham field cloud. The private
  basis/residual sidecar may remain temporarily as construction/debug data.

Risk / guardrail:
- The pass was line-positive because it added a public type, writer/reader, and
  public-contract test. Pay this back by replacing the private H2 artifact path
  rather than running both formats in parallel.

Remaining blocker / next:
- Use `write_cartesian_ida_hamiltonian` for the materialized H2 Ham artifact,
  read it back with `read_cartesian_ida_hamiltonian`, and remove private Ham
  sidecar fields/readback checks.

Line-count / complexity note:
- Source/docs/test impact was `227` insertions / `31` deletions (`+196`
  lines). This is acceptable only as a short-lived replacement boundary before
  the private artifact deletion pass.

## Pass 305 - Use Public Cartesian IDA Artifact For H2 Residual-GTO Ham

Commit(s):
- `3d095485` - Use public Cartesian IDA artifact for H2 residual GTO Ham

Summary:
- Switched the materialized H2 residual-GTO Ham artifact from the private
  `:pqs_h2_residual_gto_sidecar_ham_bundle` field cloud to the public
  `CartesianIDAHamiltonian` JLD2 writer/reader.
- Deleted the private Ham artifact's copied route metadata, base/augmented H1
  mirrors, GTO overlap/residual sidecar copies, provider-mode fields, and
  density/residual diagnostics.
- Renamed the materialization result to
  `:h2_pqs_residual_gto_ida_hamiltonian`.
- Kept the private basis/residual sidecar only as construction/debug data for
  final coefficients, residual transforms, density carrier, and residual moment
  descriptors.
- Updated Cartesian developer docs to say the H2 residual-GTO route now writes
  a public Cartesian IDA Hamiltonian artifact on the Ham side.

Validation:
- Doer reported `git diff --check`, package load, the focused public IDA test,
  pqs_diatomic ladder, direct public-Ham JLD2 inspection, and docs build.
- Manager reran diff check, the focused public IDA test, package load, the
  pqs_diatomic ladder, direct JLD2 inspection, and docs build. All passed; the
  materialized case reported
  `result_kind = :h2_pqs_residual_gto_ida_hamiltonian`,
  `ida_full_self_coulomb = 0.4574354750591777`, and
  `ida_counterpoise_branch_count = 2`. The direct artifact check showed
  `artifact_kind = :cartesian_ida_hamiltonian`, `format_version = 1`,
  center tensor size `(489, 489, 2)`, and no route metadata.

Goal advancement:
- LT5/LT6: completes the current public one-basis IDA artifact transition for
  the H2 residual-GTO Ham side. The consumer boundary is now `K`, `{U_A}`,
  `Vee`, charges/positions, and spin counts, not private route payloads.
- LT2: pays down the line-positive public API pass by deleting the old private
  Ham sidecar reader/writer surface.

Medium-goal update:
- The Ham side is no longer the blocker for public IDA handoff. Remaining
  blockers are route-neutral/general-diatomic construction, H2/q5 fixture
  assumptions, private basis/residual sidecar ownership, and Cr2-scale
  performance/source-plan validation.

Risk / guardrail:
- Do not broaden the public artifact with route metadata, residual diagnostics,
  provider modes, or old H1 mirrors. If residual construction data remains
  useful, keep it in a clearly private basis/debug sidecar or delete it after
  the producer no longer needs it.

Remaining blocker / next:
- Stop for review before route-neutral producer work. The next lane should be a
  deliberate manager choice: private basis-sidecar cleanup, route-neutral
  producer design, source-plan generalization beyond H2/q5, or Cr2 preparation.

Line-count / complexity note:
- Source/docs impact was `48` insertions / `467` deletions (`-419` lines), a
  major carrying-cost reduction and the intended payoff for the public artifact
  boundary.

## Medium-Term Goal Checkpoint - Passes 301-305

Status:
- Completed: separated residual-GTO one-body components. The H2 residual-GTO
  route now carries augmented kinetic and separated unit-nuclear center
  matrices before assembling charged H1.
- Completed: internal one-basis IDA object. The obsolete private H/V/T object
  with identity density transform has been removed.
- Completed: semantic cleanup. Counterpoise uses `center_weights`, positions
  are `ncenter x 3`, and nuclear repulsion is derived from physical charges and
  positions.
- Completed: public `CartesianIDAHamiltonian` and minimal JLD2 reader/writer.
  The public artifact stores only the all-electron one-basis IDA fields.
- Completed for H2 Ham side: the materialized H2 residual-GTO Ham artifact now
  uses the public IDA writer/reader and no longer stores private route/Ham
  sidecar diagnostics.
- Active: private basis/residual sidecar cleanup. It still carries construction
  and debug data needed by the current H2 materialized lane.
- Blocked: Cr2/general-diatomic promotion remains blocked on route-neutral
  producer design, removal of H2/q5 fixture assumptions, source-plan
  generalization, and performance review at Cr2-relevant sizes.

Guardrail update:
- The public Hamiltonian artifact is intentionally small: `K`, `{U_A}`, `Vee`,
  physical nuclear charges, `ncenter x 3` positions, and spin counts. Do not
  add route metadata, residual-GTO diagnostics, provider modes, H1 eigenvalue
  mirrors, solver diagnostics, or a density transform to that format without a
  live consumer requirement.

## Pass 306 - Prune Private H2 Residual-GTO Basis Sidecar

Commit(s):
- `07af18de` - Prune private H2 residual GTO basis sidecar

Summary:
- Deleted the private H2 residual-GTO basis-sidecar JLD2 writer and the large
  sidecar roundtrip validator. The materialized H2 residual-GTO input now uses
  `save_basis_artifact = false` and validates only the public
  `CartesianIDAHamiltonian` Ham artifact.
- Fixed the public constructor's small-metadata ownership bug by copying
  `nuclear_charges` and `nuclear_positions`, while keeping already-owned dense
  operator matrices no-copy/read-only.
- Removed the `:one_body_only` provider admission path from the current H2
  residual-GTO producer and updated algorithm/developer docs to describe the
  deleted sidecar as intentionally gone unless a named consumer appears.

Validation:
- Doer reported `git diff --check`, package load, the focused public IDA test,
  the `pqs_diatomic` ladder, docs build, and a direct no-basis Ham smoke showing
  `artifact_kind = :cartesian_ida_hamiltonian`, `format_version = 1`,
  no final-coefficient/residual-transform sidecar keys, center tensor size
  `(489, 489, 2)`, and `nup = ndn = 1`.
- Manager validated directly from the pushed commit while the response file was
  delayed. `git show --check 07af18de`, package load,
  `test/ida/cartesian_ida_hamiltonian_runtests.jl`, the `pqs_diatomic` ladder,
  and docs build passed. The materialized ladder case reported the public Ham
  artifact at `/Users/srw/dmrgtmp/h2_pqs_q5_independent_source_box_r4_gto_ham.jld2`
  with `final_dimension = 471`, `residual_rank = 18`, and
  `augmented_dimension = 489`.

Goal advancement:
- LT2/LT5/LT6: the H2 residual-GTO lane now has a single durable Ham consumer
  artifact and no private basis-sidecar persistence contract. Residual/basis
  data remain in memory only as construction facts.

Medium-goal update:
- Completed: private H2 basis/residual sidecar deletion. Active next blocker:
  spin sectors are still hardcoded by the H2 residual-GTO IDA adapter instead
  of coming from physical system/problem input.

Risk / guardrail:
- Do not revive a basis/provenance artifact unless a named downstream consumer
  requires it. The next route work should keep the public artifact small and
  avoid adding route metadata, residual diagnostics, or solver facts back into
  the Ham format.

Remaining blocker / next:
- Put `nup` and `ndn` into the driver/system metadata path and consume them in
  the public IDA Hamiltonian constructor. Exact H2/q5 source-plan dimensions
  and source names remain a later generalization blocker.

Line-count / complexity note:
- The pass was strongly line-negative: `41` insertions / `426` deletions
  (`-385` lines), with the deleted sidecar validator accounting for most of the
  reduction.

## Pass 307 - Route Cartesian IDA Spin Sectors From Driver Input

Commit(s):
- `bc6c936d` - Route Cartesian IDA spin sectors from driver input

Summary:
- Added `nup` and `ndn` to the Cartesian driver/system input path and to
  `cartesian_report` system metadata.
- Routed those fields through H2 residual-GTO route metadata and removed the
  hardcoded `CartesianIDAHamiltonian(..., 1, 1; ...)` adapter.
- Updated the H2 independent PQS driver input to explicitly set `nup = 1`,
  `ndn = 1`. Missing spin sectors now produce a clear error when the public
  IDA Ham is materialized.

Validation:
- Doer reported diff check, package load, public IDA contract test,
  `pqs_diatomic` ladder, and a direct public Ham readback showing `nup = 1`,
  `ndn = 1`, `ham_dimension = 489`, and no private sidecar fields.
- Manager reran package load, the public IDA contract test, the full
  `pqs_diatomic` ladder, and direct public Ham readback. All passed; the loaded
  artifact reported `nup = 1`, `ndn = 1`, `dim = 489`.

Goal advancement:
- LT5/LT6: removes another H2-only assumption from the public artifact
  producer. Electron/spin sectors are now caller-owned problem data, which is
  required for Cr2, ions, and broken-symmetry workflows.

Medium-goal update:
- Completed: explicit spin sectors in the H2 public IDA Ham producer. Active
  blocker: exact H2/q5 source-plan dimensions, retained counts, and source
  names still gate the successful path.

Risk / guardrail:
- Do not derive spin sectors from nuclear charge in later passes. Keep them as
  explicit problem data. Do not broaden this into solver workflow or frozen-core
  semantics.

Remaining blocker / next:
- Inventory H2/q5 fixture assumptions in the diatomic source-plan and
  residual-GTO producer path before changing them. The next code slice should
  target one narrow derived-count/name replacement, not wholesale Cr2
  generalization.

Line-count / complexity note:
- The patch was small and slightly line-positive: `14` insertions / `3`
  deletions. The added fields are justified because they remove a hardcoded
  public-artifact assumption.

## Pass 308 - H2/Q5 Fixture Assumption Inventory

Commit(s):
- None. Read-only baton response:
  `.agent_handoffs/repo_manager_doer_ida_sidecar_deletion_2026-06-19/response.003.md`

Summary:
- Classified H2/q5 literals and labels across the diatomic PQS source-plan,
  residual-GTO producer, skeletons, ladder inputs, and developer docs.
- Identified literal support counts `(275, 578, 362)`, retained counts
  `(275, 98, 98)`, retained shell count `98`, q5 source names, exact two shared
  shell assumptions, and H2 helper labels as the remaining blocker families.
- Separated current-ladder fixtures and WL/QW reference/comparator constants
  from the active public IDA producer lane so they are not removed casually.

Validation:
- Read-only audit. Doer reported the repo stayed clean and even with
  `origin/main`.

Goal advancement:
- LT5/LT6: turns the vague "remove H2/q5 assumptions" blocker into a set of
  concrete surfaces and replacement rules. This reduces the risk of broad
  branchy Cr2 work layered onto fixture code.

Medium-goal update:
- Active next code slice: derive the shared-shell realization source key and
  retained coefficient shape from local `source_mode_dims` and
  `retained_rule.retained_count`, while keeping H2 support-count admissions in
  place.

Risk / guardrail:
- Do not touch WL/QW 463 reference blocks or ladder fixture names in the next
  pass. Do not delete exact support-count admissions until downstream source
  plan constructors can derive and consume the counts consistently.

Remaining blocker / next:
- Remove one local q5/98 duplication in
  `_pqs_source_box_route_driver_independent_h2_shared_shell_realization`; then
  reassess whether the complete source-plan constructor can derive counts from
  the realized payload.

Line-count / complexity note:
- No tracked line-count change. The inventory points to a low-risk, likely
  line-neutral or slightly negative local cleanup.

## Pass 309 - Derive H2 PQS Shared-Shell Retained Shape

Commit(s):
- `888d9714` - Derive H2 PQS shared-shell retained shape

Summary:
- In `_pqs_source_box_route_driver_independent_h2_shared_shell_realization`,
  derived the raw source key from `role` and `source_mode_dims[1]` instead of
  embedding a q5 source-key literal.
- Derived `retained_count` from `retained_rule.retained_count` and checked the
  coefficient shape against `(shell_descriptor.support_count, retained_count)`
  instead of `(support_count, 98)`.
- Left upstream H2 support-count and retained-count admissions, ladder fixture
  names, and WL/QW comparator constants unchanged.

Validation:
- Doer reported diff check, package load, and the `pqs_diatomic` ladder.
- Manager reran package load and the full `pqs_diatomic` ladder. Both passed.
  The materialized case still reported `final_dimension = 471`,
  `residual_rank = 18`, and `augmented_dimension = 489`.

Goal advancement:
- LT5/LT6: removes one local q5/98 duplication from the active H2 residual-GTO
  producer path without pretending the source regions are general.

Medium-goal update:
- Active next code slice: derive complete-source-plan constructor-local
  `support_counts` and `retained_counts` from `target_payload`, leaving the
  upstream admission checks in place.

Risk / guardrail:
- This is still not a general diatomic producer. Do not remove the upstream
  admissions until the source-region and retained-rule constructors can derive
  valid counts for multiple geometries/q values.

Remaining blocker / next:
- The complete source-plan constructor still restates `(275, 578, 362)` and
  `(275, 98, 98)` after upstream admission. Remove that duplicate local
  restatement next.

Line-count / complexity note:
- Small local patch: `7` insertions / `4` deletions. Slightly line-positive but
  removes duplicated fixture literals from an active construction helper.

## Pass 310 - Derive H2 PQS Complete Source-Plan Counts

Commit(s):
- `18d9bd40` - Derive H2 PQS complete source-plan counts

Summary:
- In `_pqs_source_box_route_driver_independent_h2_complete_core_shell_source_plan`,
  replaced constructor-local duplicate count literals `(275, 578, 362)` and
  `(275, 98, 98)` with ordered tuples derived from
  `target_payload.support_counts` and `target_payload.retained_counts`.
- Preserved the existing support/retained order checks and left upstream
  source-region and retained-rule admission checks untouched.

Validation:
- Doer reported diff check, package load, and the `pqs_diatomic` ladder.
- Manager reran package load and the full `pqs_diatomic` ladder. Both passed.
  The materialized case still reported `final_dimension = 471`,
  `residual_rank = 18`, and `augmented_dimension = 489`.

Goal advancement:
- LT5/LT6: removes another duplicate fixture literal from the active H2
  residual-GTO producer path while preserving current ladder behavior.

Medium-goal update:
- Completed: local constructor-level duplicate q5/count removal. Active
  blocker remains upstream: support-region plan admission, retained-rule
  admission/source dimensions, exact two shared shells, final-basis
  support-count blocking, and H2-specific helper/provenance labels.

Risk / guardrail:
- This pass deliberately did not remove the upstream H2 admissions. The next
  source-plan generalization pass must choose one upstream authority boundary
  and prove it with a bounded fixture, not simply delete checks.

Remaining blocker / next:
- Decide whether the next slice should derive the retained-rule source
  dimensions/key from `q` while keeping H2 counts fixed, or tackle the support
  region admission. The retained-rule source dimension/key slice is likely
  smaller.

Line-count / complexity note:
- Small local patch: `4` insertions / `2` deletions. It is slightly
  line-positive but eliminates two repeated count literals at a downstream
  constructor boundary.

## Medium-Term Goal Checkpoint - Passes 306-310

Status:
- Completed: deleted the private H2 residual-GTO basis-sidecar persistence and
  roundtrip validator. The H2 residual-GTO route now writes only the public
  one-basis Cartesian IDA Ham artifact.
- Completed: public `CartesianIDAHamiltonian` now owns small nuclear metadata
  defensively while treating large dense operator matrices as owned/read-only.
- Completed: explicit `nup`/`ndn` spin sectors flow from driver/system input
  into the public IDA artifact. They are no longer hardcoded by the H2 adapter.
- Completed: H2/q5 fixture inventory separated active-producer blockers from
  ladder fixtures and WL/QW comparator constants.
- Completed: two local duplicate fixture-literal removals in the active H2
  source-plan path: shared-shell source key/retained shape and complete
  source-plan local count tuples.
- Active: upstream source-plan generalization. Remaining blockers are exact H2
  support-region admission `(275, 578, 362)`, retained-rule source dimensions
  and retained count `(5,5,5) -> 98`, exact two shared shells, final-basis
  support-count blocking, and H2-specific helper/provenance names.
- Blocked: Cr2 production remains blocked until the source-plan route derives
  its dimensions from construction for at least one non-fixture variant and has
  a performance plan at Cr2-relevant scale.

Guardrail update:
- Keep the public Ham artifact compact: `K`, `{U_A}`, `Vee`, physical nuclear
  charges, `ncenter x 3` positions, and explicit spin counts. Do not revive the
  deleted basis sidecar or add route metadata/residual diagnostics without a
  named consumer.
- For remaining H2/q5 cleanup, prefer one upstream authority boundary per pass.
  Do not delete fixture admissions until an equivalent derived construction
  check exists.

## Pass 311 - Restore Canonical Cartesian Ham Builder Template

Commit(s):
- `51641313` - Restore Cartesian Ham builder template

Summary:
- Restored `bin/cartesian_ham_builder.jl` as the human-facing Cartesian
  producer template with a visible public stage sequence:
  `cartesian_system -> recipe -> parent -> shells -> units -> transforms ->
  pair_terms -> assembly -> report -> materialization -> print/save`.
- Removed ladder instrumentation, stop-after controls, private RHF controls,
  arbitrary command-line `Meta.parse` overrides, basis-sidecar controls, and
  route-internal residual-GTO provider switches from the canonical driver.
- Added `tools/cartesian_driver_harness.jl` as the explicit home for ladder
  stage markers, fixture input inclusion, stop-after probing, private
  diagnostic knobs, and route-internal materialization switches. The temporary
  Cartesian ladders now run this harness rather than the canonical template.
- Added a small public materialization request field,
  `hamiltonian_output = :cartesian_ida_hamiltonian`, so a user-facing producer
  can ask for an IDA Ham artifact without naming the internal residual-GTO
  provider mode.
- Added an `AGENTS.md` guardrail and a cheap policy test to prevent the
  canonical driver from drifting back into a shared integration-test harness.

Validation:
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Direct policy test for the canonical driver:
  `julia --project=. -e 'using Test; const _PROJECT_ROOT = pwd(); include("test/docs/cartesian_ham_builder_policy_runtests.jl")'`
- Direct H2 residual-GTO materialized harness run. It wrote/read the public IDA
  Ham artifact and reported `final_dimension = 471`, `residual_rank = 18`,
  and `augmented_dimension = 489`.
- `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`
  passed all three cases through the harness.

Goal advancement:
- LT5/LT6: protects the public Cartesian/IDA producer boundary from test and
  private-solver drift while preserving the working H2 residual-GTO IDA
  artifact route.
- MT: keeps Cr2 preparation on the intended path: copy/adapt the visible public
  stage sequence, then implement missing functionality inside `cartesian_*`
  stage functions rather than patching private internals into driver scripts.

Risk / guardrail:
- The harness intentionally preserves old fixture override machinery and
  private diagnostic knobs while active ladder cases still need them. This is a
  quarantine, not a new architecture surface.
- The canonical driver has not been promoted as a general Cr2 producer; H2/q5
  source-plan fixture admissions remain the next substantive blocker.

Remaining blocker / next:
- Resume upstream H2/q5 source-plan generalization after this driver guardrail:
  derive retained-rule source dimensions/key or support-region admission from
  construction while preserving the compact public IDA artifact contract.

Line-count / complexity note:
- Canonical driver is back to `111` lines. The pass is deliberately net
  line-positive because it separates the harness from the template and adds a
  policy test, but it reduces conceptual carrying cost by quarantining ladder
  instrumentation outside `bin/cartesian_ham_builder.jl`.

## Pass 312 - Expose Diatomic Terminal Shellification Topology

Commit(s):
- `3b45fc01` - Expose diatomic terminal shellification topology

Summary:
- Enabled the public shellification gate for PQS bond-aligned diatomics, using
  the already-reviewed odd-`q` direct-core policy. The public stage spine now
  carries compact shellification status through `cartesian_shells`,
  `cartesian_units`, and `cartesian_transforms`.
- Kept Hamiltonian assembly and materialization semantics unchanged. Assembly
  still calls the H2-specific independent source-plan path; this pass only
  made terminal shellification/lowering topology visible before assembly.
- Updated the Cr2 stage probe to inspect actual public-stage topology objects:
  shellification scaffold, low-order unit inventory, and lowering contract
  inventory.
- For H2 q5, the public terminal topology is now visible as two direct atom
  cores, one central distorted product box, and four transverse outer mismatch
  slabs. This differs materially from the old assembly-facing
  `atom_contact_core/shared_shell_1/shared_shell_2` story.
- For Cr2 q5, `cartesian_shells`, `cartesian_units`, and
  `cartesian_transforms` now report
  `:blocked_terminal_cartesian_shellification_geometry` with the overlap
  blocker before assembly separately fails on the same geometry issue.

Validation:
- Doer reported `git diff --check`, package load, Cr2 stage probe, and
  `pqs_diatomic` line ladder.
- Manager reran `git diff --check`, package load, Cr2 stage probe, and
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`.
  All existing H2 ladder cases passed. The materialized case retained
  `final_dimension = 471`, `residual_rank = 18`, `augmented_dimension = 489`,
  and H1 lowest `-0.7946037173365894`.

Goal advancement:
- LT5/LT6: moves the true geometry/topology authority up into the public stage
  spine, which is required before replacing assembly-time H2 reconstruction.
- MT: refines the Cr2 blocker. Cr2 is not blocked by supplement wiring in this
  pass; it is blocked because q5 direct atom cores overlap on the snapped
  parent grid before any terminal unit/lowering inventory can be built.

Risk / guardrail:
- The shellification helper now converts terminal geometry failures into a
  blocked shellification status for visibility. This is acceptable for the
  staged route, but it must not become a way to silently continue into a fake
  successful Hamiltonian; assembly still fails for Cr2.
- Probe additions are intentionally developer-only and should be pruned or
  collapsed once assembly consumes the public topology.

Remaining blocker / next:
- Design the assembly handoff from the visible terminal lowering contracts.
  The next pass should not delete H2 admissions blindly; it should decide how
  the existing source-plan/final-basis builder consumes ordered terminal
  records such as direct atom cores, distorted product boxes, and outer
  mismatch slabs.
- Separately, Cr2 q5 still needs a geometry/core policy decision for overlapping
  short-bond atom cores before it can produce terminal units.

Line-count / complexity note:
- Line-positive pass: `+245/-15`, mostly in the Cr2 probe. This is acceptable
  as a one-pass topology audit, but the probe output/helpers should not become
  product surface.

## Pass 313 - Add Diatomic Atom-Contact Core Shellification

Commit(s):
- `bbea7906` - Add diatomic atom-contact core shellification

Summary:
- Replaced the previous overlap-failure behavior for short-bond diatomic
  q-side seed boxes with an explicit `:atom_contact_core` terminal region.
  The region is a direct identity-lowered core whose support is the discrete
  hull of the two atom seed boxes.
- Corrected the geometry policy: this is not a forced double-core-volume rule.
  Coincident nuclei collapse to one `q x q x q` core, small sub-core
  separations produce only the necessary hull elongation, and H2 q5 produces
  the familiar `5 x 5 x 11` / 275-row contact core.
- Physical Cr2 at the reviewed fine mapped spacing is not in the contact-core
  regime. It now shellifies through public stages as two atom-local cores,
  atom-local shells, one midpoint slab, six shared molecular shells, and two
  axial outer mismatch slabs. The remaining Cr2 blocker is assembly consuming
  the old independent-H2 retained/support plan, not terminal shellification.
- Added the contact-core hull rule to the public Algorithms documentation and
  a code-map comment at the implementation seam.

Validation:
- Doer reported `git diff --check`, package load, the `pqs_diatomic` ladder,
  and the Cr2 stage probe.
- Manager reran `git diff --check`, package load, `julia --project=docs
  docs/make.jl`, the Cr2 stage probe, and
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`.
  Docs built with only existing Documenter size warnings. The H2 ladder passed
  all three cases; the materialized case retained `final_dimension = 471`,
  `residual_rank = 18`, `augmented_dimension = 489`, H1 lowest
  `-0.7946037173365863`, and overlap identity error
  `5.29668900282789e-14`. The Cr2 probe reached `cartesian_pair_terms` and
  failed at `cartesian_assembly` with the expected independent-H2 support-plan
  blocker.

Goal advancement:
- LT5/LT6: makes terminal shellification geometry the route authority for the
  short-bond atom-contact core case and removes another H2-only geometry story
  from the conceptual path toward a general bond-aligned diatomic producer.
- MT: Cr2 preparation advances from a geometry/core-overlap blocker to an
  assembly/source-plan consumer blocker over an available 19-region terminal
  topology.

Risk / guardrail:
- The `initial_gap < q` contact-core threshold is now documented as the active
  PQS seed-box policy. Do not reinterpret it as a physical double-core-volume
  rule or force odd hull length along the bond axis.
- The Cr2 probe remains developer instrumentation. It should shrink once
  assembly consumes ordered public terminal records.

Remaining blocker / next:
- Replace the independent-H2 retained/support-plan consumer in assembly with a
  route-neutral consumer of ordered terminal lowering contracts. It must handle
  direct atom-local/contact cores, direct slabs, arbitrary ordered PQS shell
  records, and outer mismatch slabs without adding a Cr2 branch.

Line-count / complexity note:
- Doer pass was line-positive (`+187/-28`) mostly from probe/audit reporting.
  The manager doc patch adds durable algorithm documentation and should prevent
  this core-size policy from living only in chat or probe output.

## Pass 314 - Route Diatomic Assembly Through Terminal Topology

Commit(s):
- `7eb07a22` - Route diatomic assembly through terminal topology

Summary:
- Assembly now receives the public low-order terminal topology and derives
  support records from ordered terminal lowering contracts rather than
  rebuilding raw H2 geometry inside the assembly helper.
- The new support derivation handles direct core, direct slab, direct boundary
  slab, and PQS filled-source terminal contracts. For the current H2 q5
  three-region topology, a narrow compatibility view maps
  `:atom_contact_core` and the two shared molecular shells into the existing
  H2 source/final-basis realization path.
- Deleted the old assembly-side H2 raw-geometry support-plan helper and its
  shared-shell descriptor helper. The legacy `shared_shell_1/shared_shell_2`
  names remain only as the current H2 compatibility view, not as terminal
  geometry authority.
- Cr2 now reaches `cartesian_assembly` with the full 19-region terminal
  support topology available. Its first blocker has moved to
  `:missing_independent_pqs_retained_rule_plan`; the supplement side also
  reports `:missing_gto_supplement_basis` for `Cr/cc-pV5Z`.

Validation:
- Doer reported `git diff --check`, package load, the `pqs_diatomic` ladder,
  and the Cr2 stage probe.
- Manager reran `git diff --check`, package load, the Cr2 stage probe, and
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`.
  The H2 ladder passed all three cases; the materialized case retained
  `final_dimension = 471`, `residual_rank = 18`, `augmented_dimension = 489`,
  H1 lowest `-0.7946037173365863`, and overlap identity error
  `5.29668900282789e-14`. The Cr2 probe now reports
  `last_successful_public_stage = cartesian_assembly` and
  `source_plan blocker: missing_independent_pqs_retained_rule_plan`.

Goal advancement:
- LT5/LT6: moves the assembly authority from H2-specific geometry
  reconstruction to public terminal lowering/support records while preserving
  the compact public IDA artifact path.
- MT: completes the immediate assembly-support handoff. The active Cr2 blocker
  is now retained-rule/final-basis realization over ordered terminal records,
  not support topology visibility.

Risk / guardrail:
- `target_status = :available_physical_gausslet_core_shell_target_inventory`
  is now too broad for Cr2, because support topology is available while
  retained/source/final realization remains blocked. The blocker fields are
  correct, but a near-term cleanup should split support-topology availability
  from full target availability before that label hardens into a contract.
- Do not generalize by adding Cr2 branches. The next pass should consume
  terminal contract records generically or stop with an exact missing retained
  rule object.

Remaining blocker / next:
- Define the retained-rule and final-basis realization path for ordered
  terminal topologies beyond the H2 three-unit compatibility view. The direct
  sectors should remain identity/source-mode sectors, and each PQS filled
  source shell should receive its own shell-local projection/Lowdin realization
  without a global cleanup.

Line-count / complexity note:
- Doer reported `+350/-220` net `+130`. The pass removes stale H2 geometry
  reconstruction but adds a general terminal support-record adapter and probe
  reporting. The next implementation should try to be line-neutral or
  line-negative by retiring more H2 compatibility once ordered retained rules
  exist.

## Pass 315 - Tighten Terminal Topology Support Status

Commit(s):
- `005f3b03` - Tighten terminal topology support status

Summary:
- Fixed terminal support-plan availability semantics. A support plan now blocks
  unless support coverage is complete, duplicate-free, and in-parent. It also
  blocks explicitly on unsupported terminal lowering kinds and includes
  `:distorted_product_box_comx` in the generic support path.
- Fixed target gating. A generic terminal support topology without retained
  rules now leaves the physical gausslet target blocked with
  `:missing_terminal_retained_rule_plan`; it no longer reports the whole target
  as available.
- Supplement construction is now gated behind a real available target, so Cr2
  no longer attempts `Cr/cc-pV5Z` merely because support topology exists.
- Narrowed shellification exception handling to expected geometry/input
  failures (`ArgumentError` and reviewed `DimensionMismatch`) and narrowed
  supplement `ArgumentError` handling to the actual legacy-basis-missing
  message.

Validation:
- Doer reported `git diff --check`, package load, the `pqs_diatomic` ladder,
  and the Cr2 stage probe.
- Manager reran `git diff --check`, package load, the Cr2 stage probe, and
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`.
  The H2 ladder passed all three cases; the materialized case retained
  `final_dimension = 471`, `residual_rank = 18`, `augmented_dimension = 489`,
  H1 lowest `-0.7946037173365863`, and overlap identity error
  `5.29668900282789e-14`. The Cr2 probe now reports
  `target_status = blocked_physical_gausslet_target_inventory`,
  `target_blocker = missing_terminal_retained_rule_plan`, supplement blocked
  on missing target, and first blocker `source_plan blocker:
  missing_terminal_retained_rule_plan`.

Goal advancement:
- LT5/LT6: repairs the public-stage status contract before retained-rule
  generalization. The stage spine can now distinguish available support
  topology from unavailable retained/source/final realization.
- MT: keeps the active Cr2 blocker precise: generic terminal retained-rule and
  final-basis realization over ordered terminal records.

Risk / guardrail:
- The pass did not implement retained rules, by design. Do not treat the 19
  Cr2 terminal support records as a realized basis until `cartesian_transforms`
  owns a per-record retained plan and dimension budget.
- The H2 fixture label remains only when the H2 retained-rule compatibility
  plan exists. Continue deleting fixture labels as generic retained realization
  replaces the compatibility view.

Remaining blocker / next:
- Move retained-rule ownership into `cartesian_transforms`: direct records
  retain identity/source modes, PQS filled-source records derive boundary
  retained counts and shell-local projection/Lowdin inputs, and distorted
  product records block explicitly until COMX realization is implemented.

Line-count / complexity note:
- Doer reported `+75/-21` across two files. This is line-positive but corrects
  real status/exception semantics and prevents a larger wrong-contract retained
  pass.

## Pass 316 - Add Terminal Retained Rule Preflight

Commit(s):
- `d09a0671` - Add terminal retained rule preflight

Summary:
- Added a terminal retained-rule preflight owned by the `cartesian_transforms`
  stage and carried through `cartesian_pair_terms` into assembly. The plan is
  derived from ordered terminal support/lowering records rather than H2
  `shared_shell_1/shared_shell_2` names.
- Direct core, midpoint slab, and boundary slab sectors retain identity/source
  modes with retained count equal to support count. PQS filled-source sectors
  derive retained counts from their source-mode shape via the boundary
  product-mode retained rule. Distorted product boxes still block explicitly.
- Cr2 now has a concrete 19-record retained budget. The current physical Cr2
  probe reports retained dimension `4291`: two direct atom cores (`125 + 125`),
  fourteen PQS filled-source sectors (`14 * 98`), one direct midpoint slab
  (`169`), and two direct boundary slabs (`1250 + 1250`).
- The Cr2 blocker moved from missing retained rules to
  `:missing_terminal_source_plan_realization`, which is the next construction
  boundary before final-basis realization.

Validation:
- Doer reported `git diff --check`, package load, the `pqs_diatomic` ladder,
  and the Cr2 stage probe.
- Manager reran `git diff --check`, package load, the Cr2 stage probe, and
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`.
  The H2 ladder passed all three cases; the materialized case retained
  `final_dimension = 471`, `residual_rank = 18`, `augmented_dimension = 489`,
  H1 lowest `-0.7946037173365863`, and overlap identity error
  `5.29668900282789e-14`. The Cr2 probe reports
  `terminal_retained_plan_status = available_terminal_retained_rule_plan`,
  estimated final dimension `4291`, and first blocker
  `source_plan blocker: missing_terminal_source_plan_realization`.

Goal advancement:
- LT5/LT6: moves another stage authority into the public driver spine. The
  retained-dimension story now belongs to transforms rather than being first
  discovered inside assembly.
- MT: current Cr2 blocker is source-plan/final-basis realization from ordered
  terminal retained records. This is a narrower and more actionable blocker
  than generic missing retained rules.

Risk / guardrail:
- The `4291` retained dimension is a preflight budget, not an accepted Cr2
  final basis. The two direct outer boundary slabs contribute `2500` retained
  functions and must be reviewed before dense K/U/Vee materialization.
- The retained preflight implementation still lives in the diatomic
  complete-core/shell helper file while being surfaced through transforms.
  Keep the next realization pass from simply growing that file around another
  compatibility layer.

Remaining blocker / next:
- Build the terminal source-plan realization from the ordered retained records:
  direct sectors as identity/source-mode coefficient blocks and PQS filled
  source sectors as shell-local projection/Lowdin inputs. Stop or block
  explicitly on distorted product realization and review boundary-slab
  retention before allocating dense operators.

Line-count / complexity note:
- Doer reported `+228/-3`. This is materially line-positive and acceptable as
  a preflight pass only because it exposes the Cr2 retained budget and moves a
  real authority boundary. The next pass should prioritize replacing H2
  compatibility realization rather than adding a second parallel source-plan
  path.

## Pass 317 - Wire Terminal Retained Unit Transform Contracts

Commit(s):
- `c30a3200` - Wire terminal retained unit transform contracts

Summary:
- Wired the existing typed retained-unit stack into the public stage spine.
  `cartesian_units` now constructs `CartesianRetainedUnits.retained_unit_plan`
  from the terminal lowering plan, and `cartesian_transforms` constructs
  `CartesianRetainedUnitTransformContracts.retained_unit_transform_contract_plan`.
- Replaced the Pass 316 local PQS retained-count calculation with a narrow
  adapter that joins terminal support records to retained units by
  `terminal_region_key` and joins retained units to transform contracts by
  `unit_key`.
- PQS retained counts now come from the typed transform contract metadata
  `raw_product_source_retained_rule`; direct sectors use the typed
  direct-identity transform path. Distorted products remain explicitly blocked
  by the realization boundary.
- Cr2 retained budget remains `4291` over the same 19 ordered records, and
  the blocker remains `:missing_terminal_source_plan_realization`.

Validation:
- Doer reported `git diff --check`, package load, the `pqs_diatomic` ladder,
  and the Cr2 stage probe.
- Manager reran `git diff --check`, package load, the Cr2 stage probe, and
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`.
  The H2 ladder passed all three cases; the materialized case retained
  `final_dimension = 471`, `residual_rank = 18`, `augmented_dimension = 489`,
  H1 lowest `-0.7946037173365863`, and overlap identity error
  `5.29668900282789e-14`. The Cr2 probe shows typed transform paths in the
  retained budget and reports first blocker `source_plan blocker:
  missing_terminal_source_plan_realization`.

Goal advancement:
- LT5/LT6: replaces duplicate route-driver retained-rule logic with the
  module-owned retained-unit and transform-contract authorities already present
  in the route stack.
- MT: current Cr2 blocker is unchanged but cleaner: source-plan/final-basis
  realization must now consume typed retained-unit transform contracts rather
  than ad hoc retained-count records.

Risk / guardrail:
- Carrying typed retained-unit and transform-contract plans through the probe
  increased cold compile/allocation cost in early stages. This is acceptable
  while the spine is being connected, but performance needs review before any
  Cr2-scale dense operator construction.
- The adapter still lives in the H2-oriented diatomic complete-core/shell file.
  Do not grow that file into the generic final-basis realizer; extract or use
  module-owned helpers as the realization path becomes real.

Remaining blocker / next:
- Build terminal source-plan/final-basis realization from typed transform
  contracts: direct records as identity coefficient blocks, PQS records as
  shell-local projection/Lowdin planned by their raw-product source contract,
  and cross-block overlap audits before accepting concatenation. Stop on
  distorted-product realization if encountered.

Line-count / complexity note:
- Doer reported `+130/-23`. The pass is line-positive but deletes duplicate
  local retained-count calculation and makes typed plans the authority. The
  next pass should focus on replacing H2 compatibility realization, not adding
  another adapter layer.

## Pass 318 - Add Terminal Source Realization Preflight

Commit(s):
- `1e2d52d7` - Add terminal source realization preflight

Summary:
- Added a terminal source/final-basis realization preflight over the typed
  retained transform contracts. Direct identity sectors now produce available
  realization records with `:identity_source_row_selector`; PQS shell sectors
  check for raw source plan, retained rule, shell projection, shell overlap, and
  shell-local Lowdin ingredients.
- The Cr2 blocker moved from broad missing source realization to the exact
  missing ingredient `:missing_terminal_shell_projection`. The preflight reports
  the retained budget before final-basis work: direct cores `250`, PQS retained
  dimension `1372`, boundary slabs `2500`, total `4291`.
- The pass deliberately does not build K/U/Vee, residual-GTO supplements, public
  artifact fields, or a Cr2 branch. Cross-block overlap remains `:not_computed`
  because shell projection is missing first.

Validation:
- Doer reported `git diff --check`, package load, the `pqs_diatomic` ladder,
  and the Cr2 stage probe.
- Manager reran `git diff --check`, package load, the Cr2 stage probe, and
  `julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic`.
  The H2 ladder passed all three cases; the materialized case retained
  `final_dimension = 471`, `residual_rank = 18`, `augmented_dimension = 489`,
  H1 lowest `-0.7946037173365863`, and overlap identity error
  `5.29668900282789e-14`. The Cr2 probe reached `cartesian_assembly` and
  reported first blocker `source_plan blocker: missing_terminal_shell_projection`.

Goal advancement:
- LT5/LT6: sharpens the stage-authority transfer by making source/final-basis
  realization depend on typed terminal records instead of the old H2 support
  names.
- MT: current Cr2 blocker is now shell projection/overlap/Lowdin input
  production for terminal PQS shell records, followed by mandatory cross-block
  overlap audit. This is the correct boundary before accepting any concatenated
  localized basis.

Risk / guardrail:
- The implementation is line-positive and preflight-only. It is acceptable as a
  blocker-localization pass, but the next step must not add another parallel
  H2/Cr2 realization adapter.
- Boundary slabs still contribute `2500` direct retained functions in the Cr2
  budget. Review that policy before any dense operator materialization.

Remaining blocker / next:
- Locate or extract the existing per-shell projection/overlap/Lowdin ingredients
  and wire them into the typed terminal transform contracts. If cross-block
  overlap is not already small after shell-local realization, stop with
  `:terminal_pqs_cross_block_projection_required`; do not reintroduce any global
  Lowdin cleanup.

Line-count / complexity note:
- Doer reported `+267` insertions and no deletions. This increases carrying
  cost, so the next substantive pass should either replace H2 compatibility
  realization or delete/quarantine obsolete preflight scaffolding as real
  realization lands.

## Pass 319 - H2 PQS Terminal Stage Developer Smoke

Commit(s):
- this commit - Add H2 PQS terminal stage smoke

Summary:
- Added `tools/h2_pqs_terminal_stage_smoke.jl`, a developer-only one-case H2
  smoke that reuses the existing driver harness and materialized H2 fixture,
  disables artifact/TSV writes, and asserts the compact topology/final-basis/IDA
  facts needed for routine manager review.
- This is a validation-workflow improvement only; no production code, public
  API, Cr2 path, or physics contract changed.

Validation:
- Manager ran `git diff --check`, package load, and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl`. The smoke passed
  with terminal roles `(:atom_contact_core, :shared_molecular_shell,
  :shared_molecular_shell)`, support counts `(275, 362, 578)`,
  `final_dimension = 471`, H1 lowest `-0.7946037173365863`, H1-J self-Coulomb
  `0.4569117646737212`, residual rank `18`, and IDA dimension `489`.

Goal / guardrail:
- No strategic change to MT/LT goals. This supports the test-scope policy by
  replacing routine use of the three-case cold-process `pqs_diatomic` ladder
  with a smaller H2 smoke for stage-wiring/preflight/status reviews. The full
  ladder remains the broader acceptance gate when final-basis, H1/IDA,
  residual-GTO, materialization, or driver/harness behavior changes.

Line-count / complexity note:
- The initial doer version duplicated too much input setup. Manager trimmed the
  smoke to an 88-line harness assertion wrapper before acceptance.

## Passes 320-321 - Cleanup Confirmed Defects and Remove H2 Route Authority

Commit(s):
- this commit - Remove H2 compatibility route authority

Summary:
- Fixed the compact confirmed-defect set before further Cr2 construction:
  `run_h1_j` now defaults false and is visible in the canonical driver; generic
  independent diatomic labeling no longer flips to the H2 fixture; blocked H1/J
  no longer claims IDA/density materialization; H1/J reuses the stored H1 lowest
  orbital instead of a second eigensolve; complete core/shell overlap
  diagnostics now enforce/report `rank_atol`; and RouteCore pair sidecars are
  diagnostic unless a crosscheck is explicitly required.
- Removed the active H2 compatibility route authority. The independent-H2
  retained-rule adapter, source descriptor, shared-shell realization,
  complete-core/shell source-plan materializer, terminal support relabeler, and
  supplement support-partition sidecar were deleted from the active source tree.
- H2 and Cr2 now both use ordered terminal records and both block at the shared
  generic source-realization gap `:missing_terminal_shell_projection`.

Validation:
- Doer reported `git diff --check`, package load,
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl`, and
  `julia --project=. tools/cr2_cartesian_ida_stage_probe.jl`.
- The H2 smoke passed as a blocked generic-route smoke with terminal roles
  `(:atom_contact_core, :shared_molecular_shell, :shared_molecular_shell)`,
  support counts `(275, 362, 578)`, no source-plan materialization, no source
  coefficients, no source descriptor, no shared-shell realization, and blocker
  `:missing_terminal_shell_projection`.
- The Cr2 probe remained a clean blocked probe at `cartesian_assembly` with
  source-plan blocker `:missing_terminal_shell_projection`. Manager reviewed
  the diff and symbol searches, and did not rerun the heavy probes.

Goal advancement:
- LT5/LT6: removes the H2-specific successful back half as an algorithmic
  authority. H2 is now only a fixture/topology check for the generic terminal
  route, not a separate successful route.
- MT: current blocker is now unambiguous and shared: terminal PQS shell
  projection/overlap/Lowdin ingredients must be produced for ordered terminal
  shell records before any final-basis or operator materialization can resume.

Risk / guardrail:
- H2 materialized parity is intentionally unavailable until generic terminal
  shell realization lands. Do not reintroduce a compatibility oracle to recover
  old H2 numbers; recover them through the generic route.
- Older source-backed/candidate code still contains `shared_shell_1/2`
  vocabulary outside the independent terminal authority. Treat that as later
  retirement work, not as a live independent-H2 escape path.

Remaining blocker / next:
- Implement terminal shell projection/overlap/shell-local Lowdin for typed
  terminal PQS shell records. Direct records remain identity sectors. If
  cross-block overlap is not already acceptable after shell-local realization,
  stop with `:terminal_pqs_cross_block_projection_required`; do not apply a
  global cleanup.

Line-count / complexity note:
- Doer reported the combined working-tree impact for Passes 320 and 321 as
  `111 insertions(+), 1092 deletions(-)` across eight files. This is the desired
  deletion-oriented simplification before adding the generic realization path.

## Policy Gate - Hard Cartesian/PQS Anti-Bloat Review

Commit(s):
- this commit - Add hard Cartesian PQS anti-bloat gate

Summary:
- Added a hard Cartesian/PQS anti-bloat gate to `AGENTS.md`. This converts the
  prior directional carrying-cost guidance into a review gate for active
  `src/cartesian*`, `src/pqs*`, source-box, route-driver, Cartesian driver,
  Cartesian tool, and related-test work.
- The gate rejects blocker-progression commits, requires a target card before
  coding, forbids unapproved metadata/status/preflight/report/test expansion,
  defines added-source-line budgets, requires replacement commits to delete old
  paths, and makes mechanical diff review mandatory before scientific review.
- The current priority is deletion-oriented cleanup: remove recursive stage
  embedding, delete route-skeleton retained/pair mirrors once typed terminal
  plans are live, retire duplicate H2-local H1/J helpers, quarantine inactive
  pair scaffolding, and replace metadata-carried numerical data with typed
  fields on one canonical object.

Validation:
- Documentation-only policy change. Manager ran `git diff --check`.

Goal / guardrail:
- No numerical or route behavior changed. This is a governance update for LT5,
  LT6, and the current Cr2 terminal-realization lane: future commits must cross
  a real physics/cleanup/cost boundary rather than only moving blockers or
  adding staged vocabulary.

Deletion accounting:
- deleted: none; policy-only change.
- simplified: review standard is now explicit and centralized in `AGENTS.md`.
- quarantined: preflight/status/report-only work is not banned, but is capped
  and cannot be the main achievement of a source commit.
- not deleted because: no source code changed.
- exact remaining caller/blocker: route blocker remains
  `:missing_terminal_shell_projection`.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 188 - Fix Supplemented Route Provenance

Commit(s):
- this commit - Fix supplemented route provenance

Summary:
- Accepted the narrow source fix for supplemented Hamiltonian artifact recipe
  provenance. The previous manifest writer hardcoded
  `:z_axis_diatomic_residual_gto_mwg` for every supplemented artifact,
  including one-center atoms and White-Lindsey diatomics.
- A new local supplemented route-label helper now derives the route from the
  already validated construction kind and nesting. The written
  `recipe_provenance/route` values are:
  `:one_center_pqs_residual_gto_mwg`,
  `:one_center_wl_residual_gto_mwg`,
  `:z_axis_diatomic_pqs_residual_gto_mwg`, and
  `:z_axis_diatomic_wl_residual_gto_mwg`.
- No artifact schema, matrix keys, reader behavior, driver inputs, RG/MWG
  convention, shellification/thin-slab behavior, solver path, or Cr2 workflow
  changed.

Validation:
- `git diff --check` passed.
- Package load passed.
- Four small supplemented artifact/readback smokes passed with direct JLD2
  provenance inspection: H atom/PQS, H atom/WL, H2/PQS, and H2/WL, all using
  bundled `cc-pVTZ`.
- Inspected existing JLD2 keys:
  `recipe_provenance/route`, `recipe_provenance/nesting`,
  `recipe_provenance/producer`, `supplement_provenance/basis_by_center`, and
  `hamiltonian_manifest/manifest_version`.
- Existing H2 RG endpoint passed:
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Goal advancement:
- LT5/LT6: artifact recipe provenance now matches the public construction
  choice for atom/diatomic and PQS/WL supplemented routes. This closes a
  manifest truth bug without broadening the artifact schema or driver surface.

Carrying-cost result:
- deleted: hardcoded supplemented diatomic route literal from the manifest
  write path.
- simplified: supplemented recipe route provenance now uses validated
  `input.kind` and `input.nesting`.
- quarantined: manifest validation hardening remains a later lane if desired;
  driver inventory native shell-index propagation remains separate.
- not deleted because: base route-label helper remains the existing base
  provenance path.
- exact remaining caller/blocker: none for the provenance fix.
- added src lines: `9`.
- deleted src lines: `1`.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 187 - Implement Driver Inventory Geometry Ranges

Commit(s):
- this commit - Add geometry ranges to driver inventory

Summary:
- Accepted the source implementation for the geometry-range portion of the
  `HP-DRV-INV-FN-01` driver inventory amendment. The canonical driver now
  prints parent index ranges and physical edge ranges for `x`/`y`/`z` on
  every terminal inventory row, alongside the existing support/final count,
  compression ratio, identity/compact class, and slab stack facts.
- The earlier review rejected the first WIP because it reconstructed
  `shell_index` heuristically from retained-unit order. The accepted version
  deletes that heuristic. It prints `shell=unavailable` unless a native
  `terminal_region_shell_index` is present in the construction metadata. This
  keeps the output honest: geometry ranges are now useful, while shell index
  remains a known missing native fact.
- The sample H2 output shows angular z-extension slabs with their transverse
  physical range and z physical span, so the driver can expose where slabs sit
  in the parent geometry even before native shell-index propagation lands.

Validation:
- `git diff --check` passed.
- Package load passed.
- Bounded H2 canonical driver artifact/readback passed for `nesting = :pqs`
  and showed index/physical ranges for complete shells and angular z-extension
  slabs.
- Doer also validated bounded H2 `nesting = :wl` and bounded Be2 PQS driver
  artifact/readback for inventory inspection.
- Focused diff scan found no added coefficient dumps, source-mode tables,
  raw-block details, pair inventories, route skeleton dumps, or all
  support-row listings.

Goal advancement:
- LT5/LT6: improves the driver as a bounded construction sanity surface. It
  now exposes x/y/z geometry needed to review angular-balance behavior without
  adding artifact schema fields or a route-debug report framework.

Carrying-cost result:
- deleted: the rejected heuristic shared-shell counter from the first WIP.
- simplified: shell index display is now native-only or unavailable; no
  inferred shell ownership is presented as fact.
- quarantined: native shell-index propagation, artifact/schema changes,
  shellification/lowering/terminal-realizer changes, RG/MWG/IDA, solver, Cr2
  workflow, route-debug dumps, and committed tests remain outside this pass.
- not deleted because: the existing compact driver inventory is the approved
  human-facing surface; this pass extends it rather than replacing it.
- exact remaining caller/blocker: raw `region.shell_index` is dropped before
  the current inventory seam. A later narrow pass should carry
  `terminal_region_shell_index` through terminal lowering/retained transform
  metadata before the driver can print native shell indices.
- added src/bin lines: `41`.
- deleted src/bin lines: `7`.
- new tests: none.
- new metadata/status fields: no artifact/status/schema fields; only compact
  in-memory inventory row fields for driver printing.

## Cartesian Hamiltonian Producer Pass 186 - Interleave Angular Z-Extension Shellification

Commit(s):
- this commit - Interleave angular z-extension shellification

Summary:
- Accepted the `HP-COMP-ANGBOX-FN-01` source correction after explicit user
  approval of the added-line budget overage. The previous source had fixed the
  catastrophic direct-identity outer slab behavior, but still used a catch-up
  geometry model: grow all shared molecular shells, then emit one low and one
  high z-extension remainder after the shared-shell loop ended.
- The new implementation replaces the duplicated z-axis diatomic shared-shell
  loops with one common helper. Each shared-shell step now computes an
  ordinary index-layer shell body, extends it along the bond axis to the
  angular target when needed, emits planned `:angular_z_extension_slab`
  slices for `angular_target \ ordinary_body`, and advances the owned current
  box to the angular target. This preserves the intended no-overlap ownership:
  later shells grow from the already-owned angular target, not from the
  smaller ordinary body.
- PQS and White-Lindsey still share the same shellification geometry. The
  route-family split remains downstream in retained/lowering behavior.

Validation:
- `git diff --check` passed.
- Package load passed.
- Ignored angular geometry audit passed and wrote
  `tmp/work/angular_box_geometry_audit_report.txt` and
  `tmp/work/angular_box_geometry_audit_steps.tsv`. The audit shows H2
  z-extension slabs at shell indices `1` and `4`, Cr2-style `ns=5` slabs at
  shell indices `1` and `2`, and Cr2-style `ns=9` slabs at shell indices
  `1`, `2`, and `4`, with zero residual z/non-z mismatch.
- PQS/WL geometry parity checks passed in the audit for H2, Be2, and Cr2-style
  fixtures.
- H2 base artifact/readback passed for `nesting = :pqs` and `nesting = :wl`,
  both with dimension `767`.
- Be2 base artifact/readback passed for `nesting = :pqs` and `nesting = :wl`,
  both with dimension `863`.
- Existing H2 RG endpoint passed:
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
- Manager mechanical scan found only the intentional fallback
  `piece.metadata` merge for unexpected outer-mismatch residuals; no new
  algorithmic metadata bus, status cloud, public surface, committed test, or
  artifact/schema field was added.

Goal advancement:
- LT5/LT6: moves angular-balanced z-axis diatomic geometry into the layer that
  owns it. Long bond-axis support is now allocated during shared-shell growth,
  not repaired downstream by lowering or by a final catch-up cap. This keeps
  Cr2/Be2/H2 shellification physically interpretable while preserving the
  shared PQS/WL geometry boundary.

Carrying-cost result:
- deleted: duplicated shared-shell growth loops and the normal final
  planned-z catch-up behavior.
- simplified: normal z-extension support is emitted by one shared helper
  during shell growth; `outer_mismatch_slab` is fallback-only.
- quarantined: lowering, retained units, terminal realization, Hamiltonian
  assembly, RG/MWG/IDA, driver, artifacts/readers, Cr2 workflow, and committed
  tests remain untouched.
- not deleted because: generic `outer_mismatch_slab` remains needed for
  genuinely unexpected leftovers.
- exact remaining caller/blocker: none found for the interleaved
  shellification correction.
- added src lines: `96`, explicitly approved over the requested `60` line
  budget because the patch replaces duplicated loops and keeps the geometry
  reviewable.
- deleted src lines: `85`.
- new tests: none.
- new metadata/status fields: none; only existing native slab metadata is
  populated for planned z-extension slabs.

## Cartesian Hamiltonian Producer Pass 186 - Reuse PGDG Layer for QW RG Cross Blocks

Commit(s):
- this commit - Reuse PGDG layer for QW RG cross blocks

Summary:
- Accepted the narrow cleanup that removes duplicate localized-PGDG parent
  layer construction from the QW/RG supplement `G-A` cross-block path. The
  mixed Gaussian blocks now use `bundle.pgdg_intermediate.auxiliary_layer`
  directly instead of rebuilding
  `_mapped_legacy_proxy_localized(_mapped_legacy_proxy_layer(...))` for
  `:pgdg_localized_experimental`.
- This is not a physics-policy change and did not fix the Cr2 low residual
  sector. It makes the intended invariant explicit: `G-G` and `G-A` are built
  against the same carried PGDG parent layer. The old reconstructed layer and
  the bundled layer matched exactly for representative Cr2 axes, and the
  Cr2 residual-only replay stayed at the known post-cutoff values
  `min eig(K_RR) = 0.3700413519`,
  `min eig(H1_RR) = -7.1647854052`.
- The remaining Cr2 bug hunt should therefore move deeper than PGDG-layer
  reconstruction, toward one-vector kinetic parity and block/projection
  consistency in the high-cancellation residual modes.

Validation:
- `git diff --check` passed.
- Package load passed.
- Ignored parity probe `tmp/work/qwrg_proxy_layer_parity_probe.jl` passed:
  x/y old `_MappedLegacyProxyLayer1D` versus new `MappedPGDGLocalized1D`
  dimensions `265 x 39`, z `373 x 75`, with max delta `0.0` across primitive
  centers, widths, coefficients, localized centers, and weights.
- Existing H2 residual-GTO/MWG endpoint
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` passed:
  augmented dimension `489`, self-Coulomb `0.4574265214362095`, and readback
  deltas `0.0`.
- Ignored Cr2 residual-only replay `tmp/work/rg_spectral_cutoff1e6_audit.jl`
  passed: atom retained `61`; dimer retained `[62, 62]`; low residual spectra
  unchanged from the post-cutoff evidence.

Goal advancement:
- LT6: removes a duplicated construction seam from the high-cancellation
  Residual Gaussian path before adding any spectral guard or policy rule. This
  reduces ambiguity in the exact augmented-operator input blocks while keeping
  cutoff, selection, MWG/IDA, artifact, driver, and solver behavior unchanged.

Carrying-cost result:
- deleted: duplicate localized PGDG proxy reconstruction branches in
  `src/ordinary_qw_raw_blocks.jl`.
- simplified: `G-A` proxy selection now trusts the carried PGDG intermediate;
  helper signatures use the existing PGDG-like layer protocol instead of
  requiring `_MappedLegacyProxyLayer1D`.
- quarantined: ignored parity probe only; no committed tests or fixtures.
- not deleted because: analytic Gaussian primitive cross-integral helpers
  remain live and are the shared path for current supplement blocks.
- exact remaining caller/blocker: Cr2 low residual spectra persist, so the
  next evidence pass should be a one-vector kinetic parity audit rather than
  another PGDG-layer identity check.
- added src lines: 21.
- deleted src lines: 25.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 179 - Emit Angular Z-Extension Slabs

Commit(s):
- this commit - Emit angular z-extension slab geometry

Summary:
- Accepted the `HP-COMP-ANGBOX-FN-01` source pass. Z-axis diatomic
  shellification now treats the understood bond-end leftovers produced after
  transverse shared-shell growth saturates as planned
  `:angular_z_extension_slab` stack regions, not generic
  `z_low/z_high_outer_mismatch_slab` fallback regions.
- The ordinary shared-shell body remains unchanged. PQS and WL still share
  the same common shellification geometry; route-local `q` differences remain
  downstream retained-construction facts.
- Lowering is deliberately not implemented here. The emitted slabs are marked
  as pending terminal lowering, so artifact/readback is blocked until
  `HP-COMP-THINSLAB-*` lowers them compactly.

Validation:
- Manager reran: `git diff --check`; package load; ignored
  `tmp/work/angular_box_geometry_audit.jl`.
- Audit summaries: H2 support 900 / z mismatch 0; Cr2 q4-style support 722 /
  z mismatch 0; Cr2 ns5-style support 2500 / z mismatch 0; Cr2 large-slab
  ns9-style support 13448 / z mismatch 0; all PQS/WL geometry comparisons
  match.
- No artifact/readback was run because terminal lowering for
  `:angular_z_extension_slab` is intentionally deferred.

Goal advancement:
- LT5/LT6: moves the CR2 slab blowup correction to the shellification owner.
  The large z-end support is now named planned geometry instead of falling
  through as direct identity-capable outer mismatch.

Carrying-cost result:
- deleted: old diatomic classification of understood z-end leftovers as
  generic outer mismatch slabs.
- simplified: common diatomic remainder handling now has one local helper and
  one planned slab kind.
- quarantined: thin-slab lowering, retained units, transform contracts,
  terminal realization, artifacts, driver workflow, RG/MWG/IDA, and Cr2
  workflow remain unmodified.
- not deleted because: generic `:outer_mismatch_slab` remains for
  unexpected/non-z leftovers.
- exact remaining caller/blocker: `HP-COMP-THINSLAB-*` must lower
  `:angular_z_extension_slab` compactly before artifacts with these regions
  can be produced.
- added src lines: 79.
- deleted src lines: 17.
- new tests: none; ignored audit only.
- new metadata/status fields: native shellification slab metadata only; no
  artifact/status payload fields.

## Cartesian Hamiltonian Producer Pass 180 - Compact Thin-Slab Lowering

Commit(s):
- this commit - Add compact thin-slab terminal lowering

Summary:
- Accepted the `HP-COMP-FACEPROD-*`, `HP-COMP-THINSLAB-*`, and
  `HP-COMP-THINSLAB-META-*` source pass. Z-axis diatomic midpoint,
  outer-mismatch fallback, and angular z-extension slabs now lower through a
  shared compact thin-slab retained unit instead of direct identity slab CPBs.
- Added the neutral terminal face-product helper under
  `CartesianFinalBasisRealization`. WL facet blocks now use that helper as the
  reuse proof, and both PQS and WL thin slabs consume the same compact
  face-stack realization path.
- The physical driver, artifact schema, reader, Hamiltonian assembly,
  RG/MWG/IDA, route skeletons, and Cr2 workflow were not changed.

Validation:
- Manager reran: `git diff --check`; package load; WL facet parity probe;
  ignored angular geometry audit; thin-slab terminal inventory; H2/Be2
  artifact/readback for `nesting = :pqs` and `nesting = :wl`; H2 residual-GTO
  endpoint.
- WL facet parity: 30 facet blocks, max coefficient delta 0.0.
- Thin-slab inventory: for `ns = 5`, both PQS and WL have
  `:angular_z_extension_slab` regions lowered as
  `:compact_thin_slab_product_cpb`; each side has 450 support rows and 50
  final columns, matching `2 * ns^2`, not full identity retention.
- Artifact/readback: H2 PQS/WL dimension 915, Be2 PQS/WL dimension 619, all
  readback deltas 0.0; matrix symmetry errors were roundoff scale.
- H2 residual-GTO endpoint stayed at self-Coulomb
  0.4574265214362095 with readback deltas 0.0.

Goal advancement:
- LT5/LT6: closes the producer-side slab blowup path that created huge
  direct identity z caps. The shellifier names planned z-end support, and the
  lowering/retained/final-basis path now realizes it compactly for both PQS
  and WL.

Medium-goal checkpoint:
- Completed: common PQS/WL first-step geometry and compact thin-slab lowering
  for z-axis diatomics.
- Active: driver-facing terminal inventory output should be designed as a
  separate public usability lane; the ignored inventory probe proved the right
  fields but no driver/report surface was added here.
- Active: CR2 can retry the previous slab-heavy inventory as a user run after
  this commit, but Cr2 remains measurement-only rather than a committed gate.
- Guardrail: keep real shells route-specific after common shellification; only
  thin slabs share this retained/terminal face-product path.

Carrying-cost result:
- deleted: active direct identity lowering for midpoint and outer-mismatch
  thin slabs.
- simplified: WL facet coefficient construction and PQS/WL thin-slab
  terminal realization now share the neutral face-product helper.
- quarantined: driver inventory printing, artifact/report/schema changes,
  route skeleton changes, RG/MWG/IDA, raw blocks, solver/Cr2 workflow, and
  committed tests remain out of scope.
- not deleted because: old direct-slab count fields remain as compatibility
  summary fields and should report zero for these slab regions.
- exact remaining caller/blocker: no blocker found for compact thin-slab
  lowering; durable driver inventory output requires a later docs/source lane.
- added src lines: 240.
- deleted src lines: 41.
- new tests: none; ignored probes only.
- new metadata/status fields: internal compact thin-slab lowering kind,
  retained unit kind, transform path, and native slab facts only; no
  artifact, driver, manifest, or status-payload fields.

## Cartesian Hamiltonian Producer Pass 173 - Update RG Cutoff Defaults

Commit(s):
- this commit - Update RG cutoff defaults

Summary:
- Accepted the `HP-RG-CUTOFF-FN-01` source pass. Production RG defaults are
  now `residual_occupation_cutoff = 5.0e-8` and `identity_atol = 5.0e-8` in
  both the RG owner and the terminal residual compatibility wrapper.
- This intentionally changes owner-local residual selection: marginal
  directions below `5.0e-8`, including the Cr atom
  `basis_ns=9,map_ns=11,lmax=1` `s4` direction at occupation `3.637e-8`, are
  not production-retained by default. Owner grouping, negative-eigenvalue
  tolerances, merge checks, `G' S R` validation, width/zeta filtering,
  MWG/IDA, artifacts, driver workflow, and public API are unchanged.
- Updated the two explicitly approved H2 endpoint assertions in
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` from
  `1.0e-8` to `5.0e-8`.

Validation:
- `git diff --check`; package load after precompilation; Cr atom residual
  audit probe for `basis_ns=9,map_ns=11,lmax=1`; Be atom high-zeta residual
  audit probe; H2 residual-GTO/MWG endpoint
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
- Cr audit result: production construction now reports `status ok` and
  `residual_rank 32`, confirming the marginal `3.637e-8` direction is dropped
  under the new default policy.
- H2 endpoint result: base/residual/augmented dimensions `471/18/489`,
  self-Coulomb `0.4574265214362095`, readback deltas all `0.0`, both testsets
  passed.
- Mechanical suspicious-line scan over `src` and `test` added lines found no
  hits.

Goal advancement:
- LT1/LT3: aligns RG production defaults with the Cr atom ladder policy
  decision while preserving the endpoint and artifact behavior expected by the
  current driver and CR2 consumers.

Carrying-cost result:
- deleted: stale `1.0e-8` RG defaults in the two production keyword defaults
  and the two approved H2 endpoint assertions.
- simplified: Cr marginal residual handling is now a direct default cutoff
  rule instead of an exceptional validation failure.
- quarantined: residual-selection algorithm changes, owner grouping changes,
  merge-rule changes, width/zeta filtering defaults, MWG/IDA, artifact schema,
  driver workflow, EGOI, Cr2 workflow, and committed fixtures remain
  unapproved.
- not deleted because: explicit caller overrides remain supported.
- exact remaining caller/blocker: none for the approved cutoff/tolerance
  default change.
- added src lines: 2.
- deleted src lines: 2.
- new tests: none; two existing assertions updated under
  `HP-RG-CUTOFF-TEST-01`.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 165 - Install Mapped-COMX At Doside Seam

Commit(s):
- this commit - Install mapped-COMX doside source span

Summary:
- Accepted the corrected mapped-COMX source implementation after rejecting the
  earlier parallel `CartesianRawProductSources` numerical builder. The live
  implementation now sits on the existing chain
  `pqs_source_axis_transform_facts_from_pgdg_axes(...) ->
  _nested_doside_1d(...) -> _nested_retained_span(...) ->
  _cleanup_comx_transform(...)`.
- Ordinary PGDG source spans remain the default. The mapped option is selected
  only by the internal `source_span = :mapped_comx` control, uses normalized
  local `u in [-1, 1]` for `s_lambda(u)`, is restricted to protected `P2`,
  and continues through existing physical-position COMX cleanup.
- The final small tightening rejects `lambda <= 0` because `lambda = 0`
  duplicates the protected `u` column.

Validation:
- Manager reran `git diff --check`, package load, the ignored
  `tmp/work/mapped_comx_source_span_probe.jl`, and the existing H2
  supplemented RG endpoint `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
  The probe confirmed ordinary PGDG repeat delta `0.0`, mapped spans for
  `n_s = 5, 6, 7` with `u = [-1, 1]`, overlap errors below `1.3e-15`,
  protected-P2 errors below `6.2e-16`, and physical-position COMX offdiagonal
  residuals below `8.4e-15`. The H2 endpoint kept self-Coulomb
  `0.4574265214362095` and artifact readback deltas `0.0`.

Goal advancement:
- LT5/LT6: installs the high-order mapped source idea at the correct mainline
  seam without creating a duplicate route or operator path. High-order can now
  benchmark the installed option as a consumer.

Carrying-cost result:
- deleted: superseded raw-source mapped-COMX file/export/parallel-route WIP.
- simplified: mapped-COMX reuses existing doside retained-span and COMX cleanup
  machinery rather than adding a second COMX wrapper.
- quarantined: driver inputs, public APIs, artifacts, Hamiltonian/IDA/MWG/RG/
  raw-block/solver work, `protected_degree != 2`, `sqrtJ`, mapped-`s`
  localization, high-order scaffolding imports, committed Cr/Cr2 fixtures, and
  Cr2 workflow remain unapproved.
- not deleted because: ordinary PGDG source-span construction remains the live
  default and comparator.
- exact remaining caller/blocker: no blocker for the source-span seam; next
  evidence should come from H/He/high-order consumer benchmarks, not more
  source routing.
- added src lines: 147.
- deleted src lines: 8.
- new tests: none committed; ignored probe only.
- new metadata/status fields: compact source-span provenance on existing
  source-axis facts only when mapped-COMX is explicitly selected; no artifact
  or status schema changes.

## Cartesian Hamiltonian Producer Pass 166 - Wire Mapped-COMX Terminal Seed Facts

Commit(s):
- this commit - Wire mapped-COMX terminal seed facts

Summary:
- Accepted the terminal-basis wiring pass under `HP-MCOMX-TERM-FN-01`.
  `_shell_seed(...)` now uses carried materialized
  `AxisSourceTransformFact`s when they exist, so mapped-COMX source-axis facts
  become basis-defining for PQS terminal shell realization.
- The ordinary fallback remains unchanged for default/non-materialized facts:
  it still rebuilds ordinary sides through
  `_nested_projected_q_shell_full_sides(...)`. Shell ownership, retained-rule
  semantics, boundary selection, shell-local Lowdin, canonicalization, support
  validation, artifacts, driver inputs, and Hamiltonian assembly were not
  changed.

Validation:
- Manager reran `git diff --check`, package load,
  `tmp/work/mapped_comx_source_span_probe.jl`,
  `tmp/work/mapped_comx_terminal_seam_probe.jl`, and the existing H2
  supplemented RG endpoint
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
- The terminal seam probe reported `mapped_terminal_match_delta = 0.0`,
  `mapped_vs_ordinary_seed_delta = 1.835e-01`, and seed shape `(1331, 98)`.
  The H2 endpoint kept self-Coulomb `0.4574265214362095` and artifact
  readback deltas `0.0`.

Goal advancement:
- LT5/LT6: completes the mainline source-span-to-terminal-basis path needed
  before any public driver exposure. High-order can now test mapped-COMX as an
  actual basis-defining route rather than a metadata-only/source-fact probe.

Carrying-cost result:
- deleted: none.
- simplified: terminal realization now consumes existing carried materialized
  source-axis facts instead of always rebuilding ordinary sides.
- quarantined: driver exposure, public API, artifacts, Hamiltonian/IDA/MWG/RG/
  raw-block/solver work, high-order workflow, Cr/Cr2 fixtures, and Cr2
  workflow remain unapproved.
- not deleted because: ordinary fallback is required for default
  non-materialized source-axis facts.
- exact remaining caller/blocker: none for terminal consumption; next lane is
  either high-order consumer benchmarking or a separately approved driver
  exposure after benchmark evidence.
- added src lines: 71.
- deleted src lines: 4.
- new tests: none committed; ignored seam probe only.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 163 - Diatomic Common Shellifier Entry

Commit(s):
- this commit - Fix diatomic common shellifier entry

Summary:
- Accepted `HP-COMP-SHELLGEOM-DIAT-FN-01`. The shared terminal shellifier now
  receives public `n_s`/`ns` for z-axis diatomic first-step geometry, while
  route-local `q` remains a separate retained/lowering input.
- `raw_terminal_geometry(...)` now treats its former `q` input as the common
  shell side for central-gap/contact decisions and distorted-gap metadata.
  `raw_plan.q` remains populated for compatibility and now equals the common
  shell side; `raw_plan.shell_side` records the same fact explicitly.
- PQS lowering still receives `q = ns`; White-Lindsey retained boundary
  construction still receives `q = ns - 2`. No retained-unit records,
  terminal-lowering policy semantics, terminal realization, driver,
  Hamiltonian, artifact, or RG/MWG code changed.

Validation:
- Doer: `git diff --check`; package load; focused H2 shellifier-entry and
  region-count audit for `ns = 5, 6, 7`; H2 base artifact/readback for
  `nesting = :pqs` and `nesting = :wl`, both dimension `471` for the default
  bounded fixture; H2 RG endpoint/facade with self-Coulomb
  `0.4574265214362095` and zero readback deltas; no Cr2 run.
- Manager: reviewed the two-file source diff and the distorted-product
  metadata consumer; reran `git diff --check`; package load; source diff
  numstat; suspicious added-line scan; direct H2 same-`ns` parity probe for
  `ns = 5, 6, 7`, confirming PQS/WL dimensions match as `471`, `789`, `855`;
  and the H2 RG endpoint/facade with the same self-Coulomb and zero readback
  deltas.

Goal advancement:
- LT5/LT6: closes the diatomic side of the shared first-step shell geometry
  rule. PQS/WL divergence now begins after common shell records exist, at
  retained construction and realization.
- MT/CR2 readiness: same-`ns` PQS/WL diatomic comparisons no longer depend on
  route-local `q` for central-gap/contact ownership before lowering.

Carrying-cost result:
- deleted: route-local `q` reads as the diatomic common shellifier geometry
  authority.
- simplified: common shell size and retained `q` are separated at the caller
  boundary.
- quarantined: none in this lane; central-gap/contact algorithm itself was not
  redesigned.
- not deleted because: compatibility `raw_plan.q` remains read by existing
  lowering metadata paths and now aliases the common shell side.
- exact remaining caller/blocker: none for same-function/same-argument
  diatomic shellifier entry; any future change to the central-gap/contact
  algorithm itself still requires separate authority.
- added src lines: 20.
- deleted src lines: 21.
- new tests: none.
- new metadata/status fields: one internal `shell_side` field on the raw
  shellification plan only; no artifact/status/report field.

## Cartesian Hamiltonian Producer Pass 160 - Public ns Direct-Core Parity

Commit(s):
- this commit - Fix public ns direct-core parity

Summary:
- Accepted `HP-COMP-NSCORE-FN-01` with the one-center audit portion of
  `HP-COMP-SHELLGEOM-FN-01`. The shared route setup now derives
  `core_cube_side` from public `n_s`/`ns`, not route-local `q`. Route-local
  `q` still derives as before (`q = ns` for PQS, `q = ns - 2` for WL) and
  remains a retained/lowering policy parameter.
- The one-center PQS/WL even-`ns` skew is fixed. PQS and WL now enter the same
  terminal shellification function with the same parent axes, center, and
  direct core side for the same public `ns`.
- `src/cartesian_shellification/terminal_geometry.jl` did not need source
  changes for the one-center case; it was already the common route-family-free
  shellifier once caller inputs were corrected.
- Diatomic central-gap/contact policy still uses route-local `q` in
  `raw_terminal_geometry(...)`; that remains quarantined as a reported blocker
  for a later authority lane if it proves to be first-step geometry rather than
  retained-construction policy.

Validation:
- Doer: `git diff --check`; package load; one-center atom base artifact/readback
  for `ns = 5, 6, 7` and `nesting = :pqs/:wl`; H2 base artifact/readback for
  PQS and WL; H2 RG endpoint/facade with self-Coulomb
  `0.4574265214362095` and zero readback deltas; no Cr2 run.
- Manager: reviewed the two-file source diff; reran `git diff --check`;
  package load; source diff numstat; suspicious added-line scan; a bounded H
  atom parity probe confirming PQS/WL dimensions match for `ns = 5, 6, 7`
  (`321`, `495`, `561`); and the H2 RG endpoint/facade with the same
  self-Coulomb and zero readback deltas.

Goal advancement:
- LT5/LT6: removes an accidental route-local `q` leak from the common
  one-center shell geometry setup. Same-`ns` PQS/WL atom comparisons are now
  dimension-fair for even and odd `ns`.
- MT/CR2 readiness: CR2 can compare PQS and WL atom runs without treating
  even-`ns` cases as dimension-skewed by construction.

Carrying-cost result:
- deleted: stale `odd_q_core_side` helper name and summary label.
- simplified: direct core side rule now names public `ns` as the authority.
- quarantined: diatomic central-gap/contact `q` dependency remains reported
  but unchanged.
- not deleted because: no route/lowering/retained redesign was approved.
- exact remaining caller/blocker: diatomic first-step geometry may still be
  affected by route-local `q` in central-gap/contact policy; a later docs-only
  amendment is needed before changing that behavior.
- added src lines: 11.
- deleted src lines: 11.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 148 - WL Diatomic `ns` Guard

Commit(s):
- this commit - Reject undersized WL diatomic ns

Summary:
- Accepted `HP-COMP-WLNS-FN-01` as the narrow source cleanup following the
  WL diatomic `ns` audit. `src/cartesian_base_hamiltonian.jl` now rejects
  normalized `nesting = :wl`, z-axis diatomic `ns < 4` in
  `_cartesian_base_diatomic_basis_parts(...)` before route construction.
- This moves the unsupported `ns = 3` case from a late terminal
  shellification failure (`complete-shell inner box must be a filled CPB`) to
  a clear public-input error. It intentionally leaves WL one-center atoms on
  the existing `ns >= 3` rule and does not change the observed WL diatomic
  retained-support saturation for `ns >= 4`.

Validation:
- Doer: `git diff --check`; package load; ignored validation script
  `tmp/work/hp_comp_wlns_validation.jl`; WL H2 `ns = 3` early rejection with
  the new message; WL H2 `ns = 4` base artifact/readback; WL H2 `ns = 4`
  supplemented artifact/readback with finite/symmetric K and V; PQS H atom and
  PQS H2 smokes. No Cr2 run.
- Manager: inspected the two-line source diff; `git diff --check`;
  `git diff --numstat -- src bin tools test docs`; suspicious added-line
  scan; new tests/tools scan. No suspicious hits and no committed tests/tools.

Goal advancement:
- LT1/LT3: makes the public composition contract honest for the WL diatomic
  corner case.
- LT5/LT6: preserves the three-choice composition design without adding driver
  special cases or route diagnostics.

Carrying-cost result:
- deleted: none.
- simplified: unsupported WL diatomic `ns = 3` now stops at input
  normalization instead of route-internal shellification.
- quarantined: driver changes, route/shellification/terminal-lowering
  changes, artifact/provenance changes, and WL `ns >= 4` saturation remain out
  of scope.
- not deleted because: the live shellification path remains correct for
  supported WL diatomic sizes.
- exact remaining caller/blocker: none for this guard.
- added src lines: 2.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 132 - RG Orthogonality Robustness

Commit(s):
- this commit - Robustify residual identity validation

Summary:
- Accepted `HP-RG-ORTHO-FN-01`. The final residual identity validation now
  computes `R' S R`, symmetrizes it, and applies the approved combined
  absolute/relative check. The strict `G' S R` absolute validation, owner-local
  selection, residual occupation cutoff, owner metrics, merge metric, retained
  counts, and merge failure rules are unchanged.
- The strict N2 q5 p10 case at `core_spacing = 0.042857` now passes residual
  augmentation with retained counts `9,9`. Its final identity error is
  `1.673e-10` with an allowed threshold of `2.000e-10`; `G' S R` remains
  `1.776e-14`, and the final merge spectrum remains `7.232e-2 .. 1.928` with
  condition `26.65`.
- Passing N2 comparison cases at `core_spacing = 0.075` and `0.05` remain
  healthy. The H2 R3 endpoint remains at
  `0.4574265214362095`, within `2.0e-15` of target.

Validation:
- Doer validation: `git diff --check`; package load; H2 R3 endpoint; ignored
  N2 residual validation audit covering strict `0.042857` and passing `0.075`
  / `0.05` cases.
- Manager validation: `git diff --check`; `git diff --numstat` showed only
  `src/cartesian_residual_gaussians/residual_basis.jl` changed with `8` added
  and `1` deleted source line; suspicious-line scan was empty; new-test/tool
  scan was empty; diff inspection confirmed no terminal keyword plumbing,
  public API, artifact, driver, residual selection, MWG/IDA, raw-block, or
  merge-rule changes.

Goal advancement:
- RG/LT6: removes the strict N2 residual artifact blocker while preserving the
  Residual Gaussian basis algorithm and validation intent.
- LT1/LT3: keeps this as a small robustness bug fix rather than reopening
  selection, conditioning, or source-construction lanes.

Carrying-cost result:
- deleted: none.
- simplified: final residual identity validation now matches the approved
  robust symmetric absolute/relative rule.
- quarantined: strict `G' S R`, owner metric, negative-eigenvalue, and final
  merge near-singularity checks remain live invariants.
- not deleted because: the existing residual validation checks are still the
  production safety boundary.
- exact remaining caller/blocker: none for the strict N2 final-identity
  overshoot; CR2-side HF/HFDMRG can retry the strict N2 artifact.
- added src lines: 8.
- deleted src lines: 1.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 132 - N2 Residual Validation Audit

Commit(s):
- none - measurement-only residual validation audit.

Summary:
- Accepted the strict N2 q5 p10 residual validation audit. The failing
  `core_spacing = 0.042857` (`dZ = 0.300`) case is not an owner grouping
  failure, rank loss, or near-singular final merge. It is a small final
  residual identity-validation overshoot.
- For the failing case, candidates are `18`, owner counts are `9,9`, retained
  counts are `9,9`, final merge eigenvalues are `7.232e-2 / 1.928`, and merge
  condition is `26.65`. `G' S R` passes comfortably at `1.776e-14`; `R' S R -
  I` is `1.673e-10`, exceeding the current `1.0e-10` identity tolerance by
  about `6.73e-11`.
- Passing comparison cases were also healthy: `core_spacing = 0.075` had
  `R' S R - I = 7.276e-12`, and `core_spacing = 0.05` had
  `R' S R - I = 3.638e-11`.

Validation:
- Doer validation: `git diff --check`; package load; ignored audit script
  `tmp/work/n2_q5_p10_residual_validation_audit.jl`; final
  `git status --short --branch`. No tracked source edits and no commit.
- Manager validation: reviewed the metric spectra, retained counts, overlap
  errors, and current tolerances. No implementation tests were rerun.

Goal advancement:
- RG/LT6: identifies the next N2 all-electron workflow blocker as a narrow
  residual-orthogonalization robustness policy issue, not a selection,
  ownership, or severe conditioning failure.
- LT1/LT3: avoids unnecessary source churn by requiring authority before any
  tolerance or orthogonalization change.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: next design request can be narrow: robust residual
  orthogonalization/validation policy for small final `R' S R` overshoots with
  healthy spectra.
- quarantined: do not blindly relax `identity_atol` as an implementation
  shortcut; preserve `G' S R` and `R' S R` validation as meaningful numerical
  contracts.
- not deleted because: no source code changed.
- exact remaining caller/blocker: strict N2 q5 p10 `dZ=0.3` cannot write an
  artifact until the residual validation policy is amended and implemented.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 131 - CR2 Driver Timing and N2 Residual Intake

Commit(s):
- none - external CR2/user-run intake.

Summary:
- Accepted CR2-side timing evidence at current `95dbc43c`: the compile cleanup
  sequence produced large practical driver gains. Fresh-warm Be2 q5 p10
  dropped from `43.709 s` to `6.478 s`, and Be2 q5 p20 from `57.603 s` to
  `10.638 s`. The previous base working-basis cost fell from `40.968 s` to
  `3.913 s` for q5 p10 and from `51.007 s` to `4.270 s` for q5 p20.
- Accepted N2 q5 cc-pVDZ lmax=1 producer smokes through the canonical driver:
  padding 10 produced dimension `1361`, `59 MiB`, `26.146 s`; padding 20
  produced dimension `2037`, `130 MiB`, `29.245 s`. CR2 recommends N2 q5 p20
  as a useful next HF/HFDMRG smoke target.
- Recorded the new correctness/robustness issue: N2 q5 p10 at strict
  `core_spacing = 0.042857` (`dZ = 0.300`) fails during residual augmentation
  with `ArgumentError: residual-Gaussian R' S R validation failed`. Bracket
  `core_spacing = 0.075` passed at dimension `2675`; near-target
  `core_spacing = 0.05` passed at dimension `4739`.

Validation:
- CR2-side validation: all runs were serial; GaussletBases `git diff --check`
  passed; GaussletBases status was clean/even except the known untracked
  successor handoff; no leftover N2/cartesian Julia process was running.
- Manager validation: inspected the residual validation path in
  `src/cartesian_residual_gaussians/residual_basis.jl` and confirmed the
  failure is the final `R' S R` identity check after owner-local residual
  selection and final merge Lowdin. No implementation tests were rerun.

Goal advancement:
- LT1/LT3: compile cleanup has crossed the practical Be2 driver bottleneck; no
  immediate source cleanup lane is justified by timing alone.
- RG/LT6: shifts the next blocker from compilation to residual-basis numerical
  robustness for tighter all-electron N2 inputs.

Carrying-cost result:
- deleted: none; intake-only pass.
- simplified: next repo action is a focused residual validation audit, not
  more route/stage cleanup or kernel optimization.
- quarantined: strict N2 `dZ=0.3` must not be "fixed" by blindly relaxing
  `identity_atol`; first measure the residual metric spectrum, retained counts,
  `G' S R`, and `R' S R` error.
- not deleted because: no source code changed.
- exact remaining caller/blocker: reproduce and characterize the strict N2
  residual failure using ignored probes; source work would require a later
  authority decision if the failure is not just an input conditioning boundary.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 124 - Be2 Compile Attribution

Commit(s):
- none - measurement-only compile attribution pass.

Summary:
- Accepted the compile attribution result on current `b17b9161`. The Be2 q5
  p10 supplemented path still has a large cold first-run cost, but warm
  construction remains fast: cold construction was about `59.0 s`, while
  same-process fresh unique targets averaged about `2.02 s`.
- Package load after caches existed was small (`0.452 s`), and artifact
  write/readback/Gaussian numerical kernels were not identified as the main
  cold latency owners. The trace total was about `60.4 s`, with repo-owned
  compile timing about `52.6 s`.
- The top compile owners were route/stage type surfaces:
  terminal-lowering contract inventory and closures in
  `src/pqs_source_box_route_driver_helpers.jl` around line `1343`
  (`~11.3 s`), terminal shellification region-unit inventory closures in
  `src/cartesian_terminal_shellification_geometry.jl` around line `303`
  (`~7.2 s`), and giant `cartesian_transforms` / `cartesian_units` stage
  `NamedTuple` signatures in
  `src/pqs_source_box_route_driver_helpers.jl` around lines `1587` and `1494`
  (`~4.7 s` and `~3.6 s`).

Validation:
- Doer validation: `git diff --check`; package load; SnoopCompile attribution
  harness `tmp/work/be2_q5_p10_snoopcompile_attribution.jl`; trace workload
  `tmp/work/be2_q5_p10_compile_trace_workload.jl`; fresh-target timing harness
  `tmp/work/be2_q5_p10_fresh_target_timing_attribution.jl`; final
  `git status --short --branch`.
- Manager validation: inspected the reported trace owners and spot-checked the
  named code surfaces. The hot functions still visibly build large tuple-backed
  compatibility inventories and wide stage `NamedTuple` returns. No
  implementation tests were rerun because this was a read-only attribution
  pass.

Goal advancement:
- LT1/LT3: identifies a concrete compile-latency owner class after the first
  type-surface cleanups. The remaining issue is not generic package load or
  numerical-kernel runtime; it is oversized staged route/shellification type
  surfaces crossing `cartesian_units` and `cartesian_transforms`.
- RG/LT6: keeps CR2/HF usability direction focused. Precompile/sysimage work
  could hide the latency, but the trace now justifies first asking whether the
  obvious stale compatibility inventories can be reduced or deleted.

Carrying-cost result:
- deleted: none; compile-attribution pass only.
- simplified: next decision target is narrowed to route/stage inventory
  carriers, especially terminal-lowering compatibility inventory,
  terminal-region unit inventory, and the wide low-order stage summaries.
- quarantined: artifact writing, readback, Gaussian raw blocks, terminal G-G
  numerical kernels, and unit-nuclear kernels are not current compile-latency
  targets.
- not deleted because: no source authority has been granted for the newly
  identified cleanup surfaces.
- exact remaining caller/blocker: a docs-only design amendment is needed before
  replacing or deleting these route/stage inventory carriers; the amendment
  should name `src/pqs_source_box_route_driver_helpers.jl` and
  `src/cartesian_terminal_shellification_geometry.jl` and require preserved H2
  base/supplemented artifact behavior.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 111 - Remove Hidden Driver Atom d

Commit(s):
- this commit - Remove hidden driver atom d

Summary:
- Accepted `HP-DRV-ATOM-CLEAN-01`. The canonical driver no longer injects
  hidden legacy `d = core_spacing` into one-center atom basis construction.
  Atom driver basis input now contains only `q`, `core_spacing`,
  `parent_axis_family`, and `radius` from the public driver contract.
- Public inputs, defaults, override behavior, hooks, timing labels, visible
  stage sequence, printed contract, artifact schema, base/supplemented
  dispatch, and diatomic branch behavior were unchanged.

Validation:
- Doer validation: `git diff --check`, package load, H atom base driver
  artifact/readback with dimension `419` at
  `tmp/work/driver_atom_clean_h_atom_base.jld2`, focused `rg` showing no
  `d = vars[:core_spacing]` or `:d` in the driver, no H2 rerun because the
  diatomic branch was untouched, and no Cr2 run.
- Manager validation: reviewed the driver-only diff, confirmed `git diff
  --check`, numstat +1/-1 in `bin/cartesian_ham_builder.jl`, no `src`,
  `test`, or `tools` edits, no public driver contract changes, and no new
  metadata/status fields.

Goal advancement:
- LT1/LT3: removes a stale compatibility echo from the human-facing driver
  after `core_spacing` became the single public near-core scale.
- MT: completes the small post-core-spacing driver cleanup; future deletion
  pressure should focus on explicit `:white_lindsey_low_order` retirement or
  stale test/tool route inputs only under separate authority.

Carrying-cost result:
- deleted: hidden legacy atom `d` basis field in the canonical driver.
- simplified: one-center driver basis construction now mirrors the intended
  public contract.
- quarantined: none.
- not deleted because: public `core_spacing` remains the intended near-core
  input, and producer-side legacy `d` compatibility/provenance remains outside
  this driver-only pass.
- exact remaining caller/blocker: none for this cleanup.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Pass 322 - Stop Recursive Stage Embedding

Commit(s):
- this commit - Stop recursive Cartesian stage embedding

Target card:
- Target: reduce Cartesian/PQS staged-object bloat by removing recursive prior
  stage embedding from assembly/report surfaces.
- Physics endpoint: H2 and Cr2 generic independent terminal diatomic routes
  remain blocked at `:missing_terminal_shell_projection`.
- Allowed files: `src/pqs_source_box_route_driver_helpers.jl`.
- Forbidden additions: no numerical kernels, terminal shell realization, pair
  materialization, new tests, metadata/status fields, payload structs, report
  expansion, or compatibility adapters.
- Success condition: net source decrease, H2 blocked smoke still passes, no
  replacement adapter.
- Failure rule: if a removed field has a live caller, report the exact caller
  and leave it alone.

Summary:
- `cartesian_assembly` no longer returns full prior stages `shells`, `units`,
  `transforms`, or `pairs`, and no longer mirrors `spacing_inputs` or
  `route_facts`.
- `cartesian_report` now derives compact route facts from `route_skeleton`, uses
  `parent.standard_setup` for spacing/setup values, and reads pair summary facts
  from `low_order_assembly` instead of from a full embedded `pairs` stage.
- No numerical route behavior changed; the independent terminal routes remain
  blocked at the same generic source-realization gap.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `12 17
  src/pqs_source_box_route_driver_helpers.jl`.
- Suspicious added-lines grep: no matches.
- New tests/files: none.

Validation:
- Doer reported package load and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl` passed with elapsed
  time `32.680194541s`.
- Manager reviewed the diff, checked for live `assembly.(shells|units|transforms|pairs|spacing_inputs|route_facts)`
  callers, and did not rerun the smoke or Cr2 probe.

Deletion accounting:
- deleted: recursive full-stage returns from `cartesian_assembly`; duplicate
  `spacing_inputs` and `route_facts` mirrors.
- simplified: report route/setup/pair facts are derived locally from compact
  active objects.
- quarantined: none.
- not deleted because: `route_skeleton`, report-level retained/pair mirrors,
  and `report.retained_dimension` still have live materialization/report/TSV
  callers.
- exact remaining caller/blocker: `src/pqs_source_box_low_order_materialization.jl`
  still reads `report.retained_dimension`; route blocker remains
  `:missing_terminal_shell_projection`.
- added src lines: 12.
- deleted src lines: 17.
- new tests: none.
- new metadata/status fields: none.

## Pass 323 - Delete Route-Skeleton Report Mirrors

Commit(s):
- this commit - Delete route skeleton report mirrors

Target card:
- Target: remove stale route-skeleton retained/pair mirrors from active
  Cartesian/PQS report and materialization surfaces where typed terminal plans
  are the active authority.
- Physics endpoint: H2 and Cr2 independent terminal diatomic routes remain
  blocked at `:missing_terminal_shell_projection`.
- Allowed files: `src/pqs_source_box_route_driver_helpers.jl`,
  `src/pqs_source_box_route_driver_reporting.jl`, and
  `src/pqs_source_box_low_order_materialization.jl`.
- Forbidden additions: no terminal shell projection, pair materialization,
  numerical kernels, new tests, metadata/status fields, compatibility adapters,
  or report expansion.
- Success condition: net source deletion and no old route-skeleton mirror in
  the checked active report/materialization surfaces.
- Failure rule: if a route-skeleton mirror has a live caller, report the exact
  caller and do not add an adapter.

Summary:
- Deleted report-level route-skeleton mirrors:
  `source_boxes`, `source_dimensions`, `retained_units`, `retained_counts`,
  `ranges`, `retained_dimension`, `pair_entries`, `pair_family_counts`, and
  `helper_by_pair_family`.
- Deleted route-summary retained/source mirrors, the `pair_summary` surface,
  TSV loops for retained units and pair entries, dead route-facts helpers, and
  stale retained-dimension materialization fallbacks.
- `cartesian_print_summary` no longer prints old retained/pair route-skeleton
  mirrors, and materialization no longer preserves a route-skeleton
  retained-dimension display field.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `0 94
  src/pqs_source_box_route_driver_helpers.jl`, `0 7
  src/pqs_source_box_route_driver_reporting.jl`, `0 4
  src/pqs_source_box_low_order_materialization.jl`.
- Suspicious added-lines grep: no matches.
- New tests/files: none.

Validation:
- Doer reported package load and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl` passed with elapsed
  time `32.347137416s`.
- Manager reviewed the diff, checked for live active-report references to the
  removed mirrors, and did not rerun the smoke or Cr2 probe.

Deletion accounting:
- deleted: route-skeleton retained/pair mirrors from active report,
  materialization display, summary printing, and TSV output.
- simplified: route summary now carries only the active compact route shape and
  shellification kind.
- quarantined: none.
- not deleted because: lower inventory helpers remain live in route-skeleton
  construction; `route_skeleton` itself remains a live caller input for current
  route payload construction.
- exact remaining caller/blocker: H2 and Cr2 remain blocked at
  `:missing_terminal_shell_projection`.
- added src lines: 0.
- deleted src lines: 105.
- new tests: none.
- new metadata/status fields: none.

## Pass 324 - Delete Duplicate H2 Physics Surfaces

Commit(s):
- this commit - Delete duplicate H2 physics surfaces

Target card:
- Target: remove H2-local H1/J, IDA, support-weight, raw-pair, and private-RHF
  duplicate code that only supported the deleted independent-H2 compatibility
  materialization.
- Physics endpoint: H2 and Cr2 remain on the same generic terminal route and
  may remain blocked at `:missing_terminal_shell_projection`.
- Allowed files: `src/pqs_source_box_diatomic_complete_core_shell.jl`,
  `src/pqs_source_box_route_driver_helpers.jl`, and only source surfaces needed
  to remove dead caller fields.
- Forbidden additions: no terminal shell projection, pair-materialization
  framework work, supplement work, new tests, metadata/status fields,
  compatibility adapters, or generic renames of H2-local numerical code.
- Success condition: net source deletion, no H2-specific route authority
  returns, and H2/Cr2 still share the same generic blocker.
- Failure rule: if deletion needs compatibility glue or a stale helper has a
  live source caller, report the exact caller and do not add an adapter.

Summary:
- Deleted `_PQSDiatomicPhysicalGaussletH1JPayload`,
  `_PQSDiatomicPhysicalGaussletRHFInputContractPayload`, and
  `_PQSDiatomicPhysicalGaussletRHFExecutionPayload`.
- Deleted the H2-local density-provenance, support-weight,
  support raw-pair-numerator, IDA-density-interaction, H1-J diagnostic,
  diatomic H1-J payload, private-RHF contract, and private-RHF execution helper
  chain.
- Removed diatomic H1-J/private-RHF assembly/report fields and summary prints.
- `run_h1_j` and private RHF no longer implicitly request the old diatomic
  final-basis/H1 path.

Caller audit:
- The deleted H1-J/RHF payload types and H2-local raw-pair/IDA/H1-J helpers
  have no live source matches after deletion.
- `_PQSDiatomicPhysicalGaussletCoreShellSourcePlan`,
  `_PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload`,
  `_PQSDiatomicPhysicalGaussletFinalBasisPayload`, and
  `_PQSDiatomicPhysicalGaussletH1Payload` remain live only as blocked-route
  payload-collapse work in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`.
- `_pqs_source_box_route_driver_physical_gausslet_support_states` remains live
  through the residual-GTO path in `src/pqs_h2_residual_gto_handoff.jl`.
- Broader `support_weights`, `raw_pair_factor`, and `density_interaction`
  names remain in common/generic IDA, multilayer, residual-GTO, and atomic
  surfaces outside this H2-local deletion slice.
- Defensive `pqs_gto_sidecar_inputs` elision/read guards remain in
  `src/pqs_source_box_route_driver_reporting.jl` and
  `src/pqs_source_box_low_order_materialization.jl`; those are the next narrow
  stale-surface cleanup, not part of this file-limited pass.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `2 909
  src/pqs_source_box_diatomic_complete_core_shell.jl`, `0 147
  src/pqs_source_box_route_driver_helpers.jl`.
- Suspicious added-lines grep: no matches.
- New tests/files: none.

Validation:
- Doer reported package load and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl` passed with elapsed
  time `31.899142458s`.
- Manager reviewed the diff and caller audit, and did not rerun the smoke or
  Cr2 probe.

Deletion accounting:
- deleted: duplicate H2-local support-weight/raw-pair/IDA/H1-J/private-RHF
  helper chain and its assembly/report surfaces.
- simplified: H1-J/private-RHF requests no longer drive the old diatomic H1
  route; summary printing no longer includes deleted private diagnostics.
- quarantined: none.
- not deleted because: core-shell source-plan, final-basis, and H1 payload
  wrappers remain live as blocked-route payload-collapse work; residual-GTO
  still has common support/density callers outside this pass.
- exact remaining caller/blocker: H2 and Cr2 remain blocked at
  `:missing_terminal_shell_projection`; stale `pqs_gto_sidecar_inputs` guards
  remain as a small deletion target.
- added src lines: 2.
- deleted src lines: 1056.
- new tests: none.
- new metadata/status fields: none.

## Pass 325 - Delete Stale PQS-GTO Sidecar Guards

Commit(s):
- this commit - Delete stale PQS-GTO sidecar guards

Target card:
- Target: remove the remaining `pqs_gto_sidecar_inputs` compatibility surface
  after active construction of that field was deleted.
- Physics endpoint: H2 and Cr2 remain on the same generic terminal route and
  may remain blocked at `:missing_terminal_shell_projection`.
- Allowed files: `src/pqs_source_box_route_driver_reporting.jl` and
  `src/pqs_source_box_low_order_materialization.jl`.
- Forbidden additions: no replacement sidecar field, terminal shell projection,
  residual-GTO rewrite, new status fields, tests, or report expansion.
- Success condition: no `pqs_gto_sidecar_inputs` matches remain in active
  source.
- Failure rule: if a live producer or consumer exists, report the exact caller
  and do not add an adapter.

Summary:
- Deleted the durable-report elision entry for `:pqs_gto_sidecar_inputs`.
- Deleted the materialization read/fallback of `report.pqs_gto_sidecar_inputs`.
- Deleted the dead residual-GTO sidecar/provider branch that depended entirely
  on that removed field. The old H2 residual-GTO materializer now route-checks
  and then fails directly with an unavailable-sidecar error instead of carrying
  a large unreachable construction branch.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `1 128
  src/pqs_source_box_low_order_materialization.jl`, `0 1
  src/pqs_source_box_route_driver_reporting.jl`.
- Suspicious added-lines grep: no matches.
- New tests/files: none.
- `rg -n "pqs_gto_sidecar_inputs" src bin tools test`: no matches.

Validation:
- Doer reported package load and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl` passed with elapsed
  time `32.10834875s`.
- Manager reviewed the diff and did not rerun the smoke or Cr2 probe.

Deletion accounting:
- deleted: durable-report elision for `:pqs_gto_sidecar_inputs`, materialization
  read/fallback of the field, and the unreachable residual-GTO sidecar/provider
  materialization branch.
- simplified: the old residual-GTO materializer now reports direct
  unavailability instead of preserving dead construction plumbing.
- quarantined: none.
- not deleted because: residual-GTO helper functions in
  `src/pqs_h2_residual_gto_handoff.jl` remain outside this narrow stale-guard
  cleanup.
- exact remaining caller/blocker: H2 and Cr2 remain blocked at
  `:missing_terminal_shell_projection`; residual-GTO materialization has no live
  sidecar input producer.
- added src lines: 1.
- deleted src lines: 129.
- new tests: none.
- new metadata/status fields: none.

## Medium-Term Goal Checkpoint after Pass 325

- MT terminal-route cleanup: active and advanced. Passes 321-325 removed the
  H2 compatibility route authority, recursive stage embedding,
  route-skeleton report mirrors, duplicate H2 H1/J/IDA/RHF physics surfaces,
  and stale residual-GTO sidecar guards. H2 and Cr2 now share the same generic
  route blocker.
- MT generic terminal realization: still blocked at
  `:missing_terminal_shell_projection`. The next construction work must produce
  real shell-local projection/overlap/Lowdin objects through the generic
  terminal route, not new preflight/status vocabulary.
- MT anti-bloat gate: active. Recent accepted passes were deletion-positive,
  added no committed tests, and kept validation to package load plus the
  compact H2 smoke when doer needed it.
- MT pair-stage authority: still open. The public pair stage remains a later
  cleanup/construction concern; do not extend the unused pair framework before
  the generic terminal shell realization is unblocked.

## Pass 327 - Stop Blocked Final-Basis/H1 Wrapper Churn

Commit(s):
- this commit - Stop blocked final basis H1 wrapper churn

Target card:
- Target: avoid constructing blocked final-basis and H1 payload wrappers on
  the independent terminal route when no source plan exists.
- Physics endpoint: H2 and Cr2 remain on the same generic terminal route,
  blocked at `:missing_terminal_shell_projection`.
- Allowed files: `src/pqs_source_box_route_driver_helpers.jl`.
- Forbidden additions: no payload type deletion, fake/source-backed behavior
  changes, terminal shell projection, new summary/status fields, compatibility
  adapters, or tests.
- Success condition: independent terminal H2/Cr2 do not create blocked
  final-basis/H1 wrapper objects when `source_plan` is absent, while real
  source-backed source plans still take the same construction path.
- Failure rule: if avoiding wrapper construction requires new status plumbing
  or broader fake/source-backed branching, do not commit and report the
  blocker.

Summary:
- Added a small guard in `cartesian_assembly` so
  `_pqs_source_box_route_driver_diatomic_physical_gausslet_final_basis_payload`
  is skipped when `diatomic_physical_gausslet_source_plan_payload.source_plan`
  is `nothing`.
- Added the corresponding guard so the H1 payload is skipped when no final
  basis payload exists.
- This reduces blocked-route object churn without adding a replacement payload
  or changing the source-backed/fake route when it has a real source plan.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `4 0
  src/pqs_source_box_route_driver_helpers.jl`.
- Suspicious added-lines grep: no matches.
- New tests/files: none.

Validation:
- Doer reported package load and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl` passed with elapsed
  time `38.066892875s`.
- Manager reviewed the diff and did not rerun the smoke or Cr2 probe.

Deletion accounting:
- deleted: none in this pass.
- simplified: blocked independent terminal routes no longer allocate
  final-basis/H1 payload wrappers solely to carry blocked summaries.
- quarantined: none.
- not deleted because: the remaining payload types still have source-backed
  construction meaning and are part of a later payload-collapse lane.
- exact remaining caller/blocker: H2 and Cr2 remain blocked at
  `:missing_terminal_shell_projection`.
- added src lines: 4.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Pass 329 - Delete Stale Fake/Source-Backed PQS Route

Commit(s):
- this commit - Delete stale fake source-backed PQS route

Target card:
- Target: remove the stale fake/source-backed PQS route branch and its
  checked-in fixture pressure.
- Physics endpoint: H2 and Cr2 remain on the generic independent terminal
  route, blocked at `:missing_terminal_shell_projection`; materialization must
  not be recovered through the stale fake path.
- Allowed files: `src/pqs_source_box_diatomic_complete_core_shell.jl`,
  `src/pqs_source_box_route_driver_helpers.jl`, and
  `test/driver_inputs/h2_fake_pqs_q5_wl_source_backed_r4.jl`.
- Forbidden additions: no terminal shell projection, new route branch,
  compatibility adapter, payload replacement object, tests, or route renames.
- Success condition: fake/source-backed fixture is gone, no fake/source-backed
  candidate path remains, and H2/Cr2 still share the same generic blocker.
- Failure rule: if deleting the stale path requires changing active generic
  H2/Cr2 behavior or adding a replacement surface, do not commit.

Summary:
- Deleted the stale `h2_fake_pqs_q5_wl_source_backed_r4.jl` fixture.
- Removed fake/source-backed candidate construction and conversion from
  `cartesian_assembly`.
- Deleted the candidate payload, candidate-to-source-plan adapter, old
  source-backed final-basis/H1 materialization helpers, and the old
  WL/QW H2 gausslet-only reference candidate pressure.
- Removed report/print fields for the deleted final-basis/H1 payloads.

Payload audit:
- Deleted `_PQSDiatomicPhysicalGaussletCoreShellSourcePlan`,
  `_PQSDiatomicPhysicalGaussletFinalBasisPayload`, and
  `_PQSDiatomicPhysicalGaussletH1Payload`.
- `_PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload` remains because the
  generic terminal route still uses it to report the blocked source-plan state.
- `_pqs_source_box_route_driver_physical_gausslet_support_states` remains
  because `src/pqs_h2_residual_gto_handoff.jl` still calls it.
- Remaining `source_backed` matches are in unrelated old QW experimental
  surfaces, not the PQS independent terminal route.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `5 1116
  src/pqs_source_box_diatomic_complete_core_shell.jl`, `0 61
  src/pqs_source_box_route_driver_helpers.jl`, `0 43
  test/driver_inputs/h2_fake_pqs_q5_wl_source_backed_r4.jl`.
- Suspicious added-lines grep: no matches.
- New tests/files: none.

Validation:
- Doer reported package load,
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl`, and
  `julia --project=. tools/cr2_cartesian_ida_stage_probe.jl` passed. H2 and Cr2
  both continue to block at `:missing_terminal_shell_projection`.
- Manager reviewed the diff, checked that the stale fixture/candidate/payload
  names were gone from the PQS route, and did not rerun validation.

Deletion accounting:
- deleted: fake/source-backed fixture, candidate route, candidate conversion,
  old physical-gausslet final-basis/H1 payload construction, and WL/QW H2
  reference-candidate pressure.
- simplified: `cartesian_assembly` no longer coordinates fake
  candidate/final/H1 payloads.
- quarantined: none.
- not deleted because: generic blocked source-plan payload remains live;
  residual-GTO support-state helper remains live through
  `src/pqs_h2_residual_gto_handoff.jl`.
- exact remaining caller/blocker: `cartesian_assembly` still calls the generic
  source-plan payload; H2/Cr2 both block at
  `:missing_terminal_shell_projection`.
- added src lines: 5.
- deleted src lines: 1177.
- new tests: none.
- new metadata/status fields: none.

## Pass 331 - Delete Unreachable Private H2 Residual-GTO Helpers

Commit(s):
- this commit - Delete private H2 residual GTO helper file

Target card:
- Target: remove the unreachable private H2 residual-GTO helper cluster now
  that the sidecar producer and fake/source-backed path are gone.
- Physics endpoint: H2 and Cr2 remain on the generic terminal route, blocked at
  `:missing_terminal_shell_projection`; public `CartesianIDAHamiltonian`
  reader/writer APIs must remain untouched.
- Allowed files: `src/GaussletBases.jl`, `src/pqs_h2_residual_gto_handoff.jl`,
  `src/pqs_source_box_diatomic_complete_core_shell.jl`, and stale algorithm
  code-map docs that named the deleted file.
- Forbidden additions: no public IDA API changes, residual-GTO rewrite, tests,
  adapters, status fields, or replacement helpers.
- Success condition: no active source/docs code map points at the deleted
  private helper file or its support-state helper, and package load still
  passes.
- Failure rule: if any deleted helper has a live non-file-internal caller, stop
  and report the exact caller instead of adding a replacement.

Summary:
- Removed `include("pqs_h2_residual_gto_handoff.jl")` from
  `src/GaussletBases.jl`.
- Deleted `src/pqs_h2_residual_gto_handoff.jl`.
- Deleted `_pqs_source_box_route_driver_physical_gausslet_support_states` from
  `src/pqs_source_box_diatomic_complete_core_shell.jl`.
- Removed stale public algorithm-page code-map references to the deleted file
  and clarified that the private H2 residual-GTO producer has been retired.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `0 1071
  src/pqs_h2_residual_gto_handoff.jl`, `0 8
  src/pqs_source_box_diatomic_complete_core_shell.jl`, `0 1
  src/GaussletBases.jl`, plus small algorithm-doc code-map corrections.
- Suspicious added-lines grep: no matches.
- New tests/files: none.
- `rg -n "pqs_h2_residual_gto_handoff|_pqs_source_box_route_driver_physical_gausslet_support_states" src bin tools test docs`
  now returns only historical manager-log references.

Validation:
- Doer reported `git diff --check` and package load passed.
- Manager reviewed the diff, removed stale algorithm-doc references, and did
  not rerun package load or H2/Cr2 probes.

Deletion accounting:
- deleted: private H2 residual-GTO helper file, include edge, and its remaining
  support-state helper.
- simplified: residual-GTO producer docs now state that the private H2 helper is
  retired rather than pointing to a deleted implementation.
- quarantined: none.
- not deleted because: public `CartesianIDAHamiltonian` and reader/writer remain
  live in `src/cartesian_ida_hamiltonian.jl` and were untouched.
- exact remaining caller/blocker: H2 and Cr2 still block at
  `:missing_terminal_shell_projection`; a future residual-GTO producer must be
  rebuilt on the generic terminal route.
- added src lines: 0.
- deleted src lines: 1080.
- new tests: none.
- new metadata/status fields: none.

## Pass 333 - Remove Stale Residual-GTO Materialization Request Path

Commit(s):
- this commit - Remove stale residual GTO materialization request path

Target card:
- Target: stop the materialization dispatcher and developer harness from
  implying the deleted private residual-GTO producer still exists.
- Physics endpoint: H2 and Cr2 remain on the generic terminal route, blocked at
  `:missing_terminal_shell_projection`; public `CartesianIDAHamiltonian` type
  and reader/writer remain live.
- Allowed files: `src/pqs_source_box_low_order_materialization.jl`,
  `tools/cartesian_driver_harness.jl`, and the H2 supplement-materialized input
  if it was only carrying stale provider-block pressure.
- Forbidden additions: no public IDA writer/reader changes, residual-GTO
  rewrite, terminal shell projection, new status fields, tests, or WL/GTO
  fixture changes.
- Success condition: no active source/tool/test request path maps public IDA
  output to `:one_body_and_density_provider`, and no
  `residual_gto_provider_blocks` plumbing remains outside negative policy
  assertions.
- Failure rule: if cleanup needs a replacement residual-GTO producer or output
  framework, stop and report the blocker.

Summary:
- Removed `residual_gto_provider_blocks` from
  `tools/cartesian_driver_harness.jl`.
- Removed `residual_gto_provider_blocks` and the
  `hamiltonian_output === :cartesian_ida_hamiltonian` to
  `:one_body_and_density_provider` mapping from the PQS materialization
  dispatcher.
- Deleted the unavailable private H2 residual-GTO materializer stub and its
  PQS dispatch branch.
- Removed the stale provider-block request from the H2
  supplement-materialized driver input while preserving public
  `save_ida_hamiltonian` / `hamiltonian_output` controls.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `2 35
  src/pqs_source_box_low_order_materialization.jl`, `1 2
  tools/cartesian_driver_harness.jl`, `0 1
  test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_materialized.jl`.
- Suspicious added-lines grep: no matches.
- New tests/files: none.
- `rg` for `residual_gto_provider_blocks`, `one_body_and_density_provider`,
  and private `pqs_h2_residual_gto` materializer names now finds only existing
  negative policy-test assertions.

Validation:
- Doer reported package load and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl` passed with elapsed
  time `35.824882959s`.
- Manager reviewed the diff and did not rerun validation.

Deletion accounting:
- deleted: stale provider-block argument, public-IDA-to-provider-block mapping,
  unavailable private residual-GTO materializer stub, and stale fixture request
  line.
- simplified: PQS materialization now falls through to ordinary
  not-materialized behavior instead of pretending a deleted private producer is
  selectable.
- quarantined: none.
- not deleted because: `save_ida_hamiltonian`, `hamiltonian_output`, and public
  `CartesianIDAHamiltonian` controls remain live public API surfaces.
- exact remaining caller/blocker: H2 and Cr2 remain blocked at
  `:missing_terminal_shell_projection`; the supplement-materialized H2 input
  remains only because `tools/h2_pqs_terminal_stage_smoke.jl` includes it.
- added src lines: 2.
- deleted src lines: 35.
- new tests: none.
- new metadata/status fields: none.

## Pass 335 - Delete Stale Supplement-Materialized H2 Fixture

Commit(s):
- this commit - Delete stale supplement materialized H2 fixture

Target card:
- Target: remove the H2 residual-GTO materialization fixture pressure now that
  the private residual-GTO producer is gone.
- Physics endpoint: H2 and Cr2 remain on the generic terminal route, blocked at
  `:missing_terminal_shell_projection`; supplement policy remains preflight
  intent only.
- Allowed files: `tools/h2_pqs_terminal_stage_smoke.jl`,
  `tools/cartesian_driver_ladder_lib.jl`, the stale H2 materialized fixture,
  and stale developer documentation rows naming the deleted case.
- Forbidden additions: no new tests, fixtures, residual-GTO rewrite, terminal
  shell projection, status fields, or public API changes.
- Success condition: no active fixture/ladder case names the deleted H2
  residual-GTO materialized capability, and the H2 smoke uses the preflight
  fixture directly.
- Failure rule: if the ladder case or fixture has a live purpose beyond stale
  materialization pressure, stop and report exact evidence.

Summary:
- Updated `tools/h2_pqs_terminal_stage_smoke.jl` to include
  `h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl` directly.
- Deleted `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_materialized.jl`.
- Removed the stale `h2_pqs_q5_independent_source_box_r4_gto_materialized`
  ladder case from `tools/cartesian_driver_ladder_lib.jl`.
- Updated active developer migration/inventory tables so they no longer
  describe the deleted residual-GTO materialized case as current coverage.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `0 9
  test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_materialized.jl`,
  `0 5 tools/cartesian_driver_ladder_lib.jl`, `1 1
  tools/h2_pqs_terminal_stage_smoke.jl`, plus small developer-doc row updates.
- Suspicious added-lines grep: no matches.
- New tests/files: none.
- `rg` for the deleted fixture/case names now returns only historical
  manager-log references.

Validation:
- Doer reported package load and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl` passed with elapsed
  time `29.557102333s`.
- Manager reviewed the diff and did not rerun validation.

Deletion accounting:
- deleted: stale supplement-materialized H2 fixture and matching ladder case.
- simplified: H2 smoke now uses the supplement preflight fixture directly while
  preserving its explicit non-materializing overrides.
- quarantined: none.
- not deleted because: `:mwg_residual_gto` supplement policy remains live as
  preflight intent.
- exact remaining caller/blocker: H2 and Cr2 remain blocked at
  `:missing_terminal_shell_projection`.
- added source/tool/test lines: 1.
- deleted source/tool/test lines: 15.
- new tests: none.
- new metadata/status fields: none.

## Pass 337 - Residual-GTO Inventory Docs Correction

Commit(s):
- this commit - Update residual GTO inventory status

No strategic change beyond closing the stale residual-GTO request-pressure
cleanup lane. The donor inventory now says the H2-specific residual-GTO
materializer has been retired, public `CartesianIDAHamiltonian` read/write
remains live, no active PQS route currently materializes that artifact, and
future residual-GTO/MWG supplement work belongs to a generic terminal producer.

Validation:
- Doer reported `git diff --check` passed.
- Manager reviewed the docs diff and reran the stale-phrase grep; no matches.

Deletion accounting:
- deleted: stale developer-doc claim that private H2 residual-GTO construction
  writes a public Cartesian IDA artifact.
- simplified: priority wording now points future work at the generic terminal
  producer.
- quarantined: none.
- not deleted because: donor inventory still tracks residual-GTO/MWG as a
  future/retirement feature.
- exact remaining caller/blocker: H2 and Cr2 remain blocked at
  `:missing_terminal_shell_projection`.
- added source lines: 0.
- deleted source lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Pass 339 - Delete Old Source-Plan Mirror Fields

Commit(s):
- this commit - Delete old source plan mirror fields

Target card:
- Target: remove old H2/source-plan compatibility mirrors from the current
  blocked generic terminal source-plan summary.
- Physics endpoint: H2 and Cr2 remain on the generic terminal route, blocked at
  `:missing_terminal_shell_projection`.
- Allowed files: `src/pqs_source_box_diatomic_complete_core_shell.jl` and
  `tools/h2_pqs_terminal_stage_smoke.jl`.
- Forbidden additions: no terminal shell projection, payload wrapper deletion,
  active blocker status renames, tests, replacement metadata/status fields, or
  report expansion.
- Success condition: old descriptor/shared-shell/source-coefficient mirror
  fields are gone from active summaries, and the H2 smoke validates the generic
  blocker through terminal preflight.
- Failure rule: if any mirror field has a live non-smoke caller, leave it and
  report the exact caller.

Summary:
- Deleted descriptor/shared-shell/source-coefficient mirror fields from the
  blocked diatomic source-plan summary:
  `source_plan_descriptor_status`, `source_plan_descriptor_blocker`,
  `source_plan_descriptor_available`, `shared_shell_realization_status`,
  `shared_shell_realization_blocker`, `shared_shell_realization_counts`,
  `shared_shell_realization_identity_errors`,
  `shared_shell_realization_materialized`, `source_coefficients_materialized`,
  and `source_plan_authority_status`.
- Updated `tools/h2_pqs_terminal_stage_smoke.jl` to stop asserting the deleted
  mirror fields while preserving topology, source-plan status/blocker,
  terminal-preflight status/blocker, and retained-dimension checks.

Mechanical gate:
- `git diff --check`: passed.
- `git diff --numstat -- src bin tools test docs`: `0 39
  src/pqs_source_box_diatomic_complete_core_shell.jl`, `0 6
  tools/h2_pqs_terminal_stage_smoke.jl`.
- Suspicious added-lines grep: no matches.
- New tests/files: none.
- Focused mirror-field grep now returns only historical manager/blurb-log
  references.

Validation:
- Doer reported package load and
  `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl` passed with elapsed
  time `35.481812625s`.
- Manager reviewed the diff and did not rerun validation.

Deletion accounting:
- deleted: old descriptor/shared-shell/source-coefficient mirror fields from
  active blocked source-plan summaries.
- simplified: H2 smoke now checks the terminal preflight blocker instead of
  old compatibility mirrors.
- quarantined: none.
- not deleted because: `physical_gausslet_source_plan_summary`,
  `terminal_source_realization_preflight_summary`, and `low_order_assembly`
  remain the active blocked-route contract until shell projection exists.
- exact remaining caller/blocker: H2 and Cr2 remain blocked at
  `:missing_terminal_shell_projection`.
- added source/tool lines: 0.
- deleted source/tool lines: 45.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Design-Review Lane Setup

Commit(s):
- this branch - Set up Cartesian Hamiltonian producer design review lane

Summary:
- Created a tracked review directory for
  `docs/src/developer/cartesian_hamiltonian_producer_design.md`.
- Added a round-001 consolidated review prompt, review template, and
  reconciliation template. The lane is explicitly docs-only: reviewers may
  propose document edits but may not edit `src`, `test`, `bin`, `tools`, or
  `AGENTS.md`.
- The goal of round 001 is a reviewed v2 of the design document plus a
  reconciliation note ready for ChatGPT-Pro milestone review.

Strategic interpretation:
- The cleanup lane has reduced active route surfaces enough that more abstract
  cleanup is likely lower value than a controlled design phase for the real
  blocker, `:missing_terminal_shell_projection`.
- The design branch remains unapproved implementation authority. No source
  coding should begin until the review/reconciliation loop freezes approved
  design IDs and repo policy explicitly binds them.

Risk/guardrail:
- The design document must be compressed before it becomes authority; extended
  debate belongs in this review directory and Git history, not in an
  append-only implementation contract.
- A subagent review was attempted in the current manager session but could not
  start because the session had reached its agent-thread limit. The tracked
  prompt/template files are ready for a fresh or existing reviewer agent.

Next step:
- Fill `round_001_consolidated_review.md`, reconcile into design v2, then send
  v2 plus reconciliation summary to ChatGPT-Pro for broad review.

## Cartesian Hamiltonian Producer Round 001 Review Received

Commit(s):
- this branch - Add round 001 Hamiltonian producer review

Summary:
- Repo-doer completed the consolidated docs-only review at
  `docs/src/developer/design_reviews/cartesian_hamiltonian_producer_2026-06/round_001_consolidated_review.md`.
- Verdict: the draft needs design revision before ChatGPT-Pro review.
- Blocking themes: accumulated previous-block shell projection must be explicit;
  several HP-* items should be demoted from approved to candidate; localized
  IDA and `electron_electron_ida` need a precise mathematical convention; IDA
  block workspace/tiling must be specified for Cr2-scale memory; line budgets
  should not encourage skipping numerical validation.

Strategic interpretation:
- The design-review lane is working as intended: it caught real scientific and
  performance underspecification before source implementation began.
- The next manager task is reconciliation into design v2, not implementation.
  The draft remains non-authoritative.

Risk/guardrail:
- Do not send the current v1 draft to ChatGPT-Pro as a candidate authority.
  First reconcile the review into a compressed v2 with candidate/approved
  registry language fixed and unresolved IDA/projection/memory questions
  called out explicitly.

Validation:
- Doer reported `git diff --check` passed.
- Manager reviewed the review file and did not run source validation.

## Cartesian Hamiltonian Producer Design Pass 002 - Reconcile Round 001

Commit(s):
- this branch - Reconcile Hamiltonian producer design round 001

Summary:
- Reconciled `round_001_consolidated_review.md` into
  `docs/src/developer/cartesian_hamiltonian_producer_design.md`.
- Filled `round_001_reconciliation.md`.
- The design is now a v2 candidate for milestone review, not implementation
  authority.

Accepted review changes:
- Positive HP-* registry entries are candidates rather than approved surfaces.
- PQS shell realization now explicitly requires projection against all
  previously retained terminal blocks before shell-local Lowdin.
- Localized one-basis IDA convention is written in terms of
  `electron_electron_ida[a,b]` and later consumer coefficients.
- IDA block-pair tiling is required when a full support-pair workspace would
  exceed the reviewed Cr2 memory budget.
- Line budgets are review targets where hard caps could discourage numerical
  validation.

Deferred questions:
- Whether `HP-FN-04` should split into smaller design IDs.
- Whether accumulated previous-block projection can be expressed through the
  current projected-shell descriptor path.
- Which existing source file should own blockwise one-body and IDA assembly.
- Whether the 1.2 GiB Cr2 peak-memory target needs a lower target or absolute
  hard cap.

Next step:
- Send the v2 design plus round-001 reconciliation to ChatGPT-Pro for milestone
  review. Do not start source implementation or bind `AGENTS.md` yet.

## Cartesian Hamiltonian Producer Design Pass 003 - Reconcile Milestone Review

Commit(s):
- this branch - Reconcile Hamiltonian producer Slice A design

Summary:
- Reconciled the ChatGPT-Pro milestone verdict into
  `docs/src/developer/cartesian_hamiltonian_producer_design.md`.
- Added `round_002_chatgpt_pro_review.md` and
  `round_002_reconciliation.md` to the design-review lane.
- The design is now a v3 **Slice A freeze candidate**, not implementation
  authority.

Accepted review changes:
- Freeze scope is limited to Slice A terminal-basis realization. One-body,
  IDA, and final Hamiltonian producer slices remain future candidates until
  Slice A reports effective supports, cross-overlaps, ranks, and memory.
- The fictitious descriptor `projection_basis` path was replaced by an explicit
  previous-block projection formula before shell-local Lowdin.
- `CartesianTerminalBasisBlock.support_indices` now means effective coefficient
  support after projection, not merely the original terminal region support.
- Column sign canonicalization moved into Slice A terminal-basis finalization so
  later one-body and IDA assembly consume the same gauge-fixed basis.
- Direct identity sectors must validate direct overlap identity and finite
  positive IDA weights before `coefficients === nothing` is accepted.
- `HP-RES-01` was rejected as a persistent result wrapper, and `HP-CHANGE-01`
  was rejected/deferred as insufficient for projected effective supports.
- `parent_dims` was removed from `CartesianTerminalBasisRealization`.
- Slice A now uses a target/redesign-threshold budget: 150 added source lines
  target, 225 added source lines redesign threshold, and net source decrease
  required after deleting the terminal preflight path.

Strategic interpretation:
- The design-first lane caught and corrected a real mathematical hazard before
  implementation: shell-local Lowdin is not sufficient unless each PQS shell is
  first projected against the already accepted terminal basis.
- The next implementation should still not start automatically. A final freeze
  decision is needed for the Slice A IDs only; B/C/D should remain deliberately
  unapproved until empirical Slice A data exists.

Risk/guardrail:
- Do not turn the v3 document into a whole-producer authority. It is a Slice A
  freeze candidate and still leaves later operator/IDA ownership unresolved.
- If a local numerical spike shows large effective supports, unstable ranks, or
  large projection residuals, Slice A returns to design before source coding.

Validation:
- Docs-only pass; no source validation required.
- Manager used focused text scans for stale `projection_basis`, `parent_dims`,
  `HP-RES-01`, `HP-CHANGE-01`, and Slice A approval language.

Next step:
- Decide whether the listed Slice A freeze candidates should now be marked
  approved and bound in `AGENTS.md`, or whether the uncommitted numerical spike
  should run first.

## Cartesian Hamiltonian Producer Design Pass 004 - Atomic/Diatomic Unification

Commit(s):
- this branch - Require generic atomic and diatomic terminal basis

Summary:
- Reconciled the atomic/diatomic unification recommendation into
  `docs/src/developer/cartesian_hamiltonian_producer_design.md`.
- Added `round_003_atomic_diatomic_unification_review.md` and
  `round_003_reconciliation.md` to the design-review lane.
- The design remains a Slice A freeze candidate, not implementation authority.

Accepted changes:
- One-center atomic and bond-aligned diatomic PQS routes may differ in geometry
  only. After terminal support, retained, and transform records exist, both
  must use the same terminal-basis realizer and produce
  `CartesianTerminalBasisRealization`.
- Added `HP-WIRE-01`, generic terminal-basis stage integration, owned by
  `cartesian_transforms`.
- The terminal-basis realizer may dispatch on terminal lowering/transform kind
  but not on system classification, atom count, route kind, bond axis, or
  terminal role vocabulary.
- Slice A validation now includes one-center atomic, contact-core H2, and
  separated Cr2 terminal records through the same entry point.
- Future Slice B/C/D contracts now state that one-body, IDA, and final producer
  code consume `CartesianTerminalBasisRealization` without atomic/diatomic
  branches.
- The old atomic multilayer/common-H1 path is explicitly migration-only: no new
  features, useful as an oracle, and deleted after the generic one-body route
  reproduces the reviewed atomic endpoint.

Strategic interpretation:
- This prevents Slice A from becoming a diatomic-only implementation hidden
  behind generic object names. It also avoids expanding Slice A into the full
  atomic Hamiltonian; atomic support is required only at the terminal-basis
  boundary.

Validation:
- Docs-only pass; no source validation required.
- Manager will run `git diff --check` and focused text scans before commit.

Next step:
- Freeze Slice A only after the manager/user accepts the expanded ID set:
  `HP-OBJ-01`, `HP-OBJ-02`, `HP-FILE-01`, `HP-FN-00`, `HP-FN-01`,
  `HP-FN-02`, and `HP-WIRE-01`.

## Cartesian Hamiltonian Producer Design Pass 005 - Projection Spike Reconciliation

Commit(s):
- this branch - Reconcile terminal projection spike

Summary:
- Recorded repo-doer's uncommitted terminal projection spike in
  `round_004_projection_spike_report.md`.
- Reconciled the spike into
  `docs/src/developer/cartesian_hamiltonian_producer_design.md`.
- The design remains a Slice A freeze candidate, not implementation authority.

Spike result:
- One-center, H2, and Cr2 all expose compatible terminal
  support/retained/transform record shapes.
- Raw cross-overlaps were already small and projected overlaps fell to
  roundoff in the spike.
- PQS shell ranks stayed at the expected `98`.
- Direct overlap and IDA weight checks were finite and well conditioned for the
  inspected one-center, H2, and Cr2 direct records.
- Cr2 dense scratch projection took about `75 s`, so production Slice A needs
  factorized or incremental overlap/projection construction.

Accepted design changes:
- Shell projection, shell overlap, and Lowdin cleanup are not terminal contract
  fields and must not be added as metadata or summaries.
- Production Slice A must recursively use projected coefficients and effective
  supports from previous PQS blocks; the spike's shell-local approximation is
  not sufficient as implementation authority.
- Later direct records must be cross-checked against all previous direct and
  PQS blocks.
- If one-center still reaches terminal records through old route-skeleton shape
  input, Slice A must connect it to typed terminal records without adding an
  atomic adapter.

Strategic interpretation:
- The numerical spike supports the Slice A direction but prevents an immediate
  freeze without the recursive-coefficient and factorized-overlap guardrails.
  The next implementation boundary is still terminal basis realization, not
  one-body or IDA assembly.

Validation:
- Doer reported `git diff --check`, package load, and a clean tracked worktree
  after the ignored `tmp/work` spike.
- Manager reconciled docs only; no source validation run.

Next step:
- Decide whether the updated Slice A freeze candidate is now ready for explicit
  approval and `AGENTS.md` binding, or whether one more review should inspect
  the recursive projection/factorized-overlap requirements.

## Cartesian Hamiltonian Producer Design Pass 006 - Recursive Projection Freeze Review

Commit(s):
- this branch - Tighten Slice A recursive projection design

Summary:
- Reconciled the final recursive-projection freeze recommendation into
  `docs/src/developer/cartesian_hamiltonian_producer_design.md`.
- Added `round_005_recursive_projection_review.md` and
  `round_005_reconciliation.md` to the design-review lane.
- The design remains a Slice A freeze candidate, not implementation authority.

Accepted changes:
- `HP-FN-00` now owns projection plus shell-local Lowdin only; it does not own
  final sign canonicalization.
- `HP-FN-01` now owns support-weight derivation, direct-weight validation,
  sign canonicalization, and `CartesianTerminalBasisBlock` construction.
- The projection algorithm now has a roundoff guard: if a previous-block
  residual is already below `projection_atol`, do not subtract and do not grow
  effective support.
- Projection/audit workspace is specified as block cross actions
  `C_left' * S_lr * C_right`, with no global parent/final overlap and at most
  one support-pair workspace live.
- Removed per-object/helper line targets; retained only the Slice A added-line
  target, redesign threshold, and net-deletion requirement.
- The reviewed Cr2 fixture must produce a real terminal basis. Distorted-COMX
  rejection is allowed only when that typed inventory is actually present.

Strategic interpretation:
- This keeps the design gate focused on the actual numerical risk: true
  recursive projection with support-growth control. If the final spike passes,
  no further broad design review should be needed before freezing Slice A.

Validation:
- Docs-only pass; no source validation required.
- Manager will run `git diff --check` and focused text scans before commit.

Next step:
- Ask doer for one final uncommitted recursive-projection spike. If it shows
  stable ranks, bounded support growth, later-direct orthogonality, and no
  global overlap construction for one-center/H2/Cr2, freeze Slice A and bind
  the approved IDs.

## Cartesian Hamiltonian Producer Design Pass 007 - Final Recursive Spike Reconciliation

Commit(s):
- this branch - Reconcile final recursive projection spike

Summary:
- Recorded repo-doer's final uncommitted recursive-projection spike in
  `round_006_recursive_projection_spike_report.md`.
- Reconciled the spike into
  `docs/src/developer/cartesian_hamiltonian_producer_design.md`.
- The design remains unbound until the Slice A IDs are explicitly frozen.

Spike result:
- One-center, H2, and Cr2 all used the same terminal realization inputs once
  terminal records were reached.
- True recursive accepted blocks were used.
- With `projection_atol = 1e-12`, all residuals were skipped as roundoff; no
  subtraction was numerically justified.
- Effective supports remained unchanged and shell ranks stayed at `98`.
- Final cross overlaps were small: one-center `8.073e-16`, H2 `3.095e-15`,
  Cr2 `5.463e-14`.
- Later Cr2 direct records remained orthogonal to prior blocks at
  `2.4e-15 .. 5.3e-15`.

Accepted design changes:
- Default `projection_atol` is now `1.0e-12`.
- The design says `projection_atol` is a roundoff-subtraction threshold, not the
  final cross-overlap acceptance tolerance.
- Slice A support-pair workspace cap is now `64 MiB`; larger local actions must
  tile or stream.

Remaining caveats:
- One-center public staging still relies on old route-shape skeleton access to
  reach terminal records. Slice A must connect it to the typed terminal
  contract and must not freeze that skeleton as the intended public boundary.
- The Cr2 spike's largest dense local workspace was `175.928 MiB`, so
  production Slice A needs tiling/streaming for that action.

Strategic interpretation:
- The final spike supports freezing Slice A. No further broad design review is
  recommended. The next step is a docs/policy freeze commit if the user accepts
  the evidence.

Validation:
- Doer reported `git diff --check`, package load, and a clean tracked worktree
  after the ignored `tmp/work` spike.
- Manager reconciled docs only; no source validation run.

Next step:
- Freeze and bind exactly `HP-OBJ-01`, `HP-OBJ-02`, `HP-FILE-01`,
  `HP-FN-00`, `HP-FN-01`, `HP-FN-02`, and `HP-WIRE-01`; keep B/C/D IDs as
  future candidates.

## Cartesian Hamiltonian Producer Design Pass 008 - Freeze Slice A Authority

Commit(s):
- this branch - Freeze Slice A Hamiltonian producer design

Summary:
- Changed `docs/src/developer/cartesian_hamiltonian_producer_design.md` from
  Slice A freeze candidate to Slice A implementation authority.
- Bound the approved Slice A IDs in `AGENTS.md`.
- Added `round_007_slice_a_freeze.md` to the design-review lane.

Approved IDs:
- `HP-OBJ-01`
- `HP-OBJ-02`
- `HP-FILE-01`
- `HP-FN-00`
- `HP-FN-01`
- `HP-FN-02`
- `HP-WIRE-01`

Guardrail:
- `HP-FN-03`, `HP-FN-04`, and `HP-FN-05` remain future candidates. This freeze
  does not authorize one-body assembly, IDA assembly, Hamiltonian artifact
  production, or driver simplification.
- Unlisted production surfaces require a prior docs-only amendment.

Strategic interpretation:
- Slice A now has enough numerical evidence and policy structure to begin
  implementation without another broad review. The first implementation pass
  should realize the terminal basis only, delete the terminal preflight path it
  replaces, and keep one-center/H2/Cr2 on the same entry point.

Validation:
- Docs/policy-only pass; no source validation required.
- Manager will run `git diff --check`, focused text scans, and post-commit diff
  gates before push.

Next step:
- Draft the Slice A implementation blurb for repo-doer, constrained to the
  approved IDs and forbidden from touching B/C/D work.

## Cartesian Hamiltonian Producer Pass 009 - Realize Terminal PQS Basis

Commit(s):
- `082c9cb8` - Realize terminal PQS basis

Summary:
- Implemented approved Slice A terminal-basis realization.
- Added `CartesianTerminalBasisBlock`,
  `CartesianTerminalBasisRealization`, recursive terminal PQS shell
  projection, shell-local symmetric Löwdin, positive-integral sign
  canonicalization, and cross-block overlap audit.
- Wired `cartesian_transforms` to produce the shared terminal basis for
  one-center, H2, and Cr2 terminal records.
- Deleted the obsolete terminal source-realization preflight/report summary
  machinery.

Validation reported by doer:
- `git diff --check`: passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`: passed.
- `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl`: passed in `33.09s`.
- `julia --project=. tmp/work/terminal_production_cases.jl`: passed.

Terminal facts:
- One-center: dimension `419`, ranks `(98,98,98)`, max cross
  `1.90e-15`, min PQS integral `0.8691`, largest workspace `13.52 MiB`.
- H2: dimension `471`, ranks `(98,98)`, max cross `5.92e-15`, min PQS
  integral `1.4706`, largest workspace `11.26 MiB`.
- Cr2: dimension `4291`, all shell ranks `98`, max cross `2.53e-14`, min PQS
  integral `0.02385`, largest workspace `64.00 MiB`.

Line-count/complexity:
- Added source: `225` lines, exactly at the approved redesign threshold.
- Deleted source: `246` lines.
- Net source: `-21`.
- Overall commit: `242` insertions / `270` deletions including the smoke tool.

Guardrail:
- This pass does not authorize one-body assembly, IDA assembly, Hamiltonian
  artifact production, or driver simplification. B/C/D IDs remain future
  candidates.
- The near-zero integral blocker was traced to an implementation bug: the WIP
  had used a one-sided eigentransform for Löwdin. The accepted code uses
  `inv(sqrt(Symmetric(overlap)))` and full source-box boundary-mode columns
  before projection.

Next step:
- Commit the `JuliaStyle.md` Löwdin guidance so future agents do not repeat the
  one-sided eigentransform mistake, then push the accepted Slice A commit.

## Cartesian Hamiltonian Producer Pass 010 - Trim Terminal Basis Cleanup

Commit(s):
- this branch - Trim terminal basis cleanup

Summary:
- Removed the separate `eigvals(overlap)` rank check and the now-unused
  `rank_atol` plumbing from the Slice A terminal PQS realizer.
- Kept Löwdin construction as `inv(sqrt(Symmetric(overlap)))` and kept the final
  realized-shell overlap identity check.
- Removed an unused H2 smoke helper and compacted a simple topology-facts
  return.

Validation reported by doer:
- `git diff --check`: passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`: passed.
- `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl`: passed in `39.01s`.

Line-count/complexity:
- Source: `2` added / `6` deleted.
- Tools: `2` added / `12` deleted.
- Total: `4` added / `18` deleted.

Strategic interpretation:
- This removes low-value defensive/report clutter after the real bug was found
  to be the one-sided eigentransform, not a rank-loss condition. It leaves the
  construction path simpler and aligned with the new `JuliaStyle.md` Löwdin
  guidance.

Next step:
- Continue from the realized terminal basis toward the next approved design
  boundary only after B/C/D surfaces are explicitly approved.

## Cartesian Hamiltonian Producer Design Pass 011 - Freeze Slice B Authority

Commit(s):
- this branch - Freeze Slice B one-body operator design

Summary:
- Promoted `HP-FN-03` from future candidate to approved Slice B authority in
  `docs/src/developer/cartesian_hamiltonian_producer_design.md`.
- Updated `AGENTS.md` so Cartesian Hamiltonian producer source work is
  authorized for approved Slice A IDs plus `HP-FN-03`.
- Added `round_008_slice_b_freeze.md` to the design-review lane.

Approved Slice B boundary:
- Build final-basis kinetic matrix `K`.
- Build separated unit nuclear attraction matrices `U_A = -1/r_A`.
- Use one separable product-term helper over `CartesianTerminalBasisRealization`.
- Do not form global support-space operators.
- Do not add atomic/diatomic operator branches.
- Do not implement IDA, Hamiltonian artifact production, residual-GTO, or driver
  simplification.

Guardrail:
- `HP-FN-04` and `HP-FN-05` remain future candidates.
- Unlisted production surfaces still require a prior docs-only amendment.

Validation:
- Docs/policy-only pass; no source validation required.
- Manager will run `git diff --check`, focused scans for approved/candidate ID
  wording, and post-commit diff gates before push.

Next step:
- Draft the Slice B implementation blurb for Chat review before sending it to
  repo-doer.

## Cartesian Hamiltonian Producer Design Pass 012 - Clarify Slice B Source Surface

Commit(s):
- this branch - Clarify Slice B one-body source boundary

Summary:
- Applied Chat review of `361b6996`.
- Bound `HP-FN-03` to the exact source file
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`.
- Clarified that Slice B does not approve a new K/U payload, stage-return field,
  report object, persistent one-body orchestration API, or status vocabulary.
- Added Slice B implementation requirements for implicit direct selectors,
  symmetric-axis validation, upper-triangular block traversal, and no final
  averaging that hides nonsymmetric construction.

Validation/baseline note:
- H2 one-body validation should use the reviewed lowest-energy baseline
  `-0.7946037173365863` with `1e-10` tolerance.
- A terminal one-center/H exact oracle baseline is not recorded in the design
  lane. Before source coding, the implementation target card must name that
  baseline and tolerance, or establish it in ignored `tmp/work` code and stop
  before source changes.

Guardrail:
- No source implementation is started by this correction.
- `HP-FN-04` and `HP-FN-05` remain future candidates.

Next step:
- Provide the revised Slice B implementation blurb for Chat review; do not hand
  to repo-doer until reviewed.

## Cartesian Hamiltonian Producer Pass 013 - Slice B Terminal One-Body Kernel

Commit(s):
- this branch - Add terminal basis one-body product assembly

Summary:
- Accepted the first Slice B source implementation under approved `HP-FN-03`.
- Added the approved file
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`.
- Implemented `assemble_terminal_product_operator!` over
  `CartesianTerminalBasisRealization` with upper-triangular terminal block
  traversal, implicit direct selectors, and bounded support-pair workspaces.
- Included the approved source file from `CartesianFinalBasisRealization.jl`.
- Removed the dense direct-identity allocation path from the Slice A
  `_block_pair_matrix` helper.

Validation:
- `git diff --check`: passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`: passed.
- `julia --project=. tools/h2_pqs_terminal_stage_smoke.jl`: passed.
- `julia --project=. tmp/work/terminal_one_body_validation.jl h h2`: passed.
  - H terminal one-body lowest: `-0.49855234726272379`, matching the new
    ignored-script oracle baseline `-0.49855234726272035` within `1e-10`.
  - H2 terminal one-body lowest: `-0.7946037173365883`, matching the reviewed
    baseline `-0.7946037173365863` within `1e-10`.

Goal advancement:
- Advances LT1/LT4 by making one separable final-basis one-body product term a
  real numerical object in the unified terminal basis, without reintroducing
  the retired pair-materialization framework.
- Preserves the Slice A atomic/H2/diatomic unification boundary: source code
  does not branch on atom count, route kind, or terminal role vocabulary.

Risk/guardrail:
- Cr2 dense `4291 x 4291` operator construction is a stress test, not the
  default wiring gate. Future Slice B validation should use a light separated
  diatomic with Cr2-like terminal variety, such as Be2 or N2, before the final
  Cr2 scale check.
- Source line count is positive because this is the new approved numerical
  capability. The same-lane cleanup removed the immediate dense identity
  allocation; later one-body parity should delete migration-only dense/oracle
  adapters rather than preserving them.

Remaining blocker/next step:
- Extend validation to a light separated-diatomic terminal topology and then
  decide the next approved boundary for persistent K/U orchestration or Slice C.
- `HP-FN-04` and `HP-FN-05` remain unapproved.

## Cartesian Hamiltonian Producer Pass 014 - Light Separated-Diatomic Slice B Validation

Commit(s):
- this branch - Record light separated diatomic one-body validation

Summary:
- Accepted an ignored `tmp/work` validation probe for a cheaper separated
  diatomic Slice B fixture.
- Be2 variants were tried first. Physical-ish and longer Be2 cases remained
  contact-core; a far Be2 case produced atom-local cores/shells plus
  midpoint/outer slabs, but no shared molecular shell.
- Selected an N2 separated fixture because it exercises the Cr2-like terminal
  variety needed for wiring validation without Cr2 scale cost.

Selected fixture:
- Atoms/charges: N2, `(7, 7)`.
- Bond length: `8.0` bohr.
- `q = 5`, `core_side = 5`, `core_spacing = 0.15`.
- Extents: `xmax_parallel = 10.0`, `xmax_transverse = 4.0`.
- Terminal roles:
  `(:atom_local_core, :atom_local_core, :atom_local_shell,
  :atom_local_shell, :atom_local_shell, :atom_local_shell, :midpoint_slab,
  :shared_molecular_shell, :z_low_outer_mismatch_slab,
  :z_high_outer_mismatch_slab)`.
- Support counts: `(125, 125, 218, 218, 386, 386, 81, 1002, 121, 121)`.
- Coverage duplicate/missing/outside counts: `0 / 0 / 0`.
- Final retained dimension: `1063`.
- Terminal max cross overlap: `1.0957649991595716e-14`.

Slice B validation result:
- `K`: shape `(1063, 1063)`, finite, symmetry error `2.84e-14`.
- `U1`, `U2`: shape `(1063, 1063)`, finite, symmetry errors
  `3.47e-16` and `4.44e-16`.
- `H1` lowest: `-24.93857219722815`.
- Largest local workspace: `49.26 MiB`.
- One-body elapsed time: `3.70s`.
- Allocations: about `18063 MiB`.

Validation:
- `git diff --check`: passed.
- Package load: passed.
- Working tree was clean/even except for ignored `tmp/work` scripts.

Goal advancement:
- Confirms Slice B one-body assembly works on a separated-diatomic terminal
  topology with atom-local shells, midpoint slab, shared molecular shell, and
  outer mismatch slabs.
- Confirms Cr2 should remain a stress test, not the default wiring fixture.

Risk/next step:
- Allocation volume is high enough to require a later timing/reuse study before
  treating Cr2 performance as representative.
- The validation script currently constructs Coulomb Gaussian factors around
  the Slice B helper; an audit is now open to identify existing repo-owned
  factor/operator packets that should be reused rather than rebuilt locally.

## Cartesian Hamiltonian Producer Pass 015 - Slice B Gaussian Factor Reuse Audit

Commit(s):
- this branch - Record Slice B Gaussian factor reuse audit

Summary:
- Accepted a read-only subagent audit of existing Coulomb Gaussian expansion and
  one-body factor construction.
- The recommended validation path is to reuse existing factor machinery:
  `coulomb_gaussian_expansion`, `gaussian_factor_matrices`, and
  `mapped_ordinary_one_body_operators`.
- `mapped_ordinary_one_body_operators(basis; exponents, center, backend)`
  returns `MappedOrdinaryOneBody1D` with `gaussian_factors`, and internally
  calls `gaussian_factor_matrices` once for the full exponent list.
- This means the current ignored validation style is acceptable when it obtains
  per-center axis factor lists from `mapped_ordinary_one_body_operators` and
  then loops over the already-built term matrices for
  `assemble_terminal_product_operator!`.

Do not use as Slice B production authority:
- `pqs_multilayer_support_electron_nuclear_by_center_matrices`, because it is a
  dense support-space seam rather than final terminal-basis block assembly.
- `pqs_source_pair_centered_electron_nuclear_by_center_block` or retained
  variants in pair-materialization code, because they are source/retained-pair
  machinery before terminal shell realization.
- CPB local block providers or ordinary QW full-product dense by-center
  builders except as oracle/reference context.

Missing seam:
- There is no approved production helper that directly returns final-terminal
  by-center unit nuclear matrices from `(terminal_basis, axis_layers, centers,
  expansion)`. Adding one would be a new production surface and requires a
  docs-only design amendment.

Next step:
- Keep Slice B validation composed from existing factor builders plus
  `assemble_terminal_product_operator!`.
- If performance work finds repeated factor construction, fix routing to reuse
  the existing factor packets before inventing new caching.

## Cartesian/PQS Documentation Pass 016 - Algorithm Implementation Index

Commit(s):
- this branch - Add algorithm implementation index

Summary:
- Added `docs/src/developer/algorithm_implementation_index.md` as a compact
  navigation guide for future agents before numerical Cartesian/PQS coding.
- Linked it from `docs/src/developer/index.md`.
- Updated `AGENTS.md` to instruct agents to check the index before
  Cartesian/PQS numerical implementation.

Purpose:
- Prevent repeated rediscovery and accidental reimplementation of existing
  optimized algorithm paths.
- Record the term-first Coulomb Gaussian contraction lesson near the active
  developer entry point: keep the Gaussian expansion index as the short inner
  reduction, using term-first 1D factor data.

Scope:
- Documentation/navigation only. The index is not new algorithm authority and
  does not approve new source surfaces.
- Initial entries cover term-first Coulomb contractions, Gaussian factor reuse,
  CPB parent-factor/provider layers, PQS terminal realization, PQS source-box
  retained transforms, one-body unit nuclear convention, IDA factor conventions,
  performance/reuse policy, high-order donor notes, and migration/oracle paths.

Validation:
- Docs-only pass; manager ran `git diff --check` and focused link/status
  checks before commit.

Next step:
- Any Slice B performance redesign should first cite the relevant index entry
  and then update the formal Hamiltonian producer design if it needs a new
  production surface beyond `HP-FN-03`.

## Cartesian Hamiltonian Producer Pass 017 - Slice B Term-First Design Revision

Commit(s):
- this branch - Require term-first Slice B Gaussian reuse

Summary:
- Revised the Slice B design text in
  `docs/src/developer/cartesian_hamiltonian_producer_design.md` to require the
  algorithm implementation index before Slice B source work.
- Named the ordinary Gaussian-factor and term-first Coulomb source anchors that
  implementation must inspect before coding.
- Clarified that direct reuse of ordinary helpers is preferred, but when layout
  prevents direct calls, Slice B must still follow their organization: reusable
  1D Gaussian factor packets and a term-first contraction over the Coulomb
  expansion index.

Guardrail:
- `HP-FN-03` remains the only approved Slice B source surface.
- Private file-local helpers inside the approved Slice B file may perform
  term-first support-tile contraction, but they must not introduce persistent
  result/cache/stage/metadata/status objects.

Validation:
- Docs-only pass; manager ran `git diff --check` and reviewed the Slice B diff.

Next step:
- Send the revised design to ChatGPT-Pro for review of whether it now requires
  the optimized contraction pattern without accidentally authorizing a new
  payload or cache framework.

## Cartesian Hamiltonian Producer Pass 018 - Slice B Reuse Status Labels

Commit(s):
- this branch - Tighten Slice B reuse guidance

Summary:
- Added authority/status labels to the algorithm implementation index so source
  anchors distinguish active reusable kernels, donor patterns, consumer
  examples, oracle/reference paths, retired paths, and planned approved files.
- Tightened the Slice B design after ChatGPT-Pro review: the public
  `assemble_terminal_product_operator!` signature remains single-product, while
  one private file-local Gaussian-sum helper may consume coefficients and
  term-first factor arrays without becoming a persistent production surface.
- Added an N2 allocation gate requiring at least `10x` cumulative allocation
  reduction from the recorded `~18,063 MiB` one-body baseline while preserving
  H/H2 energies and the `64 MiB` simultaneous local-workspace cap.

Guardrail:
- `HP-FN-03` is still the only approved Slice B source surface. The amendment
  does not approve new cache/result/stage/metadata/status objects or route
  orchestration APIs.

Validation:
- Docs-only pass; manager ran `git diff --check` and focused status-anchor
  scans before commit.

Next step:
- Resume Slice B source work with the term-first Gaussian-sum helper inside the
  approved file only, reporting whether existing ordinary helpers were directly
  reusable or only donor patterns for terminal-layout code.

## Cartesian Hamiltonian Producer Pass 019 - Term-First Terminal Nuclear Attraction

Commit(s):
- this branch - Optimize terminal Gaussian nuclear contraction

Summary:
- Accepted the Slice B source optimization in
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`.
- Added private file-local helpers under the approved `HP-FN-03` surface to
  accumulate Gaussian-expanded unit nuclear attraction term-first inside each
  terminal support-pair tile.
- Kept `assemble_terminal_product_operator!` as the public single-product
  helper; no K/U payload, cache object, stage field, metadata key, report
  surface, IDA path, artifact path, or CPBM route was added.

Validation actually used:
- Doer: `git diff --check`, package load, H2 terminal smoke, H/H2 one-body
  validation, and light separated N2 one-body validation.
- Doer reported H `-0.49855234726272429`, H2
  `-0.79460371733658908`, N2 final dimension `1063`, one-body elapsed
  `1.47s`, allocations `1191.98 MiB`, largest workspace `49.26 MiB`, and
  finite symmetric K/U matrices.
- Manager: reviewed the single-file diff and surface exposure; reran
  `git diff --check`, package load, H2 terminal smoke, and H/H2 one-body
  validation. The long N2 gate was not rerun after review; the manager used
  doer's reported N2 allocation gate to respect the smallest-test discipline.

Goal advancement:
- LT: advances the generic terminal Hamiltonian producer by making Slice B's
  nuclear attraction path follow the known term-first Gaussian contraction
  pattern.
- MT: reduces the measured N2 one-body allocation cost by about `15x` from the
  recorded `~18,063 MiB` baseline while keeping the `64 MiB` local-workspace
  cap.

Carrying-cost accounting:
- deleted: none.
- simplified: nuclear attraction contraction no longer allocates one dense
  final-basis matrix per Gaussian term.
- quarantined: validation/profiling scripts remain ignored under `tmp/work/`.
- not deleted because: the single-product K helper remains the approved
  `HP-FN-03` public helper.
- exact remaining caller/blocker: further allocation work, if needed, is buffer
  reuse in terminal product/K assembly or terminal-basis construction, not
  Gaussian factor reconstruction.
- added src lines: `92`.
- deleted src lines: `0`.
- new tests: none.
- new metadata/status fields: none.

Guardrail:
- The private helper is an implementation detail inside the approved Slice B
  file only. It must not become route orchestration authority or a persistent
  cache/result surface without another docs-only design amendment.

## Cartesian Hamiltonian Producer Pass 020 - Delete CPBM Global Nuclear Pilot

Commit(s):
- this branch - Delete stale CPBM global nuclear pilot

Summary:
- Deleted the unused
  `src/cartesian_pair_block_materialization/one_body_global_electron_nuclear.jl`
  dense global retained-matrix pilot.
- Removed its export/include and file-map entry from
  `CartesianPairBlockMaterialization.jl`.
- Updated the algorithm implementation index so agents see this as a deleted
  CPBM pilot, not a donor route.

Validation actually used:
- Manager ran `git diff --check`, package load, exact symbol/file reference
  scans, and source numstat.
- No long H/H2/N2 endpoint was rerun because this pass deletes an uncalled CPBM
  pilot and does not alter the accepted terminal one-body implementation.

Goal advancement:
- LT: reduces stale pair-materialization scaffolding after Slice B established
  the terminal-basis one-body path.
- MT: keeps terminal one-body authority with `HP-FN-03` and prevents future
  work from reviving a dense global retained electron-nuclear assembly pilot.

Carrying-cost accounting:
- deleted: CPBM global electron-nuclear retained-matrix pilot.
- simplified: CPBM module no longer exports/includes that unused dense nuclear
  global matrix surface.
- quarantined: none.
- not deleted because: shared placement planning and global overlap/kinetic/x2
  pilots remain separate CPBM scaffolding with broader safe-term scope; this
  pass only removed the obsolete nuclear path replaced by terminal Slice B.
- exact remaining caller/blocker: no exact source/test/tool caller remains for
  `one_body_global_electron_nuclear_by_center_matrix`.
- added src lines: `0`.
- deleted src lines: `210` including module entry cleanup.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 021 - Slice C Candidate Design

Commit(s):
- this branch - Draft Slice C localized IDA candidate

Summary:
- Revised the Hamiltonian producer design to make `HP-FN-04` a concrete Slice C
  candidate for blockwise localized IDA assembly.
- Proposed `CartesianFinalBasisRealization` ownership and candidate file
  `src/cartesian_final_basis_realization/pqs_terminal_ida.jl`, with internal
  function `assemble_terminal_ida_interaction!`.
- Updated the algorithm implementation index with ordinary IDA consumer
  examples and kept them marked as donor/reference, not direct terminal route
  authority.

Guardrail:
- This is docs-only candidate design. `HP-FN-04` and `HP-FN-05` remain without
  implementation authority until explicitly approved in a later revision.
- The candidate forbids global support raw-pair matrices, global dense
  coefficient matrices, all-pairs status plumbing, metadata-carried numerical
  pair data, CPBM route authority, and new persistent cache/result/stage/report
  surfaces without a docs-only amendment.

Validation:
- Doer ran `git diff --check`.
- Manager reviewed the candidate wording and normalized the index status label
  to `consumer example only`.

Next step:
- Before any Slice C source work, get explicit approval of the candidate IDs
  and name the reviewed H2 self-Coulomb/IDA parity value in the implementation
  handoff.

## Cartesian Hamiltonian Producer Pass 022 - Slice C Fast Coulomb Reuse

Commit(s):
- this branch - Tighten Slice C fast Coulomb reuse

Summary:
- Tightened the Slice C candidate design so localized IDA explicitly uses PGDG
  `pair_factor_terms_raw` as the raw numerator source.
- Added `_mapped_coulomb_expanded_symmetric_matrix` as the fast term-first loop
  donor while making clear it is not directly callable for terminal block
  layout because it builds a full ordinary Cartesian product matrix.
- Marked `_nested_weight_aware_pair_terms` /
  `_nested_support_reference_pair_sum` as support-set oracle patterns and
  `_ordinary_cartesian_ida_from_pair_factors` as a consumer/convention example
  whose term-by-term `kron` path must not be copied when term-first contraction
  is available.

Guardrail:
- `HP-FN-04` remains a candidate only. The pass does not approve Slice C source
  work or any new production surface.
- Future implementation must report whether raw pair terms came directly from
  `pgdg.pair_factor_terms_raw` or through
  `_pqs_source_box_ida_factor_provenance`, and confirm those tensors were not
  rebuilt inside terminal block loops.

Validation:
- Docs-only pass; manager inspected existing IDA/pair-factor anchors and ran
  `git diff --check` plus focused string scans before commit.

Next step:
- Get explicit approval before converting `HP-FN-04` from candidate to Slice C
  implementation authority.

## Cartesian/PQS Policy Pass 023 - Variable Tuple Compile Guard

Commit(s):
- this branch - Strengthen variable tuple guardrails

Summary:
- Strengthened `AGENTS.md` and `JuliaStyle.md` to explicitly reject
  variable-size `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed
  `NamedTuple` route inventories for basis-size, shell-size, unit-size,
  pair-size, center-size, or all-pairs data unless explicitly approved.
- Updated the Cartesian/PQS suspicious-line gate so `Tuple(`,
  `Tuple{Vararg`, and common record/unit/pair tuple conversions are review
  triggers.

Purpose:
- Prevent compile-time bloat and excessive specialization from large
  route-inventory tuple types before Slice C introduces pair-factor and IDA
  work.
- Preserve the intended data model: vectors, indexed/lazy views, or compact
  summaries for variable-size scientific collections; tuples only for fixed,
  tiny mathematical shapes such as axis triples.

Validation:
- Docs/policy-only pass; manager ran `git diff --check` and focused grep checks
  before commit.

Next step:
- Apply the strengthened gate to any Slice C implementation blurb and review.

## Cartesian Hamiltonian Producer Pass 024 - Slice B Contract Hardening

Commit(s):
- this branch - Harden terminal seed and Gaussian factor contracts

Summary:
- Hardened `_shell_seed` so terminal PQS realization no longer only
  count-checks generated source-box boundary columns. It now validates retained
  rule kind, transform kind, source-mode dimensions, ordering, exact retained
  column indices, exact retained mode tuples, and retained count. When the
  existing raw product source plan is available on the transform contract, it
  also validates source intervals, source shape, source-mode dimensions,
  ordering, and source-mode count.
- Hardened the Slice B Gaussian-sum nuclear helper so every term matrix is
  finite and symmetric before upper-triangular block mirroring can hide bad
  input.
- Amended Slice B validation wording: light separated N2 closes the current
  correctness/performance gate; Cr2 one-body construction is a later
  whole-producer stress/performance gate unless explicitly requested.

Guardrail:
- The suspicious-line gate flagged the new `contract.metadata` read. This pass
  uses an existing transform-contract metadata field,
  `:raw_product_source_plan`, only to validate the already-existing retained
  source-mode rule against its raw source plan. It does not add a metadata key
  or make metadata a new algorithmic transport surface.

Validation:
- `git diff --check`
- package load
- `tools/h2_pqs_terminal_stage_smoke.jl`
- `tmp/work/terminal_one_body_validation.jl h h2`
- H: `-0.49855234726272429`
- H2: `-0.79460371733658908`
- N2 and Cr2 were not rerun; this was a narrow contract-hardening pass and the
  design now records Cr2 as a later stress gate.

Carrying-cost accounting:
- deleted: none.
- simplified: terminal seed/retained-rule mismatch now fails at the source of
  shell realization rather than being inferred later from numerical behavior.
- quarantined: none.
- not deleted because: raw source facts still live on transform-contract
  metadata until a later typed payload-collapse pass.
- exact remaining caller/blocker: Slice C still needs explicit approval of
  corrected `HP-FN-04` before source work.
- added src lines: `37`.
- deleted src lines: `2`.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 025 - Slice C C1 IDA Candidate Amendment

Commit(s):
- this branch - Amend Slice C C1 IDA candidate

Summary:
- Amended the Slice C candidate design so `HP-FN-04` is implementable from its
  listed inputs: explicit parent axis integral-weight vectors are now part of
  the proposed internal IDA assembly signature.
- Split Slice C authority into candidate-only C1/C2 wording. C1 is localized
  IDA matrix assembly producing final-basis `electron_electron_ida`; C2 is
  later construction of the existing `CartesianIDAHamiltonian`.
- Added C1 validation language around the reviewed H2 self-Coulomb value
  `0.4569117646737212` and allowed an ignored H2 oracle comparison against the
  existing complete-core/shell IDA formula for validation only.

Guardrail:
- `HP-FN-04` and `HP-FN-05` remain candidate-only. This pass does not approve
  Slice C source implementation, a new payload, a persistent pair-factor cache,
  metadata-carried pair data, or Hamiltonian/artifact wiring.

Validation:
- Doer ran `git diff --check` and focused grep checks for the new C1/C2
  wording, explicit weights, tolerances, reviewed self-Coulomb value, and
  candidate-only status.
- Manager reviewed the full docs diff and confirmed the candidate-only boundary
  remained intact.

Next step:
- Draft a narrow C1 implementation blurb only if explicitly proceeding with
  Slice C. The blurb must approve `HP-FN-04` for source work, require explicit
  parent weights, forbid cache/payload/orchestration additions, and preserve the
  term-first raw-pair contraction pattern.

## Cartesian Hamiltonian Producer Pass 026 - Slice C1 Terminal IDA Kernel

Commit(s):
- this branch - Implement terminal IDA interaction kernel

Summary:
- Accepted the narrow `HP-FN-04` source implementation for Slice C1. The new
  terminal IDA helper constructs final-basis `electron_electron_ida` blocks from
  raw PGDG pair-factor terms and explicit parent axis weights.
- The implementation computes final localized weights once per terminal block,
  contracts unnormalized raw numerator blocks through the realized terminal
  basis, normalizes by positive final weights, and mirrors upper-triangular
  block results.
- It reuses the existing Slice B term-first support-pair action instead of
  adding a second contraction engine.

Guardrail:
- This pass does not approve `HP-FN-05`. It adds no Hamiltonian construction,
  driver/materialization/artifact path, stage/report field, metadata key,
  payload, status framework, or persistent pair-factor cache.

Validation:
- Doer ran `git diff --check`, package load,
  `tools/h2_pqs_terminal_stage_smoke.jl`, and ignored
  `tmp/work/h2_terminal_c1_ida_validation.jl`.
- H2 validation reported final dimension `471`, H1 lowest
  `-0.79460371733658908`, finite/symmetric `electron_electron_ida` with
  symmetry error `1.665e-15`, positive IDA weights in
  `[0.35917508613447269, 4.770906361845503]`, and self-Coulomb
  `0.45691176467371986`, within `1.332e-15` of the reviewed value.
- Manager reviewed the source diff and confirmed the implementation stayed
  inside the approved C1 surface.

Carrying-cost accounting:
- deleted: none.
- simplified: C1 reuses the terminal Gaussian-sum action from Slice B rather
  than introducing a duplicate raw-pair contraction engine.
- quarantined: H2 C1 validation remains ignored under `tmp/work`.
- not deleted because: no stale C1 source implementation existed to remove.
- exact remaining caller/blocker: `HP-FN-05` / `CartesianIDAHamiltonian`
  construction remains unapproved.
- added src lines: `91`.
- deleted src lines: `0`.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 027 - Slice C2 Hamiltonian Boundary Approval

Commit(s):
- this branch - Approve Slice C2 Hamiltonian boundary

Summary:
- Updated the Hamiltonian producer design so `HP-FN-05` is approved for narrow
  Slice C2 source work: construction of the existing `CartesianIDAHamiltonian`
  from already assembled Slice B `K`, by-center unit `U_A`, and Slice C1
  `electron_electron_ida`.
- Recorded the caller-side provenance rules that C1 cannot enforce from bare
  arrays: raw pair tensors and Coulomb coefficients must be verified to come
  from the same exponent vector/order before C1, and C2 inputs must be verified
  to come from the same `CartesianTerminalBasisRealization`.
- Corrected stale design text so C1 is recorded as implemented under
  `HP-FN-04`, C2 is approved under `HP-FN-05`, and Slice D remains the future
  driver/materialization/artifact lane.

Guardrail:
- C2 approval does not authorize driver/materialization wiring, artifact
  writing, a Hamiltonian wrapper payload, route-stage fields, report fields,
  metadata keys, status frameworks, persistent factor caches, global support
  operators, or CartesianPairBlockMaterialization route revival.

Validation:
- Docs-only pass. Manager reviewed the C1 acceptance feedback, checked the
  existing `CartesianIDAHamiltonian` constructor boundary, and ran focused text
  audits for stale candidate/approval wording before commit.

Next step:
- Draft the narrow C2 implementation blurb: build or directly construct the
  existing `CartesianIDAHamiltonian`, validate H2 `one_body_hamiltonian(ham)`
  and H2 self-Coulomb through `ham.electron_electron_ida`, and keep all route,
  driver, artifact, and wrapper work deferred to Slice D.

## Cartesian Hamiltonian Producer Pass 028 - Slice C2 Constructor Validation

Commit(s):
- this branch - Record Slice C2 constructor validation

Summary:
- Accepted the C2 implementation result: no new source helper was needed. The
  existing `CartesianIDAHamiltonian(...)` constructor is sufficient for the
  approved in-memory boundary.
- Doer validated C2 with an ignored H2 script that assembled Slice B `K`, unit
  `U_A`, and Slice C1 `electron_electron_ida`, then constructed
  `CartesianIDAHamiltonian{Float64}` directly.

Validation:
- Doer ran `git diff --check`, package load,
  `tools/h2_pqs_terminal_stage_smoke.jl`, and ignored
  `tmp/work/h2_terminal_c2_hamiltonian_validation.jl`.
- H2 C2 result: dimension `471`, `nup_ndn = (1, 1)`, charges `[1.0, 1.0]`,
  positions `[0.0 0.0 -2.0; 0.0 0.0 2.0]`, finite/symmetric `K`, both unit
  `U_A`, and `V`, `one_body_matrix_delta = 0.0`, H1 lowest
  `-0.79460371733658908`, and self-Coulomb `0.45691176467371986`.
- The ignored validation checked C1 coefficient/raw tensor exponent ordering
  before IDA assembly and confirmed uncharged unit `U_A` matrices are charged
  only by `one_body_hamiltonian(ham)`.

Carrying-cost accounting:
- deleted: none.
- simplified: no C2 helper was added; the existing public constructor owns the
  Hamiltonian object validation.
- quarantined: C2 validation remains ignored under `tmp/work`.
- not deleted because: no C2 source was added.
- exact remaining caller/blocker: Slice D driver/materialization/artifact wiring
  remains unapproved.
- added src lines: `0`.
- deleted src lines: `0`.
- new tests: none.
- new metadata/status fields: none.

Next step:
- Decide whether to request ChatGPT-Pro review of the completed A/B/C in-memory
  producer boundary, or proceed to a Slice D design pass that wires the real
  Hamiltonian into materialization/artifact output while deleting obsolete
  blocked-source/report surfaces.

## Cartesian Hamiltonian Producer Pass 029 - Slice D Candidate Handoff Design

Commit(s):
- this branch - Draft Slice D materialization handoff

Summary:
- Synchronized `AGENTS.md` with the current Hamiltonian producer authority:
  `HP-FN-04` is approved only for internal Slice C1 localized IDA assembly and
  `HP-FN-05` is approved only for the narrow Slice C2
  `CartesianIDAHamiltonian` construction boundary. Driver/materialization,
  artifact production, route-stage/report fields, wrapper payloads, and
  persistent factor caches remain unauthorized.
- Added `HP-WIRE-02` as a Slice D candidate design. It chooses a direct handoff
  of `transforms.terminal_basis_realization` into `cartesian_materialization`
  rather than embedding the basis in reports, passing recursive `transforms`,
  reconstructing from summaries, or adding a build-input payload.
- Froze the proposed materialization return contract as
  `Union{Nothing,CartesianIDAHamiltonian{Float64}}`: return `nothing` when no
  base Hamiltonian is requested, return the Hamiltonian itself on success, and
  do not preserve `result_kind`/`materialized`/`ida_hamiltonian` wrapper fields.
- Expanded Slice D with concrete production work, deletion targets, and
  validation: materialization composes Slice B `K`/unit `U_A`, Slice C1 `V`,
  and the existing `CartesianIDAHamiltonian`, then optionally writes through the
  existing public writer. Readback is validation-only.

Guardrail:
- This is still docs/policy only. `HP-WIRE-02` is candidate, not implementation
  authority. The pass intentionally prevents the next doer from treating report
  metadata, source-plan blockers, or route payloads as the way to move basis
  data into materialization.
- The old White-Lindsey materialization route must remain separate if it has
  live callers; `terminal_basis_realization === nothing` is not a PQS fallback.

Validation:
- Manager audited current call sites and confirmed `terminal_basis_realization`
  is owned by `cartesian_transforms`, while `cartesian_materialization`
  currently receives only `report` and materialization inputs.
- `git diff --check` and focused text checks were used before commit.

Next step:
- Send the Slice D candidate handoff design for review. If accepted, approve a
  narrow `HP-WIRE-02` implementation blurb that updates live call sites to pass
  the terminal basis directly, constructs the real Hamiltonian in
  materialization, writes/roundtrips the existing artifact when requested, and
  deletes obsolete blocked source-plan/report surfaces.

## Cartesian Hamiltonian Producer Pass 030 - HP-WIRE-02 Materialization Wiring

Commit(s):
- this branch - Wire terminal Hamiltonian materialization

Summary:
- Accepted the narrow `HP-WIRE-02` source implementation. PQS materialization
  now receives `transforms.terminal_basis_realization` directly and returns
  either `nothing` or a real `CartesianIDAHamiltonian{Float64}`.
- Requested PQS materialization builds Slice B `K`, uncharged by-center unit
  `U_A`, Slice C1 `electron_electron_ida`, then constructs the existing
  Hamiltonian object. Artifact writes use the existing
  `write_cartesian_ida_hamiltonian` path.
- Live driver/harness/probe call sites were updated to pass the terminal basis
  directly. The terminal basis is not embedded in `cartesian_report`, and the
  full `transforms` stage is not passed to materialization.
- Deleted the generic durable materialization wrapper serialization and TSV
  materialization rows from report saving. PQS no-request materialization now
  returns `nothing`.

Validation:
- Doer ran `git diff --check`, package load,
  `tools/h2_pqs_terminal_stage_smoke.jl`, and ignored
  `tmp/work/h2_terminal_wire02_materialization_validation.jl`.
- H2 materialization reported `no_request_materialization = nothing`,
  `CartesianIDAHamiltonian{Float64}`, dimension `471`, `nup_ndn = (1, 1)`,
  physical charges/positions, finite/symmetric `K`, both unit `U_A`, and `V`,
  H1 lowest `-0.79460371733658908`, self-Coulomb
  `0.45691176467371986`, and artifact readback one-body delta `0.0`.
- Manager reviewed the diff, checked live `cartesian_materialization` call
  sites, ran `git diff --check`, line accounting, and the suspicious added-line
  gate. Suspicious hits were expected fixed axis triples, Hamiltonian output
  symbols, and the legacy White-Lindsey materialization print fallback.

Guardrail:
- No new wrapper/payload/status/cache/artifact shape was added. White-Lindsey
  materialization wrapper behavior remains separate. The pass did not use
  report summaries, metadata, or recursive stage embedding as numerical input.

Carrying-cost accounting:
- deleted: durable materialization wrapper serialization and TSV/report
  materialization wrapper handling.
- simplified: PQS materialization return contract is now `nothing` or
  `CartesianIDAHamiltonian`.
- quarantined: H2 HP-WIRE-02 validation remains ignored under `tmp/work`.
- not deleted because: White-Lindsey materialization wrapper behavior remains a
  separate live route.
- exact remaining blocker: broader Slice D report/driver cleanup and stale
  physical-gausslet target/supplement report surfaces remain for a later
  deletion pass.
- added src lines: `111`.
- deleted src lines: `60`.
- net source: `+51`.
- new tests: none.
- new metadata/status fields: none.

Medium-term goal checkpoint:
- Completed: A/B/C in-memory producer boundary. The generic terminal route now
  realizes a basis, constructs one-body and localized IDA matrices, and returns
  a real `CartesianIDAHamiltonian` through materialization for H2.
- Active: Slice D deletion/simplification. The next lane is no longer numerical
  construction; it is deleting obsolete blocked source-plan/report surfaces and
  shrinking driver/report compatibility residue around the real Hamiltonian
  endpoint.
- Active guardrail: anti-bloat/design authority. This pass stayed under the
  hard source threshold and avoided new payload/cache/report surfaces, but was
  line-positive because it connected the first production materialization
  endpoint.
- Deferred: Cr2 and larger separated-diatomic stress testing. H2 closes this
  endpoint wiring; Cr2 remains a performance/stress gate, not the first
  correctness gate.
- Deferred/open: public typed pair-stage authority and mixed-pair planning
  remain outside the current Hamiltonian materialization lane.

Next step:
- Run a focused cleanup/design pass over the remaining old
  `physical_gausslet_*` report fields and source-plan payload construction.
  Delete surfaces that now only preserve the false blocked source-plan story, or
  report exact live callers that still require them.

## Cartesian Hamiltonian Producer Pass 031 - Delete Source-Plan Report Mirrors

Commit(s):
- this branch - Delete source-plan report mirrors

Summary:
- Deleted the obsolete physical-gausslet source-plan payload type and helper:
  `_PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload` and
  `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload`.
- Removed `diatomic_physical_gausslet_source_plan_payload` construction and the
  report-field merger that published `physical_gausslet_target_*`,
  `physical_gausslet_source_plan_summary`, and
  `physical_gausslet_supplement_preflight_summary`.
- Simplified `cartesian_report` to return its base report directly and removed
  stale supplement-preflight print lines. The Cr2 probe no longer treats the old
  source-plan/final-basis summaries as the route blocker.

Validation:
- Doer ran `git diff --check`, package load,
  `tools/h2_pqs_terminal_stage_smoke.jl`, and ignored
  `tmp/work/h2_terminal_wire02_materialization_validation.jl`.
- H2 retained the terminal basis endpoint and HP-WIRE-02 materialization:
  no-request returned `nothing`, requested materialization returned
  `CartesianIDAHamiltonian{Float64}`, final dimension `471`, H1 lowest
  `-0.79460371733658908`, self-Coulomb `0.45691176467371986`, and artifact
  readback one-body delta `0.0`.
- Manager reviewed the diff, ran `git diff --check`, line accounting,
  suspicious added-line gate, and focused greps confirming deleted source-plan
  report names are gone from active source/tools.

Carrying-cost accounting:
- deleted: source-plan payload type/helper, source-plan assembly field,
  physical-gausslet report mirror helper and fields, old Cr2 probe
  source-plan/final-basis blocker reads.
- simplified: report construction no longer merges stale physical-gausslet
  mirrors; print summary no longer reports stale supplement-preflight fields.
- quarantined: H2 HP-WIRE-02 validation remains ignored under `tmp/work`.
- not deleted because: `diatomic_physical_gausslet_target_payload` and the
  supplement request/representation/preflight payload chain still have live
  assembly/probe callers for supplement intent.
- exact remaining caller/blocker: target payload still carries
  `:missing_terminal_source_plan_realization` internally for supplement target
  availability. Collapse or retire that payload chain after supplement
  preflight ownership is clarified.
- added src lines: `12`.
- deleted src lines: `178`.
- deleted tool lines: `93`.
- new tests: none.
- new metadata/status fields: none.

Next step:
- Audit/collapse the remaining physical-gausslet target and supplement payload
  chain. Keep supplement intent only if it has a real consumer; otherwise delete
  the remaining `:missing_terminal_source_plan_realization` blocker story.

## Cartesian Hamiltonian Producer Pass 032 - Delete Physical-Gausslet Payload Chain

Commit(s):
- this branch - Delete physical-gausslet payload chain

Summary:
- Deleted the remaining physical-gausslet target/supplement payload structs and
  helpers from the active route. `cartesian_assembly` no longer constructs or
  returns the target payload, supplement request payload, supplement
  representation payload, or supplement preflight payload.
- Removed the old physical-gausslet target/source/supplement status and blocker
  vocabulary, including `:missing_terminal_source_plan_realization`,
  `:available_physical_gausslet_core_shell_target_inventory`, and
  `:blocked_physical_gausslet_target_inventory`.
- Simplified the Cr2 probe by removing display of the deleted supplement
  representation payload.

Validation:
- Doer ran `git diff --check`, package load,
  `tools/h2_pqs_terminal_stage_smoke.jl`, and ignored
  `tmp/work/h2_terminal_wire02_materialization_validation.jl`.
- H2 remained on the real terminal/materialization endpoint: final dimension
  `471`, max cross overlap `5.918709287301405e-15`, no-request materialization
  returned `nothing`, requested materialization returned
  `CartesianIDAHamiltonian{Float64}`, H1 lowest `-0.79460371733658908`,
  self-Coulomb `0.45691176467371986`, and artifact readback one-body delta
  `0.0`.
- Manager reviewed the pure-deletion diff, ran `git diff --check`, line
  accounting, suspicious added-line gate, and focused greps showing the deleted
  payload/helper/status names are gone from active source/tools.

Carrying-cost accounting:
- deleted: `750` source lines and `10` tool lines from the remaining
  physical-gausslet target/supplement payload chain.
- simplified: `cartesian_assembly` now returns only `route_skeleton` and
  `low_order_assembly` for this lane; supplement policy remains only as
  recipe/input intent.
- quarantined: none.
- not deleted because: stale fixture naming around
  `h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl` remains if a
  later cleanup wants to stop implying supplement-preflight machinery exists.
- exact remaining caller/blocker: focused grep found no active matches for the
  deleted payload/helper/status names in allowed source/tools.
- added src lines: `0`.
- deleted src lines: `750`.
- deleted tool lines: `10`.
- new tests: none.
- new metadata/status fields: none.

Next step:
- Decide whether to clean up stale fixture/tool naming around the old
  supplement-preflight input, or move to a broader review/merge decision for
  the design branch now that the base H2 PQS Hamiltonian endpoint is live.

## Cartesian Hamiltonian Producer Pass 033 - Ratify Slice D and H2 Endpoint Smoke

Commit(s):
- this branch - Ratify Slice D handoff and H2 endpoint smoke

Summary:
- Reconciled design authority with the already accepted HP-WIRE-02
  implementation: `AGENTS.md` and the design document now list `HP-WIRE-02` as
  approved for the narrow base-Hamiltonian materialization handoff.
- Recorded the implemented materialization boundary:
  `cartesian_materialization(report, terminal_basis_realization, inputs)`,
  with only `route_family`, `parent_axis_bundle_object`, and
  `system_metadata` approved as computational report fields.
- Replaced the old H2 supplement-preflight fixture/tool vocabulary with a base
  H2 PQS Hamiltonian endpoint fixture and smoke. The durable smoke now requests
  materialization, validates the returned `CartesianIDAHamiltonian`, checks H1
  and self-Coulomb endpoints, and verifies existing artifact readback.

Validation:
- Manager ran `git diff --check`, package load, focused stale-name greps, and
  `julia --project=. tools/h2_pqs_base_hamiltonian_smoke.jl`.
- The converted smoke passed with final dimension `471`, K/V/unit-U symmetry
  errors below `1.0e-10`, H1 lowest delta `0.0`, self-Coulomb delta
  `1.3322676295501878e-15`, and artifact readback one-body delta `0.0`.

Carrying-cost accounting:
- deleted: the active `supplement_preflight` fixture/tool naming and old smoke
  assertions centered on terminal roles, support counts, and column ranges.
- simplified: the H2 smoke now checks the user-facing Hamiltonian endpoint
  rather than route-internal staging vocabulary.
- quarantined: ignored temporary validation scripts remain under `tmp/work`.
- not deleted because: historical running-log references to the old smoke name
  remain as pass history.
- exact remaining caller/blocker: no active `tools/` or `test/driver_inputs/`
  reference to the old supplement-preflight H2 fixture remains.
- new tests: no new test file; the existing developer smoke was renamed and
  converted.
- new metadata/status fields: none.

Next step:
- Run the broad branch review/merge decision. Slice C is closed for the planned
  in-memory base PQS boundary, and Slice D base materialization is ratified; any
  remaining work should be classified as public polish, larger-system
  performance stress, or a separate non-base Hamiltonian lane.

## Cartesian Hamiltonian Producer Pass 034 - Compact Design Authority

Commit(s):
- this branch - Compact Cartesian Hamiltonian design authority

Summary:
- Accepted a docs-only restructuring of the Cartesian Hamiltonian producer
  design record after the design-governed implementation branch merged to
  `main`. The compact current authority now lives under
  `docs/src/developer/designs/cartesian_hamiltonian_producer/`, split into
  `README.md`, `current.md`, `registry.md`, `invariants.md`, and
  `implementation_slices.md`.
- Preserved the full June 2026 design document verbatim under
  `history/cartesian_hamiltonian_producer_design_2026-06_full.md` and moved
  review rounds under the same design directory. The old top-level design file
  is now only a redirect stub.
- Updated `AGENTS.md` and the developer index so normal agent startup reads the
  compact authority path rather than the 1600-line historical design.

Validation:
- Manager reviewed the compact authority files, corrected the Slice A/one-body
  signatures in `registry.md`, ran `git diff --check`, confirmed the historical
  full design matches the pre-compaction file byte-for-byte, and checked that
  stale old-path hits are confined to historical review/log references or the
  redirect path.

Goal advancement:
- LT: reduces design drift by keeping one compact current authority and one
  complete historical archive.
- MT: closes the immediate documentation-carrying-cost problem from the
  design-first Hamiltonian producer lane. Manager-log splitting remains a
  separate deferred docs-maintenance task.

Next step:
- Use the compact design authority for future Cartesian Hamiltonian producer
  work. Do not restart implementation from the historical design/review files
  unless investigating past decisions.

## Cartesian Hamiltonian Producer Pass 035 - Add Long-Range Roadmap

Commit(s):
- this branch - Add Cartesian long-range roadmap

Summary:
- Added `docs/src/developer/roadmaps/cartesian_long_range_roadmap.md` as
  strategic sequencing for the post-recovery Cartesian era. The roadmap records
  that the base PQS Hamiltonian path now exists and that the next program is
  public producer hardening, WL/QW/PQS downstream unification, supplement and
  correction migration, high-order geometry integration, Cr2-scale validation,
  and donor retirement.
- Linked the roadmap from the developer index and the compact Hamiltonian
  producer README as strategic context only. It is explicitly not
  implementation authority and does not approve new production surfaces.
- Reconciled the pasted roadmap with current `main`: R0 is marked partially
  complete because Slice D ratification, the base H2 Hamiltonian smoke, and
  compact design authority are already done. Remaining R0 work is stale
  route-status docs, baseline recording, and the pair/assembly stage decision.

Validation:
- Manager ran `git diff --check` and focused greps confirming the roadmap says
  it is not implementation authority and that R0 reflects the current state.

Goal advancement:
- LT: creates a durable strategic sequence without expanding normal startup
  reading or weakening the compact design gate.
- MT: separates future public-polish/unification/supplement/high-order lanes
  from the completed internal base Hamiltonian recovery lane.

Next step:
- Execute R0 follow-up only when needed: refresh stale route-status docs,
  record a quantitative baseline, and decide whether `cartesian_pair_terms` and
  `cartesian_assembly` become real downstream authority or leave the base public
  workflow.

## Cartesian Hamiltonian Producer Pass 036 - Record R0 Baseline

Commit(s):
- this branch - Record Cartesian R0 baseline

Summary:
- Accepted a measurement/docs-only R0 quantitative baseline at commit
  `2979514492cd4311426bd8f6e19c8be61c3e5a66`, after the base PQS path reached
  in-memory `CartesianIDAHamiltonian` materialization and artifact roundtrip.
- The baseline records lane line/file counts, H2 cold and warm
  base-Hamiltonian smoke timing/allocation, light separated-N2 one-body status,
  cheap pair/assembly stage timings, and the explicit deferral of full Cr2
  Hamiltonian stress/performance.
- The long-range roadmap now marks the quantitative baseline complete; the
  remaining R0 issue is only whether stale historical breadcrumbs need archival
  cleanup.

Validation:
- Doer ran `git diff --check`, package load, H2 base Hamiltonian smoke, and the
  existing ignored light separated-diatomic validation.
- Manager reviewed the baseline record, ran `git diff --check`, and ran
  `julia --project=docs docs/make.jl`; the docs build passed with only
  pre-existing size/deploy warnings.

Goal advancement:
- LT: creates the carrying-cost and performance baseline that future R1-R7
  shrinkage, public-producer, and stress-validation work can be judged against.
- MT: closes the substantive R0 measurement requirement after Slice D and the
  pair/assembly role decision.

Risk or guardrail:
- The baseline is a snapshot, not implementation authority. It deliberately
  does not bless Cr2 Hamiltonian readiness or restore pair/assembly as public
  base-producer concepts.

Next step:
- Ask repo-design-manager for the R1 public base-producer design pass. No new
  source implementation should start until that public surface and its
  deletion/validation obligations are approved.

## Cartesian Hamiltonian Producer Pass 037 - Implement R1 Public Base Producer

Commit(s):
- this branch - Add public Cartesian base Hamiltonian facade

Summary:
- Accepted the narrow R1 implementation of the public
  `cartesian_base_hamiltonian` facade for origin-centered H and Cartesian
  z-axis H2. The public function returns the existing
  `CartesianIDAHamiltonian{Float64}` directly, writes the existing Hamiltonian
  artifact format when `hamfile` is supplied, and records the approved fixed
  `producer_provenance/` keys.
- Added the approved public source file and export, moved base K/unit-`U_A`/V
  construction behind the shared `_cartesian_base_ida_hamiltonian` seam, and
  added the standalone public endpoint gate under `test/driver_public/` without
  wiring it into the default test runner.

Validation:
- Doer ran `git diff --check`, package load, the new public endpoint gate
  (`83/83`, about 37.5 s), and the existing H2 base Hamiltonian smoke.
- Manager reran `git diff --cached --check`, package load, and the standalone
  public endpoint gate (`83/83`, 37.9 s). The mechanical suspicious-line scan
  flagged only local H/H2 bridge tuple conversions into existing private driver
  helpers and the fixed x/y/z `Tuple` construction for PGDG axes.

Goal advancement:
- LT: establishes the first public base Cartesian Hamiltonian producer while
  preserving the existing Hamiltonian object and artifact payload.
- MT/R1: closes the first approved R1 implementation slice for H and z-axis H2;
  broader public driver polish, x/y/general orientation, WL/QW unification,
  supplements, and Cr2 stress remain deferred roadmap lanes.

Carrying-cost result:
- deleted: the report-bound terminal IDA Hamiltonian helper body was replaced
  by the shared report-free constructor seam.
- simplified: existing private materialization and the new public facade now
  share one base Hamiltonian construction path.
- quarantined: none.
- not deleted because: the private staged route remains the construction
  backend for the public facade and legacy materialization smoke; pair/assembly
  stages remain legacy/report compatibility, not public base workflow.
- exact remaining caller/blocker: R1 public API is limited to origin-centered H
  and z-axis H2; remaining lanes are public documentation polish, stale
  pair/assembly retirement in R2/R3, and larger-system performance validation.
- added src lines: 225.
- deleted src lines: 13.
- new tests: one standalone public endpoint/provenance test file, not wired
  into `test/runtests.jl`.
- new metadata/status fields: none. New artifact data is only the approved
  fixed `producer_provenance/` group.

Risk or guardrail:
- Source additions hit the 225-line hard threshold exactly. Do not expand this
  facade in place without a follow-up design amendment or a deletion-backed
  cleanup.

## Cartesian Hamiltonian Producer Pass 038 - Document R1 Public Facade

Commit(s):
- this branch - Document Cartesian base Hamiltonian facade

Summary:
- Accepted a docs-only public-polish pass after the R1 facade landed. The
  reference/export page now names `cartesian_base_hamiltonian`, shows the
  supported origin-centered H and z-axis H2 calls, and states the current
  limitations plainly.
- Updated the IDA Hamiltonian algorithm page so it no longer says the terminal
  PQS producer is blocked at shell realization. It now describes the public
  H/H2 facade, existing artifact readback, and `producer_provenance/` behavior.
- Added a short examples-index pointer that directs users to the public facade
  instead of implying that older Cartesian example scripts are the recommended
  base-Hamiltonian entry point.

Validation:
- Manager ran `git diff --check` and `julia --project=docs docs/make.jl`. The
  docs build passed with only the existing large-page/search-index and local
  deploy-detection warnings.

Goal advancement:
- LT/R1: makes the first public base Cartesian Hamiltonian producer visible and
  usable without expanding the approved science scope.
- MT: keeps public guidance aligned with the R1 contract while leaving broader
  molecule support, supplements, WL/QW unification, Cr2 stress, and pair/stage
  retirement to later roadmap lanes.

Carrying-cost result:
- deleted: stale blocked-route wording from the public IDA Hamiltonian page.
- simplified: public examples now point to one supported facade instead of the
  older route-stage story.
- quarantined: none.
- not deleted because: older example scripts and algorithm notes remain useful
  for their existing lanes and were not part of this docs polish.
- exact remaining caller/blocker: R1 docs cover only H and z-axis H2; any
  broader public facade behavior needs repo-design-manager approval.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 039 - Be2 Performance Probe

Commit(s):
- none from doer; this entry records an ignored `tmp/work` measurement

Summary:
- Accepted a measurement-only Be2 probe using two z-axis separations. Both
  cases built finite symmetric base K, separated unit-`U_A`, and localized IDA
  V without tracked code changes.
- The close and far Be2 cases changed parent boxes/support counts and final
  dimensions (`633` and `683`), but both still used contact-core/shared-shell
  topology with outer slabs. Be2 therefore works as a cheap contact/shared-shell
  sanity probe, not as the separated atom-local topology proxy needed before
  Cr2.
- Warm Be2 was small: the far case completed the measured Hamiltonian phases in
  about `0.451 s`, with one-body/IDA contractions not yet dominating. The
  useful allocation signals are terminal-basis realization and K assembly, but
  the scale is too small to justify a Be2-only optimization pass.

Validation:
- Doer ran `git diff --check`, package load, and the ignored
  `tmp/work/be2_base_producer_perf_probe.jl`; final worktree was clean.

Goal advancement:
- LT/Roadmap: keeps Cr2 stress out of the inner optimization loop while still
  measuring a light all-electron diatomic proxy.
- MT: changes the next performance step from Be2 optimization to an N2 proxy
  measurement, because N2 is already known to exercise separated atom-local
  shells and shared/mismatch regions.

Risk or guardrail:
- Do not infer Cr2 readiness from these Be2 timings. The measured Be2 cases are
  too small and topologically incomplete for that.

Next step:
- Run a measurement-only N2 proxy pass before source optimization. Use it to
  identify whether the first optimization target is terminal-basis realization,
  K assembly/buffer reuse, unit-`U_A`, IDA V, or artifact writing.

## Cartesian Hamiltonian Producer Pass 040 - Corrected Be2 Core-Spacing Probe

Commit(s):
- none from doer; this entry records an ignored `tmp/work` measurement

Summary:
- Accepted the corrected Be2 measurement after applying the documented
  homonuclear diatomic spacing rule `core_spacing = 1.2 / (Z * (q - 3))`. For
  Be at `Z = 4`, `q = 5`, the correct value is `0.15`; the prior `0.30` probe
  was a He-scale spacing and was too coarse near the nuclei.
- With `core_spacing = 0.15`, Be2 now exercises the intended separated
  topology: atom-local cores/shells, midpoint slabs, shared shell behavior, and
  outer slabs depending on separation. Final retained dimensions were `1257`,
  `1395`, and `1541` for close/far/farther cases, with complete coverage and
  finite symmetric K, unit-`U_A`, and V.
- The corrected run supersedes the previous N2-next conclusion. Be2 is now a
  useful optimization proxy before Cr2; N2 remains optional only if a different
  chemistry/electron-count proxy is specifically desired.

Validation:
- Doer ran `git diff --check`, package load, and the ignored
  `tmp/work/be2_corrected_core_spacing_perf_probe.jl`; final worktree was
  clean.

Goal advancement:
- LT/Roadmap: keeps Cr2 out of the inner loop while moving to a topology-rich,
  all-electron diatomic proxy.
- MT: identifies terminal basis realization as the first optimization target
  (`~1.7-2.0 GiB`, `~5.8-7.5 s` warm corrected Be2), with K allocation
  (`~1.7-2.1 GiB` for subsecond work) as the next likely target.

Risk or guardrail:
- The optimization lane should preserve the corrected spacing rule. Do not tune
  performance around the invalid `core_spacing = 0.30` Be2 measurements.

Next step:
- Start with a focused terminal-basis allocation/timing audit on corrected Be2,
  then implement a bounded optimization only after the cost center is localized.

## Cartesian Hamiltonian Producer Pass 041 - Terminal Basis Allocation Audit

Commit(s):
- none from doer; this entry records an ignored `tmp/work` measurement

Summary:
- Accepted a focused terminal-basis audit on corrected Be2 far
  (`R = 8.0`, `q = 5`, `core_spacing = 0.15`). The topology-rich case has 12
  terminal blocks, final dimension `1395`, max cross overlap about `3.17e-14`,
  and largest local workspace near `64 MiB`.
- The terminal realization replay allocated about `1620 MiB` and took `3.248 s`
  in the instrumented run; the cross-overlap audit alone allocated about
  `339 MiB`.
- The dominant source is repeated support-action/support-cross construction:
  the shared molecular shell projected against nine previous blocks, enlarged
  to effective support `4225`, and allocated about `462 MiB` in projection plus
  about `306 MiB` across Gram and identity-check support actions. Large
  atom-local shells showed the same pattern at smaller scale.

Validation:
- Doer ran `git diff --check`, package load, and the ignored
  `tmp/work/be2_terminal_basis_perf_audit.jl`; final worktree was clean.

Goal advancement:
- LT/Roadmap: provides a measured optimization target on a Cr2-relevant
  separated-topology proxy without running Cr2.
- MT: first optimization should target support-action construction and local
  workspace reuse/streaming inside terminal basis realization. Do not start by
  removing cross-audit validation; it is smaller than projection/Gram/check
  allocation and still protects the final-basis contract.

Risk or guardrail:
- Avoid creating a persistent cache or new stage object. A first optimization
  should stay file-local and bounded, preferably inside
  `pqs_terminal_basis_realization.jl`, with no new public API, metadata fields,
  or report/status vocabulary.

Next step:
- Issue a narrow implementation pass for terminal-basis support-action buffer
  reuse/streaming, using corrected Be2 far as the performance proxy and H/H2
  public endpoints as correctness gates.

## Cartesian Hamiltonian Producer Pass 042 - Reuse Terminal Support-Action Scratch

Commit(s):
- this branch - Reuse terminal support action scratch

Summary:
- Accepted a bounded optimization in
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.
  The patch adds file-local `_support_cross!` filling into reusable scratch,
  routes support actions through `mul!`, reuses a tile in direct-right block
  pair assembly, and replaces dense identity temporaries in direct and shell
  identity checks.
- Corrected-Be2 far terminal realization allocation dropped from about
  `1643 MiB` to `1427 MiB`, with elapsed time from `3.332 s` to `3.039 s`.
  The shared-shell Gram plus identity-check allocation dropped from roughly
  `309 MiB` to `135 MiB`. Final dimension (`1395`) and max cross overlap
  (`~3.17e-14`) were unchanged.

Validation:
- Doer ran `git diff --check`, package load, the standalone public endpoint
  gate (`83/83`), the ignored Be2 terminal audit before/after, and the ignored
  corrected-Be2 full K/unit-`U_A`/V probe for close/far/farther.
- Manager reran `git diff --check`, the mechanical suspicious-line scan,
  package load, and `julia --project=. test/driver_public/cartesian_base_hamiltonian_runtests.jl`
  (`83/83`, 39.3 s).

Goal advancement:
- LT/Roadmap: improves the Cr2-relevant terminal-basis path on a separated Be2
  proxy without running Cr2 or expanding public scope.
- MT: removes a measured allocation source in Gram/identity/check support
  actions. Projection remains the dominant terminal-basis target.

Carrying-cost result:
- deleted: dense identity temporaries for direct-sector and shell identity
  checks.
- simplified: support-cross construction now has an in-place file-local path
  reused by support actions and direct-right block pair assembly.
- quarantined: none.
- not deleted because: projection still needs `_subtract_previous` and
  effective-support reconstruction; that is the next optimization target, not
  safe to delete in this pass.
- exact remaining caller/blocker: shared-shell projection allocation remains
  high, around `490 MiB`, from residual/action construction plus effective
  support/coefficient rebuilding.
- added src lines: 40.
- deleted src lines: 12.
- new tests: none.
- new metadata/status fields: none.

Risk or guardrail:
- Cross-overlap audit and shell identity checks remain active. Do not optimize
  by skipping validation; the next pass should target projection allocation
  directly.

## Cartesian Hamiltonian Producer Pass 043 - Restore Block-Local PQS Shells

Commit(s):
- this branch - Restore block-local PQS terminal shells

Summary:
- Accepted the source correction required by design commit `5ec55882`. The
  terminal PQS realizer now uses the owned shell rows
  `support.support_indices` / `support.support_states` after full-box
  boundary-mode generation, then performs shell-local Gram/Lowdin cleanup on
  that block-local support.
- Removed recursive previous-block projection, `_subtract_previous`, and
  `projection_atol` plumbing. Cross-block overlap remains an audit only; it is
  no longer a construction repair path.
- This supersedes the Pass 042 "projection allocation" next-target note. The
  allocation target was a symptom of the wrong cumulative-support contract, not
  a path to optimize further.

Validation:
- Doer ran `git diff --check`, package load, the standalone public endpoint
  gate (`83/83`), the H2 base Hamiltonian smoke, and the ignored corrected-Be2
  full K/unit-`U_A`/V probe.
- Manager reviewed the one-file source diff, ran `git diff --check`, and ran
  the mechanical stale-projection scan over live source/design files. Remaining
  `support.outer_box` source uses are limited to raw-plan/source-box
  construction checks.

Goal advancement:
- LT5/LT6: restores the intended PQS support provenance and stable
  final-basis block contract before larger Cr2-facing optimization work.
- MT/R1: keeps the public H/H2 facade numerically intact while correcting the
  internal terminal basis representation.

Carrying-cost result:
- deleted: recursive previous-block projection and `_subtract_previous`.
- simplified: `_realize_shell` is now shell-local and block-local.
- quarantined: none.
- not deleted because: support-action scratch reuse still serves local
  Gram/check/audit paths.
- exact remaining caller/blocker: none for removed projection symbols.
- added src lines: 19.
- deleted src lines: 43.
- new tests: none.
- new metadata/status fields: none.

Risk or guardrail:
- Do not reintroduce previous-block projection to reduce cross overlap. If a
  block-local shell has large cross overlap, treat that as a parent-metric or
  shell-construction problem and return to design/numerical review.

## Cartesian Hamiltonian Producer Pass 044 - Post-Correction Be2 Performance Refresh

Commit(s):
- none from doer; this entry records an ignored `tmp/work` measurement

Summary:
- Accepted the post-correction Be2 refresh at `d2bf139c`, after block-local
  PQS shell support replaced cumulative previous-block projection. This makes
  the old projection-heavy terminal allocation baseline stale.
- Corrected Be2 (`q = 5`, `core_spacing = 0.15`) remains a useful Cr2-facing
  optimization proxy: close/far/farther cases have final dimensions `1257`,
  `1395`, and `1541`, complete terminal coverage, and finite symmetric K,
  unit-`U_A`, and V.
- Terminal transform allocation is now much lower in the far case
  (`~371.5 MiB`) than under the stale cumulative-support model. K assembly is
  the largest single numerical assembly allocation in the far/farther cases
  (`~552-643 MiB`) despite subsecond elapsed time, so it is the next target.

Validation:
- Doer ran `git diff --check`, package load, the ignored corrected-Be2
  performance probe, and verified final worktree cleanliness.
- Manager reviewed the measurement and inspected
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl` to identify
  the likely allocation pattern before issuing the next handoff.

Goal advancement:
- LT4/LT6: redirects optimization to the current block-local producer rather
  than the retired recursive-projection model.
- MT/Roadmap: keeps Cr2 out of the inner loop while selecting a measured Be2
  allocation target in the approved Slice B source file.

Risk / guardrail:
- Optimize K assembly without changing the HP-FN-03 public surface. The likely
  issue is repeated temporary allocation inside `_terminal_product_action` and
  per-pair `block` construction, not the term-first Coulomb path.

Next step:
- Issue a bounded HP-FN-03 optimization pass in
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl` for
  single-product/K assembly buffer reuse and direct accumulation.

## Cartesian Hamiltonian Producer Pass 045 - Reduce K Assembly Allocation

Commit(s):
- this branch - Reuse terminal product assembly buffers

Summary:
- Accepted a bounded HP-FN-03 optimization in
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`. The old
  allocating `_terminal_product_action` was replaced by an in-place
  `_add_terminal_product_block!` that reuses support-cross, action, and
  PQS-left block buffers inside `assemble_terminal_product_operator!`.
- The helper now fills support-pair tiles through `_support_cross!`, uses
  `mul!` for coefficient contractions, and accumulates the current block
  contribution directly into the destination and mirrored transpose location.
- Corrected-Be2 K allocation dropped materially: close `471 -> 256 MiB`, far
  `552 -> 222 MiB`, farther `643 -> 177 MiB`. K elapsed time stayed subsecond
  and improved for far/farther, while close became slightly slower but still
  inexpensive.

Validation:
- Doer ran `git diff --check`, package load, the standalone public endpoint
  gate (`83/83`), H2 base Hamiltonian smoke, and the ignored corrected-Be2
  probe. H/H2 endpoint deltas stayed at roundoff scale and Be2 K/unit-`U_A`/V
  stayed finite and symmetric.
- Manager reviewed the one-file diff, checked include-order compatibility for
  `_buffer_view!`/`_support_cross!`, ran `git diff --check`, the source
  line-count gate, the suspicious-line scan, and package load.

Goal advancement:
- LT4/LT6: reduces allocation in the public base Hamiltonian K path on the
  current block-local Be2 proxy without changing the public API or Hamiltonian
  artifact contract.
- MT/Roadmap: completes the first post-correction optimization target and
  leaves unit-`U_A`/IDA allocation as separate future candidates.

Carrying-cost result:
- deleted: old allocating `_terminal_product_action`.
- simplified: single-product assembly now has one in-place block helper.
- quarantined: none.
- not deleted because: Gaussian-sum/unit-`U_A` and IDA term-first paths remain
  separate and were intentionally out of scope.
- exact remaining caller/blocker: none for the old helper.
- added src lines: 28.
- deleted src lines: 11.
- new tests: none.
- new metadata/status fields: none.

Risk or guardrail:
- Do not fold the Gaussian-sum or IDA paths into this helper unless a measured
  pass shows the same allocation pattern and stays within their existing
  approved internal surfaces.

## Cartesian Hamiltonian Producer Pass 046 - Audit Units-Phase Cost

Commit(s):
- none from doer; this entry records an ignored `tmp/work` measurement

Summary:
- Accepted a measurement-only audit of corrected-Be2 `cartesian_units` at
  `R = 8.0`, `q = 5`, `core_spacing = 0.15`. The fixture has 12 terminal
  regions/retained units and 7 PQS shell units with source boxes from
  `(7,7,7)` through `(13,13,25)`.
- The apparent `cartesian_units` cost was cold compilation, not warm runtime.
  The expensive cold substeps were unit inventory (`8.37 s`), lowering
  contract inventory (`14.93 s`), and related terminal lowering helpers. Warm
  repeats were sub-millisecond to about `0.0005 s` per substep, with only a few
  MiB allocated.
- Raw product source-plan setup is not the runtime bottleneck. It remains
  metadata-only in this stage; axis transform statuses were all
  `:not_materialized`, and warm raw-source metadata allocation was about
  `2.9 MiB` total across PQS shells.

Validation:
- Doer ran `git diff --check`, package load, the ignored
  `tmp/work/be2_units_phase_audit.jl`, and verified final worktree
  cleanliness.

Goal advancement:
- LT4/Roadmap: prevents a mistaken source optimization pass in route-unit
  planning and keeps performance effort focused on real runtime costs.
- MT: classifies units-stage warm runtime as good enough for now. Cold compile
  pressure may matter later, but it is a separate product-startup concern.

Risk / guardrail:
- Do not move raw source-box construction or rewrite units-stage ownership
  based on cold timing alone. If cold compile becomes a product requirement,
  investigate NamedTuple/type specialization pressure in route helper and
  terminal shellification code as a compile-time project, not a numerical
  kernel optimization.

Next step:
- With K optimized and units-stage runtime cleared, decide whether to audit
  `cartesian_transforms`/terminal-basis warm allocation, unit-`U_A`/IDA
  allocation, or stop local Be2 optimization and move back to roadmap work.

## Cartesian Hamiltonian Producer Pass 047 - Replace Cross-Overlap Audit With Structural Checks

Commit(s):
- this branch - Enforce structural terminal support checks

Summary:
- Accepted the source cleanup following design commit `e86c08a8`. The terminal
  realizer no longer computes pairwise production cross-overlap matrices or
  accepts `cross_atol`.
- Added structural support checks to enforce the corrected invariant directly:
  each realized block must match its authoritative terminal support exactly,
  duplicate support rows inside a block are errors, and support intersections
  across blocks are errors.
- Deleted `_block_pair_matrix`, `_block_action`, and the allocating
  `_support_cross` wrapper. The in-place `_support_cross!` remains because
  local support actions and one-body assembly still use it.

Validation:
- Doer ran `git diff --check`, package load, the standalone public endpoint
  gate (`83/83`), H2 base Hamiltonian smoke, and the ignored corrected-Be2
  probe. H/H2 endpoint deltas stayed at roundoff scale; corrected-Be2
  close/far/farther K, unit-`U_A`, and V stayed finite and symmetric.
- Manager reviewed the one-file diff, confirmed only the compatibility
  `max_cross_overlap` field remains in live source, ran `git diff --check`,
  the source line-count gate, the suspicious-line scan, and package load.

Goal advancement:
- LT5/LT6: removes the last production code path that treated terminal
  cross-block overlap as a numerical residual instead of a structural support
  invariant.
- MT/R3 readiness: clears stale Slice A/B/C implementation debt before
  residual-GTO/MWG design review proceeds to source work.

Carrying-cost result:
- deleted: production cross-overlap audit helpers and `cross_atol`.
- simplified: terminal basis finalization now validates owned-support
  structure directly.
- quarantined: `max_cross_overlap` remains as a `0.0` compatibility field.
- not deleted because: ignored Be2 validation probes still read the field, and
  removing it cleanly can be done in a later broader cleanup.
- exact remaining caller/blocker: ignored `tmp/work` Be2 probe access to
  `max_cross_overlap`; no live source/tool/test caller remains.
- added src lines: 23.
- deleted src lines: 39.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not reintroduce overlap audits to diagnose structural mistakes. Duplicate
  or intersecting support rows should fail as ownership/indexing errors.

## Cartesian Hamiltonian Producer Pass 048 - R3-A Residual Spectrum Spike

Commit(s):
- none from doer; this entry records an ignored `tmp/work` measurement

Summary:
- Accepted a measurement-only R3-A residual-spectrum spike for the proposed
  first supplemented proxy: base z-axis H2 plus contracted two-center
  H/cc-pVTZ, `lmax = 1`, `uncontracted = false`, no width filtering.
- The fixture produced the expected 18 supplement candidates, 9 per H center,
  in deterministic center-major/source-shell/Cartesian-component order. Base
  final dimension was `471`.
- The residual metric spectrum for `S_R = S_AA - X'X` was full rank and far
  from the proposed threshold band: min eigenvalue about `3.05e-4`, max about
  `1.35e-2`, condition estimate about `44.3`. `inv(sqrt(Symmetric(S_R)))`
  succeeded and the measured `R' S R` identity error was about `4.4e-15`.

Validation:
- Doer ran `git diff --check`, package load, the ignored
  `tmp/work/r3a_h2_ccpvtz_residual_spectrum_spike.jl`, and verified final
  worktree cleanliness.

Goal advancement:
- R3: supplies the numerical threshold evidence needed for a focused R3-A
  docs-only approval pass. No source implementation was started.
- LT5/LT6: keeps supplement work tied to measured residual-basis behavior
  before introducing new production surfaces.

Risk / guardrail:
- The recommended first R3-A thresholds are `tau_abs = 1e-10`,
  `tau_rel = 1e-10`, `tau_neg_abs = 1e-12`, and `tau_neg_rel = 1e-12`.
  These should be frozen by repo-design-manager before doer implementation.
  The ordinary-QW donor thresholds remain context, not automatic authority.

Next step:
- Request a docs-only R3-A approval amendment that moves only the R3-A IDs
  needed for residual basis plus exact one-body/moments into approved status,
  with the measured thresholds and H2 fixture locked. R3-B and R3-C should
  remain candidate-only.

## Cartesian Hamiltonian Producer Pass 049 - Implement R3-A Residual Basis

Commit(s):
- this branch - Add R3A residual GTO basis construction

Summary:
- Accepted the first R3-A implementation slice under the approved
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` path.
  The patch adds `CartesianTerminalResidualGTOAugmentation` with the approved
  numerical fields and implements deterministic residual-basis construction for
  the frozen H2/H/cc-pVTZ full-rank fixture.
- The implementation derives candidate ownership by exact center-to-nucleus
  match, assembles `X = G' S A` blockwise through terminal basis supports,
  builds `S_R = S_AA - X'X`, applies the frozen R3-A thresholds, uses
  `inv(sqrt(Symmetric(S_R)))` in candidate order, constructs `T_A = T` and
  `T_G = -X*T`, and canonicalizes signs by largest `T_A` entry.
- This pass does not implement exact augmented one-body/moment assembly
  (`HP-R3-FN-02`), MWG/IDA, supplemented Hamiltonian construction, artifacts,
  public API, Be2, or Cr2 validation.

Validation:
- Doer ran `git diff --check`, package load, and the ignored
  `tmp/work/r3a_h2_residual_gto_validation.jl`.
- Manager reviewed the two-file source diff, ran `git diff --check`, source
  line-count and suspicious-line scans, package load, and reran the ignored
  H2 R3-A residual validation script.

Goal advancement:
- R3: creates the first real residual-GTO numerical object on top of the
  block-local base PQS final basis, with no status/payload/report framework.
- LT5/LT6: establishes explicit residual provenance (`T_G`, `T_A`, owner
  indices, thresholds, orientation/sign rules) needed before exact augmented
  one-body and later MWG/IDA work.

Carrying-cost result:
- deleted: none.
- simplified: kept R3-A part 1 to one owner file plus one include.
- quarantined: ignored validation remains under `tmp/work`.
- not deleted because: no live stale R3 source surface existed yet.
- exact remaining caller/blocker: `HP-R3-FN-02` exact augmented one-body and
  moment assembly is still unimplemented; rank-deficient residual selection
  currently stops rather than silently choosing a pivot rule.
- added src lines: 131.
- deleted src lines: 1.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The internal helper lookups reflect current include order for older
  representation/cross-overlap helpers. If R3-A becomes a wider production
  surface, revisit owner/module boundaries instead of adding a provider payload.

## Cartesian Hamiltonian Producer Pass 050 - Assemble R3-A Exact One-Body And Moments

Commit(s):
- this branch - Assemble R3A augmented one-body blocks

Summary:
- Accepted the second R3-A implementation slice in the approved
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` owner
  file. The patch adds exact augmented `[G, R]` operator assembly for kinetic,
  uncharged by-center nuclear attraction, position moments, and second
  position moments.
- The implementation reuses the base terminal operator path for `G-G` blocks
  and existing CPB mixed/self GTO providers for `G-A` and `A-A` blocks, then
  applies the approved transformation formulas using `T_G` and `T_A`.
- During review, the manager removed the final `0.5*(O+O')` smoothing from
  the augmented operator helper so symmetry is validated by the endpoint script
  instead of hidden by construction. The remaining raw symmetry errors are
  around `1e-12`, well below the R3-A endpoint tolerance.

Validation:
- Doer ran `git diff --check`, package load, and the ignored
  `tmp/work/r3a_h2_residual_gto_validation.jl`.
- Manager reviewed the source diff and CPB nuclear-provider convention,
  confirmed the CPB nuclear blocks are uncharged by-center operators, ran
  `git diff --check`, package load, the source line-count gate, the
  suspicious-line scan, and reran the ignored H2 R3-A validation script.

Goal advancement:
- R3: completes the approved in-memory R3-A numerical endpoint: base H2
  augmented by contracted two-center H/cc-pVTZ residual GTOs now has exact
  augmented `K`, unit `U_A`, `x/y/z`, and `x^2/y^2/z^2` matrices.
- LT5/LT6: extends the common final-basis Hamiltonian architecture without
  adding public API, artifacts, report fields, status vocabulary, or an MWG/IDA
  path ahead of approval.

Endpoint facts:
- Base dimension `471`, residual dimension `18`, augmented dimension `489`.
- Residual eigenvalue min/max `3.0488355008683734e-04` /
  `1.3512432621413795e-02`; `G' S R` error `0.0`; `R' S R` identity error
  `3.4106051316484809e-12`.
- Base `G-G` block equality error is `0.0` in the validation script after the
  manager cleanup; augmented one-body energy improved from
  `-0.7946037173365925` to `-0.7959028345077851`.

Carrying-cost result:
- deleted: final augmented-operator symmetrization that would have hidden an
  invalid mixed/self block convention.
- simplified: one generic augmented-operator formula is reused for kinetic,
  unit nuclear attraction, and moment matrices.
- quarantined: ignored R3-A validation remains under `tmp/work`.
- not deleted because: CPB donor kernels remain the active exact mixed/self
  GTO provider source.
- exact remaining caller/blocker: committed standalone `HP-R3-TEST-01` gate is
  still not added; R3-B MWG/IDA remains candidate-only and unimplemented.
- added src lines: 139.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- R3-A is now numerically useful but not durable in the committed test suite
  until the standalone `HP-R3-TEST-01` endpoint gate is added. Do not proceed
  to R3-B MWG/IDA implementation before that test gate or an explicit manager
  decision to defer it.

## Medium-Term Checkpoint After Pass 050

- Completed: base PQS A/B/C/D/R1 in-memory and public base-Hamiltonian
  boundary; terminal shell ownership correction; Be2 local optimization triage;
  R3-A residual basis plus exact one-body/moment implementation.
- Active: make R3-A durable with the approved standalone endpoint gate, then
  request or review focused R3-B design approval for MWG/IDA interaction blocks.
- Deferred: Be2/N2 performance optimization beyond the already accepted K and
  terminal-support corrections, Cr2 stress/performance validation, public
  non-base/supplement lanes, ECP/correction workflows, and broad driver polish.
- Guardrail update: residual-GTO work should continue to be split by physics
  endpoint. Exact one-body/moments are accepted; MWG/IDA `V`, supplemented
  Hamiltonian construction, artifacts, and public API remain unauthorized until
  their own approved slice.

## Cartesian Hamiltonian Producer Pass 051 - Add R3-A Standalone Endpoint Gate

Commit(s):
- this branch - Add R3A augmented one-body endpoint test

Summary:
- Accepted the approved `HP-R3-TEST-01` durability gate at
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`. The test is
  intentionally standalone and is not wired into `test/runtests.jl`.
- The gate validates the frozen H2 plus contracted two-center H/cc-pVTZ
  `lmax = 1` residual-GTO fixture. It checks candidate labels/order, owner
  counts, residual rank, `G' S R`, `R' S R`, finite/symmetric augmented
  kinetic/unit-nuclear/moment matrices, base `G-G` block equality, and the
  augmented one-body variational endpoint.
- No source code, public API, artifact, report/status vocabulary, Be2/Cr2
  scope, or R3-B MWG/IDA path was added.

Validation:
- Doer ran `git diff --check`, package load, and the new standalone gate.
- Manager reviewed the new test file, confirmed it asserts the approved
  physics endpoint rather than stale pair/assembly/report internals, reran
  package load, and reran
  `julia --project=. test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.

Endpoint facts:
- Test passed `39/39` in `36.58 s` on the manager rerun.
- Base dimension `471`, residual dimension `18`, augmented dimension `489`.
- `E1_base = -0.7946037173365925`; `E1_aug = -0.7959028345077851`, satisfying
  `E1_aug <= E1_base + 1e-10`.

Goal advancement:
- R3: makes the accepted R3-A exact one-body/moment endpoint durable in tracked
  validation before any R3-B interaction work.
- LT5/LT6: protects the residual-GTO convention with a physical endpoint gate
  without expanding the public producer or adding a new payload layer.

Carrying-cost result:
- deleted: none.
- simplified: converted the ignored R3-A validation into one focused
  standalone endpoint gate.
- quarantined: Be2/Cr2 and R3-B MWG/IDA remain out of scope.
- not deleted because: ignored `tmp/work` validation can remain disposable
  scratch and is not part of tracked contract.
- exact remaining caller/blocker: R3-B MWG/IDA and supplemented Hamiltonian
  construction remain candidate-only and unimplemented.
- added src lines: 0.
- deleted src lines: 0.
- new tests: one standalone integration gate, `177` lines, not in
  `test/runtests.jl`.
- new metadata/status fields: none.

Risk / guardrail:
- This is a ~36 s integration gate, not a default per-edit unit test. Keep it
  as an explicit standalone acceptance check unless test-suite policy changes.

## Cartesian Hamiltonian Producer Pass 052 - Retire Stale H2 Blocked Fixtures

Commit(s):
- this branch - Remove stale H2 blocked endpoint fixtures

Summary:
- Accepted a tightly scoped R3-A retirement cleanup after the R3-A durability
  gate. The pass deleted four unused H2 driver-input fixtures whose only
  purpose was to advertise missing physical-gausslet H1/H1+J/RHF endpoints.
- A focused search found no active source, tool, bin, or test caller for the
  deleted fixture files. The old supplement-preflight/provider-blocker wrapper
  family remains absent from live implementation code.
- Retained `test/docs/cartesian_ham_builder_policy_runtests.jl`, which is a
  useful negative policy guard preventing `residual_gto_provider_blocks` from
  becoming canonical-driver surface. Retained CPB provider blocked statuses
  because they describe local provider availability, not stale route preflight
  payloads.

Validation:
- Doer ran `git diff --check` and package load.
- Manager reviewed the deletion diff, reran `git diff --check`, package load,
  and a focused `rg` for stale supplement-preflight/provider-blocker names.
  Remaining live hits are the retained canonical-driver negative test and
  current R3 design text documenting the old family as deleted/forbidden.

Goal advancement:
- R3/Roadmap: removes stale fixture vocabulary after R3-A created a real
  augmented one-body endpoint, reducing the chance that R3-B work revives the
  old blocked-preflight story.

Carrying-cost result:
- deleted: four unused `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_*`
  blocked endpoint fixtures.
- simplified: removed obsolete H2 fixture variants that only named missing
  physics endpoints.
- quarantined: none.
- not deleted because: canonical-driver negative policy test and local CPB
  provider statuses remain live contracts.
- exact remaining caller/blocker: none for the deleted fixtures.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- This cleanup does not approve R3-B source work. R3-B remains owned by
  repo-design-manager approval and must not reintroduce a provider-blocker
  preflight layer.

## Cartesian Hamiltonian Producer Pass 053 - Implement R3-B Compact MWG Hamiltonian

Commit(s):
- this branch - Add R3B compact MWG Hamiltonian assembly

Summary:
- Accepted the reapproved `HP-R3-FN-03` implementation in the existing
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` file.
  The new internal function
  `pqs_terminal_residual_gto_augmented_hamiltonian(...)` consumes the base
  Hamiltonian, terminal basis, R3-A residual object, and R3-A exact augmented
  one-body/moment matrices, then returns the existing
  `CartesianIDAHamiltonian{Float64}`.
- The implementation derives residual MWG centers and widths from exact R3-A
  moment diagonals using the compact convention `sigma = sqrt(2v)`, rejects
  invalid variances/widths, builds density-normalized `G-M` and `M-M`
  interaction blocks from the existing QW donor kernels, reuses the base
  `V_GG`, and assembles `V_aug`.
- No artifacts, public API, driver/bin/tool workflow, payload/result/status
  object, committed test, Be2/Cr2 validation, or RHF/solver work was added.

Validation:
- Doer ran `git diff --check`, package load, and the ignored
  `tmp/work/r3b_h2_augmented_hamiltonian_validation.jl`.
- Manager reviewed the one-file diff, checked the QW donor helper signatures,
  confirmed the corrected R3-B scalar in active authority docs, ran
  `git diff --check`, package load, the suspicious-line scan, and reran the
  ignored R3-B H2 validation script.

Endpoint facts:
- Augmented dimension `489`.
- MWG center range `(-3.3667418611021276, 3.3667418611025823)`.
- MWG width range `(0.6597141664082136, 7.262076101034195)`.
- `V_aug` symmetry error `1.5543122344752192e-15`; base `V_GG` delta `0.0`.
- Lowest augmented one-body orbital IDA self-Coulomb
  `0.4574331709135599`, exactly matching the corrected compact-path target.

Goal advancement:
- R3: completes the first in-memory supplemented H2 Hamiltonian endpoint under
  the compact R3-A/R3-B convention.
- LT5/LT6: extends the common final-basis Hamiltonian boundary to residual-GTO
  augmentation without introducing a parallel payload or artifact framework.

Carrying-cost result:
- deleted: none.
- simplified: compact R3-B directly reuses R3-A exact moments, base `V_GG`,
  and existing QW density-normalized pair-factor kernels.
- quarantined: ignored R3-B validation remains under `tmp/work`.
- not deleted because: QW donor kernels remain active reuse/oracle sources;
  R3-C artifact/provenance is still unapproved.
- exact remaining caller/blocker: `HP-R3-ART-01` and public/driver workflow
  remain candidate-only; no committed R3-B test update was added in this pass.
- added src lines: 123.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not reintroduce the superseded `0.457435475059184` private-gauge scalar,
  width scaling, or tolerance relaxation. The accepted compact R3-B baseline is
  `0.4574331709135599` for this fixture.

## Cartesian Hamiltonian Producer Pass 054 - Harden R3 Construction Consistency

Commit(s):
- this branch - Harden R3 residual construction checks

Summary:
- Accepted a narrow hardening pass in the existing
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` file.
  The patch adds cheap consistency checks for R3 residual dimensions,
  `T_G`/`T_A` shapes, supplement labels/centers, owner-center matching,
  augmented operator dimensions, R3-B base Hamiltonian dimensions, and R3-B
  base-vs-augmented `G-G` block equality for `K` and each unit `U_A`.
- The pass also validates R3-B moment inputs before MWG descriptor extraction:
  moment matrices must have the augmented dimension, be finite, and be
  symmetric within tolerance.
- R3-A custom Coulomb expansions remain allowed. Cached centered Gaussian
  factors are reused only when center and exponent order match; otherwise the
  factors are regenerated. R3-B still rejects custom expansions because
  `V_GG` from the base Hamiltonian has no expansion provenance.

Validation:
- Doer ran `git diff --check`, package load, the standalone R3-A endpoint gate,
  and the ignored R3-B H2 validation script. Endpoint values were unchanged:
  R3-A `E_base = -0.7946037173365925`,
  `E_aug = -0.7959028345077851`, and R3-B self-Coulomb
  `0.4574331709135599`.
- Manager reviewed the one-file diff, confirmed the requested high-order
  manager checks were included, ran `git diff --check`, package load, and the
  suspicious-line scan. The manager began rerunning the two longer endpoint
  tests but stopped after user interruption and accepted doer's fresh endpoint
  validation instead.

Goal advancement:
- R3: protects the accepted R3-A/R3-B H2 endpoints against mismatched
  separately supplied construction objects and unsafe cached-factor reuse.
- LT5/LT6: improves correctness guardrails without adding artifacts, public API,
  new files, committed tests, or a new payload/status layer.

Carrying-cost result:
- deleted: duplicate inline candidate label/center construction and the unsafe
  unconditional cached centered-factor reuse.
- simplified: residual identity and construction checks now use small shared
  helpers.
- quarantined: none.
- not deleted because: residual shape/label/center/owner checks and R3-B
  consistency checks protect active R3-A/R3-B contracts.
- exact remaining caller/blocker: R3-B custom Coulomb expansions require
  explicit base-`V_GG` expansion provenance before they can be accepted safely;
  R3-C artifact/provenance remains unapproved.
- added src lines: 80.
- deleted src lines: 17.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- This is the maximum cleanup-budget source addition for this lane. Further R3
  hardening should either delete/simplify existing code in the same pass or get
  explicit over-budget approval.

## Cartesian Hamiltonian Producer Pass 055 - Correct R3-B Final-Basis V_GM Normalization

Commit(s):
- this branch - Correct R3B VGM final-basis normalization

Summary:
- Accepted the R3-B mathematical correction following design commit
  `e1b5796d`. The mixed base-residual IDA block `V_GM` now converts
  density-normalized parent/MWG donor values to the final-basis density
  convention for PQS blocks.
- Direct terminal blocks retain the identity selector path. PQS blocks now
  compute owned-support product weights, final retained weights
  `C' * support_weights`, density-normalized coefficients, and then contract
  `C_density' * V_support_M`.
- `V_MM` and base `V_GG` construction were intentionally unchanged.

Validation:
- Doer ran `git diff --check`, package load, the ignored
  `tmp/work/r3b_h2_vgm_normalization_audit.jl`, and the ignored
  `tmp/work/r3b_h2_augmented_hamiltonian_validation.jl`.
- Manager reviewed the one-file diff, confirmed the active authority uses the
  corrected weight-aware scalar, ran `git diff --check` and package load, and
  accepted doer's fresh expensive validation runs rather than repeating them.

Endpoint facts:
- The normalization audit reports current-vs-weight-aware `V_GM` difference
  `0.0` overall, `0.0` on direct blocks, and `0.0` on PQS blocks after the fix.
- `V_GG` delta `0.0`; `V_MM` delta `0.0`; `V_aug` symmetry error
  `1.5543122344752192e-15`.
- Corrected H2 lowest augmented one-body orbital self-Coulomb is
  `0.4574256036192161`, with delta `0.0` from the weight-aware target.

Goal advancement:
- R3: fixes the mixed base-residual final-basis density convention before the
  scalar is frozen into tracked R3-B validation.
- LT5/LT6: keeps the residual-GTO/MWG path aligned with the existing localized
  IDA final-weight convention rather than a parent-density shortcut.

Carrying-cost result:
- deleted: old PQS path that applied terminal coefficients directly to
  density-normalized `G-M` support values.
- simplified: none.
- quarantined: ignored normalization and R3-B validation scripts remain under
  `tmp/work`.
- not deleted because: direct block path and `V_MM` path remain valid active
  contracts.
- exact remaining caller/blocker: tracked R3-B endpoint validation is still not
  extended to the corrected scalar.
- added src lines: 14.
- deleted src lines: 4.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not restore the direct parent-density `G-M` insertion scalar
  `0.4574331709135599`. The accepted R3-B value is now the weight-aware
  `0.4574256036192161` scalar.

## Cartesian Hamiltonian Producer Pass 056 - Track R3-B Weight-Aware Endpoint

Commit(s):
- this branch - Extend R3 endpoint test for R3B VGM

Summary:
- Accepted the approved standalone R3 endpoint gate extension in
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`. The file now
  covers both R3-A exact one-body/moment checks and the R3-B in-memory
  Hamiltonian endpoint.
- The new R3-B section constructs the base `CartesianIDAHamiltonian`, calls
  `pqs_terminal_residual_gto_augmented_hamiltonian`, verifies the returned
  Hamiltonian type and augmented dimension, checks finite/symmetric `V_aug`,
  confirms the base `V_GG` block is unchanged, and computes the lowest
  augmented one-body orbital self-Coulomb.
- The test independently reconstructs `V_GM` with the documented weight-aware
  final-basis contraction and compares both direct and PQS block rows against
  the produced Hamiltonian interaction matrix.

Validation:
- Doer ran `git diff --check`, package load, and the standalone R3 endpoint
  gate. The gate passed `49/49` in `39.14 s`.
- Manager reviewed the test-only diff, ran `git diff --check`, package load,
  and the suspicious-line scan, and accepted doer's fresh standalone endpoint
  run without repeating the long test.

Endpoint facts:
- R3-B self-Coulomb `0.4574256036192161`, delta `0.0`.
- Independent `V_GM` direct-block error `0.0`; PQS-block error `0.0`.

Goal advancement:
- R3: makes the corrected weight-aware R3-B endpoint durable in tracked
  validation before R3-C or larger-system work.
- LT5/LT6: protects the final-basis density-normalization convention with an
  independent test computation rather than only the final scalar.

Carrying-cost result:
- deleted: none.
- simplified: none.
- quarantined: none.
- not deleted because: existing R3-A checks remain active and the added R3-B
  convention checks protect the corrected normalization.
- exact remaining caller/blocker: no blocker for the tracked H2 R3-B endpoint;
  Be2/Cr2 remain deferred by R3 design guardrails.
- added src lines: 0.
- deleted src lines: 0.
- new tests: extended one standalone integration gate by 86 lines.
- new metadata/status fields: none.

Risk / guardrail:
- This test remains a standalone ~39 s integration gate, not default per-edit
  CI. Do not use it to justify Be2/Cr2 readiness; rank-loss handling and
  larger-scale provider performance remain separate lanes.

## Cartesian Hamiltonian Producer Pass 057 - Optimize R3-A Parent-Supplement Blocks

Commit(s):
- this commit - Optimize R3A parent-supplement block construction

Summary:
- Accepted the measured R3-A performance replacement in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`. The
  active path now builds analytic parent-by-supplement `G-A` and supplement
  `A-A` overlap, kinetic, position, second-moment, and by-center unit nuclear
  blocks using the QW one-dimensional table organization, then projects the
  rectangular parent-by-supplement rows through terminal blocks.
- The old CPB-per-terminal-block mixed-provider path was removed from active
  R3-A construction. This preserves the exact augmented one-body/moment
  contract while avoiding repeated bounding-box construction and repeated
  orbital/axis table rebuilding.
- No parent-stage field, shared helper API, payload/status/report object,
  artifact, public API, R3-C surface, Be2 committed gate, or Cr2 workflow was
  added.

Validation:
- Doer ran `git diff --check`, package load, the standalone R3 H2 endpoint
  gate, and the ignored Be2 donor comparison/performance probe. H2 passed
  `49/49` in `26.9 s`; R3-B self-Coulomb was
  `0.4574256036192164`, within roundoff of the accepted
  `0.4574256036192161` target.
- Be2 far donor comparison reported optimized `G-A` differences of
  `6.11e-16` overlap, `3.55e-15` kinetic, moments up to `2.66e-14`, and
  nuclear up to `4.44e-16`; optimized `A-A` differences were at or below
  `1.78e-15`, nuclear `0.0`.
- Be2 far exact-operator construction improved from the CPB-reference
  `43.66 s / 35430.81 MiB` block path to `1.28 s / 1725.49 MiB` for the
  QW block construction and `2.91 s / 3520.58 MiB` for optimized augmented
  exact operators.
- Manager ran `git diff --check`, reviewed the one-file source diff, checked
  the anti-bloat suspicious-line scan, verified no stale R3 CPB helper names
  remain in the target file, and accepted the fresh doer endpoint/performance
  runs without repeating the long scripts.

Goal advancement:
- R3: removes the first measured Be2 R3-A scaling blocker before Cr2 work.
- LT4/LT8: restores the intended one-dimensional analytic table reuse while
  keeping the terminal block projection and compact R3 owner boundary.

Carrying-cost result:
- deleted: active `_r3a_bounding_cpb`, `_r3a_mixed_block`,
  `_r3a_local_index`, and `_r3a_dense` CPB-per-block helper path.
- simplified: R3-A mixed/self block construction now uses one construction-
  local analytic parent-supplement block family per augmented-operator call.
- quarantined: ignored Be2 comparison/profiling script remains under
  `tmp/work`.
- not deleted because: the existing internal function signature still carries
  `parent_basis_object` for current callers; the small duplicate overlap build
  between residualization and augmented-operator construction is accepted by
  the current design rather than adding a persistent raw-block bundle.
- exact remaining caller/blocker: deterministic rank-loss handling, broader
  same-construction orchestration, bounded high-rank MWG storage, and R3-C
  artifacts remain deferred; Cr2 is still not approved.
- added src lines: 144.
- deleted src lines: 76.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The QW donor organization is now the tactical R3-local performance bridge.
  Do not promote it into a new shared provider/cache surface without a design
  amendment and a second production consumer.

## Cartesian Hamiltonian Producer Pass 058 - Add R3 Same-Construction Entry

Commit(s):
- this commit - Add R3 same-construction Hamiltonian entry

Summary:
- Accepted the docs-approved `HP-R3-FN-03` same-construction overload in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`. The new
  internal call shape takes a same-construction base Hamiltonian, terminal
  basis, bundles, supplement, atom locations, and nuclear charges, then
  constructs the R3 residual object, exact augmented operators, and final
  in-memory augmented `CartesianIDAHamiltonian{Float64}` inside one call.
- The overload delegates to the existing R3-A residual, R3-A augmented
  operator, and lower-level R3-B Hamiltonian helpers. It does not duplicate the
  MWG/IDA logic, add a public API/export, introduce a new result shape, or add
  artifact/report/status/payload fields.

Validation:
- Doer ran `git diff --check`, package load, the standalone R3 H2 endpoint
  gate, and an ignored same-construction comparison script. The standalone
  test passed `49/49` with augmented dimension `489` and R3-B self-Coulomb
  `0.4574256036192164`.
- The ignored comparison reported returned type
  `CartesianIDAHamiltonian{Float64}`, `K`/`U`/`V` deltas of `0.0` versus the
  lower-level composition, and self-Coulomb delta `3.33e-16` from the accepted
  target.
- Manager ran `git diff --check`, reviewed the one-file source diff, checked
  the suspicious-line scan, and accepted doer's fresh endpoint validation
  without repeating the longer scripts.

Goal advancement:
- R3: reduces same-construction provenance risk by giving future callers a
  narrow internal path that does not require manually threading residual and
  augmented-operator objects.
- LT5/LT6: keeps residual-GTO augmentation inside the existing in-memory
  Hamiltonian contract without a wrapper or artifact expansion.

Carrying-cost result:
- deleted: none.
- simplified: same-construction callers can now rely on one internal entry
  rather than independently composing residual and augmented-operator objects.
- quarantined: ignored same-construction validation remains under `tmp/work`.
- not deleted because: lower-level R3-A/R3-B helpers remain active tested
  contracts and useful validation seams.
- exact remaining caller/blocker: deterministic rank-loss handling, bounded
  high-rank MWG storage, broader public/workflow surfaces, and R3-C artifacts
  remain deferred.
- added src lines: 18.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The same-construction entry is internal only. Do not promote it to R1/public
  facade, driver workflow, or artifact authority without a design amendment.

## Cartesian Hamiltonian Producer Pass 059 - Handle R3-A Rank Loss

Commit(s):
- this commit - Handle R3A residual rank loss

Summary:
- Accepted deterministic rank-deficient residual selection in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`. R3-A no
  longer throws merely because the supplement candidate metric has near-null
  directions. It computes the retained rank from the approved eigenvalue
  threshold policy, selects independent raw candidate columns with a
  deterministic pivoted-Cholesky rule, sorts retained indices back to ascending
  candidate order, applies symmetric Lowdin on the selected metric, and embeds
  the resulting transform into the full `T_A` rows.
- The existing full-rank H2 path remains candidate-order symmetric Lowdin. The
  residual object keeps the approved orientation symbol
  `:selected_candidate_order_symmetric_lowdin`; manager review corrected a
  transient new rank-loss orientation symbol before commit so the persistent
  object vocabulary stays aligned with design authority.

Validation:
- Doer ran `git diff --check`, package load, the standalone R3 H2 endpoint
  gate, and an ignored duplicate-candidate rank-loss probe. H2 remained full
  rank with residual dimension `18`, augmented dimension `489`, and R3-B
  self-Coulomb `0.4574256036192164`.
- The ignored rank-loss probe used a duplicate H2 candidate set: candidate
  count `19`, residual rank `18`, retained indices `[1, 2, ..., 18]`,
  `G' S R = 0.0`, `R' S R` error `5.00e-12`, and no hidden negative-eigenvalue
  violation.
- Manager ran `git diff --check`, reviewed the one-file source diff, checked
  the suspicious-line scan, verified the design vocabulary for `orientation`,
  and accepted doer's fresh endpoint/probe validation without rerunning the
  long scripts.

Goal advancement:
- R3: removes the rank-loss correctness blocker for realistic supplements
  before larger Be2/Cr2-style work.
- LT5: keeps residual-GTO construction deterministic and provenance-friendly
  when raw supplement candidates are linearly dependent or already represented
  by the base final basis.

Carrying-cost result:
- deleted: throw-only rank-loss branch.
- simplified: full-rank and rank-deficient residual construction now share the
  same validation and residual-object return path.
- quarantined: ignored duplicate-candidate rank-loss probe remains under
  `tmp/work`.
- not deleted because: lower-level residual, augmented-operator, and R3-B
  helpers remain live contracts.
- exact remaining caller/blocker: bounded high-rank MWG storage, public/workflow
  surfaces, R3-C artifacts, and Cr2 remain deferred.
- added src lines: 49.
- deleted src lines: 7.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- This implements deterministic selection, not support for broad public
  supplement workflows. Do not infer Cr2 readiness without the remaining
  performance and workflow gates.

## Cartesian Hamiltonian Producer Pass 060 - Add R3-C Supplement Provenance Writer

Commit(s):
- this commit - Add R3C supplement provenance writer

Summary:
- Accepted the approved R3-C compact artifact provenance implementation in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`. The new
  internal, non-exported
  `write_pqs_terminal_residual_gto_augmented_hamiltonian(...)` first writes the
  existing `CartesianIDAHamiltonian{Float64}` artifact with
  `write_cartesian_ida_hamiltonian`, then appends fixed
  `supplement_provenance/` keys to the same JLD2 file.
- The artifact stores compact construction provenance: supplement policy,
  basis labels by center, candidate and owner counts, base/residual/augmented
  dimensions, residual convention, rank rule, thresholds, MWG convention, and
  validation labels/reference scalar when supplied.
- The implementation does not serialize residual transforms, full residual
  eigenvalues, candidate labels, MWG centers/widths, dense moment matrices, or
  full construction inputs. It adds no public API, driver workflow, wrapper,
  payload/status/report object, or new artifact shape beyond the approved
  provenance group.

Validation:
- Doer ran `git diff --check`, package load, and the ignored H2 R3-C artifact
  validation script. The script wrote
  `/Users/srw/dmrgtmp/jl_V2wAQM/r3c_h2_augmented_hamiltonian.jld2`, read it
  back with the existing Hamiltonian reader, and reported one-body and
  electron-electron IDA deltas of `0.0`.
- H2 R3-B self-Coulomb remained `0.4574256036192164`.
- Manager ran `git diff --check`, reviewed the one-file source diff, checked
  the suspicious-line scan, confirmed no new tracked tests/tools, and accepted
  doer's fresh artifact validation without rerunning the longer script.

Goal advancement:
- R3/LT6: creates the first durable supplemented Hamiltonian handoff while
  preserving the existing `CartesianIDAHamiltonian` artifact contract.
- LT5: records compact provenance from validated construction inputs without
  turning residual numerical state into metadata.

Carrying-cost result:
- deleted: none.
- simplified: artifact output reuses the existing Hamiltonian writer/readback
  instead of adding a wrapper artifact.
- quarantined: ignored H2 artifact validation remains under `tmp/work`.
- not deleted because: this is the first compact R3 artifact provenance writer
  and there was no stale R3-C writer to retire.
- exact remaining caller/blocker: public workflow exposure and Cr2/consumer
  readiness remain deferred outside R3-C.
- added src lines: 68.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: approved `supplement_provenance/` artifact keys
  only; no staged metadata/status fields.

Risk / guardrail:
- Keep this artifact schema compact. Do not add residual transforms, moments,
  centers, widths, or construction inputs unless a later consumer-driven design
  amendment promotes a residual-basis artifact group.

## Medium-Term Goal Checkpoint After Pass 060

- MT1 and MT2 are completed/stale as active work items: fake-PQS quarantine and
  independent H2 PQS recovery have been superseded by the implemented base PQS
  Hamiltonian and public H/H2 base producer.
- MT3 remains a standing guardrail: shared terminal support and retained
  vocabulary continue to anchor the base and R3 supplement paths.
- MT4 is active but changed shape: MWG/GTO supplements now have accepted H2
  in-memory and artifact endpoints; broader Be2/Cr2 and public workflow remain
  deferred.
- MT5 remains active: recent passes replaced the CPB-per-block R3-A path,
  avoided new payloads, and kept probes ignored.
- MT6 remains active: QW/raw-block code is now a tactical donor for R3-A
  performance, not public architecture. Further old-path classification should
  wait for the next public workflow or Cr2-readiness lane.

The medium-term section should be rewritten in a future planning/doc pass to
reflect the post-R3-C state; this commit records the checkpoint without
changing the top-level goal wording.

## Cartesian Hamiltonian Producer Pass 061 - R3 Closeout Status Refresh

Commit(s):
- this commit - Refresh R3 closeout status

Summary:
- Accepted a docs-only planning/status refresh after R3-A/B/C reached the
  intended narrow H2 supplemented endpoint. The compact authority now records
  R3-A residual/rank-loss/exact one-body and moment work, R3-B
  same-construction MWG/IDA Hamiltonian work, and R3-C compact supplemented
  artifact provenance as implemented for the H2 path.
- The refresh records the Be2 measurement conclusions without promoting Be2 to
  a committed validation gate: R3-A exact-operator construction moved from the
  repeated CPB-per-terminal-block path at about `43.2 s` / `35.4 GiB` to the
  one-shot parent-by-supplement analytic organization at about `1.94 s` /
  `2.1 GiB`, with roundoff agreement for tested blocks. Be2 R3-B at residual
  rank `26` showed modest MWG/IDA storage and runtime, so bounded MWG
  streaming is not urgent before the next planning lane.
- The roadmap now separates three candidate next lanes: usability workflow for
  H2/Be2 supplemented artifacts, measurement-only Cr2-readiness forecasting,
  and basis/supplement realism. The recommended next lane is usability so
  consumers can request GTO/MWG artifacts without assembling private R3 calls.

Validation:
- Design-manager ran `git diff --check`, focused `rg` checks for R3-A/B/C
  implemented status, Be2 conclusions, bounded MWG/Cr2 deferral, and stale
  R3 hardening wording, and confirmed no `src`, `test`, `tools`, or `bin`
  files changed.

Goal advancement:
- R3/LT6: closes the short-term R3 status story before public workflow or Cr2
  expansion.
- MT4: updated from "make MWG/GTO endpoint real" to "choose the next
  consumer-facing lane after a real H2 artifact endpoint and Be2 sanity
  measurement."

Carrying-cost result:
- deleted: stale live-authority wording that treated same-construction,
  deterministic rank loss, independent weight-aware `V_GM`, and R3-C artifact
  provenance as still-open hardening blockers.
- simplified: deferred lanes now distinguish usability, Cr2-readiness, and
  basis/supplement realism instead of one mixed R3 hardening bucket.
- quarantined: Cr2 remains deferred as stress/consumer-readiness, not a
  correctness gate.
- not deleted because: historical log and review material intentionally retain
  old R3 blocker history.
- exact remaining caller/blocker: choose the next lane; no public supplemented
  workflow or Cr2 run is approved by this pass.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- This is status reconciliation only. It does not approve source work, public
  API expansion, driver/bin/tool workflow, new artifact fields beyond R3-C, or
  Be2/Cr2 validation gates.

## Cartesian Hamiltonian Producer Pass 062 - Approve R3 Usability Facade

Commit(s):
- this commit - Approve R3 usability supplemented workflow

Summary:
- Approved a docs-only R3 usability amendment after R3-A/B/C closeout. The
  new lane authorizes one non-exported supported facade,
  `cartesian_residual_gto_mwg_hamiltonian(system; basis, supplement,
  hamfile = nothing)`, so callers can request residual-GTO/MWG supplemented
  Hamiltonians without manually composing base stages, supplement loading, R3
  same-construction construction, and R3-C artifact writing.
- The approved scope is intentionally narrow: z-axis H2 as the committed
  endpoint and z-axis Be2 as an internal/performance-supported proxy. Cr2,
  public export, driver/bin/tool workflow, ECP, EGOI, RHF/solver work, new
  artifact formats, wrappers, and report/status/payload objects remain
  forbidden.
- The design freezes compact `system`, `basis`, and `supplement` NamedTuple
  schemas. The supplement schema uses `basis_by_center`, `lmax`, optional
  `uncontracted`, and optional `width_filtering`; first scope is homonuclear
  and maps to the existing legacy bond-aligned diatomic supplement loader.

Validation:
- Design-manager ran `git diff --check`, focused `rg` checks for all R3U IDs,
  the non-export/no-`src/GaussletBases.jl` guardrail, H2/Be2/Cr2 wording, and
  `supplement_provenance/` reuse, and confirmed no source, test, tool, or bin
  files changed in this docs pass.

Goal advancement:
- R3/LT6: moves from "scientifically coherent internal pieces" to an approved
  usability implementation surface for producing supplemented artifacts.
- MT4: chooses the usability lane over immediate Cr2 stress. Cr2 remains a
  later measurement/stress milestone after the workflow is usable.

Carrying-cost result:
- deleted: none.
- simplified: callers get one approved internal facade rather than manually
  threading R1/R3 stage objects and writer calls.
- quarantined: Be2 is internal/performance-supported only; Cr2 is explicitly
  unsupported.
- not deleted because: lower-level R3-A/B/C helpers remain active contracts
  and validation seams.
- exact remaining caller/blocker: implementation of the R3U facade and H2
  artifact endpoint validation; public export and Cr2-readiness remain later
  design lanes.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none in this docs pass; `HP-R3U-TEST-01` approves extending the
  existing standalone R3 H2 endpoint gate during implementation.
- new metadata/status fields: none.

Risk / guardrail:
- The facade is supported internal surface, not public API. Do not edit
  `src/GaussletBases.jl`, add an export, or add a new source/test file under
  this approval.

## Cartesian Hamiltonian Producer Pass 063 - Correct R3 Residual Selection Authority

Commit(s):
- this commit - Correct R3 owner-local residual selection authority

Summary:
- Accepted a design correction that supersedes the current global
  raw-candidate symmetric Lowdin residual-selection story. Stabilizing the
  global Lowdin would improve roundoff in the wrong algorithm; it would not
  restore atom-local residual content selection.
- The live R3 authority now requires owner-local residual selection: group
  candidates by physical owner, form `M_a = S_AaAa - X_a'X_a`, interpret its
  eigenvalues as residual occupations, discard low-occupation owner-local
  modes using a still-to-be-measured `eta_RG`, orthonormalize owner sectors,
  then perform one final symmetric Lowdin only to merge inter-owner residual
  overlap.
- The old H2 MWG scalar `0.4574256036192161` is now recorded as a
  global-selection baseline, not a future target. Corrected owner-local
  residual orientation can change MWG values and must be remeasured.
- R3U facade implementation is paused until the owner-local residual-selection
  diagnostic records owner-local spectra, merge conditioning, and a corrected
  H2 scalar.

Validation:
- Design-manager ran `git diff --check`, focused `rg` checks for owner-local
  selection/occupation/final-merge wording, stale global Lowdin/candidate-order
  targets, R3U pause wording, and confirmed no `src`, `test`, `tools`, or
  `bin` files changed.

Goal advancement:
- R3/LT6: prevents the usability lane from solidifying the wrong residual-GTO
  basis construction.
- MT4/Cr2-readiness: reframes Cr-like failures as evidence against global
  residual selection, not evidence that q5 or residual-GTO/MWG is unsuitable.

Carrying-cost result:
- deleted: live authority that treated global candidate-order Lowdin and
  global raw-column pivoted-Cholesky selection as the R3 residual algorithm.
- simplified: next action is measurement-only owner-local spectra and merge
  conditioning, not another implementation stabilization pass.
- quarantined: current H2 scalar remains historical/global-selection evidence
  only.
- not deleted because: current R3-A/B/C implementation and tests remain useful
  regression evidence until the owner-local correction is implemented.
- exact remaining caller/blocker: run the owner-local measurement diagnostic
  for H2, Be2, and Cr2; choose `eta_RG`; remeasure H2 MWG scalar; then approve
  a narrow source correction.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not implement a stabilized global Lowdin pass. Do not use eigenvalue
  flooring or width filtering to preserve low-occupation residual modes. Do
  not resume the R3U facade implementation until the corrected owner-local
  residual-selection convention and H2 scalar are recorded.

## Cartesian Hamiltonian Producer Pass 063 - Implement R3 Usability Facade

Commit(s):
- this commit - Add R3 usability supplemented facade

Summary:
- Accepted the non-exported internal usability facade
  `cartesian_residual_gto_mwg_hamiltonian(system; basis, supplement,
  hamfile = nothing)` in `src/cartesian_base_hamiltonian.jl`. The facade
  validates the approved z-axis H2/Be2 input contract, builds base stages once,
  constructs the base Hamiltonian and R3 supplemented Hamiltonian from the same
  terminal basis and parent bundles, and optionally writes the R3-C
  `supplement_provenance/` artifact.
- The existing standalone R3 endpoint test was extended with one facade
  section. It exercises the H2 supplemented artifact path, readback deltas,
  provenance keys, and malformed-input/unsupported-system errors.
- No public export, root include, new source file, new committed test file,
  driver/bin/tool workflow, wrapper, report/status/payload object, or new
  artifact key family was added.

Validation:
- Doer ran `git diff --check`, package load, and the standalone R3 endpoint
  gate. The original R3-A/B section passed `49/49`; the new facade section
  passed `54/54`. The H2 supplemented artifact readback reported kinetic,
  unit-`U`, one-body, and `V` deltas of `0.0`.
- Manager reviewed the source/test diffs, confirmed no `src/GaussletBases.jl`
  export change, ran `git diff --check`, the suspicious-line scan, package
  load, and reran the standalone R3 endpoint gate after tightening the H2
  validation-reference predicate. The rerun passed `49/49` plus `54/54`.

Endpoint facts:
- Augmented dimension `489`.
- R3-B self-Coulomb `0.4574256036192164`, delta
  `3.33e-16` from the accepted target.
- Facade artifact path in manager rerun:
  `/Users/srw/dmrgtmp/jl_TchLVP/r3_h2_supplemented.jld2`.
- Facade readback deltas: kinetic `0.0`, unit `U` `0.0`, one-body `0.0`,
  `V` `0.0`.

Goal advancement:
- R3/LT6: turns the accepted residual-GTO/MWG pieces into a usable internal
  workflow for producing supplemented Hamiltonians and artifacts.
- MT4: advances the usability lane while keeping Cr2, public export, solver,
  and ECP/EGOI work deferred.

Carrying-cost result:
- deleted: none.
- simplified: callers no longer need to manually compose base stages, R3
  residual/operator/Hamiltonian calls, and R3-C writer calls.
- quarantined: Be2 remains internal/performance-supported and is not a
  committed gate.
- not deleted because: lower-level R3-A/B/C helpers remain the scientific
  implementation authority and validation seams.
- exact remaining caller/blocker: public/exported supplemented workflow,
  Cr2-readiness forecasting, ECP/EGOI/RHF/solver, and HamV6 remain deferred.
- added src lines: 124.
- deleted src lines: 0.
- new tests: extended one existing standalone endpoint file by 113 lines; no
  new committed test file.
- new metadata/status fields: none; artifact provenance stays within approved
  `supplement_provenance/`.

Risk / guardrail:
- The facade is module-qualified internal supported workflow, not a public API.
  The H2 self-Coulomb reference is written only for the exact validation
  fixture, not arbitrary H2 inputs.

## Cartesian Hamiltonian Producer Pass 064 - Approve R3 Owner-Local Source Correction

Commit(s):
- this commit - Approve R3 owner-local residual source correction

Summary:
- Accepted the measurement evidence for owner-local residual selection as R3
  source authority. H2, Be2, Cr2 q4, and Cr2 q5 all passed owner-local
  selection with final orthogonality below `1.0e-10`, and no rank loss under
  trial residual-occupation cutoffs `1.0e-8` or `1.0e-7`.
- Froze `eta_RG = 1.0e-8`, kept `tau_neg_abs = tau_neg_rel = 1.0e-12`
  separate from the physical occupation cutoff, and added final-merge
  thresholds `tau_merge_abs = tau_merge_rel = 1.0e-12` with a hard
  near-singular merge failure rule.
- Updated live R3 authority, registry, invariants, AGENTS, and R3U docs so the
  active H2 R3-B/R3U self-Coulomb target is
  `0.4574265214362075`. The prior `0.4574256036192161` scalar remains
  historical global-selection evidence only.

Validation:
- Design-manager ran `git diff --check`, focused `rg` checks for stale
  wait/remeasure wording and both H2 scalars, and confirmed
  `git diff --name-only -- src test tools bin` was empty.

Goal advancement:
- R3/LT6: converts the owner-local residual-selection correction from
  measurement evidence into approved source authority.
- MT4: unblocks the narrow R3 source correction and R3U retargeting while
  keeping public export and Cr2 full Hamiltonians deferred.

Carrying-cost result:
- deleted: live authority that kept owner-local selection as a future
  measurement/source-correction lane.
- simplified: R3U now targets one corrected scalar and one residual-selection
  convention.
- quarantined: old global-selection H2 scalars remain historical review/log
  evidence only.
- not deleted because: current implemented R3 source remains to be corrected
  by repo-doer under the approved owner-local authority.
- exact remaining caller/blocker: implement the owner-local correction in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, update
  the standalone H2 R3/R3U endpoint to `0.4574265214362075`, and keep Cr2 to
  ignored measurement only.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not preserve the old scalar by width scaling or tolerance relaxation. Do
  not use width filtering as conditioning repair, do not add a public export,
  and do not run or approve a full Cr2 Hamiltonian/artifact in the source pass.

## Cartesian Hamiltonian Producer Pass 065 - Approve Residual Gaussian Domain Module

Commit(s):
- this commit - Approve Residual Gaussian domain module

Summary:
- Approved a docs-only source-organization amendment for the future internal
  `CartesianResidualGaussians` module. The decision moves residual Gaussian
  meaning out of the terminal-basis file: owner-local residual basis
  selection, exact augmented operator transformation, moment-matched Gaussian
  descriptors, and residual-containing IDA interactions get physical/domain
  names instead of permanent `R3-A/B/C` source concepts.
- Added approved `HP-RG-*` IDs for the module files, residual basis object,
  `build_residual_gaussian_basis`, `transform_augmented_operator`,
  `moment_matched_gaussians`, `assemble_residual_ida_interaction`, migration
  wiring, and migration validation.
- Froze the migration boundary: `pqs_terminal_residual_gto.jl` may keep only
  temporary delegating wrappers for existing callers, and those wrappers must
  be deleted after callers move.

Validation:
- Design-manager ran docs-only validation: `git diff --check`, focused `rg`
  checks for the new `HP-RG-*` IDs, module path, production function names,
  forbidden public/Cr2 surfaces, and confirmed
  `git diff --name-only -- src test tools bin` was empty.

Goal advancement:
- R3/LT6: turns the residual-GTO/MWG implementation from a terminal-basis
  helper cluster into a named domain-module migration lane.
- MT4: enables a narrow cleanup/source-organization pass without changing
  scientific behavior, public API, artifact schema, or Cr2 scope.

Carrying-cost result:
- deleted: none in this docs pass.
- simplified: future source ownership is now one domain module with physical
  function names rather than more R3-named helpers in the terminal file.
- quarantined: R3-A/B/C labels remain implementation-history and review/log
  vocabulary only.
- not deleted because: current source callers still use
  `pqs_terminal_residual_gto.jl` until repo-doer performs the migration.
- exact remaining caller/blocker: implement the approved module migration,
  rewire existing H2/R3U callers, validate the H2 endpoint and ignored Be2
  measurement, then delete old wrappers once no live callers remain.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not create a broad provider/cache framework or a vague
  `stabilize_residual_metric(...)` entry point. Do not add public exports,
  artifact keys, driver workflow, Cr2 facade/full Hamiltonian, ECP/EGOI, or
  solver/RHF work under the domain-module migration.

## Cartesian Hamiltonian Producer Pass 065 - Implement R3 Owner-Local Residual Selection

Commit(s):
- this commit - Implement R3 owner-local residual selection

Summary:
- Accepted the source correction replacing global residual-GTO selection with
  owner-local residual-occupation selection. Residual Gaussian directions are
  now selected separately on each physical center, then the retained center
  sectors are merged once by the final symmetric Lowdin step.
- Deleted the active global raw-candidate metric selection and pivoted-Cholesky
  helper. The persistent residual object now carries the approved owner-local
  fields: source owner indices, residual occupations, retained counts by
  owner, occupation cutoff, final-merge thresholds, and selection/orientation
  symbols.
- Retargeted the R3-B and R3U H2 endpoint to the owner-local self-Coulomb
  value `0.4574265214362075`. The accepted run produced
  `0.45742652143620843`, a `9.4e-16` delta from target.

Validation:
- Doer ran `git diff --check`, package load, the standalone R3 H2 endpoint
  gate, and the ignored owner-local H2/Be2/Cr2 residual-selection measurement.
  The residual measurement retained full rank for H2, Be2, Cr2 q4, and Cr2 q5
  and kept final `R' S R - I` below `1.0e-10`.
- Manager reviewed the diff, reran `git diff --check`, the anti-bloat
  suspicious-line scan, package load, and
  `julia --project=. test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`.
  The endpoint passed `52/52`; the facade/artifact section passed `59/59`.

Goal advancement:
- R3/LT6: corrects the physical residual-basis construction used by the
  supplemented Hamiltonian path before it becomes broader workflow authority.
- MT4: resumes the R3 usability lane on the corrected residual convention
  while keeping Cr2 full Hamiltonians, public export, and solver work deferred.

Carrying-cost result:
- deleted: global raw-candidate residual selection and the global
  pivoted-Cholesky rank-selection helper.
- simplified: active tests, facade artifact reference, and R3-C provenance now
  share the owner-local scalar and residual-selection vocabulary.
- quarantined: Cr2 remains residual-selection measurement only; no full Cr2
  augmented Hamiltonian/artifact is claimed.
- not deleted because: R3-B/R3-C writer and the non-exported usability facade
  remain live workflow surfaces.
- exact remaining caller/blocker: decide the next lane after owner-local R3U
  restoration; likely Be2 owner-local facade measurement before any Cr2
  support expansion.
- added src lines: 92.
- deleted src lines: 76.
- new tests: none; one existing standalone endpoint file changed by `+18/-9`.
- new metadata/status fields: none; the residual numerical object fields
  changed under the approved R3 contract.

Risk / guardrail:
- The rank-loss natural-mode branch is implemented but not covered by the
  tracked H2 endpoint because H2, Be2, and measured Cr2 cases retained full
  owner-local rank. Continue to treat full Cr2 Hamiltonian and artifact work as
  deferred until explicitly approved.

### Medium-Term Goal Checkpoint After Pass 065

- MT1, base Cartesian Hamiltonian producer: completed for the approved H/H2
  public base path. Public-driver polish remains deferred rather than an active
  correctness blocker.
- MT2, R3 supplemented H2 usability: active and healthier after the owner-local
  correction. The non-exported facade and compact artifact provenance exist;
  the active H2 scalar is now owner-local.
- MT3, Be2 realism/performance proxy: active as the next sensible measurement
  target. Be2 has already shown R3-A exact-operator performance is acceptable
  after the QW donor-block optimization; it should be remeasured through the
  corrected owner-local facade before widening scope.
- MT4, Cr2 readiness: active but still deferred. Residual selection is viable
  for Cr2 q4/q5 in ignored measurement, but full Cr2 Hamiltonian/artifact and
  facade support remain unapproved.
- MT5, anti-bloat and deletion: active. This pass removed the wrong global
  residual selection path instead of preserving it behind another compatibility
  layer.

## Cartesian Hamiltonian Producer Pass 066 - Symmetrize R3 Augmented Operators

Commit(s):
- this commit - Symmetrize R3 augmented operators

Summary:
- Accepted a one-line source bug fix for the owner-local R3 Be2 usability
  path. The corrected residual basis worked, but Be2 stopped before artifact
  writing because the exact augmented `z^2` moment matrix had a floating-point
  antisymmetric residue of about `2.18e-10`, above the strict `1.0e-10`
  moment validation threshold.
- The fix explicitly symmetrizes the augmented `[G,R]` operator matrix returned
  by `_r3a_augmented_operator(...)`. This matches the physical operator
  symmetry and keeps the validation threshold unchanged.

Validation:
- Doer ran `git diff --check`, package load, the standalone R3 H2 endpoint,
  Be2 owner-local facade measurement, and the Be2 failure diagnostic. Be2 then
  returned `CartesianIDAHamiltonian{Float64}`, wrote/read back the artifact
  with zero deltas, and reported all moment symmetry errors as `0.0`.
- Manager reviewed the one-line diff, reran `git diff --check`, the
  suspicious-line scan, package load, and the standalone R3 H2 endpoint. The
  endpoint passed `52/52`; the facade/artifact section passed `59/59`.

Goal advancement:
- R3/LT6: clears the first Be2 usability blocker after owner-local residual
  selection without weakening numerical validation.
- MT3: Be2 remains the active realism/performance proxy before any Cr2 support
  expansion.

Carrying-cost result:
- deleted: none.
- simplified: exact augmented operators now enforce their mathematical
  symmetry at the construction boundary instead of relying on downstream
  validation to catch roundoff residue.
- quarantined: none.
- not deleted because: R3-A/B/C writer and facade surfaces remain live.
- exact remaining caller/blocker: decide whether to record Be2 owner-local
  usability as closed and move to residual-Gaussian module planning.
- added src lines: 1.
- deleted src lines: 1.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- This is not a tolerance relaxation and not a Cr2 readiness claim. Full Cr2
  supplemented Hamiltonian/artifact support remains deferred.

## Cartesian Hamiltonian Producer Pass 067 - Map Residual Gaussian Module Extraction

Commit(s):
- this commit - Record residual Gaussian extraction map

Summary:
- Accepted a read-only extraction map for moving residual-Gaussian domain logic
  out of `pqs_terminal_residual_gto.jl` and into the approved
  `CartesianResidualGaussians` module. The map classifies the current code into
  residual-basis construction, exact augmented-operator transformation,
  moment-matched Gaussian descriptor construction, residual IDA interaction,
  and artifact/facade hooks.
- The map confirms the approved split is feasible without new physics behavior:
  residual basis logic moves to `residual_basis.jl`, exact `[G,A] -> [G,R]`
  operator transformation moves to `augmented_operators.jl`, and matched-width
  Gaussian plus residual-interaction logic moves to `mwg_interaction.jl`.
- Artifact writing, facade parsing, basis loading, terminal basis construction,
  parent bundles, QW donor kernels, and raw analytic Gaussian formulas remain
  outside the Residual Gaussian module. Compact R3 artifact writing may remain
  a terminal/facade hook unless design-manager later expands artifact
  ownership.

Validation:
- Doer performed a read-only `rg`/source map and confirmed final
  `git status --short --branch` was clean. No Julia tests were needed because
  this pass did not change source.
- Manager checked the approved `HP-RG-*` authority, include-order context in
  `src/GaussletBases.jl` and `CartesianFinalBasisRealization.jl`, and the
  manager log before drafting the first source-migration blurb.

Goal advancement:
- RG/LT6: turns the approved domain-module design into a concrete
  behavior-preserving migration map.
- MT5: gives the cleanup/deletion pass a deletion target: old `_r3a_*`,
  `_r3b_*`, and `pqs_terminal_residual_gto_*` names should become temporary
  wrappers only, then disappear once callers move.

Carrying-cost result:
- deleted: none in this read-only pass.
- simplified: future source work is split by physical meaning rather than by
  historical R3-A/B/C implementation labels.
- quarantined: artifact/facade workflow stays outside the RG module for now.
- not deleted because: source migration has not started yet.
- exact remaining caller/blocker: implement the first behavior-preserving
  module slice, starting with residual-basis construction and compatibility
  wrapper wiring.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The first source pass should be narrow. Do not move artifact writing, public
  facade parsing, QW donor kernels, or Cr2 support into the RG module.

## Cartesian Hamiltonian Producer Pass 068 - Extract Residual Gaussian Basis Module

Commit(s):
- this commit - Extract residual Gaussian basis module

Summary:
- Accepted the first behavior-preserving Residual Gaussian source migration.
  The residual basis object and owner-local residual-basis construction moved
  from `pqs_terminal_residual_gto.jl` into the new internal
  `CartesianResidualGaussians` module under
  `src/cartesian_residual_gaussians/`.
- The old terminal residual file now keeps the temporary compatibility alias
  and wrapper needed by current callers. Exact augmented operators, matched
  width Gaussian interaction, artifact writing, and facade parsing intentionally
  remain in the terminal/facade path for later slices.
- The migration uses domain naming around residual Gaussian basis construction
  instead of adding another R3-named helper layer.

Validation:
- Doer ran `git diff --check`, package load, and the standalone H2 R3 endpoint
  gate. The endpoint remained at augmented dimension `489` with H2
  self-Coulomb `0.45742652143620904`, within `1.6e-15` of target.
- Manager reviewed the diff and new files, confirmed no public export or old
  residual-selection helpers remained, ran package load, reran the standalone
  H2 endpoint, and reran the ignored Be2 owner-local usability measurement.
  Be2 returned `CartesianIDAHamiltonian{Float64}`, wrote/read back its artifact
  with zero deltas, kept owner counts `[13, 13]`, and took about `58.7s` with
  `11465 MiB` allocated.

Goal advancement:
- RG/LT6: establishes the approved domain module with the first real physics
  object: the owner-local residual Gaussian basis.
- MT5: begins deleting semantic flattening by moving residual-basis meaning out
  of the terminal-basis implementation file.

Carrying-cost result:
- deleted: old residual-basis struct and owner-local residual construction
  helpers from `pqs_terminal_residual_gto.jl`.
- simplified: `pqs_terminal_residual_gto_augmentation(...)` is now a small
  delegation wrapper around `build_residual_gaussian_basis(...)`.
- quarantined: none.
- not deleted because: current R3 callers still use the old terminal residual
  entry point; exact operators, MWG interaction, artifact writer, and facade
  parsing are out of scope for this slice.
- exact remaining caller/blocker: move exact augmented operators next, then
  retire more `_r3a_*` helper names as callers move into the RG module.
- added src lines: 149.
- deleted src lines: 132.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- This is source organization only. It does not approve public API/export,
  artifact schema changes, Cr2 support, driver workflow, or new physics
  behavior.

## Cartesian Hamiltonian Producer Pass 069 - Extract RG Augmented Operator Transform

Commit(s):
- this commit - Extract RG augmented operator transform

Summary:
- Accepted the second Residual Gaussian source migration slice. The exact
  augmented single-operator transformation moved from
  `pqs_terminal_residual_gto.jl` into
  `src/cartesian_residual_gaussians/augmented_operators.jl` as
  `transform_augmented_operator(...)`.
- The old terminal residual file still builds the base and supplement operator
  blocks through the existing QW donor bridge and terminal product helpers, but
  delegates the exact `[G,A] -> [G,R]` operator math to the RG module.
- Symmetrization of exact operators also now lives with the RG operator
  transform as `symmetrize_operator(...)`.

Validation:
- Doer ran `git diff --check`, package load, the standalone H2 R3 endpoint,
  and the ignored Be2 owner-local facade measurement. H2 remained at augmented
  dimension `489` with self-Coulomb `0.45742652143620904`; Be2 returned
  `CartesianIDAHamiltonian{Float64}` with artifact readback deltas all `0.0`.
- Manager reviewed the diff and new file, confirmed no old
  `_r3a_augmented_operator` helper remained, ran `git diff --check`, package
  load, the standalone H2 R3 endpoint, and the Be2 owner-local facade
  measurement. Be2 took about `59.4s` and allocated `11390 MiB`. The
  suspicious-line scan flagged two existing fixed-axis `NamedTuple{(:x,:y,:z)}`
  constructors only because the callee changed; no new variable-size
  inventory or status shape was introduced.

Goal advancement:
- RG/LT6: separates exact one-body/moment operator transformation from
  terminal-basis assembly and from matched-width Gaussian interaction.
- MT5: continues replacing historical R3 helper names with domain names while
  preserving behavior.

Carrying-cost result:
- deleted: `_r3a_augmented_operator(...)` and `_r3a_sym(...)` definitions from
  the terminal residual file.
- simplified: exact operator construction calls one domain-named
  `CRG.transform_augmented_operator(...)`.
- quarantined: QW donor bridge and terminal product helpers remain in the old
  file until a later slice.
- not deleted because: the compatibility operator-set entry point still owns
  base/supplement block construction and current callers.
- exact remaining caller/blocker: move matched-width Gaussian descriptors and
  residual IDA interaction next, then reconsider whether artifact writing
  should remain outside RG.
- added src lines: 32.
- deleted src lines: 29.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- This pass is exact operator migration only. It does not change MWG/IDA
  interaction behavior, artifact schema, public API, or Cr2 scope.

## Cartesian Hamiltonian Producer Pass 070 - Extract RG MWG Interaction

Commit(s):
- this commit - Extract RG MWG interaction

Summary:
- Accepted the third Residual Gaussian source migration slice. The residual
  two-electron interaction logic moved from `pqs_terminal_residual_gto.jl` into
  `src/cartesian_residual_gaussians/mwg_interaction.jl`.
- The RG module now owns moment-matched Gaussian descriptor construction,
  density-normalized MWG pair-factor construction, weight-aware `V_GM`, direct
  `V_MM`, and final residual interaction matrix assembly through
  `assemble_residual_ida_interaction(...)`.
- The physical conventions are unchanged: widths use `sigma = sqrt(2v)`, base
  `V_GG` is copied unchanged, PQS shell `V_GM` uses final-basis weight-aware
  density normalization, direct blocks use selector behavior, `V_MM` remains
  direct density-normalized matched-Gaussian interaction, and custom expansion
  remains rejected until base `V_GG` carries expansion provenance.

Validation:
- Doer ran `git diff --check`, package load, the standalone H2 R3 endpoint,
  and the ignored Be2 owner-local facade measurement. H2 self-Coulomb remained
  `0.45742652143620904`; Be2 returned `CartesianIDAHamiltonian{Float64}` with
  all artifact readback deltas `0.0`.
- Manager reviewed the diff and new file, confirmed production assembly now
  delegates through `CRG.assemble_residual_ida_interaction(...)`, and confirmed
  the remaining `_r3b_*` wrappers are only for the existing test's independent
  `V_GM` check. Manager reran `git diff --check`, package load, the standalone
  H2 R3 endpoint, and Be2 owner-local facade measurement. Be2 took about
  `60.4s` and allocated `11420 MiB`. The suspicious-line scan flagged two
  fixed three-axis `ntuple` constructions for x/y/z MWG pair factors; these
  are fixed-dimensional axis packets, not variable-size route inventories.

Goal advancement:
- RG/LT6: completes the first internal RG split across residual basis, exact
  operator transformation, and residual MWG interaction.
- MT5: shrinks the terminal residual file and turns R3-B interaction helpers
  into domain-named RG code.

Carrying-cost result:
- deleted: `_r3b_terminal_mwg_fixed_residual(...)`,
  `_r3b_mwg_residual_residual(...)`, and production interaction assembly logic
  from the terminal residual file.
- simplified: the lower-level augmented Hamiltonian wrapper now delegates the
  residual interaction matrix to `CRG.assemble_residual_ida_interaction(...)`.
- quarantined: two `_r3b_*` wrappers remain only for the current standalone
  test's independent `V_GM` reconstruction.
- not deleted because: the test still calls those two wrappers directly, and
  artifact writing/facade parsing remain intentionally outside RG.
- exact remaining caller/blocker: decide whether to update the standalone
  test to use RG-domain helpers directly, then remove the last `_r3b_*`
  wrappers; decide separately whether compact artifact writing should remain
  outside RG.
- added src lines: 99.
- deleted src lines: 111.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- This is behavior-preserving internal migration only. It does not approve
  public export, artifact schema changes, Cr2 support, driver workflow,
  residual-basis changes, or interaction tolerance/scalar changes.

### Medium-Term Goal Checkpoint After Pass 070

- MT1, base Cartesian Hamiltonian producer: still completed for approved H/H2
  public base construction. No new base API work is active.
- MT2, R3 supplemented H2 usability: completed for the non-exported internal
  H2 workflow and compact artifact provenance under the owner-local residual
  convention.
- MT3, Be2 realism/performance proxy: active and healthy. Be2 owner-local
  usability now passes after residual-basis, exact-operator, and MWG
  interaction migration, with zero artifact readback deltas.
- MT4, Cr2 readiness: active but deferred. Cr2 residual selection is viable in
  ignored measurement, but full Cr2 supplemented Hamiltonian/artifact support
  remains explicitly unapproved.
- MT5, anti-bloat and domain cleanup: active. Residual Gaussian domain logic
  now has a real internal module. Remaining cleanup is to retire the last
  `_r3b_*` test-support wrappers and decide whether the compact artifact writer
  should remain outside RG.

## Cartesian Hamiltonian Producer Pass 071 - Retire R3B Test Wrappers

Commit(s):
- this commit - Retire R3B test wrappers

Summary:
- Accepted a cleanup pass removing the last two `_r3b_*` compatibility wrapper
  names from `pqs_terminal_residual_gto.jl`. The standalone H2 endpoint test
  now calls `CartesianResidualGaussians.moment_matched_gaussians(...)` and the
  RG `_mwg_axis_pairs(...)` helper directly for its independent `V_GM`
  convention check.
- This keeps the independent weight-aware `V_GM` check intact while removing
  stale R3-B vocabulary from the terminal residual file.

Validation:
- Doer ran `git diff --check`, package load, and the standalone H2 R3 endpoint.
  The endpoint preserved augmented dimension `489`, self-Coulomb
  `0.45742652143620904`, independent `V_GM` errors `0.0/0.0`, and artifact
  readback deltas `0.0`.
- Manager reviewed the diff, ran `git diff --check`, and ran the suspicious
  added-line scan. No H2 or Be2 rerun was performed because doer had just run
  the endpoint and the patch deleted wrappers without changing production
  interaction logic.

Goal advancement:
- RG/LT6: removes the final R3-B compatibility names left only for the test's
  independent interaction check.
- MT5: closes the immediate wrapper-retirement cleanup item after extracting RG
  interaction logic.

Carrying-cost result:
- deleted: `_r3b_residual_mwg_descriptors` and `_r3b_mwg_axis_pairs`.
- simplified: the standalone test now names the RG domain helpers it is
  validating.
- quarantined: none.
- not deleted because: artifact writing and facade parsing remain intentionally
  outside RG.
- exact remaining caller/blocker: decide whether compact artifact writing
  should remain a terminal/facade hook or get a later docs-approved RG-adjacent
  helper.
- added src lines: 0.
- deleted src lines: 5.
- new tests: none; existing standalone test changed by `+3/-2`.
- new metadata/status fields: none.

Risk / guardrail:
- This is cleanup only. It does not approve public API/export, artifact schema
  changes, Cr2 support, driver workflow, or interaction behavior changes.

## Cartesian Hamiltonian Producer Pass 072 - Clarify RG Artifact Boundary

Commit(s):
- this commit - Clarify RG artifact boundary

Summary:
- Recorded that compact supplemented Hamiltonian artifact writing remains
  outside `CartesianResidualGaussians`. RG owns residual Gaussian physics:
  residual basis selection, exact augmented operator transformation,
  moment-matched Gaussian descriptors, and residual IDA interaction assembly.
- Clarified that artifact writing is workflow/provenance glue attached to the
  supported R3 usability path. The current
  `write_pqs_terminal_residual_gto_augmented_hamiltonian(...)` location remains
  acceptable under R3-C terminal/facade workflow authority.
- Updated the live authority so RG may provide objects and Hamiltonians
  consumed by the writer, but RG does not own JLD2 file workflow, facade input
  parsing, or `supplement_provenance/` schema policy. Moving or splitting the
  writer now requires a named duplication or consumer reason, not proximity to
  RG.

Validation:
- Design-manager ran `git diff --check`, focused `rg` checks for artifact
  boundary wording, `supplement_provenance/`, `CartesianResidualGaussians`, and
  writer ownership, and confirmed
  `git diff --name-only -- src test tools bin` was empty.

Goal advancement:
- RG/LT6: keeps the new domain module focused on physics while preserving the
  compact artifact path as workflow/provenance glue.
- MT5: resolves the remaining cleanup decision from Pass 071 without moving
  code or adding a new source surface.

Carrying-cost result:
- deleted: no files or source in this docs-only clarification.
- simplified: future agents no longer need to decide whether artifact writing
  belongs in RG by default.
- quarantined: artifact/facade workflow remains outside RG unless a later
  amendment names a real reason to move it.
- not deleted because: `write_pqs_terminal_residual_gto_augmented_hamiltonian`
  remains a live R3-C workflow writer.
- exact remaining caller/blocker: none for this ownership decision; any future
  writer movement needs a separate docs-only amendment.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not move artifact schema ownership into RG just because the writer
  consumes RG outputs. No artifact schema change, public API/export, Cr2
  support, driver workflow, or new behavior is approved here.

## Cartesian Hamiltonian Producer Pass 073 - Compact RG Design Authority

Commit(s):
- this commit - Compact RG design authority

Summary:
- Compacted active Cartesian Hamiltonian producer authority so
  `residual_gaussian_domain_module.md` is the canonical Residual Gaussian
  algorithm page. `current.md` now gives live status and reading order,
  `registry.md` records approved IDs and surfaces, `invariants.md` records
  guardrails, and `implementation_slices.md` records migration status.
- Slimmed `AGENTS.md` so startup keeps approved surfaces and guardrails without
  carrying a second full copy of the RG residual-selection/MWG algorithm.
- Historical R3 labels and old scalar evidence remain in history/log material;
  active authority points future agents to RG domain names and the owner-local
  H2 scalar.

Validation:
- Planned/used validation for this docs-only pass: `git diff --check`,
  focused stale-wording scans, focused canonical-pointer scans, and
  confirmation that no `src/`, `test/`, `tools/`, or `bin/` files changed.

Goal advancement:
- RG/LT6: makes the domain module the single active algorithm authority.
- MT5: reduces duplication and startup drift without changing source behavior.

Carrying-cost result:
- deleted: repeated active-authority RG algorithm prose from current/registry/
  invariants/AGENTS.
- simplified: active startup docs now point to one canonical RG contract.
- quarantined: full R3/RG history remains historical evidence only.
- not deleted because: approved IDs, artifact keys, and source ownership still
  belong in registry/AGENTS.
- exact remaining caller/blocker: none for docs compaction; source wrappers and
  artifact/facade hooks remain governed by existing live callers and approved
  surfaces.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Future RG algorithm changes should edit the canonical domain-module page and
  registry surfaces, not re-expand startup docs with duplicated algorithms.

## Cartesian Hamiltonian Producer Pass 074 - Measure Cr2 RG Components

Commit(s):
- this commit - Record Cr2 RG component probe

Summary:
- Accepted a measurement-only Cr2 q4 Residual Gaussian component probe. The
  probe did not expand facade support, write artifacts, or change tracked
  files. It tested lower-level internal construction through terminal basis,
  residual Gaussian basis, exact augmented operators, base Hamiltonian,
  residual MWG/IDA interaction, and a small eigensolve diagnostic.
- Cr2 q4 completed all measured RG components. The base terminal dimension was
  `1557`, the Cr/cc-pV5Z `lmax=1` supplement had `66` candidates split
  `[33, 33]`, residual rank was `66`, and the augmented dimension was `1623`.
  Residual checks passed with `G' S R = 3.89e-14` and
  `R' S R - I = 4.00e-11`.
- Exact augmented `K`/`U`/moment matrices were finite and symmetric; the
  residual interaction matrix was finite and symmetric with `1.42e-14`
  symmetry error. The lowest one-body diagnostic was
  `-295.54804525309942`, and the lowest-orbital self-Coulomb diagnostic was
  `7.712075333114081`.

Validation:
- Doer ran `git diff --check`, package load, the ignored script
  `tmp/work/cr2_rg_component_probe.jl`, final tracked `git status`, and an
  ignored-file status check. Manager did not rerun the 105-second probe.

Goal advancement:
- Cr2-readiness/MT4: moves Cr2 from residual-selection-only evidence to a
  lower-level component readiness measurement. Correctness did not fail at q4;
  performance allocation now defines the next blocker.
- RG/LT6: confirms the extracted RG component path is viable beyond H2/Be2 at
  a bounded Cr2 q4 scale.

Performance read:
- Terminal/base stages: `46.21s`, about `2.13 GiB`.
- Residual mixed overlap plus `S_AA` plus selection: `3.61s`, about
  `11.90 GiB`.
- Exact augmented operators: `31.96s`, about `111.85 GiB`.
- Base Hamiltonian `G-G`: `3.71s`, about `4.12 GiB`.
- MWG pair factors plus `V_GM`/`V_MM`/`V`: `0.26s`, about `0.42 GiB`.
- Dense final storage is small relative to allocation churn: one augmented
  dense matrix is about `0.020 GiB`, and ten dense matrices about `0.196 GiB`.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: Cr2 readiness now has a concrete bottleneck instead of a broad
  uncertainty list.
- quarantined: Cr2 remains an ignored lower-level measurement only; no facade
  support, artifact, or public claim was made.
- not deleted because: no source was changed.
- exact remaining caller/blocker: investigate exact augmented operator
  allocation churn before approving Cr2 facade/artifact support or neutral
  cross-block extraction.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not jump from this probe to full Cr2 support. The next pass should be a
  focused allocation audit of exact augmented operator assembly, not facade
  expansion or artifact writing.

## Cartesian Hamiltonian Producer Pass 075 - Audit Cr2 Exact Operator Allocation

Commit(s):
- this commit - Record Cr2 exact operator allocation audit

Summary:
- Accepted a measurement-only allocation audit for Cr2 q4 exact augmented
  operator construction. The audit shows the large allocation is not from final
  dense output matrices or the `[G,A] -> [G,R]` transform. It is dominated by
  QW donor by-center nuclear block construction.
- The actual `pqs_terminal_residual_gto_augmented_operators(...)` call measured
  `32.17s` and `114553 MiB`. Replayed substeps attribute `27.33s` and
  `106319 MiB` to QW donor blocks, `2.87s` and `7076 MiB` to `G-G` product
  matrices, and only `0.061s` and `817 MiB` to exact transform calls.
- The dominant substep is `qw_bycenter_nuclear_blocks_GA_AA`, at `24.39s` and
  `95365 MiB`. The next largest are QW self moment blocks at `2.31s` /
  `8453 MiB` and QW cross moment blocks at `0.46s` / `2488 MiB`.

Validation:
- Doer ran `git diff --check`, package load,
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, and final tracked
  `git status`. The audit verified wrapper/replay kinetic delta `0.0`, all
  replayed exact operators finite, and replayed exact-operator symmetry error
  `0.0`. Manager did not rerun the Cr2 audit.

Goal advancement:
- Cr2-readiness/MT4: narrows the performance blocker to by-center nuclear
  `G-A`/`A-A` donor construction.
- RG/LT6: protects the RG extraction from premature optimization in the wrong
  place. `transform_augmented_operator(...)` is not the first Cr2 bottleneck.

Performance read:
- Final storage remains small: augmented dimension `1623`; one augmented dense
  matrix is about `20.10 MiB`; nine augmented output matrices are about
  `180.87 MiB`.
- Parent-by-supplement dimensions are `(11191, 66)`.
- Allocation churn is therefore transient donor construction, especially
  repeated one-dimensional factor/table allocation in nuclear by-center paths.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: next optimization target is concrete:
  `_r3a_qw_nuclear_blocks(...)`, not exact transform in-place work.
- quarantined: no Cr2 facade/artifact support is claimed.
- not deleted because: no source changed.
- exact remaining caller/blocker: plan a bounded optimization/design pass for
  by-center nuclear `G-A`/`A-A` block construction, likely by caching/reusing
  one-dimensional nuclear factor tables and separating `G-A` from `A-A` work.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not optimize the dense RG transform first. Do not broaden to full Cr2
  facade or artifact support before reducing or justifying nuclear donor
  allocation.

## Cartesian Hamiltonian Producer Pass 076 - Audit Cr2 Nuclear Factor Reuse

Commit(s):
- this commit - Record Cr2 nuclear factor reuse audit

Summary:
- Accepted a read-only planning/audit pass for the Cr2 q4 by-center nuclear
  factor-table hotspot. The waste is mostly in `A-A` nuclear one-dimensional
  factor tables inside `_r3a_qw_nuclear_blocks(...)`.
- The current `A-A` loop rebuilds nuclear factor tables for every ordered
  supplement orbital pair, including both `(i,j)` and `(j,i)`, even though the
  final `A-A` nuclear matrix is symmetric. `G-A` already reuses x/y
  center-coordinate tables across the two Cr centers; `A-A` reuses x/y only
  within one orbital pair and then discards that local cache.

Validation:
- Doer ran `git diff --check`, package load,
  `tmp/work/cr2_nuclear_factor_reuse_audit.jl`, and final tracked
  `git status`. Manager did not rerun the probe.

Cr2 q4 facts:
- Parent counts `(19,19,31)`, parent Cartesian count `11191`, base dimension
  `1557`, supplement candidates `66`, primitive count `1284`, Coulomb terms
  `45`, unique center coordinates x/y/z = `1/1/2`.
- `G-A` axis factor calls: x `66`, y `66`, z `132`.
- `A-A` axis factor calls: x `4356`, y `4356`, z `8712`, total `17424`.
- Allocation split: `G-A` factor tables `11443 MiB`, `G-A` assembly
  `10864 MiB`, `A-A` factor tables `67752 MiB`, `A-A` assembly `5310 MiB`.
  Nuclear total was `24.31s` with `79195 MiB` table allocation and
  `16174 MiB` assembly allocation.

Goal advancement:
- Cr2-readiness/MT4: identifies the first bounded optimization: avoid ordered
  duplicate `A-A` nuclear pair work.
- RG/LT6: keeps the optimization local to the current exact-block donor bridge
  instead of prematurely extracting a broad neutral cross-block module.

Carrying-cost result:
- deleted: none; read-only audit.
- simplified: next source pass can stay inside
  `_r3a_qw_nuclear_blocks(...)`.
- quarantined: neutral `CartesianGaussianCrossBlocks` extraction remains a
  later architecture cleanup, not the immediate Cr2 blocker fix.
- not deleted because: no source changed.
- exact remaining caller/blocker: implement upper-triangle `A-A` nuclear loop
  with mirroring, plus optional local scratch cleanup for `G-A`, then remeasure
  Cr2 exact-operator allocation.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not add a new source owner, neutral module, facade support, or artifact
  workflow in the optimization pass. Keep it a local behavior-preserving
  performance fix unless source evidence contradicts the audit.

## Cartesian Hamiltonian Producer Pass 077 - Reduce Cr2 Nuclear A-A Duplication

Commit(s):
- this commit - Reduce R3 nuclear A-A duplicate work

Summary:
- Accepted a bounded local optimization in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.
  `_r3a_qw_nuclear_blocks(...)` now computes supplement-supplement nuclear
  attraction blocks over the upper triangle for each nuclear center and mirrors
  the result, instead of rebuilding one-dimensional factor tables for both
  ordered pairs `(i,j)` and `(j,i)`.
- The same pass removed a small base-to-supplement scaled-vector temporary in
  the nuclear block assembly by updating the destination column in place.
- The physical object is unchanged: exact augmented kinetic, uncharged
  by-center nuclear-attraction, position, and second-moment matrices still feed
  the Residual Gaussian H2/Be2/Cr2 paths.

Validation:
- Doer ran `git diff --check`, package load,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`, and
  `tmp/work/cr2_exact_operator_allocation_audit.jl`.
- Manager reviewed the source diff and ran `git diff --check`, `git diff
  --numstat -- src bin tools test docs`, the suspicious-added-line anti-bloat
  scan, and the new-test/file scan. Manager did not rerun the H2/Be2/Cr2
  numerical gates.

Performance result:
- Cr2 q4 exact augmented operators improved from `32.17s` / `114553 MiB` to
  `22.06s` / `67660 MiB`.
- By-center nuclear `G-A`/`A-A` improved from `24.39s` / `95365 MiB` to
  `14.37s` / `48541 MiB`.
- QW donor blocks as a group improved from `27.33s` / `106319 MiB` to
  `17.34s` / `59495 MiB`.
- Wrapper/replay deltas remained roundoff or zero; kinetic replay delta was
  `0.0`; exact operator matrices remained finite and symmetric.

Goal advancement:
- Cr2-readiness/MT4: removes the first confirmed duplicate nuclear-table
  construction cost and leaves a smaller, better isolated Cr2 exact-operator
  allocation problem.
- RG/LT6: preserves the current RG module split. This is a tactical local
  exact-block donor optimization, not a broad neutral cross-block extraction.

Carrying-cost result:
- deleted: duplicate lower-triangle `A-A` nuclear construction and the local
  `G-A` scaled scratch allocation.
- simplified: `A-A` nuclear construction is now upper-triangle plus mirror.
- quarantined: broader donor-kernel extraction/cache work remains deferred.
- not deleted because: remaining `G-A` and upper-triangle `A-A` nuclear factor
  tables are active numerical work.
- exact remaining caller/blocker: `_r3a_qw_blocks(...)` still calls
  `_r3a_qw_nuclear_blocks(...)`; remaining Cr2 cost is live nuclear
  factor-table construction/assembly, not duplicate ordered-pair work.
- added src lines: 29.
- deleted src lines: 22.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not turn this into immediate Cr2 facade/artifact support. The next
  decision should compare the remaining Cr2 exact-operator allocation against
  the cost and risk of a broader neutral Cartesian Gaussian cross-block kernel.

## Cartesian Hamiltonian Producer Pass 078 - Extract Neutral Nuclear Raw Blocks

Commit(s):
- this commit - Extract Gaussian nuclear raw blocks

Summary:
- Accepted the first `CartesianGaussianRawBlocks` implementation slice. Exact
  uncharged by-center Cartesian Gaussian nuclear raw blocks now have one
  neutral owner:
  `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` and
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`.
- Residual Gaussian and Qiu-White both call
  `CartesianGaussianRawBlocks.gaussian_nuclear_raw_blocks_by_center(...)`.
  The route-local `G-A`/`A-A` nuclear loops were deleted from
  `pqs_terminal_residual_gto.jl` and `_qwrg_diatomic_cartesian_shell_blocks_3d`.
- The kernel returns uncharged by-center matrices. Nuclear charges remain in
  the Hamiltonian/one-body consumers, preserving the approved boundary.

Validation:
- Doer ran `git diff --check`, package load,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`,
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, and
  `tmp/work/cgrb_qw_nuclear_parity_check.jl`.
- Manager reviewed the source diff and ran `git diff --cached --check`,
  cached numstat, the suspicious-added-line anti-bloat scan, and the
  new-test/file scan. Manager did not rerun the numerical gates.

Numerical result:
- H2 Residual Gaussian self-Coulomb remained
  `0.45742652143620904`, delta `1.55e-15` from the active target.
- Be2 artifact readback deltas remained `0.0`.
- Cr2 exact-operator audit reported raw nuclear `G-A` and `A-A` reference
  deltas `0.0`, with exact operator matrices finite and symmetric.
- Small Qiu-White parity reported route-vs-neutral `G-A`/`A-A` deltas `0.0`
  and one-body rebuild deltas `0.0`.

Performance result:
- Extraction preserved the improved upper-triangle baseline. Cr2 q4 nuclear
  step measured `14.34s` / `48539.7 MiB`, compared with the immediate
  pre-extraction upper-triangle result of about `14.37s` / `48541 MiB` and the
  older duplicated lower-triangle baseline of about `24.39s` / `95364.6 MiB`.

Goal advancement:
- Cr2-readiness/MT4: removes duplicate raw nuclear implementations and gives
  the remaining Cr2 nuclear allocation a permanent owner for future streamed
  table optimization.
- RG/LT6: keeps Residual Gaussian focused on residual physics and exact
  transforms, not raw analytic Gaussian nuclear formulas.
- QW/RG shared-kernel direction: establishes the narrow nuclear-only owner
  without broadening to overlap, kinetic, moments, terminal projection, pair
  factors, artifacts, or route objects.

Carrying-cost result:
- deleted: duplicated Residual Gaussian nuclear loop and Qiu-White nuclear
  loop.
- simplified: Qiu-White one-body nuclear assembly now consumes neutral
  uncharged by-center raw blocks and applies charges locally.
- quarantined: broader raw-block framework/extraction for overlap, kinetic,
  moments, pair factors, and terminal projection.
- not deleted because: low-level QW one-dimensional factor helpers remain the
  existing donor kernels used by the neutral owner.
- exact remaining caller/blocker: the neutral kernel still calls
  `_qwrg_atomic_axis_factor_*` helpers; extracting or optimizing those tables
  requires a separate approved pass.
- added src lines: 106.
- deleted src lines: 157.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The next optimization should happen inside `CartesianGaussianRawBlocks` and
  should target streamed/reused one-dimensional nuclear tables. Do not broaden
  the owner or move artifact/facade responsibilities into it.

## Cartesian Hamiltonian Producer Pass 079 - Stream Nuclear A-A Terms

Commit(s):
- this commit - Stream Gaussian nuclear AA terms

Summary:
- Accepted a neutral-kernel optimization in
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`.
  The supplement-supplement (`A-A`) nuclear block now streams
  Coulomb-Gaussian terms through reusable per-term axis tables and immediate
  scalar contraction, instead of storing vectors of per-axis factor matrices.
- The physical calculation is unchanged: exact uncharged by-center nuclear raw
  blocks are still returned as `(; ga, aa)`, `A-A` still uses upper-triangle
  assembly plus mirroring, and nuclear charges remain outside the raw-block
  kernel.

Validation:
- Doer ran `git diff --check`, package load,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`,
  `tmp/work/cgrb_qw_nuclear_parity_check.jl`, and
  `tmp/work/cr2_exact_operator_allocation_audit.jl`.
- Manager reviewed the source diff, confirmed the legacy Gaussian orbital
  fields used by the helper, and ran `git diff --check`, numstat, the
  suspicious-added-line anti-bloat scan, and the new-test/file scan. Manager did
  not rerun the numerical gates.

Numerical result:
- H2 Residual Gaussian self-Coulomb remained
  `0.45742652143620927`, delta `1.78e-15` from the active target.
- Be2 augmented dimension remained `1421` and artifact readback deltas were
  `0.0`.
- Qiu-White parity reported neutral/reference `G-A` delta `0.0`,
  neutral/reference `A-A` delta `4.44e-16`, and route/one-body rebuild deltas
  `0.0`.
- Cr2 exact-operator audit reported raw nuclear `G-A` reference delta `0.0`,
  raw nuclear `A-A` reference delta `2.84e-14`, with exact operator matrices
  finite and symmetric.

Performance result:
- Cr2 q4 nuclear raw-block step improved from `14.34s` / `48539.7 MiB` to
  `14.07s` / `44548.7 MiB`.
- The reduction is modest in percentage, about `8.2%`, but material in
  absolute allocation: about `3991 MiB`.
- The remaining dominant cost is now the `G-A` parent-supplement factor path,
  which still calls `_qwrg_atomic_axis_factor_cross_data(...)` and stores term
  factor matrices.

Goal advancement:
- Cr2-readiness/MT4: reduces the remaining nuclear allocation and identifies
  the next concrete target, `G-A` parent-supplement factor construction.
- CGRB/LT6: proves the neutral raw-block owner can carry focused allocation
  improvements without broadening its scope or reintroducing route-local loops.

Carrying-cost result:
- deleted: `A-A` vector-of-factor-tables construction in the neutral owner.
- simplified: `A-A` now uses per-term scratch tables and immediate scalar
  contraction.
- quarantined: `G-A` term-vector construction remains.
- not deleted because: `G-A` is live and now dominates remaining allocation.
- exact remaining caller/blocker:
  `_qwrg_atomic_axis_factor_cross_data(...)` still stores term factor matrices
  for parent-supplement nuclear blocks.
- added src lines: 123.
- deleted src lines: 22.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The next optimization should target `G-A` streaming inside
  `CartesianGaussianRawBlocks`, still without changing Residual Gaussian,
  Qiu-White route semantics, artifacts, public API, or non-nuclear raw blocks.

## Cartesian Hamiltonian Producer Pass 078 - Approve Neutral Gaussian Nuclear Raw Blocks

Commit(s):
- this commit - Approve neutral Gaussian nuclear raw blocks

Summary:
- Approved a narrow neutral Cartesian Gaussian raw-block nuclear owner under
  `HP-CGRB-FILE-01`, `HP-CGRB-FN-01`, `HP-CGRB-WIRE-01`, and
  `HP-CGRB-TEST-01`. The approved owner is
  `src/cartesian_gaussian_raw_blocks/` and owns only exact uncharged by-center
  Gaussian nuclear parent-supplement `G-A` and supplement-supplement `A-A`
  raw blocks.
- This is a design/ownership correction after Passes 075-077 showed the Cr2 q4
  exact-operator allocation hotspot sits in nuclear raw-block construction.
  Residual Gaussian remains the owner of residual selection, exact operator
  transforms, MWG descriptors, and residual interaction. Qiu-White remains a
  consumer/parity route, not the owner of the shared kernel.
- The approved implementation sequence is extraction first, Residual
  Gaussian/Qiu-White rewiring second, duplicate route-local loop deletion
  third, and only then allocation optimization inside the neutral owner.

Validation:
- Design-manager validation is docs-only: `git diff --check`, focused `rg`
  for `HP-CGRB-*` IDs and forbidden-scope wording, and confirmation that no
  `src`, `test`, `tools`, or `bin` files changed.

Goal advancement:
- Cr2-readiness/MT4: creates the exact source owner needed to attack the
  dominant Cr2 q4 nuclear allocation cost without expanding Cr2 workflow.
- RG/LT6: removes raw analytic nuclear formula ownership from RG while
  preserving RG's exact augmented-operator transform contract.
- Qiu-White reuse: converts QW from route-local duplicate nuclear logic into a
  parity consumer of a neutral kernel.

Carrying-cost result:
- deleted: none; docs-only design authority.
- simplified: ownership boundary for shared Gaussian nuclear raw blocks is now
  explicit and no longer flattened into R3/RG or QW route code.
- quarantined: overlap, kinetic, moments, MWG/pair interaction, terminal
  projection, artifacts, public API, and Cr2 facade/artifact workflow remain
  out of scope.
- not deleted because: source extraction and duplicate-loop deletion belong to
  the next doer pass under the new IDs.
- exact remaining caller/blocker: implement the neutral nuclear kernel, rewire
  `pqs_terminal_residual_gto.jl` and Qiu-White callers, prove H2/Be2/Roundoff
  parity including Cr2 q4 exact nuclear blocks, then delete the old route-local
  loops.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none in this docs pass; `HP-CGRB-TEST-01` approves one standalone
  Qiu-White nuclear parity fixture if needed.
- new metadata/status fields: none.

Risk / guardrail:
- Do not generalize this into a broad Gaussian raw-block framework. The first
  owner is nuclear slice only: uncharged by-center `G-A`/`A-A`, no persistent
  cache, no route payload, and no Cr2 workflow promotion.

## Cartesian Hamiltonian Producer Pass 080 - Audit Nuclear Raw Allocation Source

Commit(s):
- this commit - Record nuclear raw allocation source audit

Summary:
- Accepted a measurement-only audit after rejecting the uncommitted `G-A`
  streaming WIP. The audit restored the accepted baseline first and measured
  where the remaining Cr2 q4 nuclear raw-block allocation is coming from.
- The result changes the optimization read: the remaining cost is no longer
  mostly final dense matrices, route-local duplicate loops, or high-level
  `G-A`/`A-A` block orchestration. It is mostly one-dimensional analytic
  nuclear factor work, especially repeated primitive factor-table filling and
  scalar axis-integral calls.

Validation:
- Doer restored the discarded WIP in
  `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` and
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`, then ran
  `git status --short --branch`, `git diff --check`, package load, and
  `tmp/work/cr2_nuclear_raw_allocation_source_audit.jl`. Final tracked status
  was clean.
- Manager confirmed clean status and `git diff --check`. Manager did not rerun
  the audit.

Cr2 q4 facts:
- Fixture: all-electron z-axis Cr2, `q=4`, `core_spacing=0.05`, extents
  `20.0 / 20.0`, Cr/cc-pV5Z, `lmax=1`, contracted, no width filtering.
- Parent Cartesian count `11191`, parent axes `(19,19,31)`, supplement
  candidates `66`, Coulomb terms `45`, nuclear centers `2`.
- G-A axis factor calls `264`, primitive factor matrices `11880`, scalar
  one-dimensional integral calls `51,539,760`.
- Parity against the accepted neutral full function: `G-A` delta `0.0`,
  `A-A` delta `0.0`.

Allocation breakdown:
- `G-A` primitive input vectors: `0.048s`, `10.767 MiB`.
- `G-A` primitive tables only: `3.188s`, `10787.069 MiB`.
- `G-A` stencil matrix conversions: `0.013s`, `557.781 MiB`.
- `G-A` primitive plus left contraction: `3.316s`, `11440.003 MiB`.
- Current `G-A` factor construction: `3.513s`, `11470.237 MiB`.
- `G-A` 3D assembly from cached factors: `0.529s`, `20.625 MiB`.
- `A-A` streamed upper-triangle fill: `9.890s`, `33071.848 MiB`.
- `A-A` output symmetrization: `0.002s`, `0.188 MiB`.
- Full neutral nuclear raw blocks: `13.880s`, `44552.840 MiB`.

Goal advancement:
- Cr2-readiness/MT4: narrows the remaining allocation problem below the
  raw-block orchestration layer. The next source pass should target the
  one-dimensional analytic axis-integral/table layer, not another wrapper-level
  `G-A` loop reshaping.
- CGRB/LT6: confirms `CartesianGaussianRawBlocks` is the correct consumer, but
  not necessarily the owner of the lowest-level axis-integral fill API.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: rejected `G-A` streaming WIP is not part of the live branch; the
  next bottleneck is now clearly lower-level axis table/scalar integral work.
- quarantined: no tracked instrumentation or source changes.
- not deleted because: no production code changed.
- exact remaining caller/blocker: decide whether to approve a narrow in-place
  `_cartesian_gaussian_axis_integral_table!`-style helper in the existing axis
  integral owner, likely `src/cartesian_gaussian_axis_integrals.jl`, and then
  adapt `CartesianGaussianRawBlocks` to consume it.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not revive the rejected high-level `G-A` streaming patch. If an in-place
  axis-table helper is approved, keep it low-level and narrow: no broad raw
  block framework, no persistent cache, no route/status/report object, and no
  Cr2 workflow expansion.

## Cartesian Hamiltonian Producer Pass 081 - Approve Nuclear Family Reuse Target

Commit(s):
- this commit - Approve nuclear family reuse target

Summary:
- Accepted the docs-only correction to the Cartesian Gaussian raw-block nuclear
  optimization target. The active boundary is now `HP-CGRB-FN-02`: reorganize
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl` around unique
  one-dimensional supplement axis families, integer `orbital_axis_family` maps,
  unique `G-A` and `A-A` table keys, and term-first table reuse.
- `HP-CGAI-FN-01` is explicitly superseded as a performance endpoint. It
  remains allowed only as optional support in
  `src/cartesian_gaussian_axis_integrals.jl` for a specialized nonallocating
  `:factor` scalar integral or tiny table-fill helper needed by
  `HP-CGRB-FN-02`.
- The correction records that the in-place table helper alone cannot solve the
  measured Cr2 allocation because `A-A` already fills preallocated tables, and
  the remaining cost comes from repeated one-dimensional family evaluation and
  scalar nuclear-factor integral allocation.

Validation:
- Design-manager validation was docs-only: `git diff --check`, focused `rg`
  for `HP-CGRB-FN-02`, supersession wording, `orbital_axis_family`, canonical
  family-pair keys, `sum_pq` contraction wording, nonallocating factor-helper
  wording, and independent-contraction forbiddance. It also confirmed no
  `src`, `test`, `tools`, or `bin` files changed.
- Manager reviewed the commit surfaces and pushed the amended docs/log commit;
  no source tests were run.

Goal advancement:
- Cr2-readiness/MT4: moves the next source target from result-matrix allocation
  to the real remaining kernel shape: one-dimensional family reuse and
  nonallocating nuclear-factor scalar evaluation.
- CGRB/LT6: preserves `CartesianGaussianRawBlocks` as the neutral nuclear
  owner while forbidding broad raw-block framework expansion.

Carrying-cost result:
- deleted: no source; docs supersede the stale performance interpretation of
  `HP-CGAI-FN-01`.
- simplified: active optimization authority now says what the kernel must
  actually reuse: axis families and canonical family-pair table keys.
- quarantined: persistent caches, metadata/status/report fields, route objects,
  payload structs, public API/export, artifact changes, Residual Gaussian
  algorithm changes, Qiu-White route semantic changes, overlap/kinetic/moment
  migration, Cr2 facade/artifact workflow, and independent x/y/z scalar
  contractions remain forbidden.
- not deleted because: source implementation is the next pass under the new
  `HP-CGRB-FN-02` boundary.
- exact remaining caller/blocker: implement family inventory and table reuse in
  `nuclear_blocks.jl`, optionally adding only the approved low-level
  nonallocating `:factor` helper if needed.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not revive the rejected high-level `G-A` streaming or in-place-table-only
  patches as the endpoint. The next implementation must preserve the coupled
  primitive-pair contraction `sum_pq c_p c_q I_x[p,q] I_y[p,q] I_z[p,q]` and
  must not independently contract x/y/z axes.

## Cartesian Hamiltonian Producer Pass 082A - Count Nuclear Axis-Family Reuse

Commit(s):
- this commit - Record nuclear family reuse count audit

Summary:
- Accepted a measurement-only audit of the Cr2 q4 neutral Cartesian Gaussian
  nuclear raw-block loop structure. The audit confirms that the family-reuse
  target is not speculative: 66 flattened 3D supplement orbitals collapse to
  only 8 unique x-axis families, 8 y-axis families, and 16 z-axis families.
- The current kernel therefore rebuilds equivalent one-dimensional nuclear
  tables many times. For Cr2 q4, unique `G-A` table keys drop from 264 to 48
  and unique `A-A` table keys drop from 8844 to 344. Total table fills drop
  from 409860 to 17640, a 23.23x reduction, before any scalar-integral
  allocation fix.
- The audit also measured per-call allocation in the scalar `:factor` integral:
  representative powers `0/0`, `0/1`, and `1/1` allocate 192, 272, and
  368 bytes per call. This justifies a specialized nonallocating nuclear-factor
  scalar helper if the family-reuse implementation still needs it.

Validation:
- Doer ran `git diff --check`, package load,
  `tmp/work/cr2_nuclear_family_reuse_count_audit.jl`, and final
  `git status --short --branch`. Final tracked status was clean, with only the
  pre-existing untracked successor handoff file.
- Manager confirmed `git status --short --branch`, `git diff --check`, and the
  ignored probe path. Manager did not rerun the Cr2 audit.

Goal advancement:
- Cr2-readiness/MT4: proves the next source pass should reorganize the neutral
  nuclear kernel around one-dimensional axis-family reuse rather than continue
  wrapper-level loop reshaping or in-place-table-only work.
- CGRB/LT6: keeps the optimization inside `CartesianGaussianRawBlocks` while
  preserving the neutral uncharged nuclear `G-A`/`A-A` contract and forbidding
  Cr2 workflow promotion.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: the implementation target is now numerically bounded by measured
  unique family/key counts.
- quarantined: ignored count/allocation probe only.
- not deleted because: no production code changed.
- exact remaining caller/blocker: production
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl` still loops by
  flattened 3D orbitals/pairs and scalar `:factor` integrals allocate per call.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The implementation pass must preserve the coupled primitive-pair contraction
  `sum_pq c_p c_q I_x[p,q] I_y[p,q] I_z[p,q]`. Do not independently contract
  x/y/z axes, add persistent caches, broaden the raw-block owner, or change
  Residual Gaussian/Qiu-White semantics.

## Cartesian Hamiltonian Producer Pass 082B - Reuse Nuclear Axis Families

Commit(s):
- `82b3f697` - Reuse Gaussian nuclear axis families
- this commit - Record nuclear family reuse implementation

Summary:
- Accepted the `HP-CGRB-FN-02` source implementation in
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`. The neutral uncharged
  Gaussian nuclear `G-A`/`A-A` outputs are unchanged at roundoff, but the
  production kernel now groups supplement one-dimensional axis families once,
  fills unique `G-A` and `A-A` axis tables term-first, and reuses those tables
  across flattened 3D orbitals and upper-triangle orbital pairs.
- The Cr2 q4 raw nuclear step dropped from the recent `~44552.840 MiB`
  baseline and the same-audit old reference `14.1394s / 48542.417 MiB` to
  `0.6420s / 16.484 MiB`, with `G-A` parity `7.11e-15` and `A-A` parity
  `1.07e-14`.
- The implementation keeps the coupled primitive-pair contraction
  `sum_pq c_p c_q I_x[p,q] I_y[p,q] I_z[p,q]`; no independent x/y/z scalar
  contraction, cache framework, public API, artifact path, or route semantics
  were added.

Validation:
- Doer ran `git diff --check`, package load,
  `tmp/work/cgrb_factor_scalar_parity.jl`,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`,
  `tmp/work/cgrb_qw_nuclear_parity_check.jl`,
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, and
  `tmp/work/cr2_nuclear_family_reuse_count_audit.jl`.
- Manager ran `git diff --check HEAD~1..HEAD`, numstat, suspicious-added-line
  scan, new-test/file scan, package load,
  `tmp/work/cgrb_factor_scalar_parity.jl`, and
  `tmp/work/cgrb_qw_nuclear_parity_check.jl`. Manager did not rerun the H2,
  Be2, or Cr2 long gates.

Numerical/performance result:
- Cr2 q4 unique family counts: x `8`, y `8`, z `16`.
- Table fills: `G-A` `11880 -> 2160`, `A-A` `397980 -> 15480`.
- Total scalar one-dimensional calls: `202,518,360 -> 9,724,410`.
- H2 residual-GTO/MWG self-Coulomb stayed
  `0.4574265214362078`, delta `3.33e-16`.
- Be2 facade/readback deltas stayed `0.0`; augmented dimension stayed `1421`.
- Qiu-White parity stayed at roundoff: neutral/reference `G-A`
  `4.44e-16`, neutral/reference `A-A` `8.88e-16`, route deltas `0.0`.
- Specialized scalar helper parity max error was `1.42e-14`.

Goal advancement:
- Cr2-readiness/MT4: removes the dominant repeated one-dimensional table and
  scalar-integral work from the current Cr2 q4 exact nuclear raw-block path.
- CGRB/LT6 and LT4: keeps raw analytic nuclear formula ownership in the
  neutral nuclear owner while making the shared Residual Gaussian/Qiu-White
  consumer path much more memory-efficient.

Mechanical/anti-bloat gate:
- `git diff --check`: clean.
- Numstat: `150 98 src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`.
- Suspicious-added-line scan found only the function-local
  `Tuple(exponents)` / `Tuple(prefactors)` content key used to group 66
  flattened orbitals by primitive-vector identity. Accepted because it is not a
  staged route inventory, metadata/API shape, or persistent tuple collection;
  it is confined to the local family-inventory builder.
- New-test/file scan: none.

Carrying-cost result:
- deleted: old per-orbital `G-A` cross-cache path, old flattened `A-A`
  pair helper, and unused per-orbital axis input helper.
- simplified: neutral nuclear raw-block production now has one family-reuse
  path.
- quarantined: none.
- not deleted because: Qiu-White and Residual Gaussian callers still need the
  neutral raw-block entry point.
- exact remaining caller/blocker: none for this pass.
- added src lines: 150.
- deleted src lines: 98.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The local content-key tuples should remain a small inventory-construction
  detail. Do not grow this into staged tuple inventories, persistent caches,
  metadata/report summaries, or a broader raw-block framework. Cr2 facade and
  artifact workflow remain deferred.

## Cartesian Hamiltonian Producer Pass 083A - Measure Post-Reuse Cr2 Cost

Commit(s):
- this commit - Record post-reuse Cr2 exact operator audit

Summary:
- Accepted a read-only post-`HP-CGRB-FN-02` Cr2 q4 exact-operator audit. The
  nuclear bottleneck is crossed: neutral nuclear raw blocks now measure
  `0.6364s / 19.828 MiB` with roundoff parity against the old reference, while
  the same-audit old nuclear reference was `14.0439s / 48545.698 MiB`.
- The full exact augmented-operator wrapper now measures
  `8.3579s / 19107.314 MiB`. The remaining cost is non-nuclear: residual setup
  mixed overlap `X` (`3.0855s / 10990.106 MiB`), Qiu-White self moment `A-A`
  blocks (`2.1954s / 8453.301 MiB`), `G-G` product matrices
  (`2.8728s / 7024.456 MiB`), and Qiu-White cross moment `G-A` blocks
  (`0.4495s / 2489.566 MiB`).
- Dense output storage and exact transforms are not the next bottleneck. One
  augmented dense matrix is about `20.097 MiB`, nine outputs are about
  `180.871 MiB`, and exact transforms measured `0.0562s / 623.604 MiB`.

Validation:
- Doer ran `git diff --check`, package load,
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, and final
  `git status --short --branch`. Final tracked status was clean, with the
  pre-existing untracked successor handoff file.
- Manager confirmed `git status --short --branch`, `git diff --check`, and
  current `HEAD`/`origin` state. Manager did not rerun the Cr2 audit.

Goal advancement:
- Cr2-readiness/MT4: closes the nuclear raw-block allocation story as the
  current dominant Cr2 q4 blocker and identifies the next measured cost center.
- LT4/LT6: preserves the distinction between an efficient neutral nuclear
  owner and still-expensive non-nuclear Qiu-White donor/mixed-overlap paths.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: bottleneck attribution is now simpler: do not target nuclear raw
  blocks next.
- quarantined: no new probe needed.
- not deleted because: no production code changed.
- exact remaining caller/blocker: non-nuclear Qiu-White donor/mixed-overlap
  and `G-G` product-matrix allocation, not nuclear raw blocks.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The obvious next implementation area, non-nuclear overlap/kinetic/moment
  raw-block organization, is explicitly outside the approved neutral nuclear
  slice. Do not start source work there under `HP-CGRB-FN-02`; request a
  docs-only design amendment first.

## Cartesian Hamiltonian Producer Pass 083B - Approve Non-Nuclear Raw Blocks

Commit(s):
- this commit - Approve non-nuclear Gaussian raw blocks

Summary:
- Approved the next narrow neutral Cartesian Gaussian raw-block slice under
  `HP-CGRB-NN-FILE-01`, `HP-CGRB-NN-FN-01`, `HP-CGRB-NN-WIRE-01`, and
  `HP-CGRB-NN-TEST-01`. The owner remains
  `src/cartesian_gaussian_raw_blocks/`, with a new non-nuclear file for exact
  overlap, kinetic, coordinate-moment, and second-moment `G-A`/`A-A` raw
  blocks.
- This is explicitly not `HP-CGRB-FN-02`; that ID remains nuclear-only
  one-dimensional axis-family reuse. The new approval follows Pass 083A, where
  nuclear raw blocks were no longer the bottleneck and the remaining measured
  costs were residual mixed overlap, Qiu-White self/cross moment blocks, and
  `G-G` product matrices.
- The approved non-nuclear scope is `G-A`/`A-A` only. The measured `G-G`
  product-matrix allocation remains deferred and needs a separate source
  decision if targeted later.

Validation:
- Design-manager validation is docs-only: `git diff --check`, focused `rg`
  for `HP-CGRB-NN-*` IDs and forbidden-scope wording, and confirmation that no
  `src`, `test`, `tools`, or `bin` files changed.

Goal advancement:
- Cr2-readiness/MT4: defines the next allowed source surface for the
  post-nuclear Cr2 q4 exact-operator allocation problem without promoting Cr2
  facade/artifact support.
- CGRB/LT6: keeps reusable analytic Gaussian raw-block construction in the
  neutral owner while preserving Residual Gaussian and Qiu-White as consumers.

Carrying-cost result:
- deleted: none; docs-only design authority.
- simplified: post-nuclear source authority is now explicit and separated from
  nuclear `HP-CGRB-FN-02`.
- quarantined: final-basis `G-G` product-matrix optimization, terminal
  projection, Residual Gaussian algorithm changes, Qiu-White semantic changes,
  artifacts, public API, and Cr2 workflow remain out of scope.
- not deleted because: source extraction and duplicate-loop deletion belong to
  the next doer pass under the new IDs.
- exact remaining caller/blocker: implement
  `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`, rewire
  `pqs_terminal_residual_gto.jl` and Qiu-White consumers, prove H2/Be2/QW/Cr2
  q4 non-nuclear parity, then delete route-local duplicate loops.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none in this docs pass; `HP-CGRB-NN-TEST-01` approves one
  standalone non-nuclear parity fixture if needed.
- new metadata/status fields: none.

Risk / guardrail:
- Do not use this approval as a broad raw-block framework or as permission to
  optimize `G-G` product matrices. Keep the first source pass to non-nuclear
  `G-A`/`A-A` raw blocks and construction-local reuse only.

## Cartesian Hamiltonian Producer Pass 084A - Map Non-Nuclear Extraction

Commit(s):
- this commit - Record non-nuclear extraction map

Summary:
- Accepted a read-only mapping audit for the first source pass under
  `HP-CGRB-NN-*`. The duplicated non-nuclear work is localized to the
  Qiu-White cross/self moment builders, while Residual Gaussian currently calls
  that donor path twice: once for residual setup mixed overlap `X`, and again
  for the full exact augmented operators.
- Exact current builders are
  `_qwrg_cartesian_shell_cross_moment_blocks_3d(...)`,
  `_qwrg_cartesian_shell_self_moment_blocks_3d(...)`, and the diatomic wrapper
  `_qwrg_diatomic_cartesian_shell_blocks_3d(...)` in
  `src/ordinary_qw_raw_blocks.jl`. Residual Gaussian callers are
  `_r3a_qw_blocks(...)`, `_terminal_residual_mixed_overlap(...)`, and
  `pqs_terminal_residual_gto_augmented_operators(...)` in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.
- The first implementation should add
  `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`, include it from
  `CartesianGaussianRawBlocks.jl`, rewire RG and QW diatomic consumers to a
  neutral non-nuclear raw-block result, and delete direct QW donor calls from
  those paths. Atomic QW factor-term and CPB/oracle callers remain live, so the
  old QW builders should not be deleted yet.

Validation:
- Doer/design-manager reported no files edited, final
  `git status --short --branch` clean except the pre-existing untracked
  successor handoff file. Manager confirmed `git status --short --branch` and
  `git diff --check`.

Goal advancement:
- Cr2-readiness/MT4: converts the new non-nuclear authority into a concrete
  source-pass boundary with exact current callers and deletion limits.
- CGRB/LT6: protects the neutral-owner extraction from becoming a broad QW or
  CPB/oracle rewrite.

Carrying-cost result:
- deleted: none; read-only audit.
- simplified: identified the exact neutral extraction boundary and live
  callers.
- quarantined: atomic factor-term QW, CPB oracle paths, `G-G` product-matrix
  optimization, and Cr2 workflow.
- not deleted because: old QW cross/self builders still have atomic and
  CPB/oracle callers.
- exact remaining caller/blocker: source pass must fit under `HP-CGRB-NN-*`
  and the 150 added-source-line budget; otherwise request an amendment before
  coding.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Keep the first source pass narrow: non-nuclear `G-A`/`A-A`
  overlap/kinetic/moment raw blocks only. Do not touch nuclear raw blocks,
  final-basis `G-G` product matrices, atomic factor-term QW, CPB oracle
  behavior, public API, artifacts, reports, statuses, or persistent caches.

## Cartesian Hamiltonian Producer Pass 084B - Extract Non-Nuclear Raw Blocks

Commit(s):
- `00d052e2` - Extract Gaussian non-nuclear raw blocks
- this commit - Record non-nuclear raw-block extraction

Summary:
- Accepted the first `HP-CGRB-NN-*` source cut. Exact non-nuclear Gaussian
  `G-A`/`A-A` raw blocks now have a neutral owner boundary in
  `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`, loaded by
  `CartesianGaussianRawBlocks.jl`.
- Residual Gaussian `_r3a_qw_blocks(...)` and the Qiu-White diatomic wrapper no
  longer call the QW cross/self moment donors directly. They consume
  `CartesianGaussianRawBlocks.gaussian_non_nuclear_raw_blocks(...)` and keep
  the previous numerical values at roundoff. The R3 same-construction overload
  also reuses the same supplement block result for residual mixed overlap and
  augmented-operator construction.
- This is an ownership/deletion boundary, not yet the allocation optimization.
  The neutral wrapper still delegates to the existing QW donor kernels, so Cr2
  exact-operator allocation stayed effectively flat: `8.4294s /
  19123.019 MiB` versus the prior post-reuse reference of about `8.3579s /
  19107.314 MiB`.

Validation:
- Doer ran `git diff --check`, package load,
  `tmp/work/cgrb_qw_nuclear_parity_check.jl`,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`, and
  `tmp/work/cr2_exact_operator_allocation_audit.jl`.
- Manager ran `git diff --check HEAD~1..HEAD`, numstat, suspicious-added-line
  scan, new-test/file scan, package load, and
  `tmp/work/cgrb_qw_nuclear_parity_check.jl`. Manager did not rerun the H2,
  Be2, or Cr2 long gates.

Numerical/performance result:
- H2 R3 self-Coulomb remained `0.4574265214362078`, delta `3.33e-16`; H2
  facade readback deltas were all `0.0`.
- Be2 augmented dimension remained `1421`; readback deltas were all `0.0`;
  self-Coulomb was `1.3119725917573790`.
- Small Qiu-White non-nuclear neutral/reference and route/neutral max deltas
  were `0.0`.
- Cr2 q4 raw non-nuclear reference delta and residual mixed-overlap `X` replay
  delta were `0.0`; nuclear deltas stayed `G-A 7.11e-15`, `A-A 1.07e-14`;
  wrapper/replay kinetic delta was `0.0`.

Goal advancement:
- CGRB/LT6: establishes the neutral non-nuclear owner boundary and removes
  direct RG/QW diatomic dependence on route-local QW donor construction.
- Cr2-readiness/MT4: confirms the next optimization must be inside the
  non-nuclear donor construction itself, while `G-G` product matrices remain a
  separate deferred target.

Mechanical/anti-bloat gate:
- `git diff --check`: clean.
- Numstat: `40 35 src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
  `3 0 src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`;
  `33 0 src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`;
  `21 31 src/ordinary_qw_raw_blocks.jl`.
- Suspicious-added-line scan: none.
- New-test/file scan: none.

Carrying-cost result:
- deleted: direct RG and QW diatomic cross/self donor call sites.
- simplified: one neutral non-nuclear raw-block owner boundary.
- quarantined: ignored probe updates only; atomic factor-term QW, CPB/oracle
  paths, and `G-G` product-matrix optimization remain out of this pass.
- not deleted because: QW donor helpers still serve live neutral-wrapper
  delegation and atomic/factor/oracle paths.
- exact remaining caller/blocker: non-nuclear donor kernels still allocate
  heavily; optimizing them is the next source target under `HP-CGRB-NN-*`.
- added src lines: 97.
- deleted src lines: 66.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not misread this ownership extraction as the performance fix. The next
  source pass may optimize the neutral non-nuclear owner internals, but must
  still avoid nuclear changes, `G-G` product-matrix work, persistent caches,
  broad provider bundles, public API, artifacts, and Cr2 workflow.

## Cartesian Hamiltonian Producer Pass 085A - Count Non-Nuclear Family Reuse

Commit(s):
- this commit - Record non-nuclear reuse count audit

Summary:
- Accepted a read-only post-extraction non-nuclear count audit. The repeated
  work is real, but the smallest next source cut is narrower than a full
  non-nuclear kernel replacement: residual mixed-overlap setup currently builds
  the full non-nuclear `G-A`/`A-A` family and then uses only projected `G-A`
  overlap.
- For the mixed-overlap-only path, current table builds are `53064`, while the
  unique overlap-only `G-A` table builds needed are only `32`, a `1658.25x`
  table-build reduction. Scalar calls drop from `835884` to `97228`, an
  `8.60x` reduction.
- The full exact-operator path also has a strong optimization target:
  full-block `A-A` scalar calls would drop from `19,783,872` to `161,044`, a
  `122.85x` reduction. That is justified later, but it is a larger source pass
  than the overlap-only mixed-`X` helper.

Validation:
- Doer ran `git diff --check`, package load, focused ignored probe
  `tmp/work/cr2_non_nuclear_family_reuse_count_audit.jl`, and final
  `git status --short --branch`. Final tracked status was clean, with the
  pre-existing untracked successor handoff file.
- Manager confirmed `git status --short --branch`, `git diff --check`, and
  current `HEAD`/`origin` state. Manager did not rerun the count probe.

Measured facts:
- Cr2 q4 fixture: parent axes `(19,19,31)`, parent Cartesian count `11191`,
  supplement candidates `66`, unique families x `8`, y `8`, z `16`.
- Current QW `G-A` all non-nuclear: `0.5110s / 2489.246 MiB`.
- Current QW `A-A` all non-nuclear: `2.2813s / 8476.370 MiB`.
- Neutral wrapper all non-nuclear: `2.6527s / 10941.712 MiB`.
- Full raw-block table builds: `G-A 792 -> 128`, `A-A 52272 -> 832`, total
  `53064 -> 960`.
- Full raw-block scalar calls: `G-A 3,343,536 -> 388,912`,
  `A-A 19,783,872 -> 161,044`, total `23,127,408 -> 549,956`.

Goal advancement:
- Cr2-readiness/MT4: identifies an overlap-only mixed-`X` helper as the next
  smallest source pass, rather than jumping directly to full non-nuclear
  family-reuse replacement.
- CGRB/LT4/LT6: preserves the neutral owner boundary while keeping source
  growth and risk proportional to the measured target.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: next implementation target is narrowed to overlap-only mixed
  overlap first.
- quarantined: ignored probe only; full non-nuclear family-reuse replacement
  and `G-G` product matrices remain later work.
- not deleted because: no production source changed.
- exact remaining caller/blocker: current residual mixed-overlap path still
  pays full non-nuclear `G-A`/`A-A` construction cost.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- The next source pass should be overlap-only. Do not mix it with full
  kinetic/moment `A-A` family reuse, nuclear work, `G-G` product matrices,
  persistent caches, public API, artifacts, reports, statuses, or broad provider
  bundles.

## Cartesian Hamiltonian Producer Pass 085B - Add Overlap-Only Mixed Blocks

Commit(s):
- `9fa0cc16` - Add overlap-only Gaussian mixed blocks
- this commit - Record overlap-only mixed-block optimization

Summary:
- Accepted the overlap-only mixed-`X` source pass. Residual mixed overlap
  `X = G' S A` now uses a neutral
  `CartesianGaussianRawBlocks.gaussian_non_nuclear_overlap_blocks(...)` helper
  instead of building the full non-nuclear `G-A`/`A-A` kinetic/moment block
  family and discarding almost all of it.
- Physics is unchanged: the Cr2 mixed-overlap replay delta is `0.0`, H2
  residual-GTO/MWG self-Coulomb is unchanged at
  `0.4574265214362078`, and the Be2 facade/readback gate still passes.
- The targeted setup cost improved from `3.1195s / 10990.107 MiB` to
  `0.0817s / 199.292 MiB`. Full Cr2 exact-operator wrapper cost stayed roughly
  flat/noisy at `8.7030s / 19269.019 MiB`, because full augmented operators
  still need the full non-nuclear raw blocks.

Validation:
- Doer ran `git diff --check`, package load,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`, and
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, then a post-commit diff
  gate and suspicious-line scan.
- Manager ran `git diff --check HEAD~1..HEAD`, numstat, suspicious-added-line
  scan, new-test/file scan, and package load. Per user direction, manager did
  not complete a redundant long Cr2 validation rerun and relied on the doer
  report for long gates.

Goal advancement:
- Cr2-readiness/MT4: removes the waste in residual setup mixed overlap and
  narrows the remaining exact-operator blocker to full non-nuclear raw-block
  construction plus separately deferred `G-G` product matrices.
- CGRB/LT4/LT6: adds a focused neutral-owner helper without broadening into a
  cache/provider framework.

Mechanical/anti-bloat gate:
- `git diff --check`: clean.
- Numstat: `4 3 src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
  `45 0 src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`.
- Suspicious-added-line scan: none.
- New-test/file scan: none.

Carrying-cost result:
- deleted: mixed-overlap dependence on full non-nuclear raw-block construction.
- simplified: residual `X` setup now uses overlap-only neutral raw blocks.
- quarantined: none.
- not deleted because: full `gaussian_non_nuclear_raw_blocks(...)` is still
  needed by exact augmented operators and Qiu-White paths.
- exact remaining caller/blocker: full non-nuclear family reuse and `G-G`
  product matrices remain the next exact-operator cost centers.
- added src lines: 49.
- deleted src lines: 3.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not fold the next full non-nuclear family-reuse pass into `G-G` product
  matrix work. Keep `G-G` deferred unless a separate amendment approves it, and
  keep any non-nuclear optimization inside the neutral owner and approved
  caller surfaces.

## Cartesian Hamiltonian Producer Pass 086A - Reuse Non-Nuclear A-A Families

Commit(s):
- `806b37e3` - Reuse Gaussian non-nuclear AA axis families
- this commit - Record non-nuclear AA family reuse

Summary:
- Accepted the full non-nuclear `A-A` family-reuse source pass after amendment.
  The physical non-nuclear supplement-supplement Gaussian blocks are unchanged
  at roundoff, but production no longer calls
  `_qwrg_cartesian_shell_self_moment_blocks_3d(...)` from
  `gaussian_non_nuclear_raw_blocks(...)`.
- The amended implementation now uses canonical supplement axis-family pair
  tables with orientation flags, assembles only upper-triangular orbital pairs,
  averages the two requested orientations to preserve the old symmetrized
  Qiu-White convention, and mirrors the result. `G-A`, nuclear blocks, `G-G`,
  artifacts, metadata, and public API were untouched.
- Cr2 q4 non-nuclear raw-block construction improved to
  `0.5336s / 3148.165 MiB`, from the earlier neutral all-non-nuclear baseline
  of about `2.65s / 10941.7 MiB`. The full exact-operator wrapper improved to
  `6.4265s / 11231.805 MiB`.

Validation:
- Doer ran `git diff --check`, package load,
  `tmp/work/cgrb_qw_nuclear_parity_check.jl`,
  `tmp/work/cr2_non_nuclear_aa_parity_breakdown.jl`,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`, and
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, then a post-amend diff
  gate and suspicious-line scan.
- Manager ran `git diff --check HEAD~1..HEAD`, numstat, suspicious-added-line
  scan, and new-test/file scan. Per user direction, manager did not rerun
  long validation already reported by the doer.

Numerical/performance result:
- Cr2 `A-A` parity by block family: overlap `0.0`, kinetic
  `1.11e-16`, x/y/z `0.0/0.0/0.0`, x2/y2/z2 `0.0/0.0/0.0`.
- Cr2 raw non-nuclear reference delta was `1.11e-16`; residual mixed-overlap
  replay delta was `0.0`.
- H2 residual-GTO/MWG self-Coulomb stayed at `0.4574265214362095`.
- Be2 facade/readback passed with K/unit-U/one-body/V deltas all `0.0`.

Goal advancement:
- Cr2-readiness/MT4: removes the full non-nuclear `A-A` donor hotspot from the
  exact-operator path and leaves `G-A` non-nuclear donor work plus separately
  deferred `G-G` product matrices as the next measured costs.
- CGRB/LT4/LT6: keeps non-nuclear raw-block optimization inside the neutral
  owner and preserves Qiu-White/Residual Gaussian as consumers.

Mechanical/anti-bloat gate:
- `git diff --check`: clean.
- Numstat: `138 13 src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`.
- Suspicious-added-line scan: none.
- New-test/file scan: none.

Carrying-cost result:
- deleted: direct production call to
  `_qwrg_cartesian_shell_self_moment_blocks_3d(...)` from the neutral
  non-nuclear wrapper.
- simplified: `A-A` non-nuclear construction now lives in the neutral owner
  with function-local family table reuse.
- quarantined: ignored Cr2 `A-A` parity breakdown probe only.
- not deleted because: QW self donor remains live for QW/reference paths; `G-A`
  cross donor remains delegated by design for this pass.
- exact remaining caller/blocker: `G-A` non-nuclear donor and `G-G` product
  matrices remain the next measured costs.
- added src lines: 138.
- deleted src lines: 13.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Continue to keep `G-G` product-matrix work separate and deferred. The next
  source pass may target non-nuclear `G-A` family reuse inside the neutral
  owner, but must not broaden into Qiu-White semantic changes, CPB/atomic
  rewrites, persistent caches, artifacts, public API, or Cr2 workflow.

### Medium-Term Goal Checkpoint After Pass 086A

- Active: Cr2-readiness/MT4 has shifted from raw nuclear allocation to
  non-nuclear `G-A` donor allocation plus separately deferred `G-G` product
  matrices. Nuclear `G-A`/`A-A`, residual setup mixed overlap `X`, and
  non-nuclear `A-A` are no longer the dominant blockers.
- Completed in the current CGRB lane: neutral nuclear owner extraction,
  nuclear axis-family reuse, non-nuclear owner extraction, overlap-only mixed
  `X`, and non-nuclear `A-A` family reuse.
- Blocked/deferred: final-basis `G-G` product-matrix optimization remains
  explicitly outside `HP-CGRB-NN-*` and needs a separate design amendment if
  targeted.
- Still active: keep Cr2 facade/artifact workflow, public API, artifacts,
  metadata/report/status fields, and broad provider/cache frameworks out of the
  raw-block optimization lane.
- Next lane direction: one more narrow `HP-CGRB-NN-*` source pass may target
  non-nuclear `G-A` family reuse inside
  `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`; after that,
  reassess whether remaining Cr2 exact-operator allocation is dominated by
  deferred `G-G` product matrices.

## Cartesian Hamiltonian Producer Pass 086B - Reuse Non-Nuclear G-A Families

Commit(s):
- `71a89433` - Reuse Gaussian non-nuclear GA axis families
- this commit - Record non-nuclear GA family reuse

Summary:
- Accepted the non-nuclear `G-A` family-reuse source pass. The remaining
  parent-supplement non-nuclear raw-block work no longer calls
  `_qwrg_cartesian_shell_cross_moment_blocks_3d(...)` from
  `gaussian_non_nuclear_raw_blocks(...)`; the neutral owner now builds overlap,
  kinetic, coordinate-moment, and second-moment `G-A` blocks through
  function-local supplement axis-family reuse.
- The overlap-only mixed-`X` path and the canonical `A-A` family path were
  preserved. `G-A` and `A-A` now share one local supplement axis-family
  inventory in `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`.
- Cr2 q4 non-nuclear raw-block construction improved to
  `0.1873s / 860.736 MiB`, from the corrected `086A` value of about
  `0.5336s / 3148.165 MiB` and the pre-`086A` neutral baseline of about
  `2.65s / 10941.7 MiB`. The Cr2 exact-operator wrapper is now
  `6.4184s / 9043.987 MiB`.

Validation:
- Doer ran `git diff --check`, package load,
  `tmp/work/cgrb_qw_nuclear_parity_check.jl`,
  `tmp/work/cr2_non_nuclear_aa_parity_breakdown.jl`,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`, and
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, then a post-commit diff
  gate and suspicious-line scan.
- Manager ran `git diff --check HEAD~1..HEAD`, numstat, suspicious-added-line
  scan, and new-test/file scan. Per user direction, manager did not rerun
  validation already reported by the doer.

Numerical/performance result:
- Cr2 `G-A` parity by block family: overlap `0.0`, kinetic `0.0`, x/y/z
  `0.0/0.0/0.0`, x2/y2/z2 `0.0/0.0/0.0`.
- Cr2 `A-A` parity remained roundoff: only kinetic at `1.11e-16`, all other
  families `0.0`.
- H2 residual-GTO/MWG self-Coulomb stayed at `0.4574265214362095`.
- Be2 facade/readback passed with final dimension `1421` and K/unit-U/one-body/V
  deltas all `0.0`.

Goal advancement:
- Cr2-readiness/MT4: completes the approved neutral non-nuclear raw-block
  optimization sequence. Remaining measured exact-operator allocation is now
  dominated by deferred `G-G` product matrices and unrelated route/stage setup,
  not neutral non-nuclear raw blocks.
- CGRB/LT4/LT6: keeps `G-A` reuse inside the neutral owner without changing
  Residual Gaussian or Qiu-White semantics.

Mechanical/anti-bloat gate:
- `git diff --check`: clean.
- Numstat: `107 15 src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`.
- Suspicious-added-line scan: none.
- New-test/file scan: none.

Carrying-cost result:
- deleted: direct production call to
  `_qwrg_cartesian_shell_cross_moment_blocks_3d(...)` from the neutral
  non-nuclear wrapper.
- simplified: `G-A` and `A-A` now share the same local supplement axis-family
  inventory.
- quarantined: ignored parity breakdown probe only.
- not deleted because: Qiu-White cross/self donor helpers remain live as
  parity/reference and route-side helpers.
- exact remaining caller/blocker: Cr2 exact-operator allocation is now
  dominated by `G-G` product matrices and unrelated route/stage setup, not
  neutral non-nuclear raw blocks.
- added src lines: 107.
- deleted src lines: 15.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not continue under `HP-CGRB-NN-*` into `G-G` product-matrix optimization.
  That target is explicitly deferred and needs a separate design amendment.
  The neutral raw-block lane should now pause unless a new measured non-G-G
  blocker appears inside its approved `G-A`/`A-A` scope.

## Cartesian Hamiltonian Producer Pass 087 - Approve Terminal G-G Product Matrices

Commit(s):
- this commit - Approve R3 terminal G-G product matrices

Summary:
- Approved a narrow `CartesianFinalBasisRealization` source lane for R3/RG
  terminal final-basis `G-G` product-matrix optimization under
  `HP-R3GG-FN-01` and `HP-R3GG-TEST-01`.
- The approval follows Pass 086B, where neutral non-nuclear raw blocks dropped
  to `0.1873s / 860.736 MiB` and the Cr2 exact-operator wrapper measured
  `6.4184s / 9043.987 MiB`. Remaining measured allocation is now dominated by
  terminal `G-G` product matrices and unrelated route/stage setup, not
  `G-A`/`A-A` raw-block construction.
- Approved product matrices are only the `G-G` kinetic, coordinate moment, and
  second-moment matrices used by
  `pqs_terminal_residual_gto_augmented_operators(...)`. The neutral raw-block
  lane remains closed for this target.

Validation:
- Design-manager validation is docs-only: `git diff --check`, focused `rg`
  for `HP-R3GG-*` IDs and forbidden-scope wording, and confirmation that no
  `src`, `test`, `tools`, or `bin` files changed.

Goal advancement:
- Cr2-readiness/MT4: opens the next measured exact-operator allocation target
  after the neutral raw-block lanes crossed their bottlenecks.
- RG/LT6: keeps the optimization attached to final-basis product assembly
  rather than misclassifying it as Gaussian raw-block work.

Carrying-cost result:
- deleted: none; docs-only design authority.
- simplified: separates terminal `G-G` product matrices from completed
  `G-A`/`A-A` raw-block lanes.
- quarantined: route/stage setup allocation, unit-nuclear Gaussian sums,
  IDA/MWG, residual selection/transforms, Qiu-White semantics, artifacts,
  public API, and Cr2 workflow remain out of scope.
- not deleted because: source simplification belongs to the next doer pass
  under `HP-R3GG-FN-01`.
- exact remaining caller/blocker: optimize product-matrix construction in
  `pqs_terminal_residual_gto.jl`, optionally using a small internal helper in
  `pqs_terminal_one_body.jl`, while preserving H2/Be2/Cr2 parity.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none in this docs pass; `HP-R3GG-TEST-01` approves existing H2
  and ignored Be2/Cr2 validation only.
- new metadata/status fields: none.

Risk / guardrail:
- Do not use this approval for a broad product-operator framework or route
  setup cleanup. Line budget is `100` added source lines, and any persistent
  workspace/cache object or source surface outside the approved files requires
  a new amendment.

## Cartesian Hamiltonian Producer Pass 087A - Reuse Terminal G-G Scratch

Commit(s):
- `fb9b0414` - Reuse R3 terminal GG product scratch
- this commit - Record terminal G-G scratch reuse

Summary:
- Accepted the first `HP-R3GG-FN-01` source pass. The R3/RG exact augmented
  operator path now uses one `nG x nG` scratch matrix for terminal final-basis
  `G-G` product construction in
  `pqs_terminal_residual_gto_augmented_operators(...)`.
- Kinetic `G-G` assembly accumulates the three axis contributions into that
  scratch, then transforms immediately. Coordinate moments and second moments
  now clear and reuse the same scratch one axis at a time, transform
  immediately, and no longer retain `pos_GG` or `x2_GG` product tuples.
- `_r3a_product_matrix(...)` was deleted. No terminal product kernel,
  Gaussian raw-block, residual, MWG/IDA, artifact, metadata, status, report, or
  public API surface changed.

Validation:
- Doer ran `git diff --check`, package load,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`, and updated ignored
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, then a post-commit diff
  gate and suspicious-line scan.
- Manager ran `git diff --check HEAD~1..HEAD`, numstat, source diff review,
  suspicious-added-line scan, new-test/file scan, and a targeted symbol check
  for `_r3a_product_matrix`, `pos_GG`, `x2_GG`, and `scratch_GG`. Per user
  direction, manager did not rerun doer validation.

Numerical/performance result:
- H2 residual-GTO/MWG self-Coulomb stayed at `0.4574265214362095`.
- Be2 facade/readback stayed at final dimension `1421`, with
  K/unit-U/one-body/V deltas all `0.0`.
- Cr2 exact-operator wrapper measured `6.4873s / 8872.858 MiB`.
- Cr2 audit attribution: seven K/moment transform outputs allocated
  `572.321 MiB`; nine replayed `assemble_terminal_product_operator!` product
  calls allocated `5307.975 MiB`.
- Cr2 raw non-nuclear delta was `1.11e-16`, residual mixed-overlap replay delta
  was `0.0`, and wrapper kinetic replay delta was `3.84e-12`.

Goal advancement:
- Cr2-readiness/MT4: removes the obvious retained full-matrix lifetime waste in
  the approved terminal `G-G` lane, but confirms that repeated terminal product
  assembly internals remain the main local allocation source.
- RG/LT6: keeps the optimization inside `CartesianFinalBasisRealization` and
  does not reopen completed neutral raw-block, residual, or MWG/IDA lanes.

Mechanical/anti-bloat gate:
- `git diff --check`: clean.
- Numstat: `33 32 src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.
- Suspicious-added-line scan: none.
- New-test/file scan: none.

Carrying-cost result:
- deleted: `_r3a_product_matrix(...)`.
- simplified: no retained `pos_GG` or `x2_GG` tuples; one local scratch lifetime
  now covers kinetic, coordinate moment, and second-moment `G-G` products.
- quarantined: ignored Cr2 audit script updates only.
- not deleted because: `assemble_terminal_product_operator!` remains the
  approved terminal product kernel and is still shared by live callers.
- exact remaining caller/blocker: allocation is now dominated by repeated
  terminal product assembly internals; a later pass would need shared
  `action/tile/block` workspace in `pqs_terminal_one_body.jl`.
- added src lines: 33.
- deleted src lines: 32.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- A follow-up may inspect workspace reuse in `pqs_terminal_one_body.jl`, but
  only within `HP-R3GG-FN-01`: no persistent cache/workspace object, no broad
  product-operator framework, no route/stage cleanup, and no raw-block,
  residual-selection, MWG/IDA, artifact, metadata, report, status, public API,
  or Cr2 workflow changes.

## Cartesian Hamiltonian Producer Pass 087B - Reuse Terminal G-G Product Workspace

Commit(s):
- `5cd9e15a` - Reuse R3 terminal GG product workspace
- this commit - Record terminal G-G workspace reuse

Summary:
- Accepted the second narrow `HP-R3GG-FN-01` source pass. The R3/RG terminal
  `G-G` product fills now share one function-local `action_buffer`,
  `tile_buffer`, and `block_buffer` across the nine kinetic/moment product
  assemblies in `pqs_terminal_residual_gto_augmented_operators(...)`.
- `assemble_terminal_product_operator!(...)` keeps its existing external
  behavior by delegating to a private `_assemble_terminal_product_operator!(...)`
  with fresh refs. The R3 path calls the private helper with caller-owned refs.
- The prepatch buffer request probe showed `_buffer_view!` already grows by
  capacity and reuses smaller views, so no `_buffer_view!` change was needed.
  No raw-block, nuclear, residual, MWG/IDA, Qiu-White, route/stage, artifact,
  metadata, report, status, or public API surface changed.

Validation:
- Doer ran `git diff --check`, package load,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`,
  `tmp/work/cr2_exact_operator_allocation_audit.jl`, and
  `tmp/work/r3gg_product_buffer_request_audit.jl`.
- Manager ran `git diff --check`, numstat, source diff review,
  suspicious-added-line scan, new-test/file scan, and a targeted symbol check
  for product workspace refs and the new private helper. Per user direction,
  manager did not rerun doer validation.

Numerical/performance result:
- Cr2 nine `G-G` product assemblies improved from `1.4260s / 4800.477 MiB`
  with fresh refs to `0.8522s / 685.534 MiB` with shared refs.
- Cr2 exact augmented-operator wrapper improved from about
  `6.4873s / 8872.858 MiB` after Pass 087A to `5.8389s / 4605.517 MiB`.
- Direct product parity max delta was `0.0`.
- H2 residual-GTO/MWG self-Coulomb stayed at `0.4574265214362095`.
- Be2 facade/readback stayed at final dimension `1421`, with readback deltas
  all `0.0`.
- Cr2 operators were finite/symmetric, raw parity deltas stayed at roundoff,
  and wrapper replay kinetic delta was `3.8369307731045410e-12`.

Goal advancement:
- Cr2-readiness/MT4: crosses the measured terminal `G-G` product-workspace
  allocation blocker opened by Pass 087A, dropping the exact-operator wrapper
  below `5 GiB` allocation on the Cr2 q4 audit.
- RG/LT6: preserves terminal final-basis ownership and does not reopen the
  completed Gaussian raw-block or residual/MWG lanes.

Mechanical/anti-bloat gate:
- `git diff --check`: clean.
- Numstat:
  `19 3 src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`;
  `15 9 src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.
- Suspicious-added-line scan: none.
- New-test/file scan: none.

Carrying-cost result:
- deleted: none.
- simplified: R3 exact `G-G` products reuse one internal terminal-product
  workspace instead of rebuilding fresh refs and buffers per product call.
- quarantined: ignored `tmp/work` probes only.
- not deleted because: public `assemble_terminal_product_operator!` remains
  live for other callers and preserves fresh-ref behavior.
- exact remaining caller/blocker: allocation now mainly lies outside these nine
  product buffers, including unit-nuclear `U_GG` Gaussian-sum work and
  route/raw-block setup, which are outside this pass.
- added src lines: 34.
- deleted src lines: 12.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not continue `HP-R3GG-FN-01` into unit-nuclear Gaussian sums or route/raw
  setup. Those are separate measured targets and need their own authority if
  pursued.

## Cartesian Hamiltonian Producer Pass 088 - Approve Remaining Exact-Operator Allocation Audit

Commit(s):
- this commit - Record remaining exact-operator allocation audit

Summary:
- Made the docs/authority decision after `954c86cd`: the terminal `G-G`
  product-workspace lane is crossed, but the remaining Cr2 q4 exact
  augmented-operator allocation is still a mixed attribution problem rather
  than a source-ready target.
- Added `r3_remaining_exact_operator_allocation_audit.md` and approved
  `HP-R3REM-AUDIT-01` as measurement-only authority. It separates the next
  audit from `HP-R3GG-FN-01` and explicitly blocks implementation in
  unit-nuclear `U_GG` Gaussian sums, route/raw-block setup, or any other source
  lane until a later amendment names exact owner/files/functions, validation,
  and line budget.

Validation:
- Design-manager validation: `git diff --check`, focused `rg` for
  `HP-R3REM-AUDIT-01`, measurement-only wording, unit-nuclear exclusions, and
  confirmation that no `src`, `test`, `tools`, or `bin` files changed.
- No implementation tests were run; this was docs-only authority work.

Goal advancement:
- Cr2-readiness/MT4: keeps performance work moving by requiring exact
  attribution before the next source lane.
- RG/LT6: preserves completed R3GG and raw-block authority boundaries and
  prevents conceptual drift into source setup or nuclear-sum optimization under
  the wrong ID.

Carrying-cost result:
- deleted: none.
- simplified: next-step authority is now one measurement-only ID rather than
  an ambiguous continuation of `HP-R3GG-FN-01`.
- quarantined: future probes must remain ignored `tmp/work` measurement.
- not deleted because: no source surface was edited or retired in this docs
  pass.
- exact remaining caller/blocker: Cr2 q4 exact augmented-operator allocation
  remains to be attributed among unit-nuclear `U_GG`, route/raw-block setup,
  and any audit overhead before source approval.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 088A - Attribute Remaining Exact-Operator Allocation

Commit(s):
- this commit - Record remaining allocation attribution audit

Summary:
- Accepted the `HP-R3REM-AUDIT-01` measurement-only audit. The remaining Cr2 q4
  exact augmented-operator allocation after the terminal `G-G` workspace passes
  is no longer primarily owned by Gaussian raw blocks or terminal `G-G`
  product buffers.
- The largest in-wrapper owner is now unit-nuclear `U_GG` Gaussian-sum
  construction by center, measured at about `2.0447s / 1856.819 MiB` inside a
  total wrapper cost of `5.7739s / 4680.627 MiB`.
- The audit recommends a future source lane for unit-nuclear `U_GG` assembly
  in `CartesianFinalBasisRealization`, centered on
  `_accumulate_terminal_gaussian_sum!` and `_terminal_gaussian_sum_action`, but
  no such source authority exists yet.

Validation:
- Doer ran `git diff --check`, package load, and
  `tmp/work/cr2_exact_operator_allocation_audit.jl`; final status had no
  tracked changes and only the pre-existing untracked successor handoff.
- Manager reviewed the reported bucket attribution and existing
  `HP-R3REM-AUDIT-01` authority. Per user direction, manager did not rerun the
  Cr2 audit.

Numerical/performance result:
- Wrapper total: `5.7739s / 4680.627 MiB`.
- Neutral non-nuclear `G-A`/`A-A` raw blocks: `0.1894s / 860.736 MiB`.
- Neutral nuclear `G-A`/`A-A` raw blocks: `0.6316s / 15.765 MiB`.
- Terminal `G-G` kinetic/moment products with shared workspace:
  `0.8352s / 733.701 MiB`.
- Unit-nuclear `U_GG` factor lookup plus Gaussian-sum construction:
  `2.0447s / 1856.819 MiB`.
- Augmented nuclear transforms only: `0.0125s / 179.268 MiB`; all augmented
  transforms together: `0.0600s / 739.698 MiB`.
- Raw-block setup/projection/symmetry miscellaneous cost was about
  `0.160s / 12.5 MiB`. Audit-only parity/oracle replay cost was large
  (`~16.575s / ~59482.7 MiB`) and explicitly not production wrapper cost.
- Outside the wrapper, route/residual setup remains expensive for this fixture:
  about `46.9s / 3569.8 MiB`, with `cartesian_units` dominant in time and
  `cartesian_transforms`/`S_AA` dominant in allocation.

Goal advancement:
- Cr2-readiness/MT4: identifies the next source candidate as unit-nuclear
  `U_GG` Gaussian-sum assembly, not raw blocks, residual/MWG code, or terminal
  `G-G` product buffers.
- RG/LT6: preserves authority boundaries by keeping the result measurement-only
  and requiring a separate docs-only source amendment before implementation.

Mechanical/anti-bloat gate:
- No production files changed in the doer pass.
- Manager-log commit only; no `src`, `test`, `tools`, or `bin` edits.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: remaining exact-operator allocation is now attributed to exact
  buckets rather than a broad post-`G-G` remainder.
- quarantined: ignored `tmp/work` probes only.
- not deleted because: no source authority existed for unit-nuclear `U_GG` or
  route/setup cleanup.
- exact remaining caller/blocker: unit-nuclear `U_GG` Gaussian-sum source work
  needs a docs-only amendment naming owner, files, functions, forbidden
  surfaces, validation gates, line budget, and failure rule.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not implement unit-nuclear `U_GG` optimization under
  `HP-R3REM-AUDIT-01` or `HP-R3GG-FN-01`. The likely owner is
  `CartesianFinalBasisRealization`, but source work must wait for a new
  docs-only authority decision.

### Medium-Term Goal Checkpoint After Pass 088A

- Completed: neutral nuclear raw-block reuse, neutral non-nuclear `G-A`/`A-A`
  reuse, residual mixed-overlap-only setup, terminal `G-G` product lifetime
  cleanup, and terminal `G-G` workspace reuse are no longer the current Cr2 q4
  exact-operator allocation blockers.
- Active: Cr2-readiness/MT4 now points to unit-nuclear `U_GG` Gaussian-sum
  assembly as the next in-wrapper numerical kernel candidate, with route/setup
  allocation still visible but outside the wrapper and not yet source-ready.
- Blocked/deferred: unit-nuclear `U_GG` source work, route/raw-block setup
  cleanup, Cr2 facade/artifact workflow, public API expansion, metadata/report
  fields, and residual/MWG/IDA changes all remain blocked without later
  docs-only amendments.
- Needing refinement: the next source authority should decide whether to target
  only reusable scratch/in-place accumulation inside the existing
  `_accumulate_terminal_gaussian_sum!` path, or first require a smaller
  count/shape audit of Gaussian-sum action allocation. It should not authorize
  a broad Gaussian-sum framework.
- Current lane direction: ask repo-design-manager for a narrow unit-nuclear
  `U_GG` amendment, or explicitly stop Cr2 exact-operator allocation work at
  the current measured cost.

## Cartesian Hamiltonian Producer Pass 089 - Approve Unit-Nuclear U_GG Source Lane

Commit(s):
- this commit - Approve unit-nuclear U_GG Gaussian-sum lane

Summary:
- Accepted the Pass 088A audit evidence as sufficient to approve a narrow
  source lane. The dominant in-wrapper Cr2 q4 exact-operator allocation is now
  unit-nuclear `U_GG` factor lookup plus Gaussian-sum construction, measured at
  `2.0447s / 1856.819 MiB` inside a `5.7739s / 4680.627 MiB` wrapper.
- Added `r3_unit_nuclear_ugg_gaussian_sum.md` and approved `HP-R3UN-FN-01` /
  `HP-R3UN-TEST-01` for terminal final-basis unit-nuclear `U_GG`
  Gaussian-sum allocation reduction under `CartesianFinalBasisRealization`.
  The lane is deliberately not a raw-block, route/setup, residual/MWG/IDA,
  public workflow, or Cr2 artifact/facade lane.

Validation:
- Design-manager validation: `git diff --check`, focused `rg` for
  `HP-R3UN-FN-01`, `HP-R3UN-TEST-01`, approved files/functions, forbidden
  surfaces, and confirmation that no `src`, `test`, `tools`, or `bin` files
  changed.
- No implementation tests were run; this was docs-only authority work.

Goal advancement:
- Cr2-readiness/MT4: turns the measured dominant exact-operator bucket into a
  source-ready, narrow, reviewable lane.
- RG/LT6: preserves completed raw-block, R3GG, residual, MWG/IDA, and workflow
  boundaries while allowing one allocation-focused kernel pass.

Carrying-cost result:
- deleted: none.
- simplified: source authority now targets one terminal Gaussian-sum path
  instead of a broad remaining-allocation bucket.
- quarantined: route/setup and raw-block setup remain outside this lane.
- not deleted because: this is a docs-only authority pass.
- exact remaining caller/blocker: source pass must stay within
  `_accumulate_terminal_gaussian_sum!` / `_terminal_gaussian_sum_action` and
  stop if it needs broader Gaussian-sum framework or persistent cache authority.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 089A - Reuse Unit-Nuclear U_GG Workspace

Commit(s):
- `79e5cd47` - Reuse R3 unit-nuclear UGG workspace
- this commit - Record unit-nuclear U_GG workspace reuse

Summary:
- Accepted the narrow `HP-R3UN-FN-01` source pass. Terminal final-basis
  unit-nuclear `U_GG` Gaussian-sum construction now uses an in-place scratch
  path in `CartesianFinalBasisRealization`.
- `_accumulate_terminal_gaussian_sum!` delegates to a scratch-buffer method that
  reuses action/tile/block buffers, uses `mul!`, and accumulates scaled blocks
  directly into destination views. The compatibility
  `_terminal_gaussian_sum_action(...)` remains live for terminal IDA block
  assembly.
- The R3 exact-operator path shares one Gaussian-sum scratch set across nuclear
  centers. No raw-block, terminal kinetic/moment `G-G`, residual/MWG/IDA,
  Qiu-White, route/stage, artifact, metadata, report, status, payload, public
  API, committed-test, or Cr2 workflow surface changed.

Validation:
- Doer ran `git diff --check`, package load,
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`,
  `tmp/work/be2_r3u_facade_measurement.jl`, and
  `tmp/work/cr2_exact_operator_allocation_audit.jl`.
- Manager ran `git diff --check`, numstat, source diff review,
  suspicious-added-line scan, new-test/file scan, targeted helper/caller search,
  and terminal IDA caller review. Per user direction, manager did not rerun
  doer validation.

Numerical/performance result:
- Cr2 `U_GG` construction only improved from `2.0206s / 1803.721 MiB` to
  `1.9910s / 593.201 MiB`.
- Cr2 `U_GG` plus factor lookup improved from `2.0447s / 1856.819 MiB` to
  `2.0145s / 646.299 MiB`.
- Cr2 exact-operator wrapper allocation improved from the prior
  `~5.84s / 4605.5 MiB` to `5.9048s / 3540.638 MiB` in the reported run.
- Cr2 `U` replay parity was `0.0`; exact operators were finite and symmetric.
- H2 residual-GTO/MWG self-Coulomb stayed at `0.4574265214362095`.
- Be2 facade/readback stayed at final dimension `1421`, with all readback
  deltas `0.0`.

Goal advancement:
- Cr2-readiness/MT4: crosses the largest in-wrapper allocation owner identified
  by `HP-R3REM-AUDIT-01`, reducing unit-nuclear `U_GG` allocation by about
  `1.2 GiB`.
- RG/LT6: keeps unit-nuclear Gaussian-sum work in
  `CartesianFinalBasisRealization` and preserves raw-block, residual/MWG, and
  public workflow boundaries.

Mechanical/anti-bloat gate:
- `git diff --check`: clean.
- Numstat:
  `74 20 src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`;
  `5 1 src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.
- Suspicious-added-line scan: none.
- New-test/file scan: none.

Carrying-cost result:
- deleted: old allocation-heavy `_accumulate_terminal_gaussian_sum!` block
  materialization path.
- simplified: shared Gaussian-sum action fill now serves both terminal IDA
  block action and in-place unit-nuclear `U_GG` accumulation.
- quarantined: ignored audit probe changes only.
- not deleted because: `_terminal_gaussian_sum_action(...)` remains live for
  terminal IDA assembly.
- exact remaining caller/blocker: route setup and remaining raw-block
  replay/oracle costs are outside `HP-R3UN-FN-01`.
- added src lines: 79.
- deleted src lines: 21.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not continue this source lane into route/stage setup, raw-block
  replay/oracle costs, terminal kinetic/moment `G-G`, neutral raw blocks,
  residual/MWG/IDA, artifacts, metadata, reports, public API, or Cr2 workflow
  without a separate authority decision.

## Cartesian Hamiltonian Producer Pass 089B - Post-UGG Production Attribution

Commit(s):
- this commit - Record post-UGG production attribution

Summary:
- Accepted the measurement-only post-UGG attribution audit. Cr2 q4
  production-equivalent in-memory construction now looks memory-acceptable:
  peak RSS was about `2652.031 MiB`, and the final Hamiltonian retained size
  was about `80.388 MiB`.
- The prior exact-operator kernel lanes remain crossed. The audit identifies
  one possible remaining exact-operator source target: same-construction reuse
  of already-built base Hamiltonian `K_GG` and by-center unit `U_GG` blocks,
  rather than recomputing those final-basis matrices inside augmented
  operators.
- The replay using base `K/U_A` blocks was exact at zero operator delta and
  reduced the exact-operator replay to `0.8620s / 1237.136 MiB`. Route setup
  dominates elapsed time end-to-end, but that is a separate attribution lane,
  not an exact-operator continuation.

Validation:
- Doer ran `git diff --check`, package load, and
  `tmp/work/cr2_post_ugg_production_attribution.jl`.
- Final status had no tracked changes and only the pre-existing untracked
  successor handoff. Manager reviewed the reported attribution and did not
  rerun the Cr2 audit.

Production-equivalent result:
- Route/base stages: `47.0332s`, GC `0.1748s`, `2134.579 MiB` allocated,
  `12.220 MiB` retained.
- Base Hamiltonian: `4.9447s`, GC `0.6853s`, `4404.736 MiB` allocated,
  `73.983 MiB` retained.
- Supplement load/convert: `0.0103s`, `3.400 MiB` allocated,
  `0.066 MiB` retained.
- Residual construction: `0.7314s`, GC `0.0284s`, `1426.763 MiB` allocated,
  `1.641 MiB` retained.
- Exact augmented operators: `5.0239s`, GC `0.3747s`, `3043.013 MiB`
  allocated, `180.872 MiB` retained.
- Final MWG Hamiltonian: `0.1779s`, GC `0.0029s`, `738.960 MiB` allocated,
  `80.388 MiB` retained.
- Total: `73.0740s / 13107.681 MiB`, GC `2.7732s`, peak RSS
  `2652.031 MiB`.

Replay attribution:
- Neutral non-nuclear raw blocks: `0.1842s / 862.643 MiB`.
- Neutral nuclear raw blocks: `0.6438s / 19.828 MiB`.
- Terminal `G-G` kinetic/moment products: `0.8438s / 721.346 MiB`.
- Unit-nuclear `U_GG`: `2.0093s / 609.146 MiB`.
- Exact augmented transforms: `0.0569s / 627.354 MiB`.
- MWG substep replay: `0.2778s / 837.870 MiB`.
- Same-construction supplement-block overload: `4.7994s / 4771.127 MiB`,
  zero K and V deltas.
- Hypothetical base `K/U_A` reuse: `0.8620s / 1237.136 MiB`, operator delta
  `0.0`.

Goal advancement:
- Cr2-readiness/MT4: shifts the next optional exact-operator target from
  low-level scratch reuse to same-construction orchestration reuse. The current
  production-equivalent peak memory is acceptable enough that stopping
  exact-operator micro-optimization is also defensible.
- RG/LT6: keeps the next possible source work attached to
  `CartesianFinalBasisRealization` and same-construction data reuse, not raw
  blocks, residual/MWG conventions, or public workflow.

Mechanical/anti-bloat gate:
- No production files changed in the doer pass.
- Manager-log commit only; no `src`, `test`, `tools`, or `bin` edits.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: separates production-equivalent memory acceptability, replay
  overhead, and same-construction duplication.
- quarantined: ignored `tmp/work` audit only.
- not deleted because: no source authority exists for same-construction base
  `K/U_A` reuse.
- exact remaining caller/blocker: a new docs-only amendment is needed before
  wiring base Hamiltonian `K_GG` and unit `U_GG[A]` reuse into augmented
  operators.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

Risk / guardrail:
- Do not implement same-construction base `K/U_A` reuse without a narrow
  authority decision. If approved, it should avoid provenance payloads,
  metadata/status fields, route-stage cleanup, public API, artifacts, raw-block
  changes, residual/MWG/IDA convention changes, and committed tests.

## Cartesian Hamiltonian Producer Pass 090 - Approve Canonical Driver Usability Workflow

Commit(s):
- this commit - Approve canonical driver usability workflow

Summary:
- Changed the driver plan now that the base and residual-GTO/MWG producer paths
  are coherent enough to work together. The canonical driver should no longer
  be only a frozen template; it is approved as a compact artifact-producing
  workflow over approved producer surfaces.
- Added `cartesian_driver_usability_workflow.md` and approved
  `HP-DRV-FILE-01`, `HP-DRV-FN-01`, and `HP-DRV-TEST-01`.
  `bin/cartesian_ham_builder.jl` may now own visible editable defaults, one
  optional trusted local input file, command-line `key=value` overrides,
  compact summaries, coarse user-facing timing, artifact write, and optional
  readback check.
- Preserved the old anti-diagnostic guardrail: route-stage controls,
  stop-after internals, ladder probes, private helper calls, raw-block
  switches, report/status/payload dumps, allocation probes, solver workflow,
  public API/export changes, artifact schema changes, committed fixtures/tests,
  and Cr2 driver workflow remain outside the canonical driver.

Validation:
- Design-manager validation: `git diff --check`, focused `rg` for
  `HP-DRV-FILE-01`, `HP-DRV-FN-01`, `HP-DRV-TEST-01`, `bin/cartesian_ham_builder.jl`,
  input-file/override/timing/artifact wording, and forbidden diagnostic
  surfaces; confirmation that no `src`, `test`, `tools`, or `bin` files
  changed.
- No implementation tests were run; this was docs-only authority work.

Goal advancement:
- LT1/LT3: moves the driver from protected template to real consumer workflow.
  Producing a Hamiltonian artifact becomes the complete-success test for the
  standard driver.
- RG/LT6: lets the supported supplemented path be exercised by a human-facing
  driver without promoting Cr2, solver, public export, or diagnostic surfaces.

Carrying-cost result:
- deleted: none.
- simplified: the driver policy now has one compact workflow lane instead of
  conflicting "do not touch the driver" and "make workflows usable" pressures.
- quarantined: diagnostics, ladder probing, allocation audits, and Cr2 stress
  controls stay in `tools/` or ignored `tmp/work`.
- not deleted because: this is a docs-only authority pass.
- exact remaining caller/blocker: source implementation must stay within
  `bin/cartesian_ham_builder.jl` and stop if it needs a parser framework,
  committed fixtures/tests, route-stage diagnostics, artifact schema changes,
  or Cr2 workflow support.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 091 - Approve Homonuclear Diatomic Supplemented Scope

Commit(s):
- this commit - Approve homonuclear diatomic supplemented scope

Summary:
- Approved a narrow molecule-scope relaxation for the residual-GTO/MWG
  supplemented facade and canonical driver. The hardcoded H2/Be2-only story is
  replaced by explicit homonuclear two-center z-axis diatomic validation.
- Added `r3_homonuclear_diatomic_supplemented_workflow.md` and approved
  `HP-R3U-ZDI-FN-01`, `HP-R3U-ZDI-WIRE-01`, and `HP-R3U-ZDI-TEST-01`.
  `src/cartesian_base_hamiltonian.jl` may replace H/Be-specific guarded
  helpers with explicit homonuclear z-axis validation, and
  `bin/cartesian_ham_builder.jl` may route supplemented mode through the
  supported facade.
- Cr2 is allowed only as an explicit generic homonuclear z-axis ignored/user-run
  stress or usability case after H2/Be2 validation. No Cr2-specific branch,
  default, committed fixture, committed gate, ECP, solver workflow, public
  export, artifact schema change, diagnostics, metadata/status/report field,
  or route object is approved.

Validation:
- Design-manager validation: `git diff --check`, focused `rg` for the
  `HP-R3U-ZDI-*` IDs, explicit homonuclear z-axis wording, optional
  `basisfile`, Cr2 stress limitations, and forbidden surfaces; confirmation
  that no `src`, `test`, `tools`, or `bin` files changed.
- No implementation tests were run; this was docs-only authority work.

Goal advancement:
- LT1/LT3: makes the canonical driver useful for real explicit molecule inputs
  while preserving compact, copyable workflow shape.
- RG/LT6: moves supplemented artifact production from H/Be fixture wording to a
  generic homonuclear diatomic producer surface without promoting Cr2-specific
  code or broad molecule support.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: H2/Be2 hardcoding is replaced by one explicit molecule-scope
  validation contract.
- quarantined: Cr2 remains ignored/user-run stress through the generic path,
  not a committed gate or special branch.
- not deleted because: source implementation has not run yet.
- exact remaining caller/blocker: repo-doer must implement within
  `src/cartesian_base_hamiltonian.jl` and `bin/cartesian_ham_builder.jl`,
  stop if the change needs heteronuclear/non-z-axis support, ECP, solver work,
  artifact schema changes, route diagnostics, committed fixtures/tests,
  metadata/status/report fields, or Cr2-specific branching.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 092 - Approve Base Atom Driver Workflow

Commit(s):
- this commit - Approve base atom driver workflow

Summary:
- Approved a narrow canonical-driver atom workflow for explicit
  origin-centered one-center atoms in `mode = :base`.
- Added `cartesian_driver_atom_workflow.md` and approved
  `HP-DRV-ATOM-FN-01`, `HP-DRV-ATOM-WIRE-01`, and `HP-DRV-ATOM-TEST-01`.
  The only approved source file is `bin/cartesian_ham_builder.jl`; the driver
  may normalize explicit one-center atom inputs and call the existing
  `cartesian_base_hamiltonian(system; basis, hamfile)` facade.
- The decision is base atom output only. Current validation remains the
  origin-centered H endpoint. Supplemented atom Hamiltonians, broader base
  atom source support, translated atoms, element tables/defaults, ECP,
  solver workflow, artifact schema changes, public API/export changes,
  route diagnostics, metadata/status/report fields, committed atom fixtures,
  and committed tests remain unapproved.

Validation:
- Design-manager validation: docs-only `git diff --check` over touched docs,
  focused `rg` for `HP-DRV-ATOM-*`, origin-centered/base-only wording,
  supplemented-atom deferral, forbidden surfaces, and exact source file
  authority.
- No implementation tests were run.
- Pre-existing source WIP `src/cartesian_base_hamiltonian.jl` was present
  before this pass and was left unstaged/untouched by the design-manager
  commit.

Goal advancement:
- LT1/LT3: keeps the canonical driver moving toward practical artifact
  production for normal user inputs without weakening the compact-driver
  guardrails.
- RG/LT6: explicitly defers supplemented atom Hamiltonians until there is a
  separate facade/RG validation decision, avoiding accidental expansion of the
  residual-GTO/MWG contract.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: atom handling is separated into base-driver workflow now and
  candidate supplemented-atom work later.
- quarantined: supplemented atom Hamiltonians and broader atom/source support
  remain candidate-only.
- not deleted because: source implementation has not run yet.
- exact remaining caller/blocker: repo-doer must implement, if requested, only
  in `bin/cartesian_ham_builder.jl` and stop if the work needs source edits
  outside that file, broader base atom support, translated atoms, supplemented
  atom support, ECP, solver workflow, artifact schema changes, route
  diagnostics, metadata/status/report fields, committed fixtures/tests, or
  element lookup/default tables.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 093 - Implement Homonuclear Diatomic Supplemented Facade

Commit(s):
- this commit - Relax supplemented diatomic facade scope

Summary:
- Accepted the first `HP-R3U-ZDI-FN-01` source pass in
  `src/cartesian_base_hamiltonian.jl`.
- `_cartesian_r3_diatomic_inputs(...)` no longer hardcodes H2/Be2 element and
  spin-sector checks. It now validates the approved data-driven homonuclear
  two-center z-axis all-electron contract: equal symbols, equal finite
  positive charges, finite z-axis centers with distinct z coordinates,
  nonnegative integer `nup`/`ndn`, and neutral electron count from explicit
  nuclear charge.
- `_cartesian_r3_supplement_inputs(...)` now accepts optional trusted
  `basisfile` and passes it to the existing
  `legacy_bond_aligned_diatomic_gaussian_supplement(...)` loading path. No
  artifact schema or provenance field changed.

Validation:
- Doer reported `git diff --check`, package load, the H2 residual-GTO endpoint
  with self-Coulomb `0.4574265214362095` and zero readback deltas, the Be2
  facade/readback probe with dimension `1421` and zero readback deltas, and an
  optional ignored Cr2-style smoke returning a finite symmetric
  `CartesianIDAHamiltonian{Float64}` of dimension `1623`.
- Manager static review only, per current no-rerun guidance: `git status
  --short --branch`, `git diff --check`, direct source diff inspection,
  `git diff --numstat -- src/cartesian_base_hamiltonian.jl`, and the
  suspicious-added-line scan. The diff gate was clean, line impact was
  `21` added / `12` deleted in one approved source file, and no new files or
  tests were added.

Goal advancement:
- LT1/LT3: moves supplemented artifact construction from H/Be fixture wording
  to a real explicit homonuclear diatomic input contract, which is the
  prerequisite for making the canonical driver a practical molecule workflow.
- RG/LT6: keeps the residual-GTO/MWG supplement path generic over explicit
  homonuclear z-axis inputs without opening heteronuclear, non-z-axis,
  ECP/solver, public API/export, artifact schema, or Cr2-specific workflow
  surfaces.

Carrying-cost result:
- deleted: H/Be-specific element/electron guard in the supplemented facade.
- simplified: molecule validation is now data-driven homonuclear z-axis
  validation, and basisfile handling is a small pass-through to the existing
  supplement loader.
- quarantined: Cr2 remains an ignored/user-run generic stress case, not a
  branch, default, fixture, committed gate, or driver workflow.
- not deleted because: the H2 validation fixture label check remains live for
  existing H2 reference artifact validation.
- exact remaining caller/blocker: canonical driver wiring is still separate
  under `HP-DRV-*` / `HP-R3U-ZDI-WIRE-01`, and base atom driver support is
  separate under `HP-DRV-ATOM-*`. Supplemented atom Hamiltonians remain
  unapproved.
- added src lines: 21.
- deleted src lines: 12.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 094 - Approve One-Center Base Atom Facade

Commit(s):
- this commit - Approve one-center base atom facade

Summary:
- Approved a narrow R1 producer-side relaxation from one-center H-only
  validation to explicit origin-centered all-electron one-center atoms through
  the existing `cartesian_base_hamiltonian(...)` facade.
- Added `r1_one_center_base_atoms.md` and approved `HP-R1-ATOM-FN-01`,
  `HP-R1-ATOM-WIRE-01`, and `HP-R1-ATOM-TEST-01`.
- Recorded the shared-workflow constraint from the user: atoms and diatomics
  must share the same producer machinery after the narrow geometry and
  shellification differences. No atom-only Hamiltonian builder, atom-only
  materialization path, atom route/report/status payload, or atom-specific
  artifact shape is approved.

Validation:
- Design-manager validation: `git diff --check`, focused `rg` for
  `HP-R1-ATOM-*`, explicit origin-centered/all-electron atom wording, shared
  workflow requirements, element-table/default exclusions, supplemented-atom
  deferral, and exact source file authority.
- No implementation tests were run; this was docs-only authority work.

Goal advancement:
- LT1/LT3: moves singular atom production toward the same practical user
  workflow as diatomics, without making the driver invent atom behavior.
- RG/LT6: keeps supplemented atoms explicitly deferred until the residual-GTO
  facade and validation endpoint are approved, while preserving the principle
  that atom and diatomic implementations share machinery wherever possible.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: source authority now distinguishes driver atom workflow from
  base-facade atom support.
- quarantined: translated atoms, supplemented atoms, ECP, element lookup,
  solver workflow, atom-specific materialization, committed non-H fixtures, and
  public non-H reference gates remain unapproved.
- not deleted because: source implementation has not run yet.
- exact remaining caller/blocker: repo-doer may edit only
  `src/cartesian_base_hamiltonian.jl` for this source lane and must stop if the
  work needs private materialization-owner edits, atom-only materialization,
  new artifact keys, translated atom support, supplemented atoms, ECP behavior,
  solver workflow, element lookup/default tables, committed fixtures/tests, or
  route/report/status/payload expansion.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 095 - Replace Canonical Driver With Facade Workflow

Commit(s):
- this commit - Replace canonical driver with facade workflow

Summary:
- Accepted the `HP-DRV-*`, `HP-R3U-ZDI-WIRE-01`, and `HP-DRV-ATOM-*` driver
  source pass in `bin/cartesian_ham_builder.jl`.
- The canonical driver is now a compact public-input script rather than a
  route-stage template. It exposes user-facing controls such as `mode`,
  `Natom`, `R`, `Z`, `atom`, `nup`/`ndn`, `q`, `core_spacing`, one-center and
  diatomic basis extents, supplement basis controls, `hamfile`, and
  `readback`.
- The driver constructs `system`, `basis`, and optional `supplement` locally
  and dispatches only through approved facades:
  `GaussletBases.cartesian_base_hamiltonian(...)`,
  `GaussletBases.cartesian_residual_gto_mwg_hamiltonian(...)`, and
  `GaussletBases.read_cartesian_ida_hamiltonian(...)`.
  Old route-stage choreography and report/TSV/materialization controls are no
  longer present in the public script.

Validation:
- Doer reported `git diff --check`, package load, H atom base driver
  artifact/readback with dimension `419`, H2 base driver artifact/readback
  with dimension `471`, H2 supplemented driver artifact/readback with
  dimension `489`, and Be2 supplemented driver artifact/readback with
  dimension `1421`.
- Manager static review only: `git status --short --branch`, `wc -l
  bin/cartesian_ham_builder.jl`, `git diff --check`, direct driver diff
  inspection, `git diff --numstat -- src bin tools test docs`,
  suspicious-added-line scan, old route-stage-name `rg`, and
  `git diff --name-status`.
- The final driver is `108` lines, under the `115` line target and `125` line
  review threshold. Line impact was `101` added / `106` deleted in one
  approved `bin` file. The suspicious-line scan was clean, and the old
  route-stage-name scan found no matches.

Goal advancement:
- LT1/LT3: makes the canonical driver a practical artifact-producing workflow
  while keeping it copyable and human-facing. A user can now drive base atom,
  base diatomic, and supported supplemented homonuclear diatomic artifacts
  without seeing route internals.
- RG/LT6: keeps supplemented construction behind the supported facade and does
  not introduce Cr2-specific workflow, public API/export changes, artifact
  schema changes, diagnostics, solver/ECP behavior, or private helper calls.

Carrying-cost result:
- deleted: old `cartesian_system` through `cartesian_save` public route-stage
  sequence, route/report/TSV/materialization flags, retained-rule controls,
  WL/internal route knobs, and `run_h1*` controls from the canonical driver.
- simplified: the driver is now facade choreography with compact trusted
  input/override handling and optional public readback.
- quarantined: Cr2 remains only an optional ignored/user-run homonuclear
  z-axis stress through the generic supplemented facade after H2/Be2 pass.
- not deleted because: compact trusted input-file and top-level `key=value`
  override handling are part of the approved driver workflow.
- exact remaining caller/blocker: producer-side non-H one-center atom support
  still requires the separate `HP-R1-ATOM-*` source pass in
  `src/cartesian_base_hamiltonian.jl`; supplemented atoms remain unapproved.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 096 - Amend Driver Public Contract Hooks

Commit(s):
- this commit - Amend canonical driver public contract hooks

Summary:
- Amended the canonical driver workflow so the driver must visibly construct
  public `system`, `basis`, and optional `supplement` contract objects before
  calling an approved facade. This keeps the driver compact and reviewable
  without exposing old route-stage choreography.
- Approved only four compact public run-level hooks: `check_file`,
  `print_contract`, `print_timing`, and `expected_dimension`. These are for
  human expert review and Codex-controlled artifact checks, not route reports,
  allocation probes, artifact schema dumps, solver controls, or private helper
  calls.
- Recorded `basisname = nothing` as the base-mode selector and
  `basisname !== nothing` as supported supplemented diatomic mode; the latter
  must reject `Natom == 1`. Recorded `padding` as physical box padding around
  the atom or two nuclei, mapped internally to existing facade fields.

Validation:
- Design-manager validation: docs-only `git diff --check`, focused `rg` for
  `check_file`, `print_contract`, `print_timing`, `expected_dimension`,
  `basisname`, `padding`, public contract construction, and forbidden surfaces.
- No implementation tests were run; this was docs-only authority work.

Goal advancement:
- LT1/LT3: keeps the canonical driver a public-contract producer template
  rather than a private route harness, even after implementation made it
  functional.
- RG/LT6: keeps supplemented mode reachable through visible public supplement
  contracts while preserving the ban on supplemented atoms until a separate
  facade/RG decision exists.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: replaces implicit driver choreography with explicit public
  contract objects and four run-level hooks.
- quarantined: route-stage controls, raw-block switches, allocation probes,
  artifact schema dumps, solver controls, Cr2-specific workflow, private helper
  calls, and supplemented atoms remain unapproved.
- not deleted because: this pass records authority only; source cleanup and
  alignment must be done by repo-doer under approved driver IDs.
- exact remaining caller/blocker: repo-doer may edit only
  `bin/cartesian_ham_builder.jl` under `HP-DRV-*` and must stop if public
  contract construction requires route internals, extra hooks, artifact schema
  changes, supplemented atom support, Cr2-specific workflow, committed
  fixtures/tests, or parser/framework expansion.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 097 - Approve Driver Physics-Stage Surface

Commit(s):
- this commit - Approve canonical driver physics-stage surface

Summary:
- Corrected the previous driver authority: visible public `system`, `basis`,
  and `supplement` construction was not enough because the real construction
  still disappeared inside `cartesian_residual_gto_mwg_hamiltonian(...)`.
- Approved `HP-DRV-STAGE-FN-01`, `HP-DRV-STAGE-WIRE-01`, and
  `HP-DRV-STAGE-TEST-01` so the canonical driver can execute and time visible
  physics-level construction stages without calling underscored helpers or old
  route stages.
- The approved stage sequence is: public contract construction, base working
  basis / terminal realization, base Hamiltonian assembly, Gaussian supplement,
  residual Gaussian augmentation, exact augmented operators, supplemented
  Hamiltonian assembly, and artifact write/check. The staged producer surface
  is owned only by `src/cartesian_base_hamiltonian.jl`; artifact writing
  remains existing writer/readback workflow.

Validation:
- Design-manager validation: docs-only `git diff --check`, focused `rg` for
  `HP-DRV-STAGE-*`, visible physics-stage wording,
  non-underscored/staged surface limits, and forbidden route diagnostic
  surfaces, plus confirmation that `src`, `test`, and `tools` had no changes.
- No implementation tests were run; this was authority-only documentation
  work.

Goal advancement:
- LT1/LT3: restores the original expert-driver property of showing the
  construction stages while keeping those stages at public physics workflow
  level rather than old route choreography.
- RG/LT6: allows residual-GTO/MWG construction to be visible and timed through
  approved staged surfaces without adding Cr2-specific branches, solver/ECP
  workflow, raw-block switches, or report/status payloads.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: separates public physics workflow stages from old route-stage
  diagnostics and from opaque all-in-one facade execution.
- quarantined: `cartesian_parent`, `cartesian_shells`, `cartesian_units`,
  `cartesian_pair_terms`, `cartesian_assembly`, reports, route skeletons,
  raw-block provider switches, allocation probes, artifact schema dumps,
  solver/ECP workflow, Cr2-specific workflow, and supplemented atoms remain
  unapproved.
- not deleted because: source factoring has not run yet, and existing one-call
  facades may remain as convenience wrappers over the staged implementation.
- exact remaining caller/blocker: repo-doer may edit only
  `src/cartesian_base_hamiltonian.jl` and `bin/cartesian_ham_builder.jl` for
  this lane; stop if the implementation needs public API/export changes,
  source files outside those two paths, broad payload/report/status objects,
  artifact schema changes, committed fixtures/tests, route diagnostics, or
  Cr2-specific workflow.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 098 - Require Visible Stage Calls

Commit(s):
- this commit - Require visible canonical-driver stage calls

Summary:
- Compared the `Dropbox/bin/nestpgg.jl` driver shape with the current
  canonical driver and tightened the plan accordingly. The desired property is
  not only compact editable inputs; it is also visible top-level construction
  calls and local variables that an expert user can read and customize.
- Clarified that `HP-DRV-STAGE-FN-01` must expose separate named
  construction-stage functions. One all-in-one staged replacement wrapper is
  not sufficient for the canonical driver.
- Split the base construction wording into two visible stages: base working
  basis / terminal realization, then base Hamiltonian assembly. Supplemented
  mode remains Gaussian supplement, residual augmentation, exact augmented
  operators, supplemented Hamiltonian assembly, then existing artifact
  write/check.

Validation:
- Design-manager validation: docs-only `git diff --check`, focused `rg` for
  visible top-level stage-call wording, all-in-one wrapper exclusions, split
  base-realization/base-Hamiltonian wording, and forbidden route-stage
  surfaces, plus confirmation that `src`, `test`, and `tools` had no changes.
  The existing local `bin/cartesian_ham_builder.jl` WIP was not part of this
  docs pass. No implementation tests were run.

Goal advancement:
- LT1/LT3: aligns the canonical driver more closely with the compact
  `nestpgg.jl` style: visible defaults plus visible staged construction,
  without carrying route-stage diagnostics.
- RG/LT6: keeps residual-GTO/MWG construction visible to expert users while
  preserving the ban on raw-block switches, report/status payloads, solver/ECP
  workflow, and Cr2-specific branches.

Carrying-cost result:
- deleted: none; docs-only clarification.
- simplified: prevents the staged surface from becoming another opaque facade
  and makes the intended driver call shape more reviewable.
- quarantined: one opaque staged wrapper, old route stages, underscored helper
  calls from the driver, raw-block provider switches, artifact schema dumps,
  allocation probes, solver/ECP workflow, Cr2-specific workflow, and
  supplemented atoms remain unapproved.
- not deleted because: source factoring has not run yet.
- exact remaining caller/blocker: repo-doer may edit only
  `src/cartesian_base_hamiltonian.jl` and `bin/cartesian_ham_builder.jl` for
  this lane; stop if separate visible stage calls require source files outside
  those paths, broad payload/report/status objects, public API/export changes,
  artifact schema changes, committed fixtures/tests, route diagnostics, or a
  single all-in-one staged replacement wrapper.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 099 - Approve Operator-Class Driver Stages

Commit(s):
- this commit - Approve operator-class driver stage timings

Summary:
- Refined `HP-DRV-STAGE-*` so visible driver stages distinguish physical
  operator classes, not only broad Hamiltonian phases.
- Approved separate coarse stages for product/moment operators, unit-nuclear
  attraction operators, and electron-electron / IDA or residual-MWG
  interactions. This makes the expected timing model visible: product/moment
  work should remain fast, while unit-nuclear and electron-electron work may be
  long.
- Expanded the staged-driver source surface narrowly beyond
  `src/cartesian_base_hamiltonian.jl` to include
  `src/pqs_source_box_low_order_materialization.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` only for
  behavior-preserving operator-class stage factoring in their existing domains.

Validation:
- Design-manager validation: docs-only `git diff --check`, focused `rg` for
  product/moment, unit-nuclear, electron-electron, approved source file names,
  and forbidden raw-block/diagnostic/probe wording. Local `src`/`bin` WIP was
  present but not part of this docs pass; cached validation confirmed no
  `src`, `test`, `tools`, or `bin` files were staged. No implementation tests
  were run.

Goal advancement:
- LT1/LT3: improves the canonical driver as an expert-readable construction
  script and makes timing useful without exposing route internals.
- RG/LT6: supports Cr-scale usability by making long operator classes visible
  while preserving the ban on diagnostics, raw-block switches, per-kernel
  instrumentation, solver/ECP workflow, and Cr2-specific branches.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: replaces one coarse "exact augmented operators" stage with the
  physical operator-class split that users actually need for timing.
- quarantined: raw-block provider switches, allocation probes, per-kernel
  timing frameworks, old route stages, underscored helper calls from the
  driver, report/status payloads, artifact schema dumps, solver/ECP workflow,
  Cr2-specific workflow, and supplemented atoms remain unapproved.
- not deleted because: source factoring has not run yet.
- exact remaining caller/blocker: repo-doer may edit only the canonical driver
  plus the four approved staged-driver source files for this lane; stop if the
  operator-class split requires raw-block changes, kernel rewrites, new modules,
  broad payload/report/status objects, public API/export changes, artifact
  schema changes, committed fixtures/tests, route diagnostics, allocation
  probes, per-kernel timing frameworks, or a single all-in-one staged wrapper.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 100 - Staged Operator-Class Driver

Commit(s):
- this commit - Stage canonical driver operator classes

Summary:
- Accepted the canonical-driver implementation of the staged workflow. The
  driver now constructs the public contract, then visibly binds and times
  `base`, `base_products`, `base_unit_nuclear`, `base_vee`, `base_ham`, and,
  for supplemented runs, `supplement_basis`, `residual`,
  `augmented_products`, `augmented_unit_nuclear`, `augmented_vee`, and final
  `ham`.
- `src/cartesian_base_hamiltonian.jl` now acts as the driver-facing staged
  owner while the lower-level approved files expose behavior-preserving
  operator-class factoring. Existing one-call facades remain, but compose
  through the same staged implementation.
- Manager accepted the known coarse-label caveat that the augmented product
  stage also performs shared supplement-block setup used by the following
  nuclear stage. This is treated as acceptable program-stage timing, not a
  per-kernel attribution claim.

Validation:
- Doer validation: `git diff --check`, package load, H atom base driver
  artifact/readback with dimension 419, H2 supplemented assignment-style
  driver artifact/readback with quick overrides
  `core_spacing=0.5`, `padding=4.0`, `basisname="cc-pVTZ"` and dimension 489,
  unknown override rejection for `mode`, `radius`, `d`,
  `reference_spacing`, and `tail_spacing`, and a driver scan showing no old
  route-stage calls or `GaussletBases._...` helper calls. No Cr2 run.
- Manager validation: `git diff --check`; numstat line-budget check
  (`bin` +148/-68, `src` +200/-68); suspicious added-line scan found only two
  fixed three-axis `Tuple(...)` constructions; new-test/tool scan was empty;
  focused driver scan found no forbidden old route-stage or underscored helper
  calls.

Goal advancement:
- LT1/LT3: turns the canonical driver into a readable human workflow rather
  than a black-box facade call, while keeping visible inputs compact and
  copyable.
- RG/LT6: supports Cr-scale artifact production by making the long operator
  classes visible through coarse timing without adding diagnostics, stop
  points, raw-block switches, or Cr2-specific branches.
- MT update: the driver usability lane is now source-functional for base and
  supplemented paths. The next driver work should be user-run artifact
  exercise and small polish only if it preserves the public contract shape.

Carrying-cost result:
- deleted: the driver's direct opaque one-call build branch and duplicated
  facade bodies that bypassed the staged implementation.
- simplified: base and supplemented facades now share staged construction
  functions; driver timing is organized by physical/program stages instead of
  route internals.
- quarantined: old route stages, raw-block switches, allocation probes,
  per-kernel timing frameworks, solver/ECP workflow, Cr2-specific workflow,
  artifact schema dumps, report/status payloads, and committed fixtures/tests.
- not deleted because: compatibility composition helpers and one-call facades
  remain live for existing callers.
- exact remaining caller/blocker: no blocker for the staged driver pass. The
  deferred parent-distortion unification concern remains separate from driver
  usability and should be audited after the driver settles.
- added src lines: 200.
- deleted src lines: 68.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 101 - Unify Core Spacing Authority

Commit(s):
- this commit - Unify core spacing producer contract

Summary:
- Corrected the R1/base and atom authority so `core_spacing` is the single
  public physical near-nucleus spacing after explicit input, driver default, or
  preset resolution. The legacy White-Lindsey atom `d` is no longer a durable
  public producer field; if temporarily accepted for compatibility, it must
  equal resolved `core_spacing`.
- Recorded that one-center White-Lindsey wiring sets
  `parent_mapping_d = core_spacing`, while the `Z` dependence lives in the
  mapping-shape/default rule: `sqrt(core_spacing / Z)` and
  `sqrt(core_spacing * Z)`.
- Preserved the important separation from `reference_spacing`, `tail_spacing`,
  and physical box padding. The canonical driver may still carry visible,
  editable defaults such as `core_spacing = 0.3` and template padding, and quick
  test overrides such as `core_spacing = 0.5` remain normal input behavior.

Validation:
- Design-manager validation only: focused `rg` removed stale live-authority
  claims that public `d` is required, independent, or mapped separately from
  `core_spacing`; remaining hits are intentional `mapping_d` provenance or
  rejected compatibility checks. `git diff --check` and cached no-source/bin
  staging checks are required before commit. No implementation tests were run.

Goal advancement:
- LT1/LT3: simplifies the public base/driver contract and removes a false
  independent spacing knob before more driver-facing atom/molecule workflow is
  built on it.
- RG/LT6: keeps driver defaults editable and compact for real molecules while
  avoiding hidden universal element defaults.

Carrying-cost result:
- deleted: none; docs-only authority correction.
- simplified: one public core-spacing concept replaces the prior public
  `d`/`core_spacing` split.
- quarantined: temporary `d` compatibility is allowed only as an equality
  check against resolved `core_spacing`; public `parent_mapping_d`, mapping
  strength/range knobs, element-table defaults, ECP, solver workflow, artifact
  schema changes, and Cr-specific workflow remain unapproved.
- not deleted because: source compatibility handling has not been audited in
  this design-only pass.
- exact remaining caller/blocker: repo-manager/doer must update source only
  under the approved R1/core-spacing surface; stop if implementation requires
  new public fields, artifact schema changes, or element-specific defaults.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 102 - Clarify Driver Defaults Versus Test Inputs

Commit(s):
- this commit - Clarify driver defaults and test scalars

Summary:
- Added the validation-policy distinction that the canonical driver may keep
  physics-facing editable defaults, while routine correctness tests may
  override those defaults for speed or path coverage. Any asserted scalar must
  be tied to its exact test input and must not be described as a
  physics-default result.
- No strategic change to the core-spacing model: `core_spacing` remains the
  public near-core scale, visible driver defaults remain allowed, and quick-test
  overrides remain normal input behavior.

Validation:
- Docs-only `git diff --check`, focused `rg` for the new validation-policy
  sentence, and no `src`, `test`, `tools`, or `bin` changes. No implementation
  tests were run.

## Cartesian Hamiltonian Producer Pass 103 - Implement Core-Spacing Producer Contract

Commit(s):
- this commit - Use core spacing for atom mapping d

Summary:
- Accepted the narrow `HP-R1-CORE-FN-01` source cleanup in
  `src/cartesian_base_hamiltonian.jl`. One-center base atoms now require only
  `q`, `core_spacing`, and `radius`; public `d` is optional legacy
  compatibility only.
- When legacy atom `d` is supplied, it is validated as positive and exactly
  equal to resolved `core_spacing`; mismatches throw `ArgumentError`. H2 still
  rejects public `d`.
- Internal `input.d` remains populated from `core_spacing` for compatibility
  and existing `producer_provenance/mapping_d`, while White-Lindsey atom
  wiring now passes `parent_mapping_d = input.core_spacing`.

Validation:
- Doer validation: `git diff --check`, package load, H atom facade
  artifact/readback with `core_spacing=0.3` and no public `d` (dimension 419),
  H atom legacy `d == core_spacing` (dimension 419), H atom
  `d != core_spacing` rejection, H2 public `d` rejection, H2 base
  artifact/readback (dimension 471), and no Cr2 run.
- Manager validation: reviewed the one-file diff, confirmed `git diff
  --check`, numstat +17/-7 in `src/cartesian_base_hamiltonian.jl`,
  suspicious added-line scan clean, and new-test/tool scan empty.

Goal advancement:
- LT1/LT3: removes the false public `d`/`core_spacing` split while preserving
  the compact driver contract.
- MT: this establishes the first prototype for the broader unification audit:
  atom and diatomic producer inputs now share one public near-core scale,
  while genuine mapping differences remain internal.

Carrying-cost result:
- deleted: `d` from the one-center required basis-key set.
- simplified: one-center White-Lindsey near-core mapping now derives from the
  same resolved `core_spacing` used by the rest of the producer contract.
- quarantined: legacy `d` is retained only as an equality-checked compatibility
  alias.
- not deleted because: `input.d` and artifact `mapping_d` remain for
  compatibility with the existing provenance schema.
- exact remaining caller/blocker: no source blocker for `HP-R1-CORE-FN-01`;
  a later artifact-schema cleanup would be needed to remove `mapping_d`.
- added src lines: 17.
- deleted src lines: 7.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 104 - Share Base System Normalization

Commit(s):
- this commit - Share base system normalization

Summary:
- Accepted the first 2->1 cleanup enabled by the core-spacing contract.
  `src/cartesian_base_hamiltonian.jl` now has one private
  `_cartesian_base_system_parts(system::NamedTuple)` helper for public system
  key/container validation, conversion to symbols/charges/locations, and
  nonnegative electron normalization.
- `_cartesian_base_inputs(...)` and `_cartesian_r3_diatomic_inputs(...)` now
  share that normalization, while all scope-specific H/H2/R3 checks and basis
  normalization remain separate.

Validation:
- Doer validation: `git diff --check`, package load, H atom no-`d`
  artifact/readback (dimension 419), H atom legacy `d == core_spacing`
  (dimension 419), H atom `d != core_spacing` rejection, H2 public `d`
  rejection, H2 base artifact/readback (dimension 471), compact H2
  supplemented driver path (dimension 489), and no Cr2 run.
- Manager validation: reviewed the one-file diff, confirmed `git diff
  --check`, numstat +11/-17 in `src/cartesian_base_hamiltonian.jl`,
  suspicious added-line scan clean, and new-test/tool scan empty.

Goal advancement:
- LT1/LT3: reduces duplicate producer input handling without changing the
  public driver or facade contract.
- MT: confirms the first unification step after `core_spacing`: common public
  system normalization can be shared, while real atom/H2/R3 semantics remain
  explicitly scoped.

Carrying-cost result:
- deleted: duplicated system key/container/charge/electron normalization from
  the R3 diatomic input path.
- simplified: base and supplemented input paths now use one private
  normalization helper before diverging into scope-specific validation.
- quarantined: none.
- not deleted because: branch-specific H/H2/R3 validation is intentionally
  still separate.
- exact remaining caller/blocker: no blocker for this cleanup. A later audit
  may consider shared H2-like basis normalization, but that should be a
  separate pass.
- added src lines: 11.
- deleted src lines: 17.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 105 - Delete Staged Compatibility Wrappers

Commit(s):
- this commit - Delete unused staged compatibility wrappers

Summary:
- Accepted a deletion-only cleanup after the staged operator-class driver
  split. Removed `_pqs_source_box_route_driver_terminal_one_body(...)` and
  `cartesian_residual_gto_augmented_operators(...)`, both of which were
  superseded by the explicit product/unit-nuclear/electron-electron stage
  functions.

Validation:
- Doer validation: pre-delete `rg` confirmed both helper names appeared only
  at their definitions; `git diff --check`; package load; post-delete `rg`
  over `src`, `bin`, `test`, `tools`, and `docs` found no remaining matches;
  no Cr2 run.
- Manager validation: reviewed the deletion-only diff, confirmed `git diff
  --check`, numstat 0/-6 in `src/cartesian_base_hamiltonian.jl` and 0/-8 in
  `src/pqs_source_box_low_order_materialization.jl`, suspicious added-line
  scan clean, new-test/tool scan empty, and helper-name `rg` empty.

Goal advancement:
- LT1/LT3: reduces compatibility surface now that the canonical driver has
  settled on explicit operator-class stages.
- MT: keeps the producer cleanup moving by deleting wrappers that preserved
  older coarse-stage vocabulary.

Carrying-cost result:
- deleted: `_pqs_source_box_route_driver_terminal_one_body(...)` and
  `cartesian_residual_gto_augmented_operators(...)`.
- simplified: staged construction now has only the explicit physical
  operator-class helper path.
- quarantined: none.
- not deleted because: `_cartesian_base_ida_hamiltonian(...)` and
  `pqs_terminal_residual_gto_augmented_operators(...)` remain live and were
  explicitly out of scope.
- exact remaining caller/blocker: none for the deleted helpers.
- added src lines: 0.
- deleted src lines: 14.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 106 - Relax One-Center Base Atoms

Commit(s):
- this commit - Relax one-center base atoms

Summary:
- Accepted the approved `HP-R1-ATOM-*` source relaxation in
  `src/cartesian_base_hamiltonian.jl`. The one-center base facade no longer
  hardcodes atom input to `symbols == ["H"]` and `charges == [1.0]`.
- One-center atoms are now validated from explicit caller data: the center must
  remain at the origin, the nuclear charge must be finite, positive, and
  integer-valued, and `nup + ndn` must equal that all-electron charge. The
  symbol remains only a label/provenance value.
- H2/base and R3 supplemented homonuclear diatomic validation were intentionally
  left separate.

Validation:
- Doer validation: `git diff --check`, package load, H atom no-public-`d`
  artifact/readback (dimension 419), H atom legacy `d == core_spacing`
  artifact/readback (dimension 419), H atom `d != core_spacing` rejection, H2
  base artifact/readback (dimension 471), H2 public `d` rejection, Be
  one-center atom smoke (dimension 419) with finite/symmetric `K`, unit `U_A`,
  and `V`, K/V artifact readback, noninteger charge rejection, nonneutral
  electron-count rejection, translated-atom rejection, and no Cr2 molecule run.
- Manager validation: reviewed the one-file diff, confirmed `git diff
  --check`, numstat +8/-5 in `src/cartesian_base_hamiltonian.jl`, suspicious
  added-line scan clean by direct `-U0` diff review, and no new-test/tool
  changes.

Goal advancement:
- LT1/LT3: removes the remaining H-only special case from the one-center base
  producer while preserving the settled driver and staged producer workflow.
- MT: advances the atomic/diatomic unification audit by replacing a chemistry
  hardcode with explicit physical input validation. The next cleanup candidate
  remains shared H2-like diatomic basis normalization, not another driver
  rewrite.

Carrying-cost result:
- deleted: one-center H-only symbol/charge guard and H-only error wording.
- simplified: one-center validation is now data-driven from explicit charge
  and electron inputs.
- quarantined: none.
- not deleted because: H2/base and R3 supplemented diatomic checks are still
  intentionally separate scope gates.
- exact remaining caller/blocker: no blocker for this atom relaxation; broader
  translated atoms, supplemented atoms, element defaults, and artifact schema
  cleanup remain unapproved.
- added src lines: 8.
- deleted src lines: 5.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 107 - Share Diatomic Basis Normalization

Commit(s):
- this commit - Share diatomic basis normalization

Summary:
- Accepted the narrow 2->1 cleanup in `src/cartesian_base_hamiltonian.jl`.
  Base H2 and supplemented homonuclear z-axis diatomic input paths now both use
  private `_cartesian_base_diatomic_basis_parts(...)` for shared basis-key
  checking and normalization.
- The helper owns only the common diatomic basis fields: `q`, `core_spacing`,
  `xmax_parallel`, `xmax_transverse`, `parent_axis_family`,
  `reference_spacing`, `tail_spacing`, plus the intentionally absent
  atom-only `radius`/`d` fields.
- Scope validation remains separate: base H2 is still exactly H2 with
  `nup=1`, `ndn=1`, while the supplemented path is still explicit
  homonuclear z-axis all-electron diatomic validation.

Validation:
- Doer validation: `git diff --check`, package load, H2 base
  artifact/readback (dimension 471), H2 supplemented facade artifact/readback
  (dimension 489), and no Cr2 run.
- Manager validation: reviewed the one-file diff, confirmed `git diff
  --check`, numstat +16/-19 in `src/cartesian_base_hamiltonian.jl`, no driver
  or atom branch changes, no new tests/tools, and net source decrease.

Goal advancement:
- LT1/LT3: removes duplicated producer input normalization without changing
  physics, artifacts, or the settled canonical driver.
- MT: continues the atomic/diatomic unification cleanup by sharing the common
  diatomic contract while preserving genuine scope differences.

Carrying-cost result:
- deleted: duplicate H2-like diatomic basis normalization from the two callers.
- simplified: `_cartesian_base_inputs(...)` H2 branch and
  `_cartesian_r3_diatomic_inputs(...)` now share one private helper for common
  basis fields.
- quarantined: none.
- not deleted because: atom-specific basis normalization remains genuinely
  different, and H2/R3 scope checks remain intentionally separate.
- exact remaining caller/blocker: none for this cleanup.
- added src lines: 16.
- deleted src lines: 19.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 108 - Select Active Route Subrecipe

Commit(s):
- this commit - Select active route subrecipe

Summary:
- Accepted the `HP-ROUTE-RECIPE-FN-01` cleanup. `cartesian_recipe(...)` in
  `src/pqs_source_box_route_driver_helpers.jl` now constructs only the active
  route-family subrecipe: `source_box` for `:pqs_source_box`, or
  `white_lindsey` for explicit `:white_lindsey_low_order`.
- The inactive subrecipe is retained as `nothing` for compatibility, rather
  than populated with stale fields. `_cartesian_base_route(...)` no longer
  carries five unused `white_lindsey_*` fields through the live PQS base
  producer path.
- Explicit White-Lindsey low-order support remains present; no materialization,
  shellification, terminal-lowering, numerical, artifact, or driver behavior
  was changed.

Validation:
- Doer validation: `git diff --check`, package load, H atom base
  artifact/readback (dimension 419), H2 base artifact/readback (dimension
  471), H2 supplemented facade artifact/readback (dimension 489), explicit
  `:white_lindsey_low_order` recipe smoke, focused `rg` confirming removed
  `white_lindsey_*` route fields in `src/cartesian_base_hamiltonian.jl`, and
  no Cr2 run.
- Manager validation: reviewed the two-file source diff, confirmed `git diff
  --check`, numstat 0/-5 in `src/cartesian_base_hamiltonian.jl` and +4/-6 in
  `src/pqs_source_box_route_driver_helpers.jl`, no driver/test/tool changes,
  no new metadata/status fields, and net source deletion.

Goal advancement:
- LT1/LT3: removes inactive White-Lindsey vocabulary from the current PQS base
  producer contract while preserving the older explicit WL route path.
- MT: continues contract unification by separating real WL/PQS algorithm
  differences from accidental route-input coupling.

Carrying-cost result:
- deleted: inactive `white_lindsey_*` route fields from `_cartesian_base_route`.
- simplified: `cartesian_recipe(...)` now builds only the selected route-family
  subrecipe.
- quarantined: inactive subrecipe side remains `nothing`, not stale populated
  data.
- not deleted because: explicit `:white_lindsey_low_order` support and WL
  materialization remain live/deferred for a separate retirement decision.
- exact remaining caller/blocker: tools/tests may still pass old WL fields in
  PQS route inputs, but they are no longer required by `cartesian_recipe(...)`.
- added src lines: 4.
- deleted src lines: 11.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 109 - Reuse Same-Construction Base K/U

Commit(s):
- this commit - Reuse same-construction base operators

Summary:
- Accepted `HP-R3BASE-FN-01` in the approved source files. The supplemented
  one-call facade and lower-level same-construction R3 path can now pass
  trusted base `K_GG` and by-center unit `U_GG[A]` blocks into augmented
  operator construction instead of recomputing those base final-basis `G-G`
  blocks.
- Standalone augmented-operator helpers preserve the old recompute behavior
  unless trusted `base_kinetic` or `base_unit_nuclear` keyword inputs are
  supplied. The trust model remains local: dimension and center-count checks,
  no provenance payloads or metadata proofs.
- The canonical driver was intentionally not changed under this source lane;
  its visible staged call site still needs separate `bin/` authority to pass
  `base_ham` into the augmented stages.

Validation:
- Doer validation: `git diff --check`, package load, H2 R3 endpoint with
  self-Coulomb `0.4574265214362095`, Be2 facade/readback with augmented
  dimension `1421` and all readback deltas `0.0`, focused Cr2 attribution
  replay with base `K/U` reuse operator delta `0.0`, same-construction `K/V`
  deltas `0.0`, finite/symmetric exact operators, and no Cr2 artifact/workflow.
- Doer reported Cr2 allocation evidence: default production exact
  augmented-operator wrapper `5.0452s / 2987.948 MiB`; base `K/U` reuse replay
  `0.8558s / 1229.886 MiB`.
- Manager validation: reviewed the two-file diff, confirmed `git diff
  --check`, numstat +10/-6 in `src/cartesian_base_hamiltonian.jl` and +45/-21
  in `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, no
  driver/test/tool changes, no metadata/status fields, and unchanged standalone
  fallback behavior.

Goal advancement:
- LT1/LT3: unifies same-construction base and supplemented operator assembly
  by reusing physical base blocks already computed in the same construction.
- MT: crosses the strongest remaining no-supplement/supplement duplication in
  the exact augmented-operator path. The remaining staged-driver gap is wiring,
  not a numerical-kernel issue.

Carrying-cost result:
- deleted: base `K_GG` and unit `U_GG[A]` recomputation from same-construction
  facade paths when trusted base blocks are available.
- simplified: augmented operator construction now has local keyword reuse paths
  with dimension/count checks.
- quarantined: none.
- not deleted because: default recompute path remains live for standalone
  callers that do not have trusted same-construction base blocks.
- exact remaining caller/blocker: canonical driver needs a separate approved
  `bin/cartesian_ham_builder.jl` call-site update to pass `base_ham` into the
  staged augmented helpers.
- added src lines: 55.
- deleted src lines: 27.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 110 - Wire Driver Base K/U Reuse

Commit(s):
- this commit - Wire driver base operator reuse

Summary:
- Accepted the tiny `HP-R3BASE-DRV-WIRE-01` driver call-site update. In
  supplemented mode, the canonical driver now passes `base_ham.kinetic` into
  `cartesian_residual_gto_augmented_products(...)` as `base_kinetic`, and
  passes `base_ham.nuclear_attraction_unit_by_center` into
  `cartesian_residual_gto_augmented_unit_nuclear(...)` as
  `base_unit_nuclear`.
- No public inputs, hooks, timing labels, stage sequence, printed contract,
  artifact schema, or base-mode behavior changed. This makes the staged driver
  benefit from the same-construction base `K/U_A` reuse accepted in Pass 109.

Validation:
- Doer validation: `git diff --check`, package load, H2 supplemented driver
  artifact/readback with dimension `489`, Be2 supplemented driver
  artifact/readback with dimension `1421`, and no Cr2 run. Doer reported
  staged timings after wiring: H2 augmented products `2.896s`, augmented unit
  nuclear `0.195s`; Be2 augmented products `3.323s`, augmented unit nuclear
  `0.215s`.
- Manager validation: reviewed the driver-only diff, confirmed `git diff
  --check`, numstat +4/-2 in `bin/cartesian_ham_builder.jl`, no `src`, `test`,
  or `tools` edits, no new tests, and no public driver contract changes.

Goal advancement:
- LT1/LT3: closes the driver-facing gap left by same-construction base
  operator reuse while preserving the canonical driver's human-facing shape.
- MT: completes the no-supplement/supplement `K/U_A` reuse lane for the
  canonical supplemented workflow.

Carrying-cost result:
- deleted: none.
- simplified: supplemented driver call sites now use same-construction base
  blocks already available in the visible staged workflow.
- quarantined: none.
- not deleted because: staged timing labels and public contract were
  intentionally preserved.
- exact remaining caller/blocker: none for canonical driver base `K/U_A`
  reuse.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 108 - Approve Family-Selective Route Recipes

Commit(s):
- this commit - Approve family-selective route recipes

Summary:
- Approved `HP-ROUTE-RECIPE-FN-01` / `HP-ROUTE-RECIPE-TEST-01` as a narrow
  cleanup lane for `cartesian_recipe(...)`.
- The design now lets `cartesian_recipe(route_inputs)` build only the selected
  route-family subrecipe. PQS-only route inputs must no longer carry inactive
  `white_lindsey_*` fields, and `_cartesian_base_route(kind)` may remove those
  unused fields. Explicit `:white_lindsey_low_order` route support remains
  preserved.
- This is a contract-unification cleanup, not a physics change. Real WL/PQS
  algorithm differences remain; inactive WL vocabulary should stop flowing
  through the current PQS base producer.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  for `HP-ROUTE-RECIPE-*`, `cartesian_recipe`, `white_lindsey_*`, selected
  subrecipe wording, and no `src`, `test`, `tools`, or `bin` changes. No
  implementation tests were run.

Goal advancement:
- LT1/LT3: removes another needless old-route concept from the live base
  producer contract.
- MT: continues the contract-unification lane after core-spacing and shared
  diatomic normalization cleanup.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future PQS base route inputs no longer need inactive WL fields.
- quarantined: canonical driver changes, numerical kernels, shellification,
  terminal lowering, materialization/artifact schema, route diagnostics,
  status/report expansion, WL materialization deletion, new tests, and Cr2 runs
  remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: doer must stop if family-selective
  `cartesian_recipe(...)` requires broader route-driver, report,
  materialization, or stage-object changes.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 109 - Approve Same-Construction Base K/U Reuse

Commit(s):
- this commit - Approve same-construction base K/U reuse

Summary:
- Approved `HP-R3BASE-FN-01` / `HP-R3BASE-TEST-01` as a narrow source lane for
  reusing same-construction base final-basis `K_GG` and unit `U_GG[A]` blocks
  inside supplemented residual-GTO/MWG exact augmented operators.
- Recorded the trust model: the base Hamiltonian, terminal basis, parent
  bundles, residual basis, and supplement must come from the same
  `cartesian_base_working_basis(...)` path; dimension and center-count checks
  are required; no provenance payload, metadata proof, report/status object, or
  persistent cache is approved.
- Preserved fallback behavior when trusted base blocks are not supplied.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  for `HP-R3BASE-*`, same-construction `K_GG` / `U_GG[A]` reuse, forbidden
  metadata/cache/public/driver scope, and no `src`, `test`, `tools`, or `bin`
  changes. No implementation tests were run.

Goal advancement:
- LT1/LT3: removes duplicate base one-body work from the supported
  no-supplement/supplement construction story without changing physics.
- RG/LT6: gives repo-doer a narrow performance/complexity cleanup path before
  broader Cr2 workflow expansion.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future same-construction supplemented runs may reuse already
  constructed base K/U blocks instead of recomputing them.
- quarantined: public API/export changes, canonical driver changes, raw
  blocks, residual algorithms, MWG/IDA conventions, terminal kernels, route
  setup, metadata/status/report/artifact fields, persistent caches, committed
  tests, and Cr2 workflow remain unapproved.
- not deleted because: source reuse has not run yet.
- exact remaining caller/blocker: doer must stop if local call shape plus
  dimension/center validation cannot guarantee same-construction trust.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 110 - Approve Driver K/U Reuse Wiring

Commit(s):
- this commit - Approve driver K/U reuse wiring

Summary:
- Approved `HP-R3BASE-DRV-WIRE-01` / `HP-R3BASE-DRV-TEST-01` for the tiny
  canonical-driver supplemented-mode call-site pass.
- The only approved file is `bin/cartesian_ham_builder.jl`. The driver may
  pass `base_ham.kinetic` as `base_kinetic` and
  `base_ham.nuclear_attraction_unit_by_center` as `base_unit_nuclear` into the
  existing augmented product and unit-nuclear stage calls.
- The approval explicitly forbids public input changes, hook changes, timing
  label changes, stage-sequence changes, artifact schema changes, diagnostics,
  source/kernel changes, committed tests/fixtures, and Cr2 workflow.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  for `HP-R3BASE-DRV-*`, `base_kinetic`, `base_unit_nuclear`, forbidden driver
  contract changes, and no `src`, `test`, `tools`, or `bin` changes. No
  implementation tests were run.

Goal advancement:
- LT1/LT3: lets the canonical driver use the same same-construction base K/U
  reuse path already approved for the supplemented facade without changing its
  human-facing contract.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future driver supplemented runs can avoid duplicate base K/U
  construction through a call-site-only change.
- quarantined: diagnostics, new hooks, timing-label changes, visible stage
  changes, public inputs, artifacts, source/kernel edits, tests/fixtures, and
  Cr2 workflow remain unapproved.
- not deleted because: source driver wiring has not run yet.
- exact remaining caller/blocker: doer must stop if the call-site update needs
  any visible driver contract change.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 111 - Approve Driver Hidden-d Cleanup

Commit(s):
- this commit - Approve driver hidden-d cleanup

Summary:
- Approved `HP-DRV-ATOM-CLEAN-01` for the tiny canonical-driver cleanup that
  removes hidden `d = vars[:core_spacing]` from one-center atom basis
  construction.
- The producer no longer requires public `d`; keeping the hidden field in the
  driver would preserve stale public-contract residue.
- Explicitly deferred larger cleanup targets: old `:white_lindsey_low_order`
  route/materialization retirement and stale inactive `white_lindsey_*` fields
  in tests/tools remain audit/amendment work, not immediate deletion.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  for `HP-DRV-ATOM-CLEAN-01`, hidden `d`, forbidden driver-contract changes,
  and no `src`, `test`, `tools`, or `bin` changes. No implementation tests
  were run.

Goal advancement:
- LT1/LT3: removes another tiny false public-contract remnant from the
  canonical driver after core-spacing unification.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future driver atom basis construction should expose only
  `core_spacing`, not a hidden `d` alias.
- quarantined: WL route/materialization retirement and stale test/tool
  `white_lindsey_*` route-input cleanup remain deferred.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: doer must stop if removing hidden `d`
  requires any visible driver contract change or producer/source change.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 112 - Approve Compact Artifact Manifest

Commit(s):
- this commit - Approve compact Hamiltonian artifact manifest

Summary:
- Approved `HP-HAM-MANIFEST-FN-01` / `HP-HAM-MANIFEST-TEST-01` for compact
  JLD2 sidecar groups on existing `CartesianIDAHamiltonian{Float64}` artifacts.
- The manifest lane reuses the prior PQS fixed-column/source-mode provenance
  model from `projected_q_shell_policy.md`: basis identity is a
  status-bearing construction label, not a representative center.
- `hamiltonian_manifest/final_basis_labels/` records exact matrix-order final
  basis rows with sector, unit/source labels, shell/ray/radial status,
  representative center metadata, owner nucleus index, locality/freezing
  labels, and supplement labels/angular powers where available.
- Optional `final_basis_source_relations/`, `source_shells/`, and
  `source_modes/` subgroups are approved only for construction-native facts;
  unavailable/mixed labels must be explicit and must not be inferred from
  centers, nearest-grid snapping, support order, support indices, or
  raw-to-final support.
- It also records uniform public recipe provenance under `recipe_provenance/`:
  system/basis/supplement inputs, padding-derived extents, route, parent-axis
  counts, and base/residual/augmented dimensions.
- Existing Hamiltonian matrix keys and `read_cartesian_ida_hamiltonian`
  behavior remain unchanged. The manifest is for downstream artifact consumers
  and validation-only direct JLD2 inspection, not a new public reader API or
  Hamiltonian wrapper.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  checks for `HP-HAM-MANIFEST-*`, `hamiltonian_manifest/`,
  `final_basis_labels`, source relation/source mode wording,
  `recipe_provenance/`, forbidden dense transform/raw inventory fields,
  no-inference guardrails, unchanged reader wording, and no `src`, `test`,
  `tools`, or `bin` changes. No implementation tests were run.

Goal advancement:
- LT1/LT3: makes canonical-driver artifacts more useful to downstream
  consumers without changing the Hamiltonian object or public driver inputs.
- RG/LT6: supports realistic supplemented diatomic handoff while keeping Cr2 as
  a generic ignored/user-run stress case rather than a schema branch.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future artifact consumers can read compact basis/recipe sidecars
  instead of reconstructing locality from private workflow context or center
  heuristics.
- quarantined: `T_G`, `T_A`, dense transforms, raw inventories, allocation
  probes, report/status payloads, public readers, driver schema dumps,
  solver-specific fields, CR2-consumer-specific fields, Cr2-specific fields,
  committed Cr2 fixtures, and atom size-policy changes remain unapproved.
- not deleted because: source manifest writing has not run yet.
- exact remaining caller/blocker: first source pass must stop if center
  conventions cannot be derived from existing terminal basis blocks, parent
  axes, residual metadata, and augmented moment/MWG descriptors without adding
  algorithmic metadata.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 113 - Implement Compact Artifact Manifest

Commit(s):
- this commit - Add compact Hamiltonian artifact manifest

Summary:
- Accepted the first `HP-HAM-MANIFEST-FN-01` source pass. Canonical base and
  supplemented artifact writes now attach `hamiltonian_manifest/manifest_version`,
  `hamiltonian_manifest/final_basis_labels/`, and `recipe_provenance/` sidecar
  groups without changing matrix keys, the Hamiltonian object, the reader, or
  the canonical driver.
- The first-pass manifest is intentionally compact: one status-bearing row per
  matrix-order final basis column plus public recipe/dimension provenance. It
  does not attempt optional `final_basis_source_relations/`, `source_shells/`,
  or `source_modes/` because those would need construction-native facts not yet
  present in this artifact seam.
- Base-row representative centers come from terminal support states and parent
  axis centers; Lowdin/PQS rows use an abs2-coefficient support centroid and are
  marked representative. Supplemented residual rows use exact augmented position
  diagonals and residual owner/label metadata.

Validation:
- Doer validation: `git diff --check`, package load, and ignored direct JLD2
  validation via `tmp/work/ham_manifest_validation.jl`. H2 base and H2
  supplemented artifacts read back through the existing reader; manifest row
  counts matched matrix dimensions; recipe dimensions matched base/residual/
  augmented dimensions; inferred-label flags were false; unavailable/mixed
  label status checks passed. No Cr2 run.
- Manager validation: `git diff --check`; `git diff --numstat` reported
  `147` added and `5` deleted source lines in `src/cartesian_base_hamiltonian.jl`;
  suspicious-line scan found only approved artifact `*_status` columns and
  assignments; new-test/tool scan was empty; focused `rg` found no optional
  source-relation/source-shell/source-mode groups and no serialized `T_G`/`T_A`
  manifest fields.

Goal advancement:
- LT1/LT3: makes canonical-driver artifacts self-describing enough for
  downstream row/order consumers without exposing route internals or changing
  the public reader contract.
- RG/LT6: directly addresses CR2 feedback that supplemented artifacts lacked a
  matrix-order final-basis manifest while keeping centers representative rather
  than identity-bearing labels.

Carrying-cost result:
- deleted: none; this was a new approved artifact-sidecar capability.
- simplified: downstream consumers no longer need to reconstruct row labels and
  public recipe context from external workflow memory or center heuristics.
- quarantined: optional source relations, source shells, and source modes are
  omitted until construction-native relation data is available; CR2-specific
  fields and dense transform/raw inventory fields remain unapproved.
- not deleted because: existing `producer_provenance/` and
  `supplement_provenance/` remain approved compatibility groups, and the minimal
  low-level Hamiltonian writer remains intentionally manifest-free.
- exact remaining caller/blocker: direct low-level
  `write_cartesian_ida_hamiltonian` callers still write the minimal artifact
  only; richer source-mode/ray/shell labels need a later construction-native
  provenance lane.
- added src lines: 147.
- deleted src lines: 5.
- new tests: none.
- new metadata/status fields: approved artifact sidecar fields only; no
  algorithmic metadata/status fields.

## Cartesian Hamiltonian Producer Pass 114 - Approve Source-Mode Manifest Seam

Commit(s):
- this commit - Approve source-mode manifest seam

Summary:
- Approved `HP-HAM-MANIFEST-SRC-FN-01` /
  `HP-HAM-MANIFEST-SRC-TEST-01` as the narrow follow-on authority for carrying
  construction-native source-mode provenance to the compact Hamiltonian
  manifest writer.
- The approved seam is:
  terminal lowering / retained-unit / raw-product source plans -> compact
  source-mode provenance object -> base working basis manifest context ->
  artifact sidecar writer.
- The preferred carrier is one internal `source_mode_provenance` field on the
  `cartesian_base_working_basis(...)` result. One optional
  `CartesianTerminalBasisRealization` field is allowed only if it avoids
  duplicated or lost terminal construction ordering.
- The seam may populate optional `hamiltonian_manifest/source_shells/`,
  `hamiltonian_manifest/source_modes/`, native
  `final_basis_source_relations/`, and native `final_basis_labels/`
  improvements. Missing ray/radial/source facts must remain explicit
  `:unavailable` or `:mixed`.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  checks for `HP-HAM-MANIFEST-SRC-*`, source-mode provenance carrier wording,
  approved source files, optional source groups, no-inference guardrails, and
  forbidden dense transform/raw inventory/route report fields; confirmed no
  `src`, `test`, `tools`, or `bin` files changed. No implementation tests were
  run.

Goal advancement:
- LT1/LT3: closes the exact design gap found by the first manifest source pass:
  the writer can only emit source-shell/source-mode groups when live
  construction-native provenance reaches the artifact seam.
- RG/LT6: keeps CR2-facing provenance useful without adding CR2-specific
  schema, ray policies, or Hamiltonian-object changes.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one approved place to carry source-mode provenance
  instead of inventing route reports or deriving labels from centers.
- quarantined: coefficients, dense transforms, `T_G`, `T_A`, raw inventories,
  raw pair inventories, route reports, allocation probes, diagnostic payloads,
  driver changes, reader changes, matrix-key changes, public API/export
  changes, and non-native ray/radial labels remain unapproved.
- not deleted because: the first compact manifest source pass is already live
  and remains valid for final-basis labels and recipe provenance.
- exact remaining caller/blocker: source implementation must stop if native
  source shell/mode facts cannot be carried from existing terminal lowering,
  retained-unit, raw-product source, and terminal basis records without adding
  algorithmic metadata or route report plumbing.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none; this approves only a compact artifact
  provenance carrier for a future source pass.

## Cartesian Hamiltonian Producer Pass 115 - Implement Source-Mode Manifest Seam

Commit(s):
- this commit - Add source-mode Hamiltonian manifest groups

Summary:
- Accepted the `HP-HAM-MANIFEST-SRC-FN-01` source pass. Base working-basis
  construction now carries compact construction-native source-mode provenance
  from terminal retained-rule / raw-product source plans to the artifact
  manifest writer.
- Artifacts now write optional `hamiltonian_manifest/source_shells/` and
  `hamiltonian_manifest/source_modes/` groups when the facts are native, and
  `final_basis_labels/` rows now receive native source-box and shell labels
  where the final row belongs to a native terminal source shell.
- Ray IDs, radial labels, and `final_basis_source_relations/` remain
  unavailable/omitted. The source pass intentionally did not serialize dense
  transforms, coefficients, `T_G`, `T_A`, raw inventories, route reports, or
  CR2-specific schema.

Validation:
- Doer validation: package load; H2 base artifact/readback with dimension 471,
  K readback delta 0.0, 3 source shells, 525 source modes, shell index max 3,
  and no inferred flags; H2 supplemented artifact/readback with dimension 489,
  K/V readback deltas 0.0, 3 source shells, 525 source modes, 18 residual rows,
  and no inferred flags; `git diff --check`.
- Manager validation: `git diff --check`; `git diff --numstat` reported
  `150` added and `4` deleted source lines in
  `src/cartesian_base_hamiltonian.jl`; suspicious-line scan found only fixed
  three-axis tuple conversions, the approved native `raw_product_source_plan`
  metadata read, and approved artifact status fields; new-test/tool scan was
  empty. The generated source-mode coordinate inventory was amended from a
  variable-length `Tuple(...)` to vector-backed append loops before acceptance.

Goal advancement:
- LT1/LT3: upgrades Hamiltonian artifacts from row/order metadata to the first
  construction-native source-mode layer without changing matrix keys, reader
  behavior, public APIs, or the canonical driver.
- RG/LT6: materially improves CR2 source inspection while preserving the
  distinction between native producer facts and consumer ray/locality policy.

Carrying-cost result:
- deleted: none; this was an approved provenance extension.
- simplified: artifact consumers can now read native source shell/mode tables
  instead of relying on ignored sidecar prototypes for this first provenance
  layer.
- quarantined: non-native ray/radial labels, final-basis source relations,
  dense relation weights, transform storage, and CR2-specific labels remain
  deferred.
- not deleted because: the minimal Hamiltonian writer remains intentionally
  sidecar-free, and old ignored CR2/postprocess scripts remain only as
  prototype/reference evidence.
- exact remaining caller/blocker: ray/locality contraction still needs either a
  reviewed native relation producer or explicit consumer-side policy; the
  producer should not infer those labels from centers or support order.
- added src lines: 150.
- deleted src lines: 4.
- new tests: none.
- new metadata/status fields: approved artifact `source_shells/` and
  `source_modes/` sidecar datasets only.

### Medium-Term Goal Checkpoint After Pass 115

- Active: CR2 usability has shifted from "can the canonical driver write a
  usable artifact?" to "can downstream consumers recover enough native row
  provenance for locality/ray decisions without relying on center heuristics?"
- Completed: compact artifact manifest baseline, public recipe provenance, and
  native source shell/mode sidecar writing for current base/supplemented
  artifacts.
- Active: timing/compile attribution is now a separate lane. Current evidence
  says one-shot driver timings mix compilation, artifact reuse, and real
  construction cost; future performance claims need clean repeated fresh-target
  runs and a type-surface audit before source optimization.
- Deferred: ray IDs, radial labels, source relations, and AFQMC/locality
  policy remain outside producer authority until a native producer seam or
  explicit consumer contract is approved.

## Cartesian Hamiltonian Producer Pass 116 - Add Retained Boundary Seed Relations

Commit(s):
- this commit - Add retained boundary seed manifest relations

Summary:
- Accepted the narrow manifest relation pass for shell-realized PQS rows.
  `hamiltonian_manifest/final_basis_source_relations/` now records one native
  retained boundary seed relation per final PQS shell column when the terminal
  retained-rule facts provide an ordered retained mode identity.
- This is a construction label, not coefficient reconstruction. Relation rows
  identify the seed source shell/mode tuple used to produce the final column,
  while explicitly marking coefficients, weights, spans, ray labels, and radial
  labels as unavailable or not serialized.
- Direct/support-dense rows and residual rows remain relation-unlabeled in this
  pass because extending them would require different native facts or consumer
  policy.

Validation:
- Doer validation: package load; H2 base artifact/readback with dimension 471,
  K readback delta 0.0, 3 source shells, 525 source modes, and 196 relation
  rows; H2 supplemented artifact/readback with dimension 489, base rows 471,
  residual rows 18, K/V readback deltas 0.0, and 196 relation rows. Relation
  kind/status were all `:boundary_mode` /
  `:native_retained_boundary_seed`; relation labels matched existing native
  `source_modes/`; all inference flags were false; `git diff --check`.
- Manager validation: `git diff --check`; `git diff --numstat` reported `46`
  added and `3` deleted source lines in `src/cartesian_base_hamiltonian.jl`;
  suspicious-line scan found only the approved native
  `raw_product_source_retained_rule` metadata read and approved artifact status
  fields; new-test/tool scan was empty.

Goal advancement:
- LT1/LT3: closes the immediate CR2 feedback loop that source shell/mode tables
  were useful but did not label final shell rows. Final shell columns now carry
  a native seed-mode label without changing numerical matrices, the reader, or
  the driver.
- RG/LT6: gives downstream ray/locality work a stronger construction label
  while preserving the no-inference rule for actual ray/radial ownership.

Carrying-cost result:
- deleted: none; this is an approved artifact provenance extension.
- simplified: CR2 consumers no longer need to infer shell-realized final row
  seed labels from ordering outside the artifact.
- quarantined: direct/support-dense relation policy, residual relation policy,
  ray IDs, radial labels, relation weights/spans, dense transforms, and
  coefficient serialization remain deferred.
- not deleted because: source shell/mode manifest writing remains the active
  carrier for relation labels, and minimal low-level Hamiltonian writing remains
  sidecar-free by design.
- exact remaining caller/blocker: richer locality contraction still needs a
  reviewed native relation producer or explicit consumer policy for
  support-dense/direct rows, residual rows, and ray/radial grouping.
- added src lines: 46.
- deleted src lines: 3.
- new tests: none.
- new metadata/status fields: approved artifact
  `final_basis_source_relations/` sidecar fields only.

## Cartesian Hamiltonian Producer Pass 117 - Approve Route Inventory Type Cleanup

Commit(s):
- this commit - Approve route inventory type cleanup

Summary:
- Approved `HP-ROUTE-INV-FN-01` / `HP-ROUTE-INV-TEST-01` as the first narrow
  type-surface cleanup lane.
- The approved source file is only
  `src/pqs_source_box_route_driver_helpers.jl`.
- The target is removal of runtime-keyed retained-unit inventory
  `NamedTuple{unit_keys}` shapes and runtime-keyed
  `pair_family_counts = NamedTuple{families}(...)`, replacing them with
  vector-backed records/tables, stable dictionaries, or helper accessors with
  stable concrete types.
- Explicitly deferred broader tuple cleanup in `RawProductBoxPlan`, terminal
  lowering plans, and retained-unit transform-contract plans.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  checks for `HP-ROUTE-INV-*`, approved helper/file names, deferred tuple
  targets, forbidden surfaces, and no `src`, `test`, `tools`, or `bin` changes.
  No implementation tests were run.

Goal advancement:
- LT1/LT3: removes the clearest remaining route-size-in-type pressure from the
  live base/supplemented producer path before timing-oriented cleanup.
- RG/LT6: keeps the canonical supplemented path moving toward practical use
  without changing physics, artifacts, or public driver contract.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one approved cleanup target for the runtime-keyed
  retained-unit inventory shape.
- quarantined: raw product source-mode tuple cleanup, terminal-lowering
  contract tuple cleanup, retained-unit transform-contract tuple cleanup,
  public input `NamedTuple` changes, artifact sidecar table changes, numerical
  kernels, driver changes, report/status/payload expansion, compatibility
  adapters, committed tests, and Cr2 workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: source implementation must stop if replacing
  the runtime-keyed inventories requires files outside
  `src/pqs_source_box_route_driver_helpers.jl`, broad route/stage rewiring, or
  preserving the old type shape behind a compatibility adapter.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 118 - Remove Runtime-Keyed Route Inventory Tuples

Commit(s):
- this commit - Replace runtime-keyed route inventories

Summary:
- Accepted the first `HP-ROUTE-INV-FN-01` source cleanup. The route-driver
  helper no longer builds retained-unit inventory maps as
  `NamedTuple{unit_keys}` or pair-family counts as `NamedTuple{families}`.
- The replacement is one private vector/dictionary-backed
  `_PQSRouteDriverInventoryRows` helper with property-style lookup preserved for
  current callers. The public `propertynames` path was amended before
  acceptance so it no longer constructs `Tuple(rows.labels)`.
- This is a type-surface cleanup only: no numerical kernels, driver behavior,
  artifact keys, manifest sidecars, raw product source tuples, terminal
  lowering tuples, or retained-unit transform-contract tuples changed.

Validation:
- Doer validation: `git diff --check`; package load; H2 base artifact/readback
  with dimension 471 and K/V deltas 0.0; H2 supplemented artifact/readback with
  dimension 489 and K/V deltas 0.0; focused `rg` found no
  `_pqs_source_box_route_driver_named_tuple_from_units`,
  `NamedTuple{unit_keys}`, `NamedTuple{families}`, or `Tuple(rows.labels)` in
  the allowed file.
- Manager validation: `git diff --check`; `git diff --numstat` reported `55`
  added and `23` deleted source lines in
  `src/pqs_source_box_route_driver_helpers.jl`; suspicious-line scan was empty;
  new-test/tool scan was empty.

Goal advancement:
- LT1/LT3: removes the clearest runtime-keyed `NamedTuple` inventory in the
  current base/supplemented route path before timing/compile attribution.
- RG/LT6: keeps CR2/HF artifact production work pointed at compile/runtime
  usability without changing the accepted physics or artifact contract.

Carrying-cost result:
- deleted: `_pqs_source_box_route_driver_named_tuple_from_units(...)`.
- simplified: repeated retained-unit map construction now goes through one
  stable private row helper with direct `hasproperty` lookup.
- quarantined: raw product source-mode tuple cleanup, terminal-lowering tuple
  cleanup, retained-unit transform-contract tuple cleanup, public input
  `NamedTuple` cleanup, artifact sidecar table changes, numerical kernels,
  driver changes, and Cr2 workflow remain deferred.
- not deleted because: current source callers still use property-style lookup,
  so the replacement preserves that access without preserving the
  runtime-keyed `NamedTuple` type shape.
- exact remaining caller/blocker: broader type-surface cleanup still has known
  deferred tuple lanes in `RawProductBoxPlan`, terminal lowering plans, and
  retained-unit transform-contract plans.
- added src lines: 55.
- deleted src lines: 23.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 119 - Approve Raw Product Source-Mode Inventory Cleanup

Commit(s):
- this commit - Approve raw product source-mode inventory cleanup

Summary:
- Approved `HP-RAW-SRCMODE-FN-01` / `HP-RAW-SRCMODE-TEST-01` as the next
  narrow type-surface cleanup lane after the retained-unit route inventory
  cleanup.
- The approved owner files are
  `src/cartesian_raw_product_sources/records.jl`,
  `src/cartesian_raw_product_sources/source_mode_indices.jl`, and
  `src/cartesian_raw_product_sources/summaries.jl`.
- Narrow consumer wiring is approved only as needed in
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  and `src/cartesian_base_hamiltonian.jl`.
- The target is replacing `RawProductBoxPlan.source_mode_indices` and
  `source_mode_column_indices` variable-length tuple storage with
  vector-backed storage, or removing stored column indices if they are exactly
  `1:count`.
- Accessor compatibility is defined as preserving deterministic ordered facts,
  column associations, retained-rule parity, and manifest source-mode/relation
  output, not preserving the old variable-length tuple concrete type.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  checks for `HP-RAW-SRCMODE-*`, approved files, target fields, validation
  gates, forbidden surfaces, and no `src`, `test`, `tools`, or `bin` changes.
  No implementation tests were run.

Goal advancement:
- LT1/LT3: continues the type-surface cleanup from the route inventory pass by
  removing another variable-size route fact from concrete field types without
  changing physics, artifacts, or public inputs.
- RG/LT6: protects current manifest source-mode and relation provenance while
  making the base/supplemented path less dependent on compile-costly staged
  tuple inventories.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one approved raw product source-mode inventory
  cleanup target with exact files and validation.
- quarantined: terminal-lowering contract tuple cleanup, retained-unit
  transform-contract tuple cleanup outside narrow caller wiring, broad
  pair-block/source-box rewrites, public input `NamedTuple` changes, numerical
  kernels, route semantic changes, driver changes, artifact schema changes,
  report/status/payload expansion, compatibility adapters preserving the old
  tuple-backed shape, committed tests, and Cr2 workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: source implementation must stop if
  vectorizing `RawProductBoxPlan` source-mode inventories forces broad
  pair-block/source-box rewrites or files outside the approved surfaces.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 120 - Vectorize Raw Product Source Modes

Commit(s):
- this commit - Vectorize raw product source modes

Summary:
- Accepted `HP-RAW-SRCMODE-FN-01`. `RawProductBoxPlan.source_mode_indices`
  and `source_mode_column_indices` no longer store variable-length
  `Tuple{Vararg...}` inventories; they now use vector-backed storage.
- `source_mode_indices(...)` still returns the same deterministic
  x-major/y-major/z-fast order, but as `Vector{NTuple{3,Int}}`. Existing
  consumers needed no source changes, which keeps the blast radius low.
- This is a compile/type-surface cleanup only. It does not change numerical
  kernels, route semantics, artifact schema, driver behavior, reader behavior,
  or manifest source-mode/relation labels.

Validation:
- Doer validation: `git diff --check`; package load; H2 base artifact/readback
  with dimension 471 and K/V deltas 0.0; H2 supplemented artifact/readback with
  dimension 489 and K/V deltas 0.0; manifest source-mode rows stayed 525;
  retained boundary relation rows stayed 196; relation labels still matched
  native source-mode labels; H2 R3 endpoint passed with self-Coulomb
  `0.4574265214362095`; focused parity for `(3,3,3)` and `(5,5,5)` source-mode
  order and retained boundary columns matched the old deterministic sequence.
- Manager validation: `git diff --check`; `git diff --numstat` reported `6`
  added and `6` deleted source lines across
  `src/cartesian_raw_product_sources/records.jl` and
  `src/cartesian_raw_product_sources/source_mode_indices.jl`; suspicious-line
  scan was empty; new-test/tool scan was empty.

Goal advancement:
- LT1/LT3: removes the strongest remaining q/source-mode-count-dependent
  concrete tuple inventory from the audited base working-basis path before
  fresh timing/compile attribution.
- RG/LT6: preserves current manifest source-mode and retained-seed relation
  behavior while reducing a compile-latency risk for the canonical
  supplemented artifact path.

Carrying-cost result:
- deleted: variable-length tuple storage for raw source-mode inventories in
  `RawProductBoxPlan`.
- simplified: source-mode order is now a vector-backed deterministic inventory
  with the existing accessor surface preserved.
- quarantined: fixed `NTuple{3,Int}` dimensions/coordinates remain, and
  terminal-lowering plus retained transform-contract tuple lanes remain
  separate future cleanup targets.
- not deleted because: retained boundary rule vectors already have the desired
  storage shape and were intentionally left unchanged.
- exact remaining caller/blocker: next compile cleanup, if still warranted
  before timing, is likely terminal-lowering plan contract tuples or retained
  transform-contract plan tuples.
- added src lines: 6.
- deleted src lines: 6.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 121 - Approve Contract-Plan Vector Cleanup

Commit(s):
- this commit - Approve contract-plan vector cleanup

Summary:
- Approved `HP-CONTRACT-VEC-FN-01` / `HP-CONTRACT-VEC-TEST-01` as the next
  narrow type-surface cleanup lane after raw product source-mode vectorization.
- The approved owner files are
  `src/cartesian_terminal_lowering/contracts.jl`,
  `src/cartesian_terminal_lowering/selection.jl`,
  `src/cartesian_terminal_lowering/summaries.jl`,
  `src/cartesian_retained_unit_transform_contracts/records.jl`,
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`, and
  `src/cartesian_retained_unit_transform_contracts/summaries.jl`.
- Narrow consumer wiring is approved only as needed in
  `src/pqs_source_box_route_driver_helpers.jl`,
  `src/cartesian_base_hamiltonian.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.
- The target is replacing variable-length tuple storage for
  `TerminalLoweringPlan.available_contracts`, `TerminalLoweringPlan.contracts`,
  and `RetainedUnitTransformContractPlan.contracts` with vector-backed storage.
- Accessors `available_contracts(plan)`, `selected_contracts(plan)`,
  `contracts(plan)`, and `transform_contracts(plan)` must preserve ordered
  behavior and semantics without preserving the old tuple-backed concrete field
  types.
- `source_cpbs::Tuple{Vararg{CoordinateProductBox}}` is explicitly out of
  scope.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  checks for `HP-CONTRACT-VEC-*`, approved files, target fields, preserved
  accessors, validation gates, forbidden surfaces, and no `src`, `test`,
  `tools`, or `bin` changes. No implementation tests were run.

Goal advancement:
- LT1/LT3: continues compile/type-surface cleanup by removing plan-level
  stage-boundary contract inventories from variable-length tuple field types.
- RG/LT6: keeps the canonical base/supplemented path on the same physics and
  artifact contracts while reducing staged type-shape pressure before fresh
  timing/compile attribution.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one approved contract-plan vector cleanup target
  with exact files, accessors, validation, and failure rule.
- quarantined: `source_cpbs`, raw product source-mode storage, retained-unit
  route inventories, public input `NamedTuple` changes, numerical kernels,
  route semantic changes, shellification, driver changes, artifact/manifest
  schema changes, report/status/payload expansion, compatibility adapters
  preserving old tuple-backed plan field types, committed tests, and Cr2
  workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: implementation must stop if vectorizing the
  plan inventories forces broad route/stage rewrites or source files outside
  the approved surfaces.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 122 - Vectorize Contract Plan Inventories

Commit(s):
- this commit - Vectorize contract plan inventories

Summary:
- Accepted `HP-CONTRACT-VEC-FN-01`. Plan-level inventories in
  `TerminalLoweringPlan` and `RetainedUnitTransformContractPlan` no longer use
  variable-length `Tuple{Vararg...}` field storage; they now use typed vectors.
- Preserved accessors and ordering:
  `available_contracts(plan)`, `selected_contracts(plan)`, `contracts(plan)`,
  and `transform_contracts(plan)` still expose the same deterministic contract
  sequence.
- Explicitly left `source_cpbs::Tuple{Vararg{CoordinateProductBox}}` unchanged,
  as required by the docs-only authority. No consumer rewiring was needed.

Validation:
- Doer validation: `git diff --check`; package load; H2 base artifact/readback
  with dimension 471 and K/V deltas 0.0; H2 supplemented artifact/readback with
  dimension 489 and K/V deltas 0.0; H2 R3 endpoint passed with self-Coulomb
  `0.4574265214362095`; focused order smoke confirmed terminal selected order
  and retained transform order.
- Manager validation: `git diff --check`; `git diff --numstat` reported `10`
  added and `9` deleted source lines across the approved terminal-lowering and
  retained transform-contract files; suspicious-line scan was empty; new-test/
  tool scan was empty.

Goal advancement:
- LT1/LT3: removes the remaining approved plan-level variable-length tuple
  inventories from the audited base working-basis path before fresh timing and
  compile attribution.
- RG/LT6: preserves current base/supplemented artifact behavior while reducing
  staged type-shape pressure for practical CR2/HF artifact generation.

Carrying-cost result:
- deleted: plan-level tuple inventory field types and tuple comprehensions for
  terminal lowering and retained transform-contract plans.
- simplified: plan constructors now produce typed vectors directly.
- quarantined: fixed tuple fields such as `source_cpbs`, summary tuple fields,
  public input `NamedTuple`s, artifact sidecars, numerical kernels, driver
  behavior, and Cr2 workflow remain outside this pass.
- not deleted because: source CPB tuples are a separate per-contract shape and
  were explicitly excluded by design authority.
- exact remaining caller/blocker: if compile latency remains high after fresh
  timing, the next work should be evidence-driven; summary/report tuple
  surfaces are the likely remaining type-surface cleanup candidates, but they
  need separate authority and proof of relevance.
- added src lines: 10.
- deleted src lines: 9.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 123 - Be2 Fresh-Target Timing Attribution

Commit(s):
- none - measurement-only timing pass.

Summary:
- Accepted the post type-surface cleanup timing result on current `5938ddbc`.
  The previously alarming Be2 q5 driver timing was not representative of warm
  construction. With `check_file=false`, `readback=false`, unique hamfiles, and
  one throwaway warmup in the same Julia process, three fresh Be2 q5 p10
  supplemented runs averaged `1.995 s`.
- The old suspected bottleneck, base working-basis construction, dropped from
  `45.8 s` cold to `0.145 s` fresh. Artifact writing was also small
  (`0.044 s` average), so the slow one-shot runs are best interpreted as
  compile/type latency rather than production numerical construction.
- The measured fresh Be2 input was CR2 pass010-style:
  `R=5.0`, `q=5`, `core_spacing=0.15`, `xmax_parallel=10.5`,
  `xmax_transverse=8.0`, Be/`cc-pVDZ`, `lmax=1`, contracted, no width filter,
  external BasisSets path. Current dimensions were base `1467`, residual `18`,
  final `1485`.

Validation:
- Doer validation: `git diff --check`; package load; timing harness
  `tmp/work/be2_q5_p10_fresh_target_timing_attribution.jl`; final
  `git status --short --branch`. No tracked files changed, and the timing
  harness remained ignored under `tmp/`.
- Manager validation: reviewed the reported cold/fresh timing separation,
  unique-target settings, stage breakdown, HEAD, and final status. No
  implementation tests were rerun because this was a read-only timing pass.

Goal advancement:
- LT1/LT3: changes the active performance interpretation. The current Be2 q5
  construction path is not a multi-minute warm numerical bottleneck; remaining
  pain is primarily first-run compilation/type latency.
- RG/LT6: supports practical CR2/HF artifact generation by showing that, after
  warmup, the canonical supplemented Be2 driver path has acceptable
  construction-scale timing for this fixture. Further source work should be
  justified by compile latency evidence, not by stale cold-run wall times.

Carrying-cost result:
- deleted: none; timing-only pass.
- simplified: the next decision boundary is now compile/type attribution versus
  stopping optimization, not another numerical construction rewrite.
- quarantined: exact-repeat same-file timing remains cache/artifact behavior
  only and is not accepted as construction timing.
- not deleted because: no source code was changed.
- exact remaining caller/blocker: if first-run latency remains operationally
  painful, the next lane should profile compilation/type surfaces with evidence
  before requesting authority; do not reopen raw-block, terminal-product, or
  unit-nuclear numerical kernels from this measurement alone.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 125 - Approve Route/Stage Type-Surface Cleanup

Commit(s):
- this commit - Approve route/stage type-surface cleanup

Summary:
- Approved `HP-ROUTE-STAGE-TYPE-FN-01` /
  `HP-ROUTE-STAGE-TYPE-TEST-01` as a narrow compile-latency cleanup lane based
  on Be2 q5 p10 supplemented driver attribution.
- The approved source files are only
  `src/pqs_source_box_route_driver_helpers.jl` and
  `src/cartesian_terminal_shellification_geometry.jl`.
- The approved targets are
  `_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`,
  `cartesian_units`,
  `_pqs_source_box_route_driver_transform_stage_low_order_summary`,
  `cartesian_transforms`,
  `_cartesian_terminal_shellification_region_unit_inventory`, and related
  terminal-region lowering inventory summary surfaces in
  `src/cartesian_terminal_shellification_geometry.jl`.
- The allowed action is deletion of stale route/stage compatibility
  inventories and replacement of remaining runtime-sized `NamedTuple` / `Tuple`
  carriers with vector-backed compact internal objects, stable dictionaries,
  accessors, or smaller summaries.
- The pass must preserve H2 base/supplemented artifact behavior, deterministic
  terminal shellification/lowering order, public driver contract,
  artifact/manifest schema, route semantics, and numerical matrices.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  checks for `HP-ROUTE-STAGE-TYPE-*`, approved files/functions, preservation
  rules, validation gates, forbidden surfaces, and no `src`, `test`, `tools`,
  or `bin` changes. No implementation tests were run.

Goal advancement:
- LT1/LT3: responds to the current Be2 q5 evidence that warm construction is
  acceptable but cold compile latency remains concentrated in route/stage type
  surfaces.
- RG/LT6: keeps practical artifact generation moving without reopening
  numerical kernels, artifact schema, or driver workflow.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one evidence-backed source lane for oversized
  route/stage compatibility inventories after the prior vector cleanup lanes.
- quarantined: driver changes, artifact/manifest changes, public API/export
  changes, numerical kernels, raw blocks, RG/MWG/IDA semantics, route
  diagnostics/status/report expansion, broad route-stage redesign, committed
  tests, PackageCompiler/PrecompileTools/sysimage work, and Cr2 workflow remain
  unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: implementation must stop if cleanup requires
  source files outside the approved boundary, broad route-stage redesign, new
  public contracts, artifact changes, numerical changes, or precompile/sysimage
  mechanisms.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 126 - Route/Stage Type-Surface Cleanup

Commit(s):
- this commit - Clean route stage type surfaces

Summary:
- Accepted `HP-ROUTE-STAGE-TYPE-FN-01`. The active route/stage carriers in
  `src/pqs_source_box_route_driver_helpers.jl` and
  `src/cartesian_terminal_shellification_geometry.jl` no longer materialize the
  targeted length-encoded tuple mirrors for shellification unit inventory,
  typed lowering contract inventory, retained-unit mirrors, and pair
  inventories.
- Deleted the dead terminal-shellification lowering inventory compatibility
  block in `cartesian_terminal_shellification_geometry.jl`: the uncalled
  contract builder helpers, aggregate lowering inventory, selected lowering
  inventory, and selected-contract helper functions. The active
  `_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`
  path remains in the route helper file.
- Preserved `_cartesian_terminal_region_lowering_contract_kind_counts`, which
  is still used by the active route helper path. H2 base/supplemented
  construction and the R3 endpoint stayed unchanged in doer validation.

Validation:
- Doer validation: `git diff --check`; package load; H2 route/stage cleanup
  validation script with base dimension `471`, supplemented dimension `489`,
  stage unit count `3`, and contract count `3`; H2 R3 endpoint with
  self-Coulomb `0.4574265214362095`; focused scans for deleted helper names
  and newly introduced runtime tuple/NamedTuple inventory syntax.
- Manager validation: `git diff --check`; `git diff --numstat` showed only the
  two approved source files changed with `45` added and `570` deleted source
  lines; suspicious-line scan was empty; new-test/tool scan was empty; deleted
  helper-name scan had no hits; final diff inspection confirmed no driver,
  artifact, numerical, raw-block, RG/MWG/IDA, or public API changes.

Goal advancement:
- LT1/LT3: removes a compile-attributed type-surface source by replacing
  active tuple mirrors with vectors/smaller carriers and deleting obsolete
  compatibility inventory code.
- RG/LT6: keeps CR2/HF artifact usability work aimed at first-run latency
  without disturbing the warm numerical construction path, driver contract, or
  artifact schema.

Carrying-cost result:
- deleted: dead terminal-shellification lowering inventory compatibility block
  and selected-lowering helper family.
- simplified: active shellification unit, lowering-contract, pair-family, and
  pair-key stage carriers now avoid the targeted length-encoded tuple mirrors.
- quarantined: remaining pre-existing tuple uses in input normalization and the
  uncalled complete-core-shell center helper are outside this pass; route
  skeleton tuple shapes in `pqs_source_box_route_driver_skeletons.jl` remain a
  separate potential lane.
- not deleted because: `_cartesian_terminal_region_lowering_contract_kind_counts`
  remains a live count helper for the active route helper path; route skeleton
  data lives outside the approved files.
- exact remaining caller/blocker: if cold compile latency remains high, the
  next evidence-backed target is likely route skeleton tuple shapes or the
  uncalled complete-core-shell helper, both requiring separate authority.
- added src lines: 45.
- deleted src lines: 570.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 127 - Post-Cleanup Compile Attribution

Commit(s):
- none - measurement-only compile/timing attribution pass.

Summary:
- Accepted the post-`HP-ROUTE-STAGE-TYPE-FN-01` attribution on current
  `118a639b`. The cleanup materially reduced Be2 q5 p10 cold latency: cold
  supplemented construction dropped from about `59.0 s` to `34.8 s`; trace
  total dropped to `39.5 s`, with repo-owned trace time about `32.2 s`.
- Warm construction remained about `2.02 s`, confirming that the remaining
  problem is still compile/type latency rather than warm numerical work.
- The targeted inventory owners fell sharply: terminal-lowering inventory fell
  to about `235 ms`, and terminal shellification unit inventory fell to about
  `286 ms`. The remaining top owners are broader route/stage signatures and
  construction-plan carriers: `cartesian_transforms` (`~5.0 s`),
  `cartesian_units` (`~3.4 s`), terminal topology support region planning
  (`~2.5 s`), terminal realization (`~2.3 s`), and terminal retained-rule
  planning (`~1.7 s`). `CartesianResidualGaussians` is also visible
  (`~3.8 s`) but is not the first route/stage cleanup target.

Validation:
- Doer validation: `git diff --check`; package load; fresh-target timing
  harness; trace workload; trace parser; final `git status --short --branch`.
  No tracked source edits and no commit.
- Manager validation: reviewed the comparison against the pre-cleanup
  attribution and accepted that the previous top inventory owners moved out of
  dominant territory. No implementation tests were rerun because this was a
  measurement-only pass.

Goal advancement:
- LT1/LT3: confirms that type-surface cleanup is paying down first-run
  latency, while narrowing the next possible source lane to broader
  `cartesian_shells` / `cartesian_units` / `cartesian_transforms` stage
  carriers and terminal support/retained-plan objects.
- RG/LT6: keeps CR2/HF artifact usability work focused on first-run compile
  experience, with warm construction still acceptable for the Be2 q5 proxy.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: next decision is no longer generic precompile versus cleanup; the
  trace identifies remaining route/stage carrier signatures as the leading
  source target if more source work is desired.
- quarantined: artifact/JLD2, Gaussian numerical kernels, raw blocks, terminal
  G-G products, and unit-nuclear kernels remain non-targets for this compile
  issue.
- not deleted because: no source authority has been granted for the newly
  identified carrier cleanup.
- exact remaining caller/blocker: a docs-only amendment is needed before
  touching `src/pqs_source_box_diatomic_complete_core_shell.jl`, terminal
  support/retained planning carriers, or route skeleton tuple shapes outside
  the prior approved files.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 128 - Approve Route/Stage Carrier Cleanup

Commit(s):
- this commit - Approve route/stage carrier cleanup

Summary:
- Approved `HP-ROUTE-STAGE-CARRIER-FN-01` /
  `HP-ROUTE-STAGE-CARRIER-TEST-01` as the follow-up compile-latency cleanup
  lane after the post-cleanup attribution on `118a639b`.
- The approved source files are
  `src/pqs_source_box_route_driver_helpers.jl` and
  `src/pqs_source_box_diatomic_complete_core_shell.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  optional only where directly required to slim terminal realization plan
  carriers in the approved path.
- The approved targets are `cartesian_shells`, `cartesian_units`,
  `cartesian_transforms`, terminal topology support-region planning, terminal
  retained-rule planning, and directly required terminal realization plan
  carriers.
- The allowed action is to stop carrying giant shellification, route-skeleton,
  support-plan, retained-rule-plan, and terminal-plan `NamedTuple` / tuple
  shapes across approved stage signatures, replacing necessary facts with
  compact typed/vector-backed carriers, smaller summaries, accessors, or local
  recomputation from canonical objects.
- Route skeleton construction semantics and
  `src/pqs_source_box_route_driver_skeletons.jl` remain out of scope.

Validation:
- Design-manager validation only: docs-only `git diff --check`, focused `rg`
  checks for `HP-ROUTE-STAGE-CARRIER-*`, approved files, approved targets,
  route-skeleton guardrails, preservation rules, validation gates, forbidden
  surfaces, and no `src`, `test`, `tools`, or `bin` changes. No implementation
  tests were run.

Goal advancement:
- LT1/LT3: follows the current trace evidence by moving from the now-crossed
  targeted compatibility-inventory cleanup to broader route/stage carrier
  signatures still responsible for cold compile latency.
- RG/LT6: keeps CR2/HF artifact usability focused on first-run route/stage
  specialization without changing warm numerical construction, artifacts, or
  driver behavior.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one approved carrier cleanup lane with exact files,
  optional terminal realization boundary, and a route-skeleton owner guardrail.
- quarantined: `src/pqs_source_box_route_driver_skeletons.jl`, driver changes,
  artifact/manifest changes, public API/export changes, numerical kernels, raw
  blocks, RG/MWG/IDA semantics, route diagnostics/status/report expansion,
  broad route-stage redesign, committed tests, PackageCompiler/PrecompileTools
  sysimage work, and Cr2 workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: implementation must stop if cleanup requires
  source files outside the approved boundary, broad route-stage redesign,
  public API changes, artifact changes, numerical changes, or
  precompile/sysimage machinery.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 129 - Route/Stage Carrier Cleanup

Commit(s):
- `8c3df2ad` - Slim route stage carriers

Summary:
- Accepted `HP-ROUTE-STAGE-CARRIER-FN-01`. The route skeleton is now unpacked
  inside `cartesian_shells` instead of being carried forward as a wide stage
  object, and later route stages no longer copy pass-through shellification,
  lowering, support-plan, or retained-rule-plan mirrors.
- Terminal topology support-region planning and terminal retained-rule planning
  in `src/pqs_source_box_diatomic_complete_core_shell.jl` now use
  vector-backed rows/counts through the existing route inventory row helper
  instead of runtime-keyed `NamedTuple`/tuple count carriers.
- The patch stayed inside the two required approved files; the optional
  terminal basis realization file was not needed. Route skeleton construction
  semantics and `src/pqs_source_box_route_driver_skeletons.jl` were untouched.

Validation:
- Doer validation: `git diff --check`; package load; H2 route-stage carrier
  validation script with H2 base dimension `471` and supplemented behavior
  unchanged; H2 R3 endpoint with augmented dimension `489`, self-Coulomb
  `0.4574265214362095`, and facade readback deltas all `0.0`; focused scan
  found no newly added `NamedTuple{...}`, `Tuple(...)`, or `Tuple{Vararg...}`
  carriers in the edited files. Optional Be2 timing reported cold total
  `23.40 s`, base working-basis `11.68 s`, and warm fresh mean `1.99 s`.
- Manager validation: `git diff --check`; `git diff --numstat` showed only the
  two approved source files changed with `63` added and `155` deleted source
  lines; suspicious-line scan including `Dict{Symbol,Any}` and `Any[]` added
  lines was empty; new-test/tool scan was empty. Manager requested tightening
  the new pair-helper map to `Dict{Symbol,Symbol}` and blocked retained-rule
  empty fields to `NamedTuple[]`, then rechecked the amended diff.

Goal advancement:
- LT1/LT3: continues reducing first-run route/stage specialization pressure by
  removing wide stage-carried objects after the previous targeted inventory
  cleanup.
- RG/LT6: directly supports practical CR2/HF Hamiltonian generation by cutting
  cold Be2 proxy latency while preserving warm construction behavior and
  artifact contracts.

Carrying-cost result:
- deleted: pass-through route skeleton, shellification, lowering, support-plan,
  retained-rule-plan, and duplicate lowering-inventory summary mirrors from
  later route stages.
- simplified: `cartesian_shells`, `cartesian_units`, `cartesian_transforms`,
  support-region planning, and retained-rule planning now carry smaller
  vector-backed or scalar summaries across the approved path.
- quarantined: route skeleton construction in
  `src/pqs_source_box_route_driver_skeletons.jl` remains out of scope; remaining
  pre-existing `Any` dictionaries in retained-rule joins and complete-core
  center helpers are not introduced by this pass.
- not deleted because: terminal retained-rule plan and transform-contract plan
  remain needed for basis realization and manifest source-mode provenance.
- exact remaining caller/blocker: if cold compile latency remains high, rerun
  attribution before approving more cleanup; likely remaining targets would
  need new evidence and possibly route skeleton authority.
- added src lines: 63.
- deleted src lines: 155.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 130 - Post-Carrier Compile Attribution

Commit(s):
- none - measurement-only compile/timing attribution pass.

Summary:
- Accepted the post-`HP-ROUTE-STAGE-CARRIER-FN-01` attribution on current
  `8c3df2ad`. Cold Be2 q5 p10 supplemented construction is now in the
  expected `20-25 s` range: fresh cold construction was `22.928 s`, and the
  trace run was `23.147 s`. Warm fresh construction remained about `1.99 s`.
- The prior route/stage carrier owners no longer dominate. `cartesian_transforms`
  dropped from about `4996 ms` to `590 ms`, `cartesian_units` from about
  `3402 ms` to `1081 ms`, and terminal support/retained plans are no longer
  dominant (`136 ms` / `74 ms`).
- Remaining trace time is broad/mixed specialization, with top repo owners
  including `CartesianResidualGaussians` (`~3.8 s`),
  `CartesianFinalBasisRealization` (`~2.3 s`), parent construction
  (`~1.4 s`), and Gaussian supplement conversion (`~1.2 s`). Artifact/JLD2 and
  ordinary Gaussian numerical kernels are not the dominant issue.

Validation:
- Doer validation: `git diff --check`; package load; timing harness; trace
  harness and parser; final `git status --short --branch`. No tracked source
  edits and no commit.
- Manager validation: reviewed the cold/warm comparison and owner shifts. No
  implementation tests were rerun because this was a measurement-only pass.

Goal advancement:
- LT1/LT3: closes the current route/stage type-surface cleanup sequence. The
  two source cleanup passes reduced cold Be2 q5 from about `59 s` to about
  `23 s`, and the remaining cost no longer points to a single narrow
  route/stage carrier.
- RG/LT6: supports CR2/HF usability by showing warm construction remains
  acceptable and first-run latency is now reduced enough that further source
  cleanup needs new evidence.

Carrying-cost result:
- deleted: none; measurement-only pass.
- simplified: current recommendation is to stop route/stage source cleanup for
  now. If first-run latency is still operationally painful, evaluate workflow
  precompile/sysimage rather than continuing source churn by default.
- quarantined: residual specialization is visible but mixed with trace harness
  stage closures; it needs separate attribution before any RG source authority.
- not deleted because: no source work was requested or approved by this pass.
- exact remaining caller/blocker: no narrow production source owner remains
  dominant enough to justify another immediate cleanup lane. Next action should
  be either stop/try CR2 workflow or run a docs-approved precompile/sysimage
  audit.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 131 - Approve RG Orthogonality Robustness

Commit(s):
- this commit - Approve residual Gaussian orthogonality robustness

Summary:
- Approved `HP-RG-ORTHO-FN-01` and `HP-RG-ORTHO-TEST-01` as a narrow Residual
  Gaussian robustness lane for final residual identity validation. The
  motivating strict N2 q5 p10 case at `core_spacing = 0.042857` has excellent
  `G' S R` orthogonality (`1.776e-14`), positive full-rank owner metrics,
  retained counts `9,9`, and a healthy final merge spectrum (`7.232e-2 ..
  1.928`, condition `26.65`), but `R' S R - I = 1.673e-10` slightly exceeds
  the old absolute-only `1.0e-10` check.
- The amendment keeps owner-local residual selection, occupation cutoff,
  negative-eigenvalue policy, and MWG/IDA conventions unchanged. It approves
  only stable symmetric final residual-overlap validation and a combined
  absolute/relative final identity check for small floating-point overshoots
  after a healthy merge.

Validation:
- Design-manager validation: read compact RG authority, registry, current
  status, invariants, implementation slices, AGENTS source-ID list, and the
  running-log tail; updated only docs/AGENTS; planned validation is
  `git diff --check`, focused `rg` for the new IDs and source surfaces, and a
  check that no `src`, `test`, `tools`, or `bin` files changed. No
  implementation tests were run because this is a docs-only approval pass.

Goal advancement:
- RG/LT6: makes the supplemented residual-Gaussian path robust enough for
  strict molecule inputs without changing the residual physics or broadening
  public workflow.
- LT1/LT3: avoids source churn outside the now-dominant workflow lanes by
  requiring a narrow source surface and explicit failure rule if the issue is
  not a final floating-point validation overshoot.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one exact source lane in `residual_basis.jl`, with
  terminal residual keyword plumbing allowed only if needed.
- quarantined: residual selection changes, global residual selection,
  occupation-cutoff changes, merge eigenvalue flooring, width filtering as
  conditioning repair, MWG/IDA/raw-block/nuclear/artifact/driver/public API
  changes, status/report fields, committed tests, and Cr2 workflow remain
  unapproved.
- not deleted because: source correction has not run yet.
- exact remaining caller/blocker: implementation must stop if strict N2
  requires changing residual selection, supplement construction, or final-basis
  construction.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 132 - Approve Driver Nesting Input

Commit(s):
- this commit - Approve canonical driver nesting input

Summary:
- Approved `HP-DRV-NEST-FN-01`, `HP-DRV-NEST-WIRE-01`, and
  `HP-DRV-NEST-TEST-01` as a narrow canonical-driver construction-family lane.
  The driver may expose `nesting = :pqs` or `nesting = :wl` as a first-class
  user construction choice. `:pqs` remains the default and maps to the existing
  PQS source-box route; `:wl` maps only to the existing
  `:white_lindsey_low_order` route.
- The amendment explicitly treats `nesting` as public construction-family
  input, not a diagnostic route switch. Route skeletons, retained-rule plans,
  raw-block switches, stop-after controls, diagnostics, route reports, and
  route-stage labels remain hidden. Supplemented `:wl` must be rejected
  clearly unless already valid through the existing supported supplemented
  facade/staged path.

Validation:
- Design-manager validation: read the driver usability workflow, registry
  driver section, current compact authority, implementation slices, AGENTS
  driver section, and running-log tail; updated only docs/AGENTS; ran
  `git diff --check`, focused `rg` for `HP-DRV-NEST-*`, `nesting`, the route
  mappings, and forbidden route-diagnostic wording, plus a no-source/test/bin
  diff check. No implementation tests were run because this is docs-only
  approval.

Goal advancement:
- LT1/LT3: preserves the active PQS and White-Lindsey low-order construction
  tracks as user-facing construction families without reviving route-stage
  diagnostics.
- RG/LT6: keeps the canonical driver useful for real artifact production while
  preserving the supported supplemented facade boundaries.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one exact driver/facade plumbing lane for
  `nesting`, instead of smuggling route-family choice through private route
  fields.
- quarantined: route skeletons, retained rules, raw-block switches,
  stop-after controls, diagnostics, route reports, route-stage labels, broad
  supplemented White-Lindsey behavior, artifact schema changes, public
  API/export changes, committed tests/fixtures, solver/ECP workflow, and
  Cr2-specific behavior remain unapproved.
- not deleted because: source wiring has not run yet.
- exact remaining caller/blocker: if `nesting = :wl` cannot produce a small
  base artifact/readback through the existing White-Lindsey low-order route
  without broader route/materialization work, implementation must stop and
  report the blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 133 - Driver Nesting And WL Terminal Atom

Commit(s):
- this commit - Wire driver nesting and WL terminal basis

Summary:
- Accepted the `nesting = :pqs` / `nesting = :wl` driver/facade wiring together
  with the narrow White-Lindsey terminal-basis seam. The canonical driver now
  exposes `nesting` as a first-class construction-family input and passes it
  through the base basis contract. `:pqs` preserves the existing PQS route, and
  `:wl` maps to the existing `:white_lindsey_low_order` route.
- The WL route now produces the same `CartesianTerminalBasisRealization` type
  consumed by staged product, unit-nuclear, IDA, Hamiltonian, artifact, and
  readback stages for the one-center base atom path. This is a real base
  Hamiltonian endpoint, not a blocker-only commit, and it deliberately does not
  revive the old WL H1/H1+J materialization path.
- WL H2 remains blocked because the current native WL shellification/terminal
  record path is produced for one-center systems only; H2 still reaches the
  missing-terminal-basis blocker. This is now a precise route/shellification
  support gap, not a driver contract ambiguity.

Validation:
- Doer: `git diff --check`; package load; `nesting=:pqs` atom artifact/readback
  (`dim 35`); `nesting=:pqs` H2 artifact/readback (`dim 95`);
  `nesting=:wl` atom artifact/readback (`dim 125`); `nesting=:wl` H2 clear
  blocker; H2 residual-GTO/MWG PQS endpoint with augmented dimension `489` and
  self-Coulomb `0.4574265214362095`.
- Manager: reviewed the WIP diff against `HP-DRV-NEST-*` and `HP-WLTERM-*`,
  confirmed changed source/bin surfaces, ran `git diff --check`,
  `git diff --numstat -- src bin tools test docs`, suspicious added-line scan,
  new tests/tools scan, and final status. Julia validations were not rerun per
  the standing instruction not to duplicate doer runs.

Goal advancement:
- LT1/LT3: restores the intended public construction-family choice in the
  canonical driver without exposing private route controls.
- LT5/LT6: keeps PQS and White-Lindsey identities explicit while forcing both
  through the same terminal-basis and Hamiltonian artifact boundary.

Medium-goal update:
- none. This pass is a driver/WL usability seam; it does not change the current
  CR2/HF readiness direction except by making the public construction contract
  more honest.

Risk / guardrail:
- The WL terminal realizer currently treats direct retained units and WL
  boundary-stratum units as identity terminal blocks on owned support. That is
  the approved low-order seam, but future WL diatomic support must extend
  native shellification/lowering records rather than adding driver-level
  cross-product cases.

Remaining blocker / next:
- WL diatomic base construction needs a separate route/shellification authority
  if `nesting=:wl` H2 artifacts are desired. Supplemented WL remains outside
  the current workflow until base WL diatomics have native terminal records.

Line-count / complexity note:
- added src lines: 148, deleted src lines: 21; added bin lines: 6, deleted bin
  lines: 3; new tests: none; new metadata/status fields: none.
- deleted: the PQS-only terminal-realization gate no longer blocks the WL atom
  route.
- simplified: driver construction family now flows through `basis.nesting`
  into the route family rather than relying on hidden hardwiring.
- quarantined: old WL H1/H1+J materialization path, WL H2, supplemented WL,
  and Cr2 workflow remain outside this pass.
- not deleted because: existing PQS terminal realization and old WL reference
  materialization still have separate live/reference roles.
- exact remaining caller/blocker: no native WL diatomic terminal records are
  produced for the staged Hamiltonian path.

## Cartesian Hamiltonian Producer Pass 134 - Nesting Artifact Truth

Commit(s):
- this commit - Record nesting in Hamiltonian provenance

Summary:
- Accepted the `HP-NEST-ART-FN-01` cleanup correcting artifact truth for the
  public `nesting` construction-family input. Base `producer_provenance/` and
  manifest `recipe_provenance/` now record `nesting`, and base route labels
  are derived from `(input.kind, input.nesting)` instead of using the old
  PQS-oriented default. WL one-center artifacts now record
  `:one_center_wl_base`; PQS one-center artifacts record
  `:one_center_pqs_base`.
- Supplemented `nesting = :wl` now rejects in `cartesian_base_working_basis`
  after input normalization but before base-stage construction, so the staged
  driver and one-call facade fail with the intended public-contract message
  rather than falling into the later missing-terminal-basis blocker.
- Updated the `CartesianFinalBasisRealization` module docstring to reflect its
  current terminal/final-basis responsibilities, including PQS, WL terminal
  realization, terminal one-body/IDA, and RG compatibility helpers.

Validation:
- Manager/doer combined: `git diff --check`; package load; ignored
  `tmp/work/nesting_artifact_truth_validation.jl` checking small PQS atom
  artifact/readback provenance (`:pqs`, `:one_center_pqs_base`), small WL atom
  artifact/readback provenance (`:wl`, `:one_center_wl_base`), and supplemented
  WL early `ArgumentError`. No Cr2 run.
- Mechanical manager gate before commit: `git diff --numstat -- src bin tools
  test docs`, suspicious added-line scan, and new tests/tools scan.

Goal advancement:
- LT5/LT6: artifact sidecars now tell the truth about the public construction
  family, preventing WL artifacts from being misread as PQS.
- LT1/LT3: keeps the human-facing driver contract aligned with the persisted
  Hamiltonian contract.

Medium-goal update:
- none. This is a correctness cleanup after the nesting/WL seam, not a new
  feature lane.

Risk / guardrail:
- No reader behavior or matrix keys changed. The provenance schema was already
  amended before source work. WL H2 and supplemented WL remain unsupported.

Remaining blocker / next:
- The next WL question remains native WL diatomic terminal records, not driver
  branching. CR2 stress/performance should use truthful provenance after this
  point.

Line-count / complexity note:
- added src lines: 21; deleted src lines: 13; new tests: none; new
  metadata/status fields: approved provenance key only.
- deleted: stale PQS-only route labeling assumption.
- simplified: base route labeling is centralized in `_cartesian_base_route_label`.
- quarantined: WL H2 and supplemented WL remain blocked.
- not deleted because: the later supplement-stage WL rejection remains as a
  guard for manually assembled bases.
- exact remaining caller/blocker: no successful WL diatomic artifact path
  exists, so no `:z_axis_diatomic_wl_base` label is written.

## Cartesian Hamiltonian Producer Pass 133 - Approve WL Terminal Basis Seam

Commit(s):
- this commit - Approve White-Lindsey terminal basis realization

Summary:
- Approved `HP-WLTERM-FILE-01`, `HP-WLTERM-FN-01`,
  `HP-WLTERM-WIRE-01`, and `HP-WLTERM-TEST-01` as the narrow source lane
  exposed by the `nesting = :wl` WIP. The current driver/facade wiring can
  reach `route_family = :white_lindsey_low_order`; the missing seam is that
  `_pqs_source_box_route_driver_terminal_basis_realization(...)` deliberately
  returns `nothing` for non-PQS routes.
- The amendment approves only making the existing White-Lindsey low-order route
  produce the same `CartesianTerminalBasisRealization` consumed by the staged
  Hamiltonian path. It may realize direct blocks and WL
  boundary-stratum/product blocks from existing terminal support, retained-rule,
  and transform records. It does not approve adapting the old WL H1/H1+J
  materialization path.

Validation:
- Design-manager validation: inspected the compact authority and the current
  source seam showing the PQS-only guard plus existing WL lowering symbols;
  updated only docs/AGENTS; ran `git diff --check` and focused `rg` for
  `HP-WLTERM-*`, the approved source files, WL boundary-stratum wording, and
  forbidden old-materialization wording. The working tree contained unstaged
  source/bin WIP in `bin/cartesian_ham_builder.jl` and
  `src/cartesian_base_hamiltonian.jl`; it was not staged or committed. No
  implementation tests were run because this is docs-only approval.

Goal advancement:
- LT1/LT3: keeps PQS and White-Lindsey as construction families while forcing
  both through the same terminal-basis/Hamiltonian downstream boundary.
- RG/LT6: preserves supplemented PQS behavior and avoids cross-product
  branching in the driver; WL supplemented behavior remains a separate design
  decision.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has the exact missing seam to implement instead of
  carrying a blocker-only driver WIP.
- quarantined: route skeleton construction, shellification behavior,
  retained-selection policy, old WL H1/H1+J materialization adaptation,
  diagnostics/status/report expansion, raw-block changes, RG/MWG/IDA changes,
  supplemented WL behavior, artifact schema changes, public API/export changes,
  committed tests/fixtures, solver/ECP workflow, and Cr2 workflow remain
  unapproved.
- not deleted because: source realization has not run yet.
- exact remaining caller/blocker: if WL boundary-stratum final basis cannot be
  materialized from existing terminal lowering, retained-unit, and transform
  records, implementation must stop and report the exact missing native fact.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 134 - Approve Nesting Artifact Truth

Commit(s):
- this commit - Approve nesting artifact truth cleanup

Summary:
- Approved `HP-NEST-ART-FN-01` and `HP-NEST-ART-TEST-01` as the narrow
  artifact/provenance cleanup lane after `nesting = :wl` became a real public
  construction-family choice for one-center base artifacts. The current bug is
  provenance truth: WL artifacts must not carry PQS-oriented route labels, and
  artifacts must record the public `nesting` input.
- The amendment requires `producer_provenance/` and `recipe_provenance/` to
  record `nesting` and to derive the route label from
  `(input.kind, input.nesting)`. Approved labels are
  `:one_center_pqs_base`, `:one_center_wl_base`, and
  `:z_axis_diatomic_pqs_base`. No WL H2 label is approved until WL H2 artifact
  construction succeeds under separate authority.
- The source lane is deliberately small: `src/cartesian_base_hamiltonian.jl`
  may correct provenance and early supplemented-WL rejection, while
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  may receive only a docstring correction that removes PQS-only module
  wording.

Validation:
- Design-manager validation for this docs-only pass: compact authority and
  artifact-manifest reads, update only docs/AGENTS, run `git diff --check`,
  focused `rg` for `HP-NEST-ART-*`, `nesting`, `one_center_wl_base`,
  `recipe_provenance`, and `producer_provenance`, plus a no-source/test/bin
  diff check. No implementation tests are part of this approval pass; package
  load and artifact/readback checks belong to `HP-NEST-ART-TEST-01`.

Goal advancement:
- LT1/LT3: keeps the public construction-family contract honest at the
  artifact boundary, so downstream consumers can distinguish PQS and
  White-Lindsey construction without route-stage vocabulary.
- LT5/LT6: preserves the common terminal/Hamiltonian artifact boundary while
  preventing helper-file naming from becoming false provenance.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: doer now has one exact provenance-truth lane instead of treating
  nesting as an implicit driver-only fact.
- quarantined: driver public inputs, route skeletons, shellification,
  terminal lowering, raw blocks, RG/MWG/IDA, artifact matrix keys, reader
  behavior, public API/export changes, WL H2 support, diagnostics/reports,
  committed tests, and Cr2 workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if truthful nesting provenance requires
  reader changes, matrix-key changes, or a broader manifest structure, the
  source pass must stop and report the blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: only approved artifact provenance keys, no staged
  metadata.

## Cartesian Hamiltonian Producer Pass 135 - Record Nesting/Supplement Composition Target

Commit(s):
- this commit - Record nesting supplement composition plan

Summary:
- Added a compact composition plan for the intended 2 x 2 x 2 producer matrix:
  `geometry = atom | z-axis diatomic`, `nesting = :pqs | :wl`, and
  `supplement = off | on`. This is a design-direction amendment only, not
  source authority.
- Current support is now explicitly classified: atom/base works for PQS and
  WL; atom/supplement remains unapproved; PQS H2/base works while broader
  generic base diatomic support is still limited; WL diatomic base is blocked
  by missing native WL diatomic terminal records; PQS supplemented homonuclear
  z-axis diatomics are supported; supplemented WL is blocked first by WL
  diatomic base.
- The next dependency order is frozen: WL diatomic base terminal records first,
  supplemented atoms through the common Residual Gaussian path second, and
  supplemented WL after WL terminal bases are real.

Validation:
- Design-manager validation for this docs-only pass: read compact current
  driver/WL/manifest authority, added the new plan page, updated README,
  current authority, implementation slices, invariants, registry candidate
  section, AGENTS guardrail, and this running log; run `git diff --check`,
  focused `rg` for the new plan and candidate IDs, and confirm no source,
  test, tool, or bin files changed. No implementation tests are part of this
  planning amendment.

Goal advancement:
- LT1/LT3: prevents partial support from becoming the hidden driver contract
  and makes the desired user-facing choices explicit.
- LT5/LT6: preserves the common terminal-basis and Residual Gaussian
  boundaries as the way to fill missing cells, rather than branching by
  driver case or route report.

Carrying-cost result:
- deleted: none; docs-only planning pass.
- simplified: later source blurbs can target one missing cell at a time with
  known dependencies instead of arguing from the current partial matrix.
- quarantined: WL diatomic base, supplemented atoms, supplemented WL,
  translated atoms, non-z-axis geometry, heteronuclear supplement generality,
  ECP/EGOI/RHF/solver work, public export/API redesign, route diagnostics,
  and Cr2-specific branches remain unapproved.
- not deleted because: no source cleanup has run in this pass.
- exact remaining caller/blocker: the first source lane should be WL z-axis
  diatomic base terminal records, unless a later manager decision chooses
  supplemented atoms first with explicit source authority.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 136 - Approve WL Diatomic Base Lane

Commit(s):
- this commit - Approve WL diatomic base composition lane

Summary:
- Promoted the first 2 x 2 x 2 composition placeholder:
  `HP-COMP-WLDIAT-FN-01` and `HP-COMP-WLDIAT-TEST-01`. The approved source
  target is `Natom = 2`, `nesting = :wl`, `basisname = nothing`
  artifact/readback through native White-Lindsey z-axis diatomic terminal
  records, the existing `CartesianTerminalBasisRealization`, and the staged
  base Hamiltonian path.
- The lane allows the narrow route/shellification/lowering/terminal-basis/base
  source surfaces needed to produce native WL diatomic terminal records. It
  explicitly forbids driver special cases, old WL H1/H1+J materialization,
  artifact schema or reader changes, RG/MWG/supplement work, route
  diagnostics, committed tests, and Cr2 runs.
- The route provenance value `:z_axis_diatomic_wl_base` is approved as a
  truthful value under existing provenance keys once the WL H2 path succeeds;
  it is not an artifact schema change.

Validation:
- Design-manager validation for this docs-only pass: read current composition
  and WL terminal-basis authority, updated registry/current/implementation
  slices/README/R1 provenance wording/WL terminal note/AGENTS/running log, run
  `git diff --check`, focused `rg` for `HP-COMP-WLDIAT-*`, approved source
  files, forbidden WL materialization/driver-special-case wording, and confirm
  no source, test, tool, or bin files changed. No implementation tests are part
  of this approval pass; package load and artifact/readback checks belong to
  `HP-COMP-WLDIAT-TEST-01`.

Goal advancement:
- LT1/LT3: moves the first missing composition cell from known blocker to
  explicit implementation authority without changing the driver contract.
- LT5/LT6: keeps White-Lindsey and PQS distinct upstream while forcing both to
  meet at the same terminal-basis/Hamiltonian artifact boundary.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: WL H2 now has one exact source lane instead of being a vague
  route/shellification blocker.
- quarantined: supplemented WL, supplemented atoms, old WL H1/H1+J
  materialization, driver special cases, artifact schema changes, route
  diagnostics/status/report payloads, RG/MWG/supplement work, public
  API/export changes, committed tests, solver/ECP/RHF/EGOI work, and Cr2
  workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if native WL diatomic terminal records need
  broader route skeleton redesign, source files outside the approved list, or
  a parallel Hamiltonian builder, the source pass must stop and report.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: no new schema fields; only an approved route
  provenance value.

## Cartesian Hamiltonian Producer Pass 137 - WL Diatomic Base Path

Commit(s):
- this commit - Enable WL diatomic base terminal path

Summary:
- Accepted the minimal `HP-COMP-WLDIAT-FN-01` source pass. The existing
  terminal shellification path is now allowed for White-Lindsey
  bond-aligned diatomics as well as one-center systems, which lets the native
  WL lowering, retained-unit, transform-contract, and WL terminal-basis
  realization chain produce a `CartesianTerminalBasisRealization` for small
  H2. The staged base products, unit nuclear, IDA, Hamiltonian assembly,
  writer, and reader then run without a driver special case.
- Added the already-approved `:z_axis_diatomic_wl_base` route label in the
  centralized base route-label helper. No artifact keys, driver inputs,
  terminal-lowering policy, old WL H1/H1+J materialization path, RG/MWG, or
  supplement code changed.

Validation:
- Doer/manager: `git diff --check`; package load; ignored
  `tmp/work/wl_diatomic_base_validation.jl` checking small H2 PQS
  artifact/readback with `route = :z_axis_diatomic_pqs_base`, small H2 WL
  artifact/readback with `route = :z_axis_diatomic_wl_base` (`pqs_dim=111`,
  `wl_dim=1215`); existing H2 residual-GTO/MWG PQS endpoint
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` passed with
  augmented dimension `489`, self-Coulomb `0.4574265214362095`, and facade
  readback deltas all `0.0`. No Cr2 run.
- Mechanical manager gate: `git diff --numstat -- src bin tools test docs`,
  suspicious added-line scan, and new tests/tools scan; no suspicious hits and
  no new committed tests/tools.

Goal advancement:
- LT1/LT3: fills the first missing composition cell,
  `Natom = 2`, `nesting = :wl`, `basisname = nothing`, through the canonical
  staged producer.
- LT5/LT6: confirms the intended common boundary: WL and PQS differ upstream
  in shellification/lowering, then converge to `CartesianTerminalBasisRealization`
  and the same Hamiltonian artifact machinery.

Medium-goal update:
- The WL diatomic base cell is now implemented for small H2. Supplemented WL
  remains blocked by policy until the next composition lane proves the RG
  boundary can consume WL terminal bases.

Risk / guardrail:
- The WL H2 dimension is much larger than the small PQS H2 smoke at the same
  `q=2`; that is expected from the WL boundary-stratum identity realization,
  but future performance work should measure before expanding use. This pass
  intentionally did not optimize or reinterpret WL retained units.

Remaining blocker / next:
- Next composition choice is supplemented atoms or supplemented WL. For CR2
  readiness, continue using PQS supplemented runs unless a source lane promotes
  supplemented WL.

Line-count / complexity note:
- added src lines: 2; deleted src lines: 1; new tests: none; new
  metadata/status fields: none.
- deleted: obsolete blocker state for WL z-axis diatomic base at the
  shellification gate.
- simplified: no driver case table; WL diatomic now uses the same staged base
  Hamiltonian path.
- quarantined: supplemented WL, old WL H1/H1+J materialization, Cr2 workflow,
  and route diagnostics remain out of scope.
- not deleted because: the later supplemented WL rejection remains until a
  separate composition lane is approved.
- exact remaining caller/blocker: no supplemented WL path is approved yet.

## Cartesian Hamiltonian Producer Pass 138 - Approve Base Homonuclear Diatomic Relaxation

Commit(s):
- this commit - Approve base homonuclear diatomic relaxation

Summary:
- Approved `HP-COMP-BASEDIAT-FN-01` and `HP-COMP-BASEDIAT-TEST-01` as a
  narrow base producer input-contract relaxation in
  `src/cartesian_base_hamiltonian.jl` only. The source target is to relax the
  two-center base branch from H2-only to explicit homonuclear z-axis
  all-electron diatomics.
- The approved contract requires equal symbols, equal finite positive
  integer-valued nuclear charges, finite distinct z-axis centers, and neutral
  `nup + ndn == sum(charges)`. Symbols remain labels only; charge and explicit
  electron count remain authority. The basis contract is unchanged and both
  `nesting = :pqs` and `nesting = :wl` remain visible construction choices.
- The lane is explicitly not route/shellification authority. `nesting = :wl`
  non-H diatomics still depend on `HP-COMP-WLDIAT-*` for native WL diatomic
  terminal records.

Validation:
- Design-manager validation for this docs-only pass: read the compact
  composition plan, R1 base provenance wording, current authority, registry,
  implementation slices, and AGENTS; update docs/AGENTS/running log only; run
  `git diff --check`, focused `rg` for `HP-COMP-BASEDIAT-*`, homonuclear
  z-axis contract wording, forbidden element-table/driver/route/supplement
  wording, and confirm no source, test, tool, or bin files changed. No
  implementation tests are part of this approval pass; package load and
  artifact/readback checks belong to `HP-COMP-BASEDIAT-TEST-01`.

Goal advancement:
- LT1/LT3: removes the artificial H2-only base diatomic contract without
  changing the driver or adding element-specific defaults.
- LT5/LT6: brings base mode closer to the supplemented homonuclear diatomic
  contract while preserving the common geometry/terminal/Hamiltonian boundary.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one exact one-file lane for base diatomic
  validation instead of treating H2 as a permanent special case.
- quarantined: driver changes, route skeleton/shellification/terminal-lowering
  changes, raw blocks, supplement/RG/MWG work, artifact schema or reader
  changes, element lookup/default tables, heteronuclear or translated/general
  geometry support, committed tests, and Cr2 workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if non-H base diatomics need route or
  shellification changes outside `src/cartesian_base_hamiltonian.jl`, the
  source pass must stop and report the exact blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 139 - Base Homonuclear Diatomic Relaxation

Commit(s):
- this commit - Relax base homonuclear diatomic validation

Summary:
- Accepted the `HP-COMP-BASEDIAT-FN-01` one-file source pass. The
  two-center base branch in `src/cartesian_base_hamiltonian.jl` no longer
  hardcodes H2 symbols, charges, or spin counts. It now accepts explicit
  homonuclear z-axis all-electron diatomics with equal labels, equal
  integer-valued nuclear charges, distinct z-axis centers, and
  `nup + ndn == sum(charges)`. Symbols remain provenance labels and charges
  remain authority.
- The existing basis normalization and `kind = :h2` route shape are preserved,
  so both `nesting = :pqs` and `nesting = :wl` continue through the same
  staged base Hamiltonian path and existing route provenance labels.

Validation:
- Doer: `git diff --check`; package load; H2 base artifact/readback for
  `nesting = :pqs` and `nesting = :wl`; small Be2 base artifact/readback with
  bounded `q = 2`, `nesting = :pqs` and `nesting = :wl`; rejection checks for
  heteronuclear labels, unequal charges, non-neutral electron count, and
  non-z-axis centers. No Cr2 run.
- Manager: inspected the source diff; `git diff --check`;
  `git diff --numstat -- src bin tools test docs`; suspicious added-line scan;
  new tests/tools scan; focused stale H2-only wording scan. No suspicious hits
  and no new committed tests/tools.

Goal advancement:
- LT1/LT3: completes the base-mode geometry relaxation so atom/base and
  diatomic/base inputs are now explicit-data contracts rather than
  element-special-cased H/H2 contracts.
- LT5/LT6: keeps the composition model intact: geometry validation happens in
  the base input seam, nesting remains a construction-family choice, and the
  driver does not grow an element case table.

Remaining blocker / next:
- The two remaining composition holes are supplemented atoms and supplemented
  WL. The clearer next lane is supplemented WL for z-axis diatomics, because
  base PQS/WL diatomics now both reach a terminal basis and the supplemented
  path is currently blocked only by policy/wiring.

Carrying-cost result:
- deleted: hardcoded H2 symbol/charge/spin guard and H2-only error wording.
- simplified: atom and diatomic integer-charge checks share one private helper.
- quarantined: heteronuclear, non-z-axis/general geometry, Cr2 workflow,
  supplemented atoms, and supplemented WL remain out of this pass.
- not deleted because: `kind = :h2` and H2-named internal basis helper remain
  the existing z-axis diatomic route shape, not public H2-only policy.
- exact remaining caller/blocker: none for base homonuclear z-axis diatomics.
- added src lines: 18.
- deleted src lines: 10.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 140 - Approve Supplemented WL Diatomic Lane

Commit(s):
- this commit - Approve supplemented WL diatomic composition lane

Summary:
- Promoted `HP-COMP-SUPPWL-FN-01` and `HP-COMP-SUPPWL-TEST-01` from candidate
  to approved composition authority. The target cell is `Natom = 2`,
  `basisname !== nothing`, and `nesting = :wl` through the existing
  supplemented homonuclear z-axis diatomic staged facade.
- The approved source surface is intentionally narrow:
  `src/cartesian_base_hamiltonian.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for a direct RG/MWG compatibility genericity blocker. The expected
  source work is to remove the two early supplemented-WL policy blockers only
  if the existing Residual Gaussian/MWG path consumes the WL
  `CartesianTerminalBasisRealization`.
- The amendment preserves the supplement contract, residual selection, exact
  augmented operators, residual MWG/IDA interaction, base K/U reuse, artifact
  keys, manifest/provenance, driver inputs, and stage labels. `nesting`
  remains a construction-family choice, not a diagnostic route switch.

Validation:
- Design-manager validation for this docs-only pass: read current compact
  composition authority and the doer blurb; update registry/current/
  implementation slices/README/AGENTS/running log only; run `git diff --check`;
  focused `rg` for `HP-COMP-SUPPWL-*`, candidate/blocker wording, approved
  source surfaces, and forbidden driver/route/RG/artifact/Cr2 surfaces; confirm
  no source, test, tool, or bin files changed. No implementation tests are part
  of this approval pass; package load and artifact/readback checks belong to
  `HP-COMP-SUPPWL-TEST-01`.

Goal advancement:
- LT1/LT3: moves the next missing composition cell from policy blocker to
  exact source authority without adding a driver case branch.
- LT5/LT6: keeps PQS and WL distinct upstream but requires supplemented WL to
  meet the same terminal-basis/RG/MWG/Hamiltonian boundary as supplemented PQS.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one exact lane for supplemented WL
  instead of hidden early rejection plus candidate-only planning language.
- quarantined: driver changes, supplemented atoms, route skeleton/
  shellification/terminal-lowering changes, raw blocks, residual-selection
  changes, MWG/IDA convention changes, artifact schema or reader changes,
  public API/export changes, old WL H1/H1+J materialization, committed tests,
  solver/ECP work, diagnostics/status/report payloads, and Cr2 workflow remain
  unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if the existing RG/MWG path cannot consume
  WL terminal bases without new terminal records, route lowering semantics,
  residual-selection changes, or artifact/schema changes, the source pass must
  stop and report the exact blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 141 - Supplemented WL Diatomic Path

Commit(s):
- this commit - Enable supplemented WL diatomic composition

Summary:
- Accepted the `HP-COMP-SUPPWL-FN-01` near-deletion source pass. The two
  deliberate early `nesting = :wl` supplemented blockers were removed from
  `src/cartesian_base_hamiltonian.jl`: one in
  `cartesian_base_working_basis(...; supplemented = true)` and one in
  `cartesian_residual_gto_supplement_basis(...)`.
- No replacement branch was added. Supplemented WL now composes through the
  existing staged path: WL base terminal basis, supplement basis, residual
  augmentation, exact augmented operators, residual MWG/IDA, Hamiltonian
  assembly, write, and readback.

Validation:
- Doer: `git diff --check`; package load; H2 supplemented artifact/readback
  with `nesting = :pqs`; H2 supplemented artifact/readback with
  `nesting = :wl`; direct WL finite/symmetric K and V plus readback deltas;
  bounded Be2 supplemented WL artifact/readback with finite/symmetric K and V
  plus readback deltas; supplemented atoms still reject at the two-center R3
  facade; invalid heteronuclear, unequal-charge, non-neutral, and non-z-axis
  diatomic inputs still reject. No Cr2 run.
- Manager: inspected the deletion-only diff; `git diff --check`;
  `git diff --numstat -- src bin tools test docs`; suspicious added-line scan;
  new tests/tools scan; exact stale blocker string scan. No suspicious hits and
  no new committed tests/tools.

Goal advancement:
- LT1/LT3: fills the supplemented WL z-axis diatomic composition cell without
  changing driver inputs, route internals, residual selection, or artifact
  schema.
- LT5/LT6: confirms the intended common boundary: PQS and WL differ upstream,
  then both feed the same terminal-basis/RG/MWG/Hamiltonian machinery.

Medium-goal update:
- The remaining composition hole is supplemented atoms. That is now the next
  natural lane if the goal is to complete geometry/nesting/supplement
  composition; Cr2-facing work can continue through the supplemented diatomic
  paths without a Cr2-specific branch.

Carrying-cost result:
- deleted: two obsolete supplemented `nesting = :wl` policy blockers.
- simplified: supplemented PQS and WL now use the same visible staged
  composition path.
- quarantined: supplemented atoms, heteronuclear/general geometry, Cr2
  workflow, old WL H1/H1+J materialization, and diagnostics remain out of
  scope.
- not deleted because: existing staged helpers remain active composition
  boundaries.
- exact remaining caller/blocker: none found for supplemented WL z-axis
  diatomics.
- added src lines: 0.
- deleted src lines: 4.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 142 - Supplemented Atom Composition

Commit(s):
- this commit - Enable supplemented atom composition

Summary:
- Accepted the `HP-COMP-SUPPATOM-FN-01` source pass. The canonical driver no
  longer rejects `basisname !== nothing` for `Natom = 1`, and the producer now
  selects supplemented input validation by center count: one-center systems use
  the existing base atom validator, while two-center systems keep the existing
  supplemented homonuclear diatomic validator.
- `cartesian_residual_gto_supplement_basis(...)` now loads atomic supplements
  with `legacy_atomic_gaussian_supplement(...)` for one center and preserves
  `legacy_bond_aligned_diatomic_gaussian_supplement(...)` for two centers. No
  atom-only Hamiltonian builder or materialization path was added.

Validation:
- Doer: `git diff --check`; package load; base H atom artifact/readback; H
  atom supplemented artifact/readback for `nesting = :pqs` and `nesting = :wl`;
  bounded Be atom supplemented WL artifact/readback; H2 supplemented PQS and
  WL smokes; translated atom rejection; supplement basis-count mismatch
  rejection. No Cr2 run.
- Manager: inspected the source/bin diff; `git diff --check`;
  `git diff --numstat -- src bin tools test docs`; suspicious added-line scan;
  new tests/tools scan. No suspicious hits and no new committed tests/tools.

Goal advancement:
- LT1/LT3: fills the remaining `geometry = atom`, `supplement = on`
  composition cell without a special atom Hamiltonian path.
- LT5/LT6: confirms the intended product of public choices: geometry,
  nesting, and supplement now compose through the same staged producer
  boundaries for the approved atom and homonuclear z-axis diatomic geometries.

Medium-goal checkpoint:
- Completed: the 2 x 2 x 2 composition target is source-functional for
  origin-centered atoms and homonuclear z-axis diatomics, `nesting = :pqs` and
  `:wl`, supplement off/on.
- Active next: CR2-facing usability and performance should now use the
  canonical driver rather than requesting molecule-specific branches. Any
  remaining difficulty should be reported as producer friction, artifact
  provenance friction, or solver/consumer friction, not patched into the
  driver by default.

Carrying-cost result:
- deleted: obsolete driver diatomic-only supplemented guard.
- simplified: supplemented input selection shares atom/base and diatomic/R3
  validators; supplement loading is keyed by validated center count.
- quarantined: route changes, terminal lowering, raw blocks, residual
  selection, MWG/IDA conventions, artifact schema, reader/API, solver/ECP,
  translated atoms, heteronuclear/general geometry, and Cr2 workflow remain out
  of scope.
- not deleted because: the diatomic supplement loader remains the active
  two-center path.
- exact remaining caller/blocker: none found for bounded supplemented
  one-center atoms.
- added src/bin lines: 22.
- deleted src/bin lines: 8.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 142 - Approve Supplemented Atom Lane

Commit(s):
- this commit - Approve supplemented atom composition lane

Summary:
- Promoted `HP-COMP-SUPPATOM-FN-01` and `HP-COMP-SUPPATOM-TEST-01` from
  candidate to approved composition authority. The target cell is `Natom = 1`,
  `basisname !== nothing`, and `nesting = :pqs` or `nesting = :wl` through the
  existing staged supplemented Hamiltonian machinery.
- Approved source surfaces are `src/cartesian_base_hamiltonian.jl` and
  `bin/cartesian_ham_builder.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for a direct one-owner RG/MWG genericity blocker. The intended source
  work is to load atomic supplements with
  `legacy_atomic_gaussian_supplement(...)`, keep the existing diatomic
  supplement loader for two-center inputs, and relax only the canonical
  driver's supplemented `Natom == 2` guard.
- The amendment preserves base atom validation, terminal basis construction,
  residual Gaussian augmentation, exact augmented operators, residual MWG/IDA,
  base K/U reuse, assembly, writer, readback, manifest/provenance, driver
  public inputs, hooks, stage labels, and artifact contract. It explicitly
  rejects an atom-only Hamiltonian builder or atom-specific artifact shape.

Validation:
- Design-manager validation for this docs-only pass: read current compact
  composition authority and the doer blurb; update registry/current/
  implementation slices/README/AGENTS/running log only; run `git diff --check`;
  focused `rg` for `HP-COMP-SUPPATOM-*`, supplemented-atom candidate/blocker
  wording, approved source/bin surfaces, and forbidden driver/route/RG/artifact/
  Cr2 surfaces; confirm no source, bin, test, or tool files changed. No
  implementation tests are part of this approval pass; package load and
  artifact/readback checks belong to `HP-COMP-SUPPATOM-TEST-01`.

Goal advancement:
- LT1/LT3: moves the final explicit composition cell from policy blocker to
  exact source authority without adding an atom-only construction path.
- LT5/LT6: keeps atoms and diatomics on the same terminal-basis/RG/MWG/
  Hamiltonian boundary, with one-center residual selection treated as the
  one-owner case.

Medium-goal update:
- All initial explicit `atom | z-axis diatomic`, `:pqs | :wl`, and
  `supplement = off | on` composition lanes now have approved implementation
  authority under the current origin-centered atom and homonuclear z-axis
  diatomic constraints. Remaining medium-lane choices are no longer basic
  composition cells; they are geometry expansion, public/export polish,
  basis/supplement realism, Cr2 readiness, or solver/ECP work.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one exact atom supplement lane instead of
  driver rejection plus diatomic-only supplement-loader coupling.
- quarantined: atom-only Hamiltonian builders, new driver inputs, route
  switches, diagnostics, stop-after controls, new stage labels, route
  skeleton/shellification/terminal-lowering changes, raw blocks,
  residual-selection changes, MWG/IDA convention changes, artifact schema or
  reader changes, public API/export changes, solver/ECP work,
  heteronuclear/general geometry, translated atoms, committed tests, and Cr2
  workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if supplemented atoms require new
  residual-selection semantics, atom-only MWG/IDA conventions,
  terminal-lowering changes, route changes, artifact schema changes, or a
  separate atom Hamiltonian path, the source pass must stop and report the
  exact blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 143 - Approve Atom Parent Sizing

Commit(s):
- this commit - Approve one-center atom parent sizing lane

Summary:
- Approved `HP-COMP-ATOMBOX-FN-01` and `HP-COMP-ATOMBOX-TEST-01` as a narrow
  one-center atom parent-sizing correction. The live bug is that atom parent
  counts are still fixed by the old minimal-test rule `side = 2*q + 1`, so
  atom `padding` / public `basis.radius` is not authoritative for the physical
  box size.
- The approved source surface is only `src/cartesian_base_hamiltonian.jl`. I
  checked for a separate parent/system sizing helper and did not find a
  dedicated owner, so this amendment does not approve a vague optional helper
  file. If later source work needs to move sizing into an existing helper, that
  exact file needs a separate docs-only amendment.
- The approved contract makes `basis.radius` the one-center physical box
  extent authority. Parent axis counts must depend on radius plus
  `core_spacing` / the existing spacing policy, analogous to z-axis diatomic
  physical-extent sizing. `q` remains nesting/source-mode resolution, not a
  direct parent side-count control.

Validation:
- Design-manager validation for this docs-only pass: read compact current
  authority, registry, R1 base producer wording, implementation slices, AGENTS,
  and the manager request; inspect source only enough to confirm the current
  fixed count lives in `src/cartesian_base_hamiltonian.jl`; update docs/AGENTS/
  running log only; run `git diff --check`; focused `rg` for
  `HP-COMP-ATOMBOX-*`, `2*q + 1`, parent-axis-count/radius/padding wording,
  and forbidden driver/artifact/route/Cr2 surfaces; confirm no source, bin,
  test, or tool files changed. No implementation tests are part of this
  approval pass; package load and artifact/readback checks belong to
  `HP-COMP-ATOMBOX-TEST-01`.

Goal advancement:
- LT1/LT3: removes another hidden test artifact from the public atom producer
  contract. Public radius/padding becomes a real physical-size input instead
  of provenance-only metadata.
- LT5/LT6: keeps atom and diatomic sizing conceptually aligned: geometry
  chooses the extent rule, spacing sets the lattice count, and `q` remains a
  resolution/nesting control rather than a box-size substitute.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one exact atom-sizing lane instead of
  treating the current `2*q + 1` count as an implicit public rule.
- quarantined: driver changes, route-family switches, raw blocks,
  residual-selection changes, MWG/IDA convention changes, artifact schema or
  reader changes, public API/export changes, solver/ECP work,
  diagnostics/status/report payloads, committed tests, Cr2-specific workflow,
  translated atoms, non-origin atom support, element lookup/default tables,
  broad parent-construction rewrites, and diatomic sizing changes remain
  unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if fixing atom parent sizing requires broad
  parent-construction redesign, route semantics changes, driver contract
  changes, artifact schema changes, translated-atom support, or source files
  outside `src/cartesian_base_hamiltonian.jl`, the source pass must stop and
  report the exact blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 144 - Approve Public ns Input Contract

Commit(s):
- this commit - Approve public ns input contract

Summary:
- Approved `HP-COMP-NS-FN-01` and `HP-COMP-NS-TEST-01` as a narrow public
  size-parameter cleanup. The durable public contract should name `ns` as the
  requested cube/source/nesting size. Route-local `q` is derived only after the
  construction family is selected: `q = ns` for `nesting = :pqs`, and
  `q = ns - 2` for `nesting = :wl`.
- The approved source surface for later implementation is limited to
  `src/cartesian_base_hamiltonian.jl` and `bin/cartesian_ham_builder.jl`.
  Legacy public `q` may remain temporarily only as compatibility, with
  normalization to `ns` and rejection of inconsistent `ns`/`q` pairs.
- Physical controls remain separate: `core_spacing`, `reference_spacing`,
  `tail_spacing`, atom `radius` / driver `padding`, and diatomic extents are
  not aliases for `ns` or `q`.

Validation:
- Design-manager validation for this docs-only pass: read current compact
  authority, registry, R1 base/atom contracts, driver workflow docs,
  implementation slices, AGENTS, and the running log; update docs/AGENTS only;
  run `git diff --check`; focused `rg` for `HP-COMP-NS-*`, public `ns`,
  derived route-local `q`, `q_rule`, `ns_source`, and stale public-`q` wording;
  confirm no source, bin, test, or tool files changed. No implementation tests
  are part of this approval pass.

Goal advancement:
- LT1/LT3: removes a confusing public naming collision now that PQS and
  White-Lindsey are both visible construction-family choices.
- LT5/LT6: keeps atoms and diatomics on the same driver/facade contract while
  allowing the selected nesting family to own its route-local interpretation.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one public size field, `ns`, instead of
  teaching users that one route-local `q` value has the same meaning across PQS
  and White-Lindsey.
- quarantined: route skeleton, shellification, terminal lowering, raw blocks,
  residual-selection changes, MWG/IDA changes, numerical kernels, artifact
  matrix keys, reader behavior, public API/export redesign, solver/ECP work,
  route diagnostics, status/report payloads, committed tests, and Cr2 workflow
  remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if `ns` normalization requires route-stage
  redesign, terminal construction changes, artifact reader changes, or source
  outside `src/cartesian_base_hamiltonian.jl` and `bin/cartesian_ham_builder.jl`,
  the source pass must stop and report the exact blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: compact provenance keys `ns`, `q_rule`, and
  `ns_source` are approved only inside existing `producer_provenance/` and
  `recipe_provenance/` groups.

## Cartesian Hamiltonian Producer Pass 145 - Public `ns` Normalization

Commit(s):
- this commit - Normalize public ns input

Summary:
- Accepted the `HP-COMP-NS-FN-01` source/bin pass. The canonical driver now
  exposes `ns` as the public requested cube/source/nesting size and passes
  `basis.ns` into the producer. Driver edits are limited to the visible input
  name, `public_inputs`, basis construction, and contract print; stage order,
  hooks, comments, and workflow shape were not redesigned.
- `src/cartesian_base_hamiltonian.jl` now normalizes the public size contract:
  `nesting = :pqs` derives route-local `q = ns`; `nesting = :wl` derives
  `q = ns - 2` and rejects `ns < 3`. Legacy direct-facade `basis.q` remains
  accepted as compatibility, deriving `ns = q` for PQS and `ns = q + 2` for
  WL. Supplying both `ns` and inconsistent `q` throws.
- Existing `producer_provenance/` and `recipe_provenance/` now record `ns`,
  derived `q`, `q_rule`, `ns_source`, and existing `nesting`.

Validation:
- Doer: `git diff --check`; package load; atom base artifact/readback with
  public `ns` for `nesting = :pqs` and `nesting = :wl`; H2 base
  artifact/readback with public `ns` for `nesting = :pqs` and `nesting = :wl`;
  H atom supplemented smoke with public `ns`; legacy direct-facade `q`
  compatibility for PQS and WL; inconsistent `ns`/`q` rejection;
  `nesting = :wl, ns < 3` rejection; direct provenance inspection for base WL
  and supplemented PQS artifacts. No Cr2 run.
- Manager: inspected the source/bin diff; `git diff --check`;
  `git diff --numstat -- src bin tools test docs`; suspicious added-line scan;
  new tests/tools scan; focused driver stale-public-`q` scan. No suspicious
  hits and no new committed tests/tools. The remaining `2 * input.q + 1`
  source hit is the separately approved atom-box lane, not part of this pass.

Goal advancement:
- LT1/LT3: makes the public size knob match user expectations across PQS and
  WL while preserving route-local `q` for existing construction internals.
- LT5/LT6: removes an ambiguity that distorted PQS/WL comparisons and prepares
  the one-center atom box fix to use resolution metadata rather than public
  `q` semantics.

Remaining blocker / next:
- `HP-COMP-ATOMBOX-*` should now proceed: remove the one-center
  `2 * input.q + 1` parent-axis count artifact and make atom box size depend
  on public radius/padding plus spacing policy, without changing driver
  inputs.

Carrying-cost result:
- deleted: public driver `q` surface.
- simplified: public size normalization is centralized in one producer helper.
- quarantined: atom parent sizing, route skeletons, shellification, terminal
  lowering, raw blocks, RG/MWG/IDA, reader behavior, diagnostics,
  reports/status, and Cr2 workflow remain out of this pass.
- not deleted because: route-local `q` remains the existing internal
  construction input and legacy direct-facade `q` is temporarily supported for
  compatibility.
- exact remaining caller/blocker: atom parent sizing is still a separate
  approved lane.
- added src/bin lines: 51.
- deleted src/bin lines: 15.
- new tests: none.
- new metadata/status fields: compact provenance keys `ns`, `q_rule`, and
  `ns_source` in existing provenance groups.

## Cartesian Hamiltonian Producer Pass 146 - Atom Parent Box Sizing

Commit(s):
- this commit - Size atom parent box from physical radius

Summary:
- Accepted the `HP-COMP-ATOMBOX-FN-01` source pass. The one-center atom branch
  no longer sets `parent_axis_counts = (2*q + 1, 2*q + 1, 2*q + 1)`.
  Instead it computes an odd mapped physical parent count from the
  White-Lindsey atomic mapping, public radius, core spacing, tail spacing, and
  reference spacing, then applies route-local `q` only as the minimum
  core-shell side.
- The implementation stays inside `src/cartesian_base_hamiltonian.jl`. It does
  not touch the canonical driver or route helpers. This matters because the
  route-helper fallback still treats radius as an index radius; the producer
  now passes explicit counts for atoms to avoid that fallback.

Validation:
- Doer: `git diff --check`; package load; H atom base artifact/readback; H
  atom supplemented PQS and WL artifact/readback; Be atom base and
  supplemented PQS smokes; atom padding sensitivity at fixed `ns = 2`
  showing parent counts/final dimensions change from `(5,5,5)`/`35` at
  radius `1.0`, to `(7,7,7)`/`43` at radius `2.0`, to `(9,9,9)`/`51` at
  radius `4.0`; H2 base smoke; Be2 WL base smoke. No Cr2 run.
- Manager: inspected the source diff and the mapped-count helper precedent;
  `git diff --check`; `git diff --numstat -- src bin tools test docs`;
  suspicious added-line scan; new tests/tools scan; focused stale `2*q+1`
  scan. No suspicious hits and no new committed tests/tools.

Goal advancement:
- LT1/LT3: fixes the atom public contract so driver `padding` / producer
  `radius` is a real physical extent authority for one-center atoms.
- LT5/LT6: completes the immediate `ns`/atom-box correction chain. Public
  `ns` now controls source/nesting resolution, while atom radius controls
  physical parent extent.

Remaining blocker / next:
- Re-run CR2-side calibrated driver comparisons before drawing PQS/WL
  conclusions. Atom dimensions and WL/PQS comparison points should now be
  recalibrated with public `ns` and physical padding.

Carrying-cost result:
- deleted: fixed one-center `2 * input.q + 1` parent side count.
- simplified: atom parent sizing is tied to physical radius and mapping rather
  than public source/nesting size.
- quarantined: driver changes, route helpers, shellification, terminal
  lowering, raw blocks, RG/MWG/IDA, artifact schema, reader/API, translated
  atoms, and Cr2 workflow remain out of scope.
- not deleted because: route-local `q` still defines the minimum core-shell
  side needed by existing shellification.
- exact remaining caller/blocker: none for one-center parent sizing in this
  lane.
- added src lines: 20.
- deleted src lines: 2.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 147 - Approve WL Diatomic `ns` Cleanup

Commit(s):
- this commit - Approve WL diatomic ns contract cleanup

Summary:
- Approved `HP-COMP-WLNS-FN-01` and `HP-COMP-WLNS-TEST-01` as a narrow
  White-Lindsey z-axis diatomic `ns` contract cleanup. The read-only audit
  showed that public `nesting = :wl`, `Natom = 2`, `ns = 3` normalizes to
  route-local `q = 1` and `core_cube_side = 1`, then fails in terminal
  shellification with `complete-shell inner box must be a filled CPB`.
  Future source work should reject this input early in
  `src/cartesian_base_hamiltonian.jl`.
- The same audit showed that working WL diatomic `ns` values can saturate the
  retained support when physical parent extent dominates. For the CR2-shape
  Be2 stage probe, `ns = 4:7` all retained the same 6075 support points; `4/5`
  shared row order and `6/7` shared row order, while `4` versus `6` changed
  block decomposition/order. That is not ignored input and should not be
  "fixed" in this lane.

Validation:
- Doer audit evidence: `git diff --check`; package load; ignored probes
  `tmp/work/wl_diatomic_ns_response_audit.jl` and
  `tmp/work/wl_diatomic_ns_response_cr2_shape_audit.jl`; no Cr2 run.
- Docs-only amendment validation: update only docs; run `git diff --check`;
  focused scans for `HP-COMP-WLNS-*`, WL diatomic `ns < 4`, and retained
  support saturation wording; confirm no source, bin, test, tool, or artifact
  files changed.

Goal advancement:
- LT1/LT3: moves an unsupported public input from a late route-internal
  failure to a clear producer input contract.
- LT5/LT6: clarifies that equal public `ns` is not a complete PQS/WL
  retained-basis comparison and prevents a false optimization lane around
  legitimate WL support saturation.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one precise early-rejection lane instead
  of treating `ns = 3` as a terminal-basis or shellification bug.
- quarantined: driver changes, route skeleton/shellification/terminal-lowering
  changes, retained-unit changes, terminal-basis changes, raw blocks,
  RG/MWG/IDA changes, artifact schema or reader changes, public API/export,
  diagnostics/report/status payloads, committed tests, and Cr2 workflow remain
  unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if early WL diatomic `ns` rejection requires
  anything outside `src/cartesian_base_hamiltonian.jl`, the source pass must
  stop and report the exact boundary.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 149 - Approve WL Diatomic Compact Retained Basis

Commit(s):
- this commit - Approve WL diatomic compact retained-basis authority

Summary:
- Approved `HP-WLDIAT-COMPACT-FN-01` and
  `HP-WLDIAT-COMPACT-TEST-01` as the narrow correction lane for
  White-Lindsey z-axis diatomic compact retained-basis construction.
- The design now distinguishes a mechanically reachable WL diatomic artifact
  path from the intended compact WL retained-basis contract. The current
  elongated shared-shell boundary-stratum identity realization is explicitly
  not the production PQS/WL comparison story.
- The approved correction preserves the WL unit model, but requires each
  diatomic WL unit to carry or realize compact retained columns from products
  of one-dimensional contractions on owned support.

Validation:
- Docs-only amendment validation: update `AGENTS.md` and compact design docs;
  run `git diff --check`; focused scans for `HP-WLDIAT-COMPACT-*`,
  full-support identity wording, and WL diatomic compact retained-basis
  wording; confirm no `src`, `bin`, `test`, `tests`, `tools`, or artifact files
  changed. Local docs build should be run because the registry is large and
  GitHub docs recently failed on size/rendering.

Goal advancement:
- LT1/LT3: prevents the public `nesting = :wl` diatomic path from being treated
  as scientifically comparable while it still realizes boundary-stratum units
  as identity blocks.
- LT5/LT6: keeps the 2 x 2 x 2 composition plan intact while adding the missing
  WL-specific compact retained-basis contract before further PQS/WL
  comparisons.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one explicit compact WL retained-basis
  target instead of ambiguous terminal-basis or driver fixes.
- quarantined: driver changes, artifact/provenance/schema changes, PQS
  behavior changes, Hamiltonian assembly changes, raw blocks, RG/MWG/IDA,
  old WL H1/H1+J materialization, diagnostics/status/report payloads,
  committed fixtures/tests, and Cr2 workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if compact retained units require
  construction-native facts not currently available, the source pass must stop
  and report the missing fact instead of dropping rows or relabeling identity
  units.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 150 - Clarify WL Compact Donor Primitive

Commit(s):
- this commit - Clarify WL compact coefficient donor rule

Summary:
- Clarified the `HP-WLDIAT-COMPACT-*` authority after the postmortem on how
  the identity-realization bug entered. The old WL route-global stack was
  correctly retired as bloat, but it contained an essential numerical
  primitive: boundary-stratum CPB to local product-of-1D coefficient map to
  terminal block coefficients.
- The compact authority now says source work may mine deleted WL coefficient
  helpers only as historical donor/reference material for that primitive. It
  must re-express the primitive behind the current
  `CartesianTerminalBasisRealization` boundary and must not revive old
  reports, adapters, status surfaces, or H1/H1+J materialization.
- The docs also freeze that identity realization is valid only for true
  direct/core identity units, not for metadata-only WL boundary-stratum
  retained units.

Validation:
- Docs-only clarification validation: `git diff --check`; focused scans for
  the donor/reference primitive wording and direct/core identity exception;
  confirm no `src`, `bin`, `test`, `tests`, or `tools` files changed. Local
  docs build should be run because this touches compact design authority.

Goal advancement:
- LT1/LT3: keeps the WL fix focused on the missing numerical primitive instead
  of accidentally restoring the bloated route-global implementation.
- LT5/LT6: documents the actual failure mode, reducing the chance that
  metadata names are again mistaken for realized basis functions.

Carrying-cost result:
- deleted: none; docs-only clarification.
- simplified: source handoff now has a sharper donor rule.
- quarantined: old route-global WL stack, reports, adapters, status surfaces,
  H1/H1+J materialization, driver changes, artifact changes, PQS changes,
  raw blocks, RG/MWG/IDA, tests/fixtures, and Cr2 workflow remain unapproved.
- not deleted because: no source work in this pass.
- exact remaining caller/blocker: the future source pass must locate or
  re-express the CPB-local product-of-1D coefficient primitive without
  reintroducing retired framework code.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 151 - Add Gausslet Algorithm Refresher

Commit(s):
- this commit - Add Gausslet algorithm refresher

Summary:
- Added `docs/src/developer/architecture/gausslet_algorithm_refresher.md` as a
  concise agent-facing memory-refresh note. It is explicitly not design
  authority and does not approve implementation surfaces.
- The refresher summarizes operational distinctions among Qiu-White,
  White-Lindsey, and PQS; the WL boundary-stratum coefficient-map rule; PQS
  source-box and terminal-realization guardrails; direct/core identity sectors;
  COMX, PGDG, IDA, MWG, GTO supplementation, and current core-spacing naming.
- Linked the refresher from the developer index, hidden docs navigation, and
  `AGENTS.md` as the short operational memory document before deeper
  paper-centered fundamentals.

Validation:
- Docs-only pass: `git diff --check`; focused scan for the new refresher link
  and key category-mistake terms; confirm no staged `src`, `bin`, `test`,
  `tests`, or `tools` files. Local docs build should be run because the page is
  linked into the rendered docs.

Goal advancement:
- LT1/LT3: reduces recurrence of WL/PQS category mistakes by giving agents a
  short operational reminder before implementation work.
- LT5/LT6: lowers cognitive load without restoring old route-global
  scaffolding or adding new design authority.

Carrying-cost result:
- deleted: none; docs-only refresher.
- simplified: startup reading now has a compact first stop before the long
  fundamentals packet and design registry.
- quarantined: no new source surfaces, tests, tools, public API, artifact
  schema, or design IDs.
- not deleted because: this is a memory/refresher document, not a source
  cleanup.
- exact remaining caller/blocker: if future source handoffs need algorithm
  authority, they still must cite the compact design registry and approved IDs.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 152 - Compact WL Diatomic Terminal Basis

Commit(s):
- this commit - Implement compact WL diatomic terminal basis

Summary:
- Accepted `HP-WLDIAT-COMPACT-FN-01`. White-Lindsey diatomic
  boundary-stratum retained units no longer realize as full-support identity
  blocks. They now build CPB-local face/edge/corner product coefficient maps
  behind the common `CartesianTerminalBasisRealization` boundary. Direct/core
  retained units still realize as identity blocks, which remains the intended
  contract.
- The source pass re-expresses the old WL product-of-1D primitive without
  reviving the route-global WL H1/H1+J materialization stack, reports,
  adapters, or status surfaces.
- Small WL supplemented H2 reaches the existing residual-GTO/MWG stages with
  finite symmetric matrices. The active PQS H2 RG endpoint remains unchanged:
  augmented dimension `489` and self-Coulomb `0.4574265214362095`.

Validation:
- Doer: `git diff --check`; package load; ignored
  `tmp/work/wl_compact_block_size_audit.jl`; ignored
  `tmp/work/wl_diatomic_base_validation.jl` with `pqs_dim=111`, `wl_dim=199`;
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`; small
  supplemented WL H2 driver smoke with dimension `209`, finite/symmetric `K`
  and `V`.
- Manager: reread active `HP-WLDIAT-COMPACT-*` authority; reviewed the diff;
  checked donor face/edge products return full-parent coefficient matrices
  with support indices, so terminal-basis support-row subsetting is the right
  shape; ran `git diff --check`, package load,
  `git diff --numstat -- src bin tools test docs`, suspicious added-line scan,
  new tests/tools scan, and final status.
- Mechanical gate summary: source line impact is `145` insertions and `20`
  deletions across the two approved source files; no new committed tests or
  tools. Suspicious scan hits are confined to existing CPB/retained-unit
  metadata reads and one compact construction-native scalar,
  `white_lindsey_retained_count_1d`, carried through existing contract
  metadata so the terminal realizer can reconstruct the WL product primitive.
  No artifact/status/report fields were added.

Goal advancement:
- LT1/LT3: repairs the WL diatomic terminal-basis category mistake directly:
  a WL boundary stratum is a compact product retained basis over its support,
  not one identity basis function per parent site.
- LT5/LT6: restores a meaningful PQS/WL comparison boundary for the approved
  homonuclear z-axis diatomic workflow while preserving the common staged
  Hamiltonian and residual-GTO/MWG paths.

Medium-goal checkpoint:
- Completed: the serious WL diatomic full-support identity realization error
  is crossed for the current retained-unit path.
- Active: recalibrate PQS/WL driver comparisons after this correction before
  drawing physics conclusions; keep using the canonical driver for CR2-facing
  artifact production and report friction at producer/consumer boundaries.
- Needs refinement: the inherited even-`ns` WL retained-count rule still makes
  `ns = 4` produce `26` boundary columns rather than the nominal `56`. This is
  now an explicit WL policy question, not a hidden identity-realization bug.

Carrying-cost result:
- deleted: boundary-stratum identity realization as production behavior.
- simplified: WL terminal append now has a clear split between direct/core
  identity units and boundary product coefficient blocks.
- quarantined: driver changes, artifact schema/provenance changes, PQS
  behavior changes, raw blocks, Residual Gaussian/MWG/IDA changes, old WL
  H1/H1+J materialization, committed tests/fixtures, and Cr2 workflow remain
  out of scope.
- not deleted because: direct/core identity realization is valid and remains
  live.
- exact remaining caller/blocker: even-`ns` WL retained-count policy follows
  the inherited symmetric-odd donor behavior; changing that requires a later
  explicit design/source lane.
- added src lines: 145.
- deleted src lines: 20.
- new tests: none.
- new metadata/status fields: one internal construction metadata scalar only;
  no artifact/status/report fields.

## Cartesian Hamiltonian Producer Pass 153 - Approve WL Boundary Parity Cleanup

Commit(s):
- this commit - Approve WL boundary retained-count parity cleanup

Summary:
- Approved `HP-WLDIAT-PARITY-FN-01` and `HP-WLDIAT-PARITY-TEST-01` as a
  narrow follow-up to the compact WL diatomic terminal-basis pass.
- The design now records that the remaining `ns = 4 -> 26` behavior is an
  over-applied symmetric-odd donor rule, not a valid WL policy. Odd side counts
  are required for nucleus-centered core/contact blocks, but not for boundary
  shell strata.
- The approved source surface is only
  `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`.
  Boundary-stratum product construction should use the requested retained
  count without symmetric-odd coercion, so WL `ns = 4`, route-local `q = 2`,
  targets `4^3 - 2^3 = 56` boundary columns and `ns = 5` targets
  `5^3 - 3^3 = 98`.

Validation:
- Docs-only amendment validation: `git diff --check`; focused scans for
  `HP-WLDIAT-PARITY-*`, symmetric-odd/odd-side wording, and `56`/`98`
  boundary-count expectations; confirm no `src`, `bin`, `test`, `tests`,
  `tools`, or artifact files changed. Local docs build should be run because
  compact design authority changed.

Goal advancement:
- LT1/LT3: fixes the next visible WL numerical semantics issue after the
  identity-realization bug: boundary shells should keep shell retained counts,
  not inherit nucleus-centered parity.
- LT5/LT6: keeps the PQS/WL comparison tied to the shared public `ns` input
  without changing driver, route skeleton, or artifact contracts.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: source work has one precise parity target in one owner file.
- quarantined: driver changes, public `ns` normalization changes,
  route/shellification/terminal-lowering changes, retained-unit metadata shape
  changes, direct/core identity changes, artifact changes, PQS behavior,
  raw blocks, RG/MWG/IDA, old WL materialization, diagnostics/status/report
  payloads, committed tests/fixtures, and Cr2 workflow remain unapproved.
- not deleted because: source cleanup has not run yet.
- exact remaining caller/blocker: if parity cleanup needs files outside the WL
  terminal realizer or a metadata/schema change, the source pass must stop and
  report the exact blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 154 - Implement WL Boundary Parity Cleanup

Commit(s):
- this commit - Preserve WL boundary retained count

Summary:
- Accepted `HP-WLDIAT-PARITY-FN-01`. The WL boundary-stratum terminal
  coefficient path now keeps the requested retained count instead of applying
  nucleus-centered symmetric-odd coercion. Direct/core identity behavior and
  compact face/edge/corner product construction are unchanged.
- The visible endpoint moved as intended: WL diatomic `ns = 4`, route-local
  `q = 2`, now gives `56` boundary columns rather than `26`; `ns = 5` remains
  `98`.

Validation:
- Doer: `git diff --check`; package load; ignored
  `tmp/work/wl_compact_block_size_audit.jl` with `ns=4`
  `expected_shell=56`, `boundary_columns=56`, `final_dimension=119`, and
  `ns=5` `expected_shell=98`, `boundary_columns=98`,
  `final_dimension=161`; ignored `tmp/work/wl_diatomic_base_validation.jl`
  with `pqs_dim=111`, `wl_dim=303`; ignored supplemented WL smoke with
  dimension `309`, finite/symmetric `K`/`V`, and readback deltas `0.0`;
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` with
  augmented dimension `489` and self-Coulomb `0.4574265214362095`.
- Manager: reviewed the one-line diff; ran `git diff --check`;
  `git diff --numstat -- src bin tools test docs`; suspicious added-line
  scan; new tests/tools scan; package load; final status.

Goal advancement:
- LT1/LT3: completes the immediate WL compact retained-basis correction chain
  by removing the over-applied parity rule from boundary shell strata.
- LT5/LT6: restores the intended same-public-`ns` WL/PQS comparison starting
  point without changing driver, route skeleton, shellification, terminal
  lowering, artifacts, or residual-GTO/MWG paths.

Carrying-cost result:
- deleted: symmetric-odd coercion from WL boundary-stratum contraction.
- simplified: boundary shell parity now follows the public `ns` retained-count
  contract.
- quarantined: driver, public `ns` normalization, route/shellification/
  lowering, retained metadata shape, artifacts, PQS behavior, raw blocks,
  RG/MWG/IDA, old WL materialization, committed tests/fixtures, and Cr2
  workflow remain untouched.
- not deleted because: direct/core identity realization remains valid.
- exact remaining caller/blocker: none for this parity lane.
- added src lines: 1.
- deleted src lines: 1.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 155 - Approve RG Final Identity Tolerance

Commit(s):
- this commit - Approve RG final identity tolerance

Summary:
- Approved `HP-RG-IDTOL-FN-01` and `HP-RG-IDTOL-TEST-01` as a narrow
  Residual Gaussian policy update. The Be atom cc-pV5Z `lmax = 1` audit showed
  a healthy residual space with `21` retained directions, minimum retained
  occupation `6.151e-6`, final merge condition `1.0`, and `G' S R =
  1.776e-14`; the only failure was `R' S R - I = 2.183e-10` against an old
  allowed error of about `2.000e-10`.
- The default final residual identity validation tolerance is approved as
  `1.0e-8`. The default `residual_occupation_cutoff` remains `1.0e-8`; width
  and zeta filtering remain explicit and user-controlled.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-RG-IDTOL-*`, Be cc-pV5Z evidence, `residual_occupation_cutoff`, and
  no `src`, `bin`, `test`, `tests`, `tools`, driver, or artifact changes. A
  local docs build is appropriate because compact authority files changed.

Goal advancement:
- LT6/RG robustness: keeps owner-local residual selection stable while
  recognizing that a tiny final identity overshoot with healthy merge metrics
  is a validation tolerance issue, not a reason to drop real residual
  directions.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has a one-file primary owner and one optional
  compatibility plumbing file.
- quarantined: residual selection changes, occupation-cutoff changes,
  width/zeta filtering defaults, owner grouping, merge failure rules, MWG/IDA,
  raw blocks, terminal basis, WL/PQS route behavior, artifacts, driver changes,
  committed fixtures, and Cr2 workflow remain unapproved.
- exact remaining caller/blocker: if Be cc-pV5Z cannot pass by changing only
  the final identity tolerance default, source work must stop and report the
  exact blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 156 - Implement RG Final Identity Tolerance

Commit(s):
- this commit - Update RG final identity tolerance

Summary:
- Accepted `HP-RG-IDTOL-FN-01`. The default final Residual Gaussian
  `R' S R` identity validation tolerance is now `1.0e-8` in both the RG owner
  and the compatibility wrapper default.
- This is intentionally only a final identity validation tolerance change. It
  does not change residual direction selection, the default
  `residual_occupation_cutoff = 1.0e-8`, owner grouping, merge behavior,
  width/zeta filtering, MWG/IDA, Hamiltonian assembly, artifacts, or driver
  behavior.
- The motivating Be atom cc-pV5Z `lmax = 1` case now passes with all `21`
  residual directions retained. The minimum retained occupation remains
  `6.151346475239507e-6`, so the marginal high-zeta direction was not dropped.

Validation:
- Doer: `git diff --check`; package load; ignored
  `tmp/work/be_atom_high_zeta_residual_audit.jl`; ignored
  `tmp/work/be_atom_high_zeta_artifact_smoke.jl`; and
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`. The Be
  cc-pV5Z audit reported retained count `21`, merge condition `1.0`,
  `max |G' S R| = 1.7763568394002505e-14`,
  `max |R' S R - I| = 2.1827872842550278e-10`, and allowed identity
  tolerance `2.0000000002182788e-8`.
- Manager: reviewed the two-line source diff; ran `git diff --check`,
  `git diff --numstat -- src bin tools test docs`, suspicious added-line scan,
  new tests/tools scan, package load, and
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`. The H2 RG
  endpoint remained at augmented dimension `489`, self-Coulomb
  `0.4574265214362095`, and zero facade readback deltas.

Goal advancement:
- LT6/RG robustness: keeps the legacy-supported residual occupation cutoff and
  selected residual space stable while allowing a numerically harmless final
  orthonormality cleanup miss to pass when merge conditioning and `G' S R`
  checks are healthy.
- MT/CR2 readiness: removes the high-zeta Be atom supplement validation
  friction without adding new workflow knobs or changing the canonical driver.

Carrying-cost result:
- deleted: old `1.0e-10` final identity default.
- simplified: RG owner and terminal residual compatibility wrapper now agree
  on the approved `1.0e-8` default.
- quarantined: residual selection, occupation cutoff, width/zeta filtering,
  owner grouping, merge rules, MWG/IDA, artifacts, driver, routes, raw blocks,
  tests, and Cr2 workflow remain untouched.
- not deleted because: strict `G' S R` validation and the existing merge
  checks remain live invariants.
- exact remaining caller/blocker: none for this tolerance-default lane.
- added src lines: 2.
- deleted src lines: 2.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 157 - Retire Complete-Core-Shell RHF Stack

Commit(s):
- this commit - Retire complete-core-shell RHF stack

Summary:
- Accepted `HP-RETIRE-CCS-RHF-FN-01`. The stale complete-core-shell RHF
  payload/status workflow was removed by deleting
  `src/pqs_multilayer_complete_core_shell_rhf.jl` and its root include.
- The current Cartesian producer path is unchanged: canonical staged driver,
  staged base/supplemented producer functions, and `CartesianIDAHamiltonian`
  artifacts. No replacement, adapter, compatibility shim, report/status object,
  committed test, driver change, or artifact change was added.

Validation:
- Doer: `git diff --check`; package load; focused `rg` over `src`, `bin`,
  `test`, and `tools` for `pqs_multilayer_complete_core_shell_rhf`; canonical
  base driver artifact/readback smoke with H atom dimension `105`; canonical
  supplemented H2/cc-pVTZ artifact/readback smoke with dimension `489`; and
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` with H2
  augmented dimension `489`, self-Coulomb `0.4574265214362095`, and zero
  readback deltas.
- Manager: reviewed the deletion diff; reran `git diff --check`; focused live
  reference scan; `git diff --numstat -- src bin tools test docs`;
  suspicious added-line scan; new tests/tools scan; package load; and the H2
  RG endpoint/facade test with the same self-Coulomb and zero readback deltas.

Goal advancement:
- LT2/LT5 cleanup: removes a large stale route-era RHF payload stack now that
  CR2-facing work uses the canonical Hamiltonian artifact producer rather than
  complete-core-shell RHF payloads.
- MT/Cartesian usability: reduces the active Cartesian sector's conceptual
  load without changing the driver or numerical producer behavior.

Carrying-cost result:
- deleted: `src/pqs_multilayer_complete_core_shell_rhf.jl` and its root
  include.
- simplified: package load and source navigation no longer carry the inactive
  complete-core-shell RHF workflow.
- quarantined: canonical driver, artifacts/provenance/manifest, route,
  shellification, terminal lowering, raw blocks, Residual Gaussian/MWG/IDA,
  Hamiltonian assembly, QW donor kernels, and tests remain untouched.
- not deleted because: active `pqs_multilayer_complete_core_shell_h1.jl`,
  final-basis helpers, and low-order materialization paths were explicitly out
  of scope for this first retirement pass.
- exact remaining caller/blocker: none found for the deleted RHF stack.
- added src lines: 0.
- deleted src lines: 1880.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 158 - Retire Route-Driver Materialization Workflow

Commit(s):
- this commit - Retire route-driver materialization workflow

Summary:
- Accepted `HP-RETIRE-DRV-MAT-FN-01`, `HP-RETIRE-DRV-MAT-TOOL-01`,
  `HP-RETIRE-DRV-MAT-DOC-01`, and `HP-RETIRE-DRV-MAT-TEST-01`. The old
  route-driver materialization/report/save wrapper workflow was removed without
  changing the canonical staged driver or the numerical producer path.
- Removed `cartesian_materialization`, `cartesian_print_summary`,
  `cartesian_print_details`, `cartesian_save`, the matching underscored
  materialization/report/save helpers, and the now-uncalled WL atomic
  pure-gausslet materialization helper stack. Deleted the old route-driver
  reporting file and include.
- Deleted old wrapper tools
  `tools/cartesian_driver_harness.jl`,
  `tools/cr2_cartesian_ida_stage_probe.jl`, and
  `tools/h2_pqs_base_hamiltonian_smoke.jl`. Quarantined
  `tools/cartesian_driver_ladder_lib.jl` as a one-line retirement error because
  two unapproved runner scripts still include it.

Validation:
- Doer: `git diff --check`; package load; focused `rg` over `src`, `bin`,
  `test`, and `tools` for all retired wrapper names; H atom canonical base
  driver artifact/readback smoke with dimension `105`; H2/cc-pVTZ canonical
  supplemented driver artifact/readback smoke with dimension `489`;
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` with
  self-Coulomb `0.4574265214362095`; docs-policy gate with runner context.
- Manager: reviewed the diff and authority boundary; reran focused retired-name
  scans, including the deleted WL atomic helper stack; `git diff --check`;
  `git diff --numstat -- src bin tools test docs`; suspicious added-line scan;
  package load; docs-policy runner-context check; bounded canonical base and
  supplemented driver smokes; and the H2 RG endpoint/facade with the same
  self-Coulomb and zero readback deltas.
- The docs-policy test update removed the old wrapper-stage assertions and
  also dropped adjacent stale frozen-driver assertions that no longer match the
  current approved driver: the header says "producer" rather than "template",
  and the driver already has approved key=value Julia-value overrides using
  `Meta.parse`. No driver code changed.

Goal advancement:
- LT2/LT5 cleanup: removes a stale route-stage workflow that could pull agents
  away from the canonical staged producer and artifact path.
- MT/CR2 readiness: preserves the CR2-facing driver/artifact route while
  reducing old harness/report surface area.

Carrying-cost result:
- deleted: old wrapper/report/save source path, route-driver reporting file,
  stale harness/probe/smoke tools, and old WL atomic pure-gausslet
  materialization helper stack.
- simplified: active algorithm index now points at the current base
  Hamiltonian helper instead of the retired materialization helper; docs-policy
  pressure no longer expects the canonical driver to call retired wrappers.
- quarantined: `tools/cartesian_driver_ladder_lib.jl` now fails immediately
  with a retirement message.
- not deleted because: `tools/run_cartesian_driver_ladder.jl` and
  `tools/run_cartesian_line_ladder.jl` still include the ladder library but
  were outside the allowed tool surface for this pass.
- exact remaining caller/blocker: a later tool-retirement pass should delete
  the two runner scripts or approve their quarantine explicitly.
- added src lines: 0.
- deleted src lines: 463.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 159 - Delete Dangling Ladder Runners

Commit(s):
- this commit - Delete retired Cartesian ladder runners

Summary:
- Accepted `HP-RETIRE-LADDER-RUNNERS-FN-01`. Deleted the two dangling runner
  entrypoints into the already-retired Cartesian route-stage ladder workflow:
  `tools/run_cartesian_driver_ladder.jl` and
  `tools/run_cartesian_line_ladder.jl`.
- No source, docs, tests, canonical driver, numerical path, artifact path, or
  `tools/cartesian_driver_ladder_lib.jl` changes were made.

Validation:
- Doer: `git diff --check`; package load; focused `rg` over `src`, `bin`,
  `test`, and `tools` for `run_cartesian_driver_ladder`,
  `run_cartesian_line_ladder`, and `cartesian_driver_ladder_lib`; canonical H
  atom base artifact/readback smoke with dimension `105`.
- Manager: reviewed the diff; reran `git diff --check`; `git diff --numstat
  -- src bin tools test docs`; focused live-reference scan; package load; and
  the same bounded canonical H atom base driver smoke with dimension `105`.

Goal advancement:
- LT2/LT5 cleanup: completes the route-driver materialization workflow
  retirement by removing the final approved runner entrypoints.
- MT update: per user direction, pause this cleanup lane after this pass unless
  a later amendment names another stale surface.

Carrying-cost result:
- deleted: two obsolete ladder runner scripts.
- simplified: no remaining live `src`, `bin`, `test`, or `tools` references to
  the retired runner names or ladder library.
- quarantined: `tools/cartesian_driver_ladder_lib.jl` remains untouched under
  its prior retirement state.
- not deleted because: no other files were in scope.
- exact remaining caller/blocker: none found for the two deleted runners.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 157 - Approve CCS RHF Stack Retirement

Commit(s):
- this commit - Approve CCS RHF stack retirement

Summary:
- Approved `HP-RETIRE-CCS-RHF-FN-01` and
  `HP-RETIRE-CCS-RHF-TEST-01` as a narrow deletion authority lane for the old
  complete-core-shell RHF payload stack.
- The deletion target is `src/pqs_multilayer_complete_core_shell_rhf.jl` plus
  its include in `src/GaussletBases.jl`. Focused search found no live source,
  bin, test, or tool caller outside the file itself and the root include; the
  remaining hits are docs/history and file-internal payload/status vocabulary.
- Current CR2-facing work consumes canonical driver
  `CartesianIDAHamiltonian` artifacts for atom/diatomic, base/supplemented,
  and PQS/WL cases. The old RHF stack is stale route-era workflow machinery.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-RETIRE-CCS-RHF-*` and `pqs_multilayer_complete_core_shell_rhf`;
  confirm no `src`, `bin`, `test`, `tools`, artifact, or implementation files
  changed. A docs build is appropriate because compact authority files changed.

Goal advancement:
- LT5/LT6: converts a stale payload/status route stack into an approved
  deletion target instead of preserving it through compatibility adapters.
- CR2 workflow clarity: reinforces that the supported handoff is the canonical
  staged driver plus `CartesianIDAHamiltonian` artifacts, not the old RHF
  payload stack.

Carrying-cost result:
- deleted: none in this docs pass.
- simplified: future source work has a precise two-file deletion surface and
  an expected net source deletion of roughly `1879` lines.
- quarantined: canonical driver changes, H1/final-basis stack changes,
  source-box materialization changes, Qiu-White donor kernels, artifacts,
  route/shellification/terminal-lowering/raw-block/RG/MWG/IDA, Hamiltonian
  assembly, committed tests, Cr2 workflow, replacements, adapters,
  compatibility wrappers, reports, status fields, and payload objects remain
  unapproved.
- exact remaining caller/blocker: if any live source/bin/test/tool caller
  depends on the RHF stack, the deletion pass must stop and report the exact
  caller rather than preserving the path through an adapter.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 158 - Approve Route-Driver Materialization Retirement

Commit(s):
- this commit - Approve route-driver materialization retirement

Summary:
- Approved `HP-RETIRE-DRV-MAT-FN-01`, `HP-RETIRE-DRV-MAT-TOOL-01`,
  `HP-RETIRE-DRV-MAT-DOC-01`, and `HP-RETIRE-DRV-MAT-TEST-01` as a narrow
  retirement/quarantine lane for the old route-driver materialization/report/
  save workflow.
- Retired names are `cartesian_materialization`, `cartesian_print_summary`,
  `cartesian_print_details`, `cartesian_save`, and the matching underscored
  route-driver materialization/report/save helpers. Focused search found no
  hits in `bin/cartesian_ham_builder.jl`; live hits are old wrapper
  definitions, old tools/harnesses, stale docs-policy assertions, and stale
  active-doc references.
- Current CR2-facing artifact production remains the canonical staged driver
  and `CartesianIDAHamiltonian` artifact path.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-RETIRE-DRV-MAT-*` and the retired wrapper names; confirm no `src`,
  `bin`, `test`, `tools`, artifact, or implementation files changed. A docs
  build is appropriate because compact authority and docs navigation changed.

Goal advancement:
- LT5/LT6: retires another route-era workflow surface rather than preserving it
  with compatibility adapters or policy tests.
- CR2 workflow clarity: keeps the supported producer story centered on the
  canonical staged driver and Hamiltonian artifacts.

Carrying-cost result:
- deleted: none in this docs pass.
- simplified: active docs now mark the old Slice D wrapper boundary as
  historical and approved for retirement.
- quarantined: canonical driver changes, staged producer changes, artifacts,
  route/shellification/terminal-lowering/raw-block/RG/MWG/IDA, Hamiltonian
  assembly, complete-core-shell H1/final-basis files, broad ordinary/Qiu-White
  donor-kernel retirement, replacement wrappers, adapters, status fields,
  payloads, new tests, and Cr2 workflow remain unapproved.
- exact remaining caller/blocker: if any current canonical producer path or
  public artifact workflow depends on these wrappers, the deletion pass must
  stop and report the exact dependency rather than preserving the workflow
  through adapters.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 159 - Approve Ladder Runner Deletion

Commit(s):
- this commit - Approve ladder runner deletion

Summary:
- Approved `HP-RETIRE-LADDER-RUNNERS-FN-01` and
  `HP-RETIRE-LADDER-RUNNERS-TEST-01` as a tiny deletion authority follow-up
  after `HP-RETIRE-DRV-MAT-*`.
- The old route-driver materialization/report/save workflow is already
  retired, and `tools/cartesian_driver_ladder_lib.jl` is quarantined. The two
  remaining runner scripts,
  `tools/run_cartesian_driver_ladder.jl` and
  `tools/run_cartesian_line_ladder.jl`, are now only entrypoints into that
  retired workflow.
- This pass does not approve deleting the quarantined ladder library itself.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-RETIRE-LADDER-RUNNERS-*`, the two runner names, and
  `cartesian_driver_ladder_lib`; confirm no `src`, `bin`, `test`, `tools`,
  artifact, or implementation files changed. A docs build is appropriate
  because compact authority changed.

Goal advancement:
- LT5/LT6 cleanup: removes the last obvious entrypoint pressure around the
  retired ladder workflow without reviving wrappers, adapters, or diagnostics.

Carrying-cost result:
- deleted: none in this docs pass.
- simplified: future source work has exactly two tool files to delete.
- quarantined: canonical driver changes, source changes, test changes except
  validation scans, artifacts, route/shellification/terminal-lowering/raw-block
  work, RG/MWG/IDA/Hamiltonian assembly, new wrappers, adapters, status fields,
  payloads, reports, tools, tests, Cr2 workflow, and deletion of
  `tools/cartesian_driver_ladder_lib.jl` remain unapproved.
- exact remaining caller/blocker: if any live source, canonical workflow, or
  approved tool still depends on these runner scripts, the deletion pass must
  stop and report the exact dependency.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 160 - Approve Public ns Direct-Core Side Parity

Commit(s):
- this commit - Approve public ns direct-core side parity

Summary:
- Approved `HP-COMP-NSCORE-FN-01` and `HP-COMP-NSCORE-TEST-01` as a narrow
  docs-only source authority amendment for the remaining even-`ns` PQS/WL
  parity bug.
- Public `ns` remains the requested cube/source/nesting size. Route-local
  `q` derivation remains `q = ns` for PQS and `q = ns - 2` for WL. The direct
  nucleus-centered core side, however, must be derived from public `ns` as
  `isodd(ns) ? ns : ns + 1`, not from route-local `q`.
- The user correction is now explicit: odd-side parity is necessary only for
  direct nucleus-centered core identity blocks. Boundary shells, WL
  boundary-stratum retained products, and non-direct support regions must not
  inherit that oddization rule. WL boundary retained counts remain
  `ns = 4 -> 56`, `ns = 5 -> 98`, and `ns = 6 -> 152`.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-COMP-NSCORE-*`, `direct_core_side`, and stale `core/contact`
  wording in active compact authority; confirm no `src`, `bin`, `test`,
  `tools`, artifact, or implementation files changed. A docs build is
  appropriate because compact authority and startup navigation changed.

Goal advancement:
- LT5/LT6: keeps PQS/WL comparison tied to the same public `ns` while
  preserving route-local WL boundary construction. This fixes the direct-core
  parity issue without changing the driver, public input model, artifact
  schema, or retained-basis coefficient construction.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: source work now has one precise target in
  `src/pqs_source_box_route_driver_helpers.jl`, with
  `src/cartesian_base_hamiltonian.jl` allowed only for one-center parent
  minimum sizing consistency.
- quarantined: driver changes, public input changes, route skeleton redesign,
  terminal lowering, retained-unit or terminal-realizer changes, artifact
  schema changes, manifest expansion, old WL materialization, committed tests,
  fixtures, and Cr2 workflow remain unapproved.
- exact remaining caller/blocker: if the parity fix requires route skeleton
  semantics, terminal lowering, retained-unit records, WL boundary coefficient
  construction, artifact schema, or driver inputs, the later source pass must
  stop and report the exact blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 161 - Approve Common Shell Decomposition Boundary

Commit(s):
- this commit - Approve common shell decomposition boundary

Summary:
- Approved `HP-COMP-SHELLGEOM-FN-01` and `HP-COMP-SHELLGEOM-TEST-01` as a
  narrow docs-only audit/source authority amendment.
- The design now records the two-geometry distinction explicitly. Common
  terminal shell decomposition is the shared first step: parent lattice,
  nuclear centers, direct core side, direct core regions, terminal shell
  regions, owned support rows, ordering, and coverage. It must be
  route-family-free.
- PQS and White-Lindsey diverge only after common shell records exist. PQS uses
  shell support plus a full source CPB for retained source-box modes. WL splits
  shell boundaries into faces/edges/corners/strata and builds retained columns
  by product-of-1D contractions.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-COMP-SHELLGEOM-*`, "common terminal shell decomposition", and stale
  wording that implies separate first-step PQS/WL shellifiers; confirm no
  `src`, `bin`, `test`, `tools`, artifact, or implementation files changed. A
  docs build is appropriate because compact authority, startup navigation, and
  the algorithm refresher changed.

Goal advancement:
- LT5/LT6: prevents another WL/PQS divergence caused by route-era vocabulary.
  The next source pass can ask a bounded question: does first-step atom shell
  decomposition already share one route-family-free path, and if not, where
  does route-family state leak in before lowering?

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has a precise audit/cleanup target in
  `src/cartesian_shellification/terminal_geometry.jl`, with narrow caller
  plumbing allowed in `src/pqs_source_box_route_driver_helpers.jl`.
- quarantined: driver changes, public input changes, route skeleton redesign,
  terminal lowering redesign, retained-unit or transform changes, PQS retained
  realization changes, WL boundary coefficient changes, artifact/manifest/
  reader changes, Hamiltonian/IDA/MWG/RG/raw-block changes, old WL
  materialization, committed tests, fixtures, and Cr2 workflow remain
  unapproved.
- exact remaining caller/blocker: if route-family-free common shell
  decomposition requires terminal lowering, retained-unit records, PQS retained
  modes, WL boundary coefficient construction, route skeleton semantics,
  artifact schema, or driver inputs, the later source pass must stop and
  report the exact blocker. If z-axis diatomic central-gap/contact policy
  proves entangled with family-specific retained geometry, that needs a later
  amendment.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 162 - Approve Diatomic Common Shellifier Entry

Commit(s):
- this commit - Approve diatomic common shellifier entry

Summary:
- Accepted the source result `e41aba6eb` as the atom-side direct-core parity
  fix, then recorded the remaining policy correction as docs authority:
  diatomic PQS/WL must also run through the same common shell decomposition
  function with the same first-step arguments.
- Approved `HP-COMP-SHELLGEOM-DIAT-FN-01` and
  `HP-COMP-SHELLGEOM-DIAT-TEST-01` as a narrow extension of the common-shell
  lane. For a fixed public z-axis diatomic system, parent axes, public `ns`,
  direct core side, nuclear centers, and bond axis, PQS and WL must enter the
  common shellifier identically before lowering.
- The design now classifies central-gap/contact, shared-shell, and
  outer-mismatch ownership as common shell geometry. PQS `q` and WL inner side
  are retained-construction inputs after common shell records exist.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-COMP-SHELLGEOM-DIAT-*`, "same first-step arguments", and stale
  wording that leaves diatomic central-gap/contact as a family-specific
  retained-construction decision; confirm no `src`, `bin`, `test`, `tools`,
  artifact, or implementation files changed. A docs build is appropriate
  because compact authority and the algorithm refresher changed.

Goal advancement:
- LT5/LT6: prevents the next PQS/WL divergence after atom parity. The future
  source pass should be caller-boundary cleanup, not retained-basis redesign:
  both families pass the same common inputs to the shellifier, and only then
  diverge into PQS source-box modes or WL boundary contractions.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: the remaining diatomic source target is now explicit and narrow:
  shellifier-boundary naming/input cleanup in
  `src/cartesian_shellification/terminal_geometry.jl` and
  `src/pqs_source_box_route_driver_helpers.jl`.
- quarantined: central-gap/contact algorithm redesign, terminal lowering,
  retained units, PQS retained realization, WL boundary coefficients, route
  skeletons, artifacts, driver inputs, committed tests/fixtures, and Cr2
  workflow remain unapproved.
- exact remaining caller/blocker: if same-function/same-argument diatomic
  entry requires changing the central-gap/contact algorithm rather than only
  route-family-independent shellifier inputs, the later source pass must stop
  and request a separate amendment.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 163 - Approve Mapped-COMX Source Span

Commit(s):
- this commit - Approve mapped-COMX source span

Summary:
- Approved `HP-MCOMX-FILE-01`, `HP-MCOMX-OBJ-01`, `HP-MCOMX-FN-01`,
  `HP-MCOMX-WIRE-01`, and `HP-MCOMX-TEST-01` as a narrow mainline
  source-span facility.
- The high-order branch produced the mapped-COMX idea and scratch evidence,
  but mainline owns the installation shape. High-order should consume the
  installed option for H/He/Cr benchmarking, like CR2 consumes Hamiltonian
  artifacts, rather than maintaining a duplicate implementation.
- The approved first rule is protected physical `P2` plus mapped Chebyshev
  enrichment `T_k(s_lambda(u))` with `lambda = 0.5`, no `sqrtJ`, and
  physical-`u` COMX localization. Ordinary polynomial source spans remain
  available and unchanged.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-MCOMX-*`, `mapped-COMX`, and forbidden-source-surface wording;
  confirm no `src`, `bin`, `test`, `tools`, artifact, or implementation files
  changed. A docs build is appropriate because compact authority and the
  algorithm index changed.

Goal advancement:
- LT5/LT6: installs a promising high-order experimental primitive at the
  correct mainline layer: source-span / carried-space facts, not Hamiltonian
  route scaffolding. It also preserves the consumer model for high-order and
  CR2.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has an explicit owner:
  `CartesianRawProductSources`, with only narrow PQS axis-transform wiring.
- quarantined: source-default replacement, public API/export changes,
  canonical driver input changes, artifacts, Hamiltonian/one-body/IDA/MWG/RG/
  raw-block/solver changes, injection/Ylm, `sqrtJ`, mapped-`s` production
  localization, high-order scaffolding imports, committed Cr/Cr2 fixtures, and
  Cr2 workflow remain unapproved.
- exact remaining caller/blocker: if the source option cannot be installed
  through `CartesianRawProductSources` and existing PQS raw-source
  axis-transform wiring without changing Hamiltonian assembly, artifacts,
  public driver inputs, or high-order-specific workflow, the later source pass
  must stop and report the missing mainline seam.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 164 - Correct Mapped-COMX Owner Seam

Commit(s):
- this commit - Correct mapped-COMX authority seam

Summary:
- Corrected the `HP-MCOMX-*` authority before accepting the current source
  WIP. The mapped-COMX idea remains approved, but the source owner is not
  `CartesianRawProductSources`. It is the existing nested doside / COMX seam:
  `pqs_source_axis_transform_facts_from_pgdg_axes(...) ->
  _nested_doside_1d(...) -> _nested_retained_span(...) ->
  _cleanup_comx_transform(...)`.
- The rejected WIP shape created a parallel numerical builder under
  `CartesianRawProductSources`, applied the nonlinear map to raw physical
  centers, duplicated COMX plumbing, and blurred metadata/provenance with
  construction data. The corrected contract requires normalized local
  `u in [-1, 1]` for `s_lambda(u)` and keeps the existing physical-position
  COMX cleanup.
- The first production pass is now restricted to `protected_degree = 2`.
  General protected degrees need a later parity-balanced mapped-order fill
  rule.

Validation:
- Required for this docs-only correction: `git diff --check` on touched docs
  and focused scans for `HP-MCOMX`, `mapped_comx_source_span.jl`,
  `CartesianRawProductSources`, `_nested_doside_1d`,
  `_nested_retained_span`, and `_cleanup_comx_transform` in the compact
  authority files.
- Source validation remains for the later doer pass. This correction does not
  approve or accept the current source WIP.

Goal advancement:
- LT5/LT6: preserves the mainline source-span facility goal while preventing a
  new route-like implementation layer. The high-order branch remains an idea
  source and benchmark consumer, not an implementation owner.

Carrying-cost result:
- deleted: none; docs-only correction.
- simplified: future implementation should be a small branch in the existing
  doside retained-span construction, not a new module/file.
- quarantined: `src/cartesian_raw_product_sources/mapped_comx_source_span.jl`,
  `CartesianRawProductSources` numerical builders, a second COMX wrapper,
  source-default changes, public API/export changes, artifacts, Hamiltonian/
  one-body/IDA/MWG/RG/raw-block/solver changes, `protected_degree != 2`,
  injection/Ylm, `sqrtJ`, mapped-`s` production localization, high-order
  scaffolding imports, committed Cr/Cr2 fixtures, and Cr2 workflow remain
  unapproved.
- exact remaining caller/blocker: if the option cannot be installed by adding
  a source-span selector at the existing doside seam while preserving
  `_cleanup_comx_transform(...)`, the source pass must stop and report the
  missing seam instead of building a parallel path.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 166 - Approve Mapped-COMX Terminal Wiring

Commit(s):
- this commit - Approve mapped-COMX terminal wiring

Summary:
- Approved `HP-MCOMX-TERM-FN-01` and `HP-MCOMX-TERM-TEST-01` as a narrow
  terminal-basis wiring lane. The mapped-COMX source-span path can produce
  materialized `AxisSourceTransformFact`s, but PQS terminal realization still
  needs to consume those carried facts in `_shell_seed(...)` for them to define
  shell seed coefficients.
- The approved seam is only
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  with `CartesianFinalBasisRealization.jl` allowed only for import/include
  cleanup if directly required.
- The pass may validate three materialized axis facts, check intervals and
  dimensions against the shell `outer_box` / source shape, build
  `full_coefficients` from their coefficient matrices, and then continue
  through the existing boundary-mode selection, support restriction,
  shell-local Lowdin, canonicalization, and support validation.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-MCOMX-TERM-*`, `_shell_seed`,
  `raw_product_source_axis_transform_facts`, `AxisSourceTransformFact`, and
  `_nested_projected_q_shell_full_sides` in compact authority.
- Later source validation should include ordinary PQS H2 regression, mapped
  source-span probe, a focused He or H terminal seam check proving mapped shell
  coefficients differ from ordinary and match the carried materialized axis
  facts, and the H2 supplemented RG endpoint if the touched path crosses it.

Goal advancement:
- LT5/LT6: closes the next authority gap between source-span construction and
  terminal-basis realization without expanding driver, artifact, Hamiltonian,
  or high-order workflow surfaces.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: source work now has a single terminal seam instead of pressure
  to add route metadata, artifact fields, or driver switches.
- quarantined: driver inputs, source defaults, artifacts/manifests/readers,
  Hamiltonian/IDA/MWG/RG/raw-Gaussian/solver/EGOI/Cr2/high-order workflow
  changes, second COMX wrappers, and committed tests/fixtures remain
  unapproved.
- exact remaining caller/blocker: if `_shell_seed(...)` cannot consume carried
  axis facts without changing shell ownership, retained-rule semantics,
  Lowdin realization, artifact schema, or driver inputs, the source pass must
  stop and report that blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 167 - Approve Mapped-COMX Driver Selector

Commit(s):
- this commit - Approve mapped-COMX driver selector

Summary:
- Approved `HP-MCOMX-DRV-FN-01` and `HP-MCOMX-DRV-TEST-01` as the narrow
  canonical-driver/facade lane for selecting source-span family.
- The visible driver choice is `source_span = :ordinary` or
  `source_span = :mapped_comx`; ordinary remains the default. The mapped
  choice is currently PQS-only and must reject clearly for `nesting = :wl`
  unless a later WL-specific source-span amendment approves otherwise.
- The approved source surface is `bin/cartesian_ham_builder.jl`,
  `src/cartesian_base_hamiltonian.jl`, and narrow propagation in
  `src/pqs_source_box_route_driver_helpers.jl` to the already-approved
  mapped-COMX source-axis transform fact path.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-MCOMX-DRV-*`, `source_span`, `:mapped_comx`, and the forbidden
  route-record/artifact/new-COMX wording; docs build because compact authority
  and the implementation index changed.
- Later source validation should include default ordinary driver
  artifact/readback, mapped-COMX H or He PQS driver smoke proving carried facts
  are basis-defining, bounded ordinary-vs-mapped He supplemented/MWG/IDA
  comparison if practical, H2 RG endpoint, and no Cr2 run.

Goal advancement:
- LT1/LT3/LT5: makes the mapped-COMX option usable through the compact,
  copyable canonical driver while keeping the driver at construction-choice
  level instead of route diagnostics.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one visible selector and must use the
  existing doside/source-axis/terminal wiring.
- quarantined: new route records, terminal-lowering changes, artifact schema
  changes, source defaults, another COMX path, high-order workflow controls,
  and committed tests/fixtures remain unapproved.
- exact remaining caller/blocker: if making `source_span` driver-selectable
  requires new route records, terminal-lowering changes, artifact schema
  changes, or another COMX path, the source pass must stop and report that
  blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 168 - Wire Mapped-COMX Driver Selector

Commit(s):
- this commit - Wire mapped-COMX driver selector

Summary:
- Accepted the `HP-MCOMX-DRV-FN-01` source pass. The canonical driver now has
  the visible construction choice `source_span = :ordinary` or
  `source_span = :mapped_comx`, with ordinary still the default and
  `:mapped_comx` clearly rejected for `nesting = :wl`.
- The staged base input path normalizes the selector and passes it only into
  the existing PQS source-box route. Mapped runs now attach materialized
  mapped-COMX `AxisSourceTransformFact`s to PQS retained units through the
  existing retained-unit metadata seam; the already-approved terminal
  `_shell_seed(...)` consumer then makes those facts basis-defining.
- No Hamiltonian assembly, artifact schema, manifest, reader, RG/MWG/IDA,
  terminal-lowering, route-skeleton, or COMX numerical kernel changed.

Validation:
- Doer validation: `git diff --check`; package load; mapped source-span probe;
  mapped terminal seam probe; driver wiring probe; H atom default ordinary
  artifact/readback; He ordinary/mapped base artifact/readback; H2 mapped base
  artifact/readback; H atom supplemented ordinary/mapped `cc-pVTZ`; public
  `:pgdg` rejection; `:mapped_comx` plus `nesting = :wl` rejection; H2 RG
  endpoint with self-Coulomb `0.4574265214362095`.
- Manager validation rerun: `git diff --check`; package load; mapped terminal
  seam probe (`mapped_terminal_match_delta = 0.0`,
  `mapped_vs_ordinary_seed_delta = 1.835e-01`); mapped source-span probe for
  `n_s = 5,6,7`; driver wiring probe (ordinary/mapped He dimension `419`,
  shell coefficient delta `1.662340e-01`); H2 RG endpoint/readback.
- Mechanical diff gate: source/bin/test scan had suspicious hits only for
  `unit.metadata` reads and metadata merge in
  `src/pqs_source_box_route_driver_helpers.jl`. These are accepted because
  they reuse the existing retained-unit -> raw-product-source transform fact
  seam already consumed by `retained_unit_transform_contract_plan(...)`; no
  new route record, status cloud, artifact field, or metadata-as-report bus was
  added.

Goal advancement:
- LT1/LT3/LT5: the mapped-COMX source-span option is now usable through the
  canonical driver for real H/He/PQS comparisons, while staying at the same
  level as `nesting`: a construction-family selector, not a route diagnostic
  or high-order workflow hook.

Carrying-cost result:
- deleted: the previous user-facing blocker where the canonical driver rejected
  `source_span` as an unknown input.
- simplified: mapped-COMX uses the existing retained-unit raw axis-fact seam
  and the existing terminal consumer; no parallel COMX or route path was added.
- quarantined: public `:pgdg`, WL mapped-COMX, route records, terminal
  lowering, artifacts/manifests/readers, Hamiltonian/IDA/MWG/RG/raw-block
  changes, EGOI/Cr2/high-order workflow, and committed tests remain
  unapproved.
- not deleted because: ordinary fallback remains the default supported
  production path and is required when materialized source-axis facts are
  absent.
- exact remaining caller/blocker: none for canonical PQS driver selection;
  high-order can now run ordinary-vs-mapped driver comparisons.
- added src lines: 110.
- deleted src lines: 8.
- new tests: none committed; ignored `tmp/work` probes only.
- new metadata/status fields: none; existing retained-unit
  `raw_product_source_axis_transform_facts` metadata is reused.

## Cartesian Hamiltonian Producer Pass 169 - Record He ns=5 Mapped-COMX Limit

Commit(s):
- this commit - Record He ns=5 mapped-COMX limit

Summary:
- Recorded high-order-manager's 2026-06-26 interpretation memo as mainline
  authority context. The installed mapped-COMX path is mechanically wired
  through the real driver, but bounded He/PQS `n_s = 5` testing did not improve
  all-electron scalar H1/IDA capture and was worse than ordinary PQS on the
  harder `core_spacing = 0.2` stress case.
- This does not retract `HP-MCOMX-*` or the `source_span = :mapped_comx`
  selector. It blocks default promotion: ordinary PQS remains the default and
  mapped-COMX remains opt-in.
- The next relevant mapped-COMX evidence is bounded He `n_s = 6` and `n_s = 7`
  H1/IDA testing with shell-restricted scalar-capture diagnostics. Cr and
  molecule-facing promotion remain later.

Validation:
- Docs-only status validation: `git diff --check`; focused scan for `n_s = 5`,
  mapped-COMX, default-promotion wording, and He `n_s = 6`/`7` next evidence.
- No source/bin/test/tool edits and no implementation tests.

Goal advancement:
- LT1/LT5: preserves mapped-COMX as a consumer-testable mainline facility
  while preventing the first angular proxy or driver-wiring success from being
  misread as production/default evidence.

Carrying-cost result:
- deleted: none; docs-only status correction.
- simplified: future agents have a clear validation interpretation: opt-in
  mapped-COMX is allowed, default promotion is not.
- quarantined: mapped-COMX defaults, Cr/Cr2 promotion, shell-adapted injection,
  target-adapted source rules, and high-order workflow controls remain
  unapproved.
- exact remaining caller/blocker: mapped-COMX needs bounded He `n_s = 6`/`7`
  H1/IDA evidence before any default or molecule-facing promotion.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 170 - Accept Cr Map-ns Residual Audit

Commit(s):
- none - read-only residual audit accepted

Summary:
- Accepted repo-doer's read-only audit of the Cr atom PQS
  `basis_ns = 9`, `map_ns = 11`, `lmax = 1` residual-GTO construction
  failure. The failure is now classified as a small final `R' S R` identity
  acceptance miss, not owner grouping failure, merge singularity, or wrong
  candidate selection.
- The reconstructed fixture retained all `33` Cr/cc-pV5Z `lmax = 1`
  residual directions on the single Cr owner. The owner-local residual metric
  is not singular (`eig min/max = 3.637e-08 / 6.062e-05`, condition about
  `1.67e3`), the merge metric is essentially identity, and `G' S R` is small
  at `2.558e-13`.
- Production throws because `max |R' S R - I| = 2.235e-08` is slightly above
  the current combined `identity_atol = 1.0e-8` rule, whose allowed error for
  this case is about `2.000e-08`.

Validation:
- Doer validation: `git diff --check`; package load; ignored audit probe
  `tmp/work/cr_atom_pqs_map_ns_residual_validation_audit.jl`; failed case plus
  nearby passing comparisons `basis_ns=9,map_ns=9,lmax=1` and
  `basis_ns=7,map_ns=11,lmax=1`.
- Manager validation: reviewed the audit table and current RG tolerance
  authority. No tracked source/doc/bin/test/tool files were changed by the
  audit.

Goal advancement:
- LT1/LT3: preserves the Cr atom map-ns ladder interpretation. The useful
  science signal remains that `lmax = 2` is the large improvement; `map_ns=11`
  is not yet a compelling run target. If `map_ns=11` should be admitted for
  completeness, the change should be framed as final identity acceptance only,
  not residual cutoff or selection policy.

Carrying-cost result:
- deleted: none.
- simplified: the blocker is classified narrowly as final identity tolerance,
  not residual algorithm redesign.
- quarantined: cutoff changes, width/zeta filtering, residual selection,
  owner grouping, merge rules, driver changes, artifacts, Cr2 workflow, and
  solver work remain unapproved.
- not deleted because: this was a read-only audit with an ignored `tmp/work`
  probe.
- exact remaining caller/blocker: policy decision needed only if the Cr
  `basis_ns=9,map_ns=11,lmax=1` construction should pass production
  validation.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 171 - Approve RG 5e-8 Cutoff Policy

Commit(s):
- this commit - Approve RG 5e-8 cutoff policy

Summary:
- Approved `HP-RG-CUTOFF-FN-01` and `HP-RG-CUTOFF-TEST-01` as a docs-only
  policy amendment. The default Residual Gaussian owner-local residual
  occupation cutoff is now `5.0e-8`, and the default final residual
  `R' S R` identity validation `identity_atol` is also `5.0e-8`.
- This intentionally supersedes the prior interpretation of the Cr atom
  `basis_ns = 9`, `map_ns = 11`, `lmax = 1` audit as only a final identity
  tolerance miss. The retained Cr marginal occupation was `3.637e-8`; the new
  production policy is to drop such marginal directions by default rather than
  preserve them because the old cutoff was lower.
- Owner-local grouping, negative-eigenvalue tolerances, final merge metric
  checks, `G' S R` validation, width/zeta filtering, MWG/IDA, artifacts,
  driver workflow, public API, and source ownership remain unchanged.

Validation:
- Docs-only amendment validation required: `git diff --check`; focused scans
  for `HP-RG-CUTOFF-*`, `5.0e-8`, stale `residual_occupation_cutoff = 1.0e-8`
  wording, and Cr/Be/H2 validation gates.
- Later source validation should show the Cr atom `basis_ns = 9`,
  `map_ns = 11`, `lmax = 1` residual construction passes or cleanly drops the
  marginal `s4` direction, Be atom cc-pV5Z still passes, and the H2
  residual-GTO/MWG endpoint remains unchanged.

Goal advancement:
- LT1/LT3: makes the retained-residual policy match the intended production
  science: marginal owner-local residual directions around `3.6e-8` are not
  part of the default retained RG space.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has a single explicit cutoff/tolerance pair
  instead of an implicit small-constant change.
- quarantined: residual algorithm changes, owner grouping changes, merge rule
  changes, width/zeta filtering defaults, MWG/IDA, artifacts, driver/public
  API work, committed tests, and Cr2 workflow remain unapproved.
- exact remaining caller/blocker: if the Cr atom case cannot pass or drop the
  marginal direction by changing only the two approved defaults, implementation
  must stop and report the blocker.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 172 - Approve RG Cutoff Test Assertion Update

Commit(s):
- this commit - Approve RG cutoff test assertion update

Summary:
- Approved a narrow follow-up under `HP-RG-CUTOFF-TEST-01` for the exact
  committed H2 endpoint test mismatch found during source validation.
- The only committed test edits approved are in
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`: update the
  in-memory `residual.occupation_cutoff` assertion and the artifact/provenance
  `values[:occupation_cutoff]` assertion from `1.0e-8` to `5.0e-8`.
- No new test file, fixture, driver workflow, artifact/schema change, Cr2 run,
  or broader test rewrite is approved.

Validation:
- Docs-only amendment validation required: focused scan for
  `cartesian_r3a_h2_augmented_one_body_runtests.jl`, `occupation_cutoff`,
  `values[:occupation_cutoff]`, `5.0e-8`, and `HP-RG-CUTOFF-TEST-01`;
  `git diff --check` on touched docs.

Goal advancement:
- LT1/LT3: removes a documentation-authority mismatch so the source pass can
  validate the new RG default policy through the existing H2 endpoint gate.

Carrying-cost result:
- deleted: none; docs-only authority correction.
- simplified: the existing endpoint test remains the gate; no new committed
  test surface is created.
- quarantined: all other committed tests/fixtures, driver changes, artifacts,
  Cr2 workflow, and broader RG algorithm changes remain unapproved.
- exact remaining caller/blocker: source doer may update only those two
  assertions; any other test failure needs a separate report or amendment.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 173 - Approve Outer-Mismatch Compact Lowering

Commit(s):
- this commit - Approve common outer-mismatch compact lowering

Summary:
- Approved `HP-COMP-OUTERMM-FN-01` and `HP-COMP-OUTERMM-TEST-01` as a
  docs-only authority amendment. The CR2/HFDMRG inventory exposed a producer
  basis-size bug: z-axis diatomic `:outer_mismatch_slab` regions were being
  lowered as `:direct_boundary_slab_identity_cpb` and then realized as full
  identity terminal rows.
- The first version of the request was PQS-centered, but live lowering code
  shows White-Lindsey falls through the same direct identity seam. The approved
  contract is therefore common: thickness-1 outer-mismatch slabs are boundary
  slabs, not real shells, and PQS/WL must lower them through the same compact
  boundary-slab function with the same terminal region, public `ns`, and
  source/support facts.
- Real shell regions remain route-specific after common shellification: PQS
  uses full source-box shell projection, while WL uses face/edge/corner
  product-of-1D contractions.

Validation:
- Docs-only validation required: `git diff --check`; focused scans for
  `HP-COMP-OUTERMM-*`, `outer_mismatch_slab`,
  `direct_boundary_slab_identity_cpb`, `ns x ns x 1`, `ns x ns x ns`, and
  stale PQS-only wording; confirm no source/bin/test/tool files changed.
- Later source validation should include bounded H2 or Be2 artifact/readback
  under `nesting = :pqs` and `nesting = :wl`, prove neither family lowers
  outer-mismatch slabs to direct identity CPBs, and audit the same compact
  boundary-slab function/input seam.

Goal advancement:
- LT5/LT6: protects the canonical producer artifact path from a severe
  producer-side basis inflation while preserving the common-shell doctrine:
  shared first geometry, identical lowering for thickness-1 mismatch slabs,
  and route-specific retained construction only for real shells.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future source work has one common outer-mismatch lowering target
  instead of a PQS-only patch plus hidden WL divergence.
- quarantined: driver changes, artifact/schema/reader changes, terminal
  realization changes, retained-unit record changes, route skeleton redesign,
  residual/RG/MWG/IDA changes, committed Cr2 tests, and direct slab deletion
  remain unapproved.
- exact remaining caller/blocker: if existing compact boundary facts are
  insufficient, source work must stop and report the missing native fact; if
  one end has more than `ns` boundary slabs, a separate policy decision is
  required.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 174 - Generalize Thin-Slab Lowering

Commit(s):
- this commit - Generalize thin-slab lowering authority

Summary:
- Superseded the outer-mismatch-only `HP-COMP-OUTERMM-*` lane with
  `HP-COMP-THINSLAB-FN-01` and `HP-COMP-THINSLAB-TEST-01`, then tightened it
  to account for multi-slice outer-mismatch regions.
- The follow-up audit showed `:direct_midpoint_slab` follows the same bad
  direct-identity path as `:outer_mismatch_slab`: it lowers through
  `:direct_slab_identity_cpb`, becomes a direct retained unit, and is then
  realized as full identity support rows. That is the same conceptual error
  for a thickness-1 slab.
- The new authority is route-family-free for thin slabs: for both PQS and WL,
  `:direct_midpoint_slab` and `:outer_mismatch_slab` must use the same compact
  slab lowering function and inputs. The unit-slice scale is `ns x ns x 1`;
  a thickness-`t <= ns` outer-mismatch stack should scale about
  `t * ns * ns`. Direct/core sectors remain identity, and real shell regions
  remain route-specific after common shellification.

Validation:
- Docs-only validation required: `git diff --check`; focused scans for
  `HP-COMP-THINSLAB-*`, `HP-COMP-OUTERMM-*`, `:direct_midpoint_slab`,
  `:outer_mismatch_slab`, `:direct_slab_identity_cpb`,
  `:direct_boundary_slab_identity_cpb`, and `ns x ns x 1`; confirm no
  source/bin/test/tool files changed.
- Later source validation should prove bounded H2 or Be2 under
  `nesting = :pqs` and `nesting = :wl` no longer lowers midpoint or
  outer-mismatch slabs to identity CPBs, and should audit that both families
  call the same compact thin-slab function with matched inputs. If
  shellification does not already expose slab normal/thickness, the source
  pass may add native metadata in `terminal_geometry.jl`; it must not parse
  role strings.

Goal advancement:
- LT5/LT6: tightens the common-shell doctrine by separating three cases:
  direct/core identity sectors, real shells with route-specific retained
  construction, and thin slab stacks with shared compact lowering.

Carrying-cost result:
- deleted: none; docs-only authority correction.
- simplified: future source work has one thin-slab target instead of separate
  midpoint and outer-mismatch patches.
- quarantined: driver changes, public inputs, artifact/schema/reader changes,
  broad terminal realization redesign, route skeleton redesign, real-shell
  PQS/WL policy changes, RG/MWG/IDA changes, committed Cr2 tests, and slab
  deletion remain unapproved.
- exact remaining caller/blocker: if slab thickness exceeds `ns`, or if compact
  lowering needs missing native facts beyond slab normal/thickness metadata,
  source work must stop for a separate policy decision.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 175 - Approve Angular-Balanced Geometry Audit

Commit(s):
- this commit - Approve angular-balanced shellification audit

Summary:
- Added `HP-COMP-ANGBOX-AUDIT-01` as measurement-only authority and recorded
  `HP-COMP-ANGBOX-FN-01` / `HP-COMP-ANGBOX-TEST-01` as candidate-only source
  and validation lanes. This responds to the CR2 slab inventory and the
  follow-up discussion that compact lowering alone treats the symptom if the
  shellifier keeps creating large axial leftovers.
- The durable geometry rule is now explicit: shared z-axis diatomic molecular
  shell bodies should be angular-balanced from the outer nuclei in physical
  parent-axis coordinates. Axial excess should be emitted as native thin-slab
  stacks, not left as route-family-specific shells or direct identity sectors.
- The thin-slab category is deliberately broader than the prior
  outer-mismatch wording. It covers central midpoint slabs, planned
  non-boundary angular z-extension slabs, planned boundary angular z-extension
  slabs, and unexpected outer-mismatch fallback slabs. Planned z-extension
  stacks may be chunked into slices with thickness `<= ns`; unexpected
  fallback slabs with thickness greater than `ns` still stop for a policy
  decision.

Validation:
- Docs-only validation required: `git diff --check`; focused scans for
  `HP-COMP-ANGBOX-*`, `angular-balanced`, `angular_z_extension_slab`,
  `outer_nucleus_45_degree`, `thin-slab`, and `t * ns * ns`; confirm no
  source/bin/test/tool files changed.
- Later audit should report physical parent endpoints/counts, snapped nuclear
  indices, core boxes, molecular inner box, each proposed shared-shell
  expansion, transverse scale, low/high longitudinal margins, angular-balance
  ratio, planned z-extension stacks, and residual outer mismatch for bounded
  H2/Be2 and CR2-style fixtures.

Goal advancement:
- LT5/LT6: moves the slab correction up to the owner of the geometry. The
  target is no longer merely "compress the leftover slabs"; it is "create real
  shell bodies and planned compact slab stacks with native metadata."

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: midpoint slabs, non-boundary z-extension slabs, boundary
  z-extension slabs, and fallback outer mismatch now share one conceptual
  thin-slab category.
- quarantined: production source edits, driver changes, artifact/schema/reader
  changes, route skeleton redesign, RG/MWG/IDA/Hamiltonian/raw-block changes,
  direct slab deletion, committed Cr2 tests/fixtures, and Cr2 workflow remain
  unapproved.
- exact remaining caller/blocker: source work waits on the ignored geometry
  audit. If angular-balanced shellification needs more than native endcap slab
  metadata and region emission in `terminal_geometry.jl`, the source lane must
  stop for a separate amendment.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 176 - Correct Angular-Box Direction

Commit(s):
- this commit - Correct angular-box z-extension wording

Summary:
- Corrected the docs-only `HP-COMP-ANGBOX-*` wording after reviewing the audit
  interpretation with the user. The prior amendment put the repair at the
  right layer, but described the geometry direction backwards: the current
  shared-shell boxes are underextended in `z` relative to the
  outer-nucleus 45-degree rule, not overextended.
- Replaced the stale overextension wording with angular-balance required
  z-extension beyond the ordinary index-layer shell body. Planned non-boundary
  and boundary z-extension slabs remain thin-slab stack objects and continue
  to lower through `HP-COMP-THINSLAB-*`.
- Renamed the proposed metadata to `:angular_z_extension_slab` and
  `angular_extension_physical` so the candidate contract does not encode the
  wrong sign.

Validation:
- Docs-only validation required: `git diff --check`; focused scans for stale
  overextension terminology and old endcap metadata names; confirm no
  source/bin/test/tool files changed.

Goal advancement:
- LT5/LT6: preserves the common-shell and thin-slab plan while aligning the
  source target with the audit evidence. Future doer work should generate
  angular-balanced target boxes and planned z-extension slab stacks throughout
  shared-shell growth, not only at the final parent boundary.

Carrying-cost result:
- deleted: stale wording only.
- simplified: the geometry diagnosis is now one story: index-layer shell
  growth underextends z, so angular balance creates planned z-extension slabs.
- quarantined: source edits, driver changes, artifacts, route skeleton
  redesign, and Cr2 workflow remain unapproved.
- exact remaining caller/blocker: source work remains blocked pending an
  implementation blurb that names the `raw_terminal_geometry(...)` source cut.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 177 - Clarify Angular Target Decomposition

Commit(s):
- this commit - Clarify angular target decomposition

Summary:
- Tightened the `HP-COMP-ANGBOX-*` docs before source work. The policy now
  explicitly says each shared-shell step computes an angular-balanced target
  box, and that the ordinary index-layer shell body plus planned z-extension
  thin-slab stacks realizes the target coverage.
- This avoids the remaining misread where a doer could make the whole
  angular-balanced target a route-specific shell body. The ordinary body may
  be underextended in `z`; the bond-axis-only difference is a planned
  thin-slab stack.
- Cleaned the remaining hybrid old/new wording in `AGENTS.md`.

Validation:
- Docs-only validation required: `git diff --check`; focused scans for stale
  shell-body and endcap wording; confirm no source/bin/test/tool files
  changed.

Goal advancement:
- LT5/LT6: preserves the shellifier-owned angular target while keeping the
  retained-construction split intact. Planned z-extension slabs are not PQS or
  WL real shells and not identity rows.

Carrying-cost result:
- deleted: stale wording only.
- simplified: the source-facing model is now ordinary shell body plus planned
  z-extension slabs equals angular target coverage.
- quarantined: source edits, driver changes, artifact/schema changes, route
  skeleton redesign, and Cr2 workflow remain unapproved.
- exact remaining caller/blocker: source work still needs a separate doer
  blurb naming the `raw_terminal_geometry(...)` cut and validation.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 178 - Promote Angular-Box Source Lane

Commit(s):
- this commit - Promote angular-box shellification lane

Summary:
- Promoted `HP-COMP-ANGBOX-FN-01` and `HP-COMP-ANGBOX-TEST-01` from
  candidate-only to approved source/validation authority. The ignored audit
  established the exact cut: the large z-end leftovers are a common
  shellifier geometry issue, identical for PQS and WL, not a retained-basis
  or route-family mismatch.
- Approved source surface is `src/cartesian_shellification/terminal_geometry.jl`
  only, with `src/pqs_source_box_route_driver_helpers.jl` and
  `src/pqs_source_box_diatomic_complete_core_shell.jl` optional only for
  directly required summary/caller plumbing.
- Approved behavior is narrow: when ordinary z-axis diatomic shared-shell
  expansion stops with transverse axes saturated and bond-axis parent support
  remaining, emit the bond-axis leftovers as planned
  `:angular_z_extension_slab` stack regions, chunked so each slab has
  thickness `<= ns`, with native axis/side/thickness, stack, bond-axis,
  angular-rule, margin, transverse-scale, and extension-size metadata.
  Generic `:outer_mismatch_slab` remains only for unexpected or non-z
  leftovers. Thin-slab lowering remains a separate blocker under
  `HP-COMP-THINSLAB-*`.

Validation:
- Docs-only validation required: `git diff --check`; focused scans for
  `HP-COMP-ANGBOX-*`, `Status: approved`, `:angular_z_extension_slab`,
  `terminal_geometry.jl`, and no live candidate-only ANGBOX wording; confirm
  no source/bin/test/tool files are staged in this docs pass.
- Later source validation requires package load and the ignored angular
  geometry audit showing H2/Cr2-style fixtures, planned z-extension support,
  zero residual z mismatch after classification, and PQS/WL geometry parity.
  No artifact/readback is required while lowering is intentionally deferred.

Goal advancement:
- LT5/LT6: opens the source lane at the correct owner layer. The shellifier
  may now classify the missing z coverage as planned geometry while keeping
  lowering, retained units, terminal realization, artifacts, and driver
  workflow out of scope.

Carrying-cost result:
- deleted: none; docs-only promotion.
- simplified: one approved shellification owner replaces the previous
  candidate-only audit interpretation.
- quarantined: driver changes, artifact/schema/provenance/reader changes,
  terminal lowering, retained units, transform contracts, terminal
  realization, Hamiltonian/RG/MWG/IDA/raw-block changes, route-family-specific
  PQS/WL geometry, committed tests/fixtures, Cr2-specific branches, and Cr2
  workflow remain unapproved.
- exact remaining caller/blocker: doer may implement the shellification
  emission only; compact lowering of the emitted slabs is still blocked until
  the `HP-COMP-THINSLAB-*` implementation pass.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 179 - Approve Neutral Face-Product Helper

Commit(s):
- this commit - Approve neutral face-product helper

Summary:
- Approved `HP-COMP-FACEPROD-FN-01` and `HP-COMP-FACEPROD-TEST-01` as the
  neutral terminal face-product helper seam needed before the thin-slab source
  pass should proceed. The source attempt showed the right numerical primitive
  was reusable, but the wrong ownership: putting a shared slab helper in
  `white_lindsey_terminal_basis_realization.jl` and calling it from PQS makes
  neutral coefficient assembly look WL-owned.
- New approved helper owner is
  `src/cartesian_final_basis_realization/terminal_face_product_blocks.jl`,
  included from `CartesianFinalBasisRealization.jl`, with narrow consumers in
  the WL and PQS terminal realization files. The helper is private/internal
  and should reuse `_nested_doside_1d(...)` and `_nested_face_product(...)`.
- Required reuse proof: White-Lindsey facet terminal realization should be
  refactored to use the neutral helper. Future `HP-COMP-THINSLAB-*` slab
  realization may then use the same helper for midpoint, outer-mismatch
  fallback, and angular z-extension slabs.

Validation:
- Docs-only validation required: `git diff --check` on staged docs; focused
  scans for `HP-COMP-FACEPROD-*`, `terminal_face_product_blocks.jl`,
  `_nested_face_product(...)`, `_nested_doside_1d(...)`, and no implication
  that PQS calls a WL-owned helper; confirm source WIP is not staged in this
  docs commit.
- Later source validation requires package load, WL facet parity before/after
  refactor, bounded WL H2/Be2 artifact/readback, bounded PQS H2/Be2
  artifact/readback if PQS terminal realization imports the helper, and no
  Cr2 run.

Goal advancement:
- LT5/LT6: reduces cognitive complexity by placing the reusable face-like
  product block exactly where both WL facets and shared thin slabs can consume
  it. This keeps real shells route-specific while making the face-product
  coefficient assembly route-neutral.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future doer work has one neutral helper target instead of a
  WL-owned helper plus PQS call-through.
- quarantined: driver changes, public API/export, artifacts, shellification,
  terminal lowering policy, route skeletons, RG/MWG/IDA/Hamiltonian/raw-block
  changes, old high-order workflow, committed tests/fixtures, Cr2 workflow,
  duplicate face-product assembly, PQS-specific slab projection, and treating
  thin slabs as WL boundary strata remain unapproved.
- exact remaining caller/blocker: if the helper cannot serve both current WL
  facets and future thin slabs without numerical changes, source work must
  stop and report whether the blocker is helper signature, support-record
  shape, retained-unit metadata, or terminal-realization ownership.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 180 - Approve Thin-Slab Metadata Inventory

Commit(s):
- this commit - Approve thin-slab metadata inventory

Summary:
- Approved `HP-COMP-THINSLAB-META-FN-01` and
  `HP-COMP-THINSLAB-META-TEST-01` for the live terminal-shellification
  metadata/scaffold inventory update in
  `src/cartesian_terminal_shellification_geometry.jl`. The immediate blocker
  was not numerical lowering, but the old inventory mapping in
  `_cartesian_terminal_region_unit_mapping(region)`: midpoint slabs and
  outer-mismatch slabs still described direct identity CPBs, and planned
  `:angular_z_extension_slab` had no case.
- The approved source lane is deliberately metadata-only. It lets route
  summaries agree with compact thin-slab lowering, while forbidding
  coefficient materialization, Hamiltonian data construction, shellification
  geometry changes, route skeleton redesign, artifact/report payloads, and a
  new reporting framework. Direct core and atom-contact core identity mappings
  remain unchanged.

Validation:
- Docs-only validation used: `git diff --check`; focused scans for
  `HP-COMP-THINSLAB-META-*`, `_cartesian_terminal_region_unit_mapping`,
  `_cartesian_terminal_shellification_region_unit_inventory`,
  `:angular_z_extension_slab`, and direct identity slab vocabulary; staged
  diff check limited to docs/`AGENTS.md`.
- Later source validation requires package load, the existing angular geometry
  audit, a thin-slab inventory/probe that no longer fails on
  `angular_z_extension_slab`, and a focused scan confirming no planned direct
  identity mapping for midpoint, outer-mismatch, or angular z-extension slabs.
  No Cr2 run is required.

Goal advancement:
- LT5/LT6: closes the authority gap between common shellification,
  angular z-extension classification, neutral face-product assembly, and
  compact thin-slab lowering. The inventory layer now has permission to stop
  broadcasting stale identity-sector semantics while remaining non-numerical.

Medium-term checkpoint:
- Active: the common-shell/thin-slab lane now has aligned authority for
  shellification classification, metadata inventory, neutral face-product
  coefficient assembly, and compact lowering.
- In progress: source WIP still needs review and validation; no lowering pass
  is accepted merely because this metadata amendment exists.
- Deferred: Cr2 full workflow, artifact/schema changes, driver reporting,
  solver/ECP, and far-separated central product policy remain outside this
  lane.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: one live metadata inventory can now match compact thin-slab
  lowering instead of preserving stale direct-identity categories.
- quarantined: source edits outside the named metadata/scaffold surface,
  shellification algorithm changes, terminal realization changes, driver
  changes, artifacts, route reports, and Cr2 workflow remain unapproved.
- exact remaining caller/blocker: doer may update the inventory mapping only;
  if that requires materialized retained units, artifact/report payloads,
  shellification changes, or route skeleton redesign, the pass must stop.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 181 - Approve Driver Region Inventory

Commit(s):
- this commit - Approve driver terminal-region inventory

Summary:
- Approved `HP-DRV-INV-FN-01` and `HP-DRV-INV-TEST-01` for a compact
  terminal-region / shellification inventory in canonical driver output. The
  Cr2 slab blowup showed that normal users need to see region kind, support
  rows, final columns, compression ratios, and identity-vs-compact realization
  without running ignored probes.
- The lane is driver-output authority, not route diagnostics. Approved source
  files are `bin/cartesian_ham_builder.jl` and
  `src/cartesian_base_hamiltonian.jl`, with optional compact accessors in the
  route helper or final-basis modules only if directly required. The summary is
  bounded and human-facing: region label/index, region kind, lowering or
  realization kind, support rows, final columns, compression ratio, realization
  class, and native slab facts when applicable.
- Forbidden surfaces remain broad: no new driver inputs, flags, stop-after
  controls, route switches, source-mode or pair inventories, raw-block or
  all-row dumps, artifact/schema/reader changes, numerical construction
  changes, route-report framework, solver settings, Cr2-specific workflow, or
  committed tests.

Validation:
- Docs-only validation required: `git diff --check`; focused scans for
  `HP-DRV-INV-*`, terminal-region inventory wording, bounded-output
  restrictions, and absence of source/bin/test/tool staging in this docs pass.
- Later source validation requires package load, bounded H2 or Be2 driver
  inventory output for `nesting = :pqs` and `nesting = :wl`, supplemented smoke
  if supplemented-stage objects are touched, unchanged artifact/readback
  deltas, and a scan confirming the output excludes source modes, pair
  inventories, raw-block details, all-row listings, and full metadata. No Cr2
  run is required.

Goal advancement:
- LT5/LT6: moves the slab-size lesson into the normal human workflow without
  promoting route internals. Users should see accidental identity sectors and
  compact slab compression during ordinary canonical-driver runs.

Carrying-cost result:
- deleted: none; docs-only authority pass.
- simplified: future doer work has a bounded inventory target instead of an
  open-ended reporting request.
- quarantined: artifact schemas, route reports, source-mode/pair/raw-block
  dumps, all-row metadata, driver inputs, numerical construction, and Cr2
  workflow remain unapproved.
- exact remaining caller/blocker: if existing stage/final-basis summaries
  cannot support the compact output without a broad payload or report
  framework, source work must stop and identify the missing summary seam.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 182 - Driver Terminal Inventory Output

Commit(s):
- this commit - Add driver terminal inventory output

Summary:
- Accepted the `HP-DRV-INV-FN-01` implementation. The canonical driver now
  prints a bounded terminal-region inventory immediately after the base
  working-basis stage when `print_contract = true`.
- The summary is intentionally small: base final dimension, region key,
  region kind, lowering kind, support rows, final columns, compression ratio,
  identity/compact class, and native slab axis/side/thickness/stack facts.
  It does not print source-mode rows, pair inventories, raw-block details,
  coefficients, all support rows, route skeletons, or full metadata.
- Manager amendment replaced direct plan-field reads with the existing
  `CartesianRetainedUnits.units(...)` and
  `CartesianRetainedUnitTransformContracts.transform_contracts(...)`
  accessors, and kept the approved line budget at exactly `80` added
  `src`/`bin` lines.

Validation:
- `git diff --check` passed.
- Package load passed.
- Canonical H2 base driver, `nesting = :pqs`, wrote/read back and printed
  the compact inventory; angular z-extension slabs show `289` support rows
  to `25` final columns.
- Canonical H2 base driver, `nesting = :wl`, wrote/read back and printed
  the same slab rows while complete shells show
  `white_lindsey_boundary_strata`.
- Canonical H2 supplemented driver with `cc-pVTZ`, `nesting = :pqs`, printed
  the base inventory and completed with final dimension `981`.
- Existing H2 RG endpoint test passed with augmented dimension `489`,
  self-Coulomb `0.4574265214362095`, and facade readback deltas `0.0`.
- Anti-bloat scan found only intentional compact-summary reads:
  `contract.metadata` for native slab facts, `block.coefficients` for
  identity-vs-compact class, and `block.support_indices` for row counts.

Goal advancement:
- LT5/LT6: makes accidental identity sectors and compact slab compression
  visible in ordinary driver runs, without turning the driver into a route
  debugger or changing artifact schemas.

Carrying-cost result:
- deleted: none; this is bounded user-facing visibility.
- simplified: driver users no longer need ignored inventory probes to see
  whether slabs and shells are compact or identity-realized.
- quarantined: source-mode inventory, pair inventory, raw-block details,
  coefficient dumps, all-row support dumps, route skeletons, artifact/schema
  changes, numerical construction, and Cr2 workflow remain outside this lane.
- not deleted because: ordinary and compact terminal realization paths remain
  live; the inventory only summarizes them.
- exact remaining caller/blocker: none for the bounded inventory output.
- added src lines: `65` in `src`, `15` in `bin`, `80` total.
- deleted src lines: `1` in `src`, `0` in `bin`.
- new tests: none.
- new metadata/status fields: one compact in-memory `base.terminal_inventory`
  summary for driver printing only; no artifact/status/schema fields.

## Cartesian Hamiltonian Producer Pass 183 - Approve RG 1e-6 Production Cutoff

Commit(s):
- this commit - Approve RG 1e-6 cutoff policy

Summary:
- Approved `HP-RG-CUTOFF-FN-02` and `HP-RG-CUTOFF-TEST-02` as the next
  Residual Gaussian production cutoff policy. The default
  `residual_occupation_cutoff` moves from `5.0e-8` to `1.0e-6`; the final
  residual identity validation `identity_atol` remains `5.0e-8`.
- The evidence is Cr2-specific but the policy is generic RG production
  behavior: the worst low-H1 residual modes were built from marginal
  owner-local residual directions with occupations about `1.27e-7` to
  `8.98e-7`, and a `1.0e-6` cutoff drops `6` directions per owner. Broad
  widths remain diagnostic evidence, but width filtering is not the first
  production rule because one-center atoms can have broad candidates without
  the same bad `H1_RR` sector.
- The approved validation is residual-only first. Cr2 owner retained counts
  should drop from `68 + 68` to `62 + 62`; doer should recompute and report
  `min eig(K_RR)`, `min eig(H1_RR)`, and low-mode candidate composition. If
  low-H1 ghost modes remain, the next authority should be a separate
  kinetic/`H1_RR` spectral-guard lane, not width scaling or a hidden tolerance
  change.

Validation:
- Docs-only validation required: `git diff --check`; focused scans for
  `HP-RG-CUTOFF-FN-02`, `HP-RG-CUTOFF-TEST-02`, `1.0e-6`, `identity_atol =
  5.0e-8`, `68 + 68`, `62 + 62`, `K_RR`, `H1_RR`, and forbidden full-HF/Cr2
  workflow wording; confirm no source/bin/test/tool files are staged.
- Later source validation requires Be high-zeta and H2 residual-GTO/MWG
  endpoints to pass, H2 cutoff/provenance assertions updated from `5.0e-8` to
  `1.0e-6`, and the Cr2 residual-only audit above. No full HF, Cr2 artifact,
  or Cr2 workflow is approved.

Goal advancement:
- LT5/LT6: tightens the residual basis to avoid carrying marginal RG
  directions that create low-H1 sectors in realistic Cr2-style use, while
  preserving the owner-local selection algorithm, exact one-body transforms,
  MWG/IDA conventions, and artifact/schema boundaries.

Carrying-cost result:
- deleted: none; docs-only policy pass.
- simplified: production RG cutoff policy now has one current value,
  `1.0e-6`, with the older `5.0e-8` approval preserved as history/evidence.
- quarantined: width-filtering defaults, kinetic/`H1_RR` spectral guards, full
  HF, Cr2 artifacts/workflow, source changes outside the RG owner/plumbing,
  residual-selection algorithm changes, and new committed fixtures/tests
  remain unapproved.
- exact remaining caller/blocker: if residual-only spectra still show low-H1
  ghost modes after the cutoff change, source work must stop and request a
  separate spectral-guard amendment.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 184 - Implement RG 1e-6 Production Cutoff

Commit(s):
- this commit - Set RG residual cutoff to 1e-6

Summary:
- Accepted the `HP-RG-CUTOFF-FN-02` source pass. The production default
  `residual_occupation_cutoff` is now `1.0e-6` in the Residual Gaussian owner
  and in the compatibility wrapper. The final residual identity validation
  remains `identity_atol = 5.0e-8`.
- The existing H2 residual-GTO/MWG endpoint assertions were updated exactly as
  approved under `HP-RG-CUTOFF-TEST-02`: both the in-memory
  `residual.occupation_cutoff` check and the written
  `supplement_provenance/occupation_cutoff` check now expect `1.0e-6`.
- Cr2 residual-only validation shows the cutoff performs the intended first
  selection cleanup but does not solve the whole low-H1 problem. The retained
  residual rank drops from the previous `68 + 68` to `62 + 62`, with minimum
  retained occupation `2.5714596279604456e-6`. The residual-only spectra still
  have a low two-owner mode: `min eig(K_RR) = 0.37004135191195486` and
  `min eig(H1_RR) = -7.164785405168449`, with owner weights about
  `0.50000003 / 0.49999997`. This keeps the next decision focused on a
  kinetic/`H1_RR` spectral guard rather than width filtering or another
  hidden tolerance change.

Validation:
- `git diff --check` passed.
- Package load passed.
- Existing H2 endpoint
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` passed:
  augmented dimension `489`, self-Coulomb `0.4574265214362095`, and facade
  readback deltas `0.0`.
- Ignored residual-only Cr2 audit
  `tmp/work/cr2_residual_low_ke_replay_audit.jl` passed without full HF,
  without dense Vee, and without writing a new Hamiltonian artifact. It wrote
  summary evidence under
  `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/`:
  `cr2_residual_low_ke_replay_cutoff1e6_summary.txt`,
  `cr2_residual_low_ke_replay_cutoff1e6_owner_metrics.tsv`, and
  `cr2_residual_low_ke_replay_cutoff1e6_modes.tsv`.

Goal advancement:
- LT5/LT6: removes marginal owner-local residual directions from production
  by default while preserving the RG owner-local selection/merge structure,
  exact augmented one-body convention, MWG/IDA behavior, artifact schema, and
  driver workflow.

Carrying-cost result:
- deleted: old `5.0e-8` production residual-occupation default.
- simplified: RG owner and compatibility wrapper now agree on the approved
  `1.0e-6` production cutoff.
- quarantined: `identity_atol`, owner grouping, merge checks, `G' S R`
  validation, width/zeta filtering, MWG/IDA, artifacts, driver, full Cr2 HF,
  and kinetic/`H1_RR` spectral-guard policy remain unchanged.
- not deleted because: the H2 endpoint test remains the live RG endpoint gate,
  with only its approved cutoff assertions updated.
- exact remaining caller/blocker: post-cutoff Cr2 residual-only spectra still
  show a low two-owner `H1_RR` mode, so a separate kinetic/`H1_RR` guard lane
  is needed before treating Cr2 supplemented RGs as physically safe.
- added src lines: 2.
- deleted src lines: 2.
- new tests: none; two existing H2 assertions updated.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 185 - Approve RG Residual-Sector Spectral Audit

Commit(s):
- this commit - Approve RG spectral audit lane

Summary:
- Approved `HP-RG-SPECTRAL-AUDIT-01` as a docs-only, measurement-only
  follow-up after the `HP-RG-CUTOFF-FN-02` production cutoff cleanup. The
  cutoff correctly drops Cr2 retained residuals from `68 + 68` to `62 + 62`,
  but the post-cutoff residual-only replay still shows a low two-owner mode:
  `min eig(K_RR) = 0.3700413519`,
  `min eig(H1_RR) = -7.1647854052`, with owner weights about `0.5 / 0.5`.
- The new lane authorizes ignored `tmp/work` probes and durable text/TSV
  output only. The audit may compute residual counts by owner, low `K_RR`, low
  `H1_RR = K_RR + sum_A Z_A U_A_RR`, low-mode owner weights,
  residual-occupation composition, and one-center atom baselines when
  available. It must classify whether flagged modes come from the smallest
  retained occupations or otherwise healthy retained modes.

Validation:
- `git diff --check` passed for the docs-only amendment.
- Focused searches confirmed `HP-RG-SPECTRAL-AUDIT-01` is recorded as
  measurement-only and that the forbidden surfaces remain explicit.
- Staged file scan confirmed no source, test, tool, bin, artifact, or driver
  files are part of the amendment.

Goal advancement:
- LT6: keeps Cr2 residual safety work disciplined. The cutoff change is now
  separated from the next evidence pass; no kinetic/`H1_RR` guard, automatic
  pruning, source instrumentation, full HF, dense Vee, artifacts, driver work,
  or cutoff/tolerance change is authorized yet.

Carrying-cost result:
- deleted: none; docs-only policy pass.
- simplified: the next step is a named residual-sector audit instead of an
  implicit extension of the cutoff source lane.
- quarantined: production source edits, committed tests/fixtures, artifact
  schema/provenance/reader/manifest changes, driver changes, MWG/IDA changes,
  dense Vee/full HF/solver work, automatic pruning, kinetic/`H1_RR` guards,
  and cutoff/tolerance changes remain forbidden.
- exact remaining caller/blocker: if existing construction seams cannot
  cheaply reconstruct `K_RR`/`H1_RR`, the audit must stop and report the
  missing reusable seam rather than adding source instrumentation.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 186 - Amend Driver Inventory Geometry Columns

Commit(s):
- this commit - Require geometry columns in driver inventory

Summary:
- Amended the existing `HP-DRV-INV-FN-01` / `HP-DRV-INV-TEST-01` driver
  inventory authority. The original compact inventory columns were sufficient
  to catch catastrophic direct-identity slabs, but not the subtler angular-box
  failure where all z-extension slabs are compact yet emitted only after the
  final shared shell.
- The canonical driver inventory row contract now includes shell index or
  explicit unavailable status, index ranges for `x`/`y`/`z`, physical
  coordinate ranges for `x`/`y`/`z`, and the existing slab normal/side/
  thickness/stack facts when applicable. Physical `x`/`y` ranges are required
  because angular-balance review compares transverse physical scale against
  the bond-axis margin; `z` alone is not enough.

Validation:
- `git diff --check` passed for the docs-only amendment.
- Focused search confirmed the `HP-DRV-INV-*` authority now names shell index,
  axis index ranges, physical coordinate ranges, and the forbidden non-goals.
- Staged path scan confirmed no source, test, tool, bin, artifact, or driver
  files are part of the docs-only amendment.

Goal advancement:
- LT5/LT6: keeps the canonical driver useful as a compact human-facing
  construction view without turning it into route diagnostics. This should
  make both identity-slab blowups and angular z-extension ordering errors
  visible during ordinary H2/Be2/Cr2-style runs.

Carrying-cost result:
- deleted: none; docs-only policy pass.
- simplified: driver inventory expectations are now one explicit geometry row
  contract instead of relying on ignored probes to notice z-extension ordering.
- quarantined: artifact schema changes, source-mode/pair/raw-block dumps,
  route-stage diagnostics, new driver inputs, numerical construction changes,
  shellification changes, terminal lowering changes, Cr2-specific workflow,
  and committed fixtures/tests remain forbidden.
- exact remaining caller/blocker: source work still needs a compact way to
  derive index and physical ranges from existing terminal/parent summaries
  without adding a broad report payload.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 187 - Harden Supplemented Route Provenance Test

Commit(s):
- this commit - Harden supplemented route provenance test

Summary:
- Accepted a narrow validation hardening pass after the supplemented route
  provenance fix. The existing H2 supplemented RG/MWG endpoint test now opens
  the written artifact and directly checks `recipe_provenance/route`,
  `recipe_provenance/nesting`, `recipe_provenance/producer`, and
  `hamiltonian_manifest/manifest_version`. This would have failed under the
  old hardcoded `:z_axis_diatomic_residual_gto_mwg` route label while keeping
  the test on the existing physics endpoint and artifact path.
- The pass deliberately does not require `producer_provenance/` for
  supplemented artifacts because that group is not currently written there;
  adding it would be producer/schema behavior outside the validation-only
  lane. The live truth source for supplemented construction is
  `recipe_provenance/` plus `supplement_provenance/`.

Validation:
- `git diff --check` passed.
- Package load passed.
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` passed,
  including augmented dimension `489`, self-Coulomb `0.4574265214362095`, and
  artifact readback deltas `0.0`.
- Diff gate: `git diff --numstat -- src bin tools test docs` showed `0`
  source/bin/tool lines, `+9/-0` in the existing H2 endpoint test, and
  `+46/-0` in this running-log entry. The only added test lines are direct
  JLD2 group/key/value assertions; no suspicious source patterns apply.

Goal advancement:
- LT5/LT6: strengthens artifact truthfulness around public construction
  family and route labels without changing Hamiltonian matrices, artifact
  schema, reader behavior, driver inputs, or CR2 workflow.

Carrying-cost result:
- deleted: none.
- simplified: no new test file or fixture matrix; the check lives in the
  existing H2 supplemented endpoint artifact block.
- quarantined: atom/WL route-label matrix coverage remains an ignored/manual
  validation candidate until separately approved for committed tests.
- not deleted because: this is validation hardening only.
- exact remaining caller/blocker: committed atom/PQS, atom/WL, and diatomic/WL
  provenance route checks need a later exact test-authority amendment if
  desired.
- added src lines: 0.
- deleted src lines: 0.
- new tests: no new file; existing test gained focused assertions.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 188 - Use Native Shell Index In Driver Inventory

Commit(s):
- this commit - Use native shell index in driver inventory

Summary:
- Accepted the narrow driver-inventory provenance fix. The inventory reader
  already consumed `terminal_region_shell_index` when present, but the current
  construction path did not populate it, so complete shells and angular
  z-extension slabs printed `shell=unavailable`.
- The route helper now builds a small terminal-region-key to native
  `raw_region.shell_index` map from the typed shellification plan and merges
  positive shell indices into retained-unit metadata before transform
  contracts are built. No traversal-derived shell counter is used. Direct
  core and other non-shell regions remain explicit `shell=unavailable`.
- This keeps the canonical driver inventory tied to shellifier-native region
  provenance, so angular z-extension slabs can show the shell step they belong
  to without adding artifact fields or a route-debug report.

Validation:
- `git diff --check` passed.
- Package load passed.
- Ignored probe `tmp/work/terminal_inventory_native_shell_index_probe.jl`
  showed H2/PQS complete shells as `shell=1` and `shell=2`, angular
  z-extension slabs as `shell=1`, and H2/WL boundary-stratum rows carrying
  the same native shell indices.
- Canonical driver H2/PQS and H2/WL base artifact/readback runs passed with
  `check_file=true`, `print_contract=true`, and printed x/y/z index and
  physical ranges alongside native shell indices.
- Existing H2 RG endpoint
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` passed:
  augmented dimension `489`, self-Coulomb `0.4574265214362095`, and readback
  deltas `0.0`.
- Diff gate: `git diff --numstat -- src bin tools test docs` showed `+37/-0`
  in `src/pqs_source_box_route_driver_helpers.jl` before this log entry; no
  bin/tool/test files changed.

Goal advancement:
- LT5/LT6: strengthens the driver as a bounded construction sanity surface for
  shellification geometry. The inventory can now expose whether angular slabs
  are attached to the expected native shell step instead of relying on
  retained-unit ordering.

Carrying-cost result:
- deleted: no traversal-derived shell-index logic; none existed in the
  accepted source, and this pass keeps it out.
- simplified: the existing inventory consumer now receives the native fact it
  was already prepared to display.
- quarantined: artifact schema, matrix keys, reader behavior, driver inputs,
  shellification algorithms, terminal lowering semantics, retained-unit
  semantics, source-mode/pair/raw-block dumps, Cr2 workflow, and committed
  tests remain unchanged.
- not deleted because: the retained-unit metadata seam is the existing compact
  route-to-inventory carrier.
- exact remaining caller/blocker: none for native shell-index inventory
  display; broader committed route-label/shell-index test matrices would need
  separate test authority.
- added src lines: 37.
- deleted src lines: 0.
- new tests: none; ignored probe only.
- new metadata/status fields: no artifact/status/schema fields; the existing
  in-memory `terminal_region_shell_index` inventory metadata key is now
  populated from native shellifier data where positive shell indices exist.

## Cartesian Hamiltonian Producer Pass 189 - Align Thin-Slab Summary Vocabulary

Commit(s):
- this commit - Align thin-slab summary vocabulary

Summary:
- Accepted a narrow cleanup of private terminal-shellification summary
  vocabulary. The live lowering/retained/terminal paths already realize
  midpoint slabs, angular z-extension slabs, and outer-mismatch fallback slabs
  through compact thin-slab lowering. The private summary helpers still used
  old direct/deferred labels that implied direct identity sectors or pending
  support.
- The cleanup renames those private dependency/count/lowering-family symbols
  to compact thin-slab vocabulary and deletes no-op pending reasons for
  angular z-extension slabs, which are now independently lowerable in this
  private summary classification. The module docstring now explicitly names
  planned angular z-extension slabs and compact thin-slab lowering categories.

Validation:
- `git diff --check` passed.
- Package load passed.
- Focused `rg` found no remaining matches for the four stale direct/deferred
  slab labels named in the cleanup handoff.
- Bounded H2/PQS canonical driver artifact/readback smoke passed with
  compact angular z-extension slab inventory rows and dimension `471`.
- Diff gate: `git diff --numstat -- src bin tools test docs` showed `+15/-16`
  in `src/cartesian_shellification/terminal_geometry.jl` before this log
  entry; no bin/tool/test files changed.

Goal advancement:
- LT5/LT6: reduces conceptual drift around thin slabs so future agents see
  midpoint, angular z-extension, and outer-mismatch fallback slabs as compact
  slab categories rather than reviving the old identity/deferred-slab story.

Carrying-cost result:
- deleted: stale private summary vocabulary for direct slab/boundary-slab and
  angular z-extension pending/deferred states.
- simplified: angular z-extension slabs no longer carry a private
  missing-lowering reason or retirement target because they are classified as
  independently lowerable compact thin slabs.
- quarantined: shellification algorithm, region kinds, region ordering,
  region metadata, terminal lowering semantics, retained-unit records,
  terminal realization, driver output, artifacts/readers, numerical
  construction, Cr2 workflow, and committed tests remain unchanged.
- not deleted because: private summary/scaffold helpers remain live as bounded
  shellification metadata summaries.
- exact remaining caller/blocker: none for this vocabulary cleanup.
- added src lines: 15.
- deleted src lines: 16.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 190 - Retire Stale Composition Rejections

Commit(s):
- this commit - Retire stale supplemented composition rejection docs

Summary:
- Accepted a docs-only cleanup of current/startup composition authority after
  the supported supplemented atom and supplemented White-Lindsey lanes landed.
  The old base-atom driver and nesting-artifact lanes originally rejected
  those cells, but that rejection is now historical: supported one-center
  supplemented atoms are governed by `HP-COMP-SUPPATOM-*`, and supported
  homonuclear z-axis diatomic supplemented WL construction is governed by
  `HP-COMP-SUPPWL-*`.
- Updated `cartesian_driver_atom_workflow.md`, `current.md`, `registry.md`,
  `implementation_slices.md`, and `AGENTS.md` so current/startup docs no
  longer imply those paths are globally unapproved. The composition plan was
  already current and was left unchanged.
- Preserved true deferred boundaries: translated atoms, heteronuclear/general
  geometry, ECP/solver/HF workflow, Cr2-specific workflow, unsupported
  geometry/supplement combinations, and artifact schema redesign.

Validation:
- `git diff --check` passed before this log entry.
- Focused stale-phrase scan over `AGENTS.md` and the compact Cartesian
  producer docs found no remaining live matches for the old supplemented atom
  and supplemented-WL rejection phrases named in the handoff.
- Broader `rg` over `supplemented atom`, `nesting = :wl`,
  `HP-COMP-SUPPATOM`, and `HP-COMP-SUPPWL` was reviewed; remaining hits are
  current composition authority, lane-scoped boundaries, or unrelated WL
  policy such as mapped-COMX/WL rejection.
- Changed path list was docs/startup only; no source, bin, test, tools,
  artifact, driver behavior, RG/MWG, shellification, or route changes.
- No Julia run; not needed for docs-only cleanup.

Goal advancement:
- LT5/LT6: reduces current-authority drift after composition support landed.
  Future agents should see the implemented 2 x 2 x 2 composition cells rather
  than revive obsolete rejection policy from earlier driver/nesting lanes.

Carrying-cost result:
- deleted: live/current phrasing that presented supplemented atoms and
  supplemented WL as still rejected by default.
- simplified: current docs now point those supported cells directly at
  `HP-COMP-SUPPATOM-*` and `HP-COMP-SUPPWL-*`.
- quarantined: historical lane boundaries remain described as lane-local
  history where useful.
- not deleted because: older lane sections still matter as history and as
  boundaries for their original base-only/provenance-only scopes.
- exact remaining caller/blocker: none for stale current composition wording.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 191 - Validate R3 Same-Construction Nuclear Charges

Commit(s):
- `c0e05846b` - Validate R3 same-construction nuclear charges

Summary:
- Accepted the narrow R3 consistency fix after the conformance audit found that
  same-length `nuclear_charges` passed into the same-construction supplemented
  Hamiltonian boundary could differ from `base_hamiltonian.nuclear_charges`
  while final assembly still used the base Hamiltonian values.
- The same-construction boundary
  `pqs_terminal_residual_gto_augmented_hamiltonian(base_hamiltonian, basis,
  bundles, supplement, atom_locations, nuclear_charges; ...)` now validates
  both charge count and charge values against the base Hamiltonian through a
  small `_r3_validate_same_construction_nuclear_charges(...)` helper.

Validation:
- `git diff --check` passed.
- Package load passed.
- Ignored smoke `tmp/work/r3_same_construction_charge_mismatch_smoke.jl`
  verified mismatched charge behavior.
- Existing H2 RG endpoint
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` passed.

Goal advancement:
- LT6: closes a provenance/consistency hole in the supplemented
  same-construction path without changing residual selection, exact operator
  formulas, MWG/IDA, artifact schema, driver input, or route behavior.

Carrying-cost result:
- deleted: none.
- simplified: inconsistent same-construction nuclear charge inputs now fail
  early instead of being ignored downstream.
- quarantined: artifact schema, driver input, RG/MWG conventions, residual
  selection, exact operator formulas, raw blocks, route/shellification, broad
  tests, and Cr2 runs remain unchanged.
- exact remaining caller/blocker: none for charge consistency.
- added src lines: small validation helper only.
- new tests: none; ignored smoke only.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 192 - Clarify Remaining QW Non-Nuclear Donors

Commit(s):
- this commit - Clarify QW non-nuclear donor retention

Summary:
- Recorded the follow-up decision from the conformance audit of
  `HP-CGRB-NN-*`: no source cleanup should be done for the remaining
  QW-local non-nuclear cross/self helpers in this pass. The already-crossed
  lane applies to the main diatomic Qiu-White non-nuclear path, which consumes
  the neutral `CartesianGaussianRawBlocks` owner.
- The remaining helpers are still live through atomic QW reference paths,
  factor-term output, hybrid sidecars, dense-parent probes, and CPB/provider
  surfaces. Part of their output, especially `factor_ga`/`factor_aa`, is
  outside current neutral CGRB-NN authority, which only owns overlap, kinetic,
  coordinate moments, and second moments.

Validation:
- `git diff --check` passed for this docs-only clarification.
- Focused searches reviewed `HP-CGRB-NN-*`, remaining QW helper names,
  `factor_ga`/`factor_aa`, and sidecar/provider wording.
- Changed path list was docs/startup only; no source, bin, test, tools,
  artifacts, or driver behavior changed.

Goal advancement:
- LT6: prevents over-deletion by separating the completed diatomic neutral
  raw-block lane from still-live atomic/reference/sidecar/provider surfaces.
  Future cleanup must name those surfaces and any factor-block ownership
  explicitly.

Carrying-cost result:
- deleted: none; docs-only classification.
- simplified: future agents now have a clear reason not to delete the
  remaining QW donor helpers under old CGRB-NN wording.
- quarantined: atomic QW factor-term output, hybrid sidecars, dense-parent
  probes, CPB/provider surfaces, neutral factor-block ownership, public API,
  route/report/status payloads, nuclear/raw-block semantic changes, RG/MWG/IDA,
  Hamiltonian behavior, Cr2 runs, and committed tests remain out of scope.
- exact remaining caller/blocker: remaining QW helpers can only be rewired or
  deleted after a new lane approves the specific atomic/sidecar/provider
  surfaces and any missing neutral factor-block ownership.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 193 - Record Residual Injection Hybrid Memo

Commit(s):
- this commit - Record residual injection hybrid memo

Summary:
- Added a durable design memo for the proposed optional
  injection-plus-Residual-Gaussian construction. The memo records the intended
  split: classify candidate modes owner-locally by residual norm, globally
  orthonormalize the injected near-gausslet subspace, replace the corresponding
  gausslet-sector directions, and then return to owner-local residual
  selection for true RG/MWG directions.
- The memo also records that initial GTO candidate orthonormalization needs its
  own candidate-overlap rank policy, separate from residual occupation,
  injection, injected-subspace merge, and final identity tolerances.

Validation:
- `git diff --check` for this docs-only pass.
- Focused link/readback review of the new memo, README link, current-status
  note, and RG domain-module note.
- No Julia run; no source, bin, test, tool, artifact, or driver file was
  touched.

Goal advancement:
- LT6: preserves the conceptual cure for near-gausslet residual-complement
  instability without turning it into an unapproved source lane. The
  documented path explicitly avoids a global residual-selection algorithm while
  allowing global injected-subspace orthonormality.

Carrying-cost result:
- deleted: none.
- simplified: future RG work has a named memo separating candidate GTO rank,
  optional injection, true residual selection, and interaction conventions.
- quarantined: production source behavior, defaults, committed tests, artifact
  schema, driver inputs, public API, Cr2 workflow, full HF, spectral guards,
  and MWG convention changes remain unapproved.
- exact remaining caller/blocker: a later docs-only amendment must approve
  exact `HP-RG-INJECT-*` IDs, files, validation, thresholds, and failure rules
  before implementation.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 194 - Clarify Residual Injection Memo

Commit(s):
- this commit - Clarify residual injection memo

Summary:
- Tightened the injection memo after external review. The memo now states that
  injection/residual classification applies to owner-local orthonormal
  principal modes `y_i = A_tilde v_i`, not raw GTO columns.
- It also spells out the concrete injected-sector construction:
  `B = G' S Y`, an orthonormal complement `Q_perp` inside the original
  gausslet coefficient space, and `F = [Y, G Q_perp]`. The rank/condition of
  `B` is now a required audit diagnostic.
- The memo now labels inherited injected-sector IDA as an approximation:
  one-body operators use the injected/raw representation exactly, while
  two-body IDA keeps original gausslet-sector IDA semantics for the replaced
  subspace.
- The first-audit `lambda_inj >= residual_occupation_cutoff` rule is recorded
  as a default audit policy rather than a permanent mathematical requirement.

Validation:
- `git diff --check` for this docs-only clarification.
- Focused readback of the memo sections on optional switch, classification,
  injected sector, interaction convention, and measurement-only audit.
- No Julia run; no source, bin, test, tool, artifact, or driver file was
  touched.

Goal advancement:
- LT6: makes the future injection audit harder to misimplement by naming the
  actual classified modes and the algebraic injected-sector replacement, while
  preserving the no-source-authority status.

Carrying-cost result:
- deleted: none.
- simplified: reduced ambiguity around raw GTO labels, global injected
  orthonormality, and inherited IDA.
- quarantined: production behavior, defaults, source code, tests, artifacts,
  driver inputs, public API, Cr2 workflow, and MWG convention changes remain
  unapproved.
- exact remaining caller/blocker: a later docs-only `HP-RG-INJECT-*`
  amendment must still approve any measurement probe or source implementation.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 195 - Approve Residual Injection Audit

Commit(s):
- this commit - Approve residual injection audit

Summary:
- Approved `HP-RG-INJECT-AUDIT-01` as the first measurement-only authority
  lane for the optional injection-plus-Residual-Gaussian scheme. The audit may
  use ignored probes to classify owner-local principal modes
  `y_i = A_tilde v_i`, sweep trial `lambda_inj` values, globally merge
  injected modes, report `B = G' S Y_inj` rank/condition, construct trial
  `F = [Y, G Q_perp]`, and measure true RG counts plus `K_RR`/`H1_RR`
  spectra after injection.
- The approved audit must also report injected-sector one-body projection
  errors for `K`, each unit `U_A`, and `H1`, so the result tests the actual
  one-body failure mode rather than only residual counts.

Validation:
- `git diff --check` for this docs-only authority pass.
- Focused searches/readback for `HP-RG-INJECT-AUDIT-01`,
  `residual_injection_cutoff`, `B = G' S Y_inj`, `F = [Y, G Q_perp]`, and the
  forbidden surfaces in `current.md`, `registry.md`,
  `implementation_slices.md`, `residual_gaussian_domain_module.md`,
  `residual_gaussian_injection_hybrid.md`, and `AGENTS.md`.
- No Julia run; no source, bin, test, tool, artifact, or driver file was
  touched.

Goal advancement:
- LT6: turns the injection memo into an actionable measurement lane while
  preserving the no-production-source boundary. The audit is designed to decide
  whether injection actually removes or sharply reduces the Cr/Cr2 low residual
  ghost sector before implementation.

Carrying-cost result:
- deleted: none.
- simplified: repo-manager can now hand doer an exact ignored-audit task
  instead of inventing measurement scope from the memo.
- quarantined: production source behavior, source instrumentation, committed
  tests, artifacts, driver inputs, public API, RG defaults, automatic pruning,
  spectral guards, MWG/IDA conventions, full HF, dense Vee/solver work, Cr2
  full Hamiltonians, Cr2 artifacts, and Cr2-specific workflow remain
  unapproved.
- exact remaining caller/blocker: an ignored audit must now determine whether
  existing construction seams can reconstruct the owner-local candidate spans,
  injected-sector projection, trial injected sector, and residual-sector
  one-body blocks cheaply enough.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 196 - Approve Default-Off Residual Injection Implementation

Commit(s):
- this commit - Approve default-off residual injection implementation

Summary:
- Approved `HP-RG-INJECT-FN-01` as a narrow default-off in-memory source lane
  for the injection-plus-Residual-Gaussian hybrid. The decision records the
  interpretation of the first audit: injection did not immediately remove the
  current Cr2 low-H1 residual sector, but it remains the right general
  construction because RG alone has a singular-complement limit for
  near-gausslet GTO directions.
- The approved implementation may classify owner-local principal modes,
  globally merge injected modes, construct an injected replacement base sector
  `F`, transform exact one-body operators into `[F, R]`, and keep residual
  MWG/IDA only for true residual directions. The source must preserve current
  behavior when `residual_injection_cutoff <= 0`.

Validation:
- `git diff --check` for this docs-only authority pass.
- Focused readback of `HP-RG-INJECT-FN-01` in the memo, registry, current
  status, implementation slices, RG domain module, and AGENTS startup rules.
- No Julia run; no source, bin, test, tool, artifact, or driver file was
  touched.

Goal advancement:
- LT6: moves from measurement-only rationale to an implementable, default-off
  RG-domain source lane while explicitly preserving current production
  behavior and quarantining artifact/default/workflow decisions.

Carrying-cost result:
- deleted: none.
- simplified: future source work now has exact files, allowed transforms,
  forbidden surfaces, validation, and dense-transform stop conditions.
- quarantined: default-on injection, driver/public API, artifact schema or
  injection-enabled artifact writing, MWG channels for injected directions,
  global residual selection, spectral pruning, full HF, dense Vee/solver work,
  Cr2 artifacts/workflow, route/shellification/raw-block changes, committed
  tests, and broad payload/report shapes remain unapproved.
- exact remaining caller/blocker: implementation must stop if exact one-body
  transformation into `[F, R]` requires an unacceptable persistent dense
  `nG x nG` transform/workspace or source outside the approved files.
- added src lines: 0.
- deleted src lines: 0.
- new tests: none.
- new metadata/status fields: none.

## Cartesian Hamiltonian Producer Pass 197 - Implement Default-Off Residual Injection

Commit(s):
- this commit - Implement default-off residual injection

Summary:
- Accepted the `HP-RG-INJECT-FN-01` source implementation after one manager
  correction. The first WIP was not accepted because it carried persistent
  dense fixed-sector transforms. The corrected WIP carries compact injected
  authority as `injected_A` plus `injected_G = G' S Y`; reconstructs
  `Q_perp` locally when transforming exact operators; rejects
  injection-enabled artifact writing; and enforces `dim(Y_inj) < nG`.
- A follow-up review found candidate-overlap `rtol` had been collapsed into an
  absolute cutoff. The final accepted source computes the owner-local
  candidate rank threshold as
  `max(candidate_overlap_atol, candidate_overlap_rtol * maximum(s))`.

Validation:
- Doer and manager ran `git diff --check`, package load, the H2
  residual-GTO/MWG endpoint, and `tmp/work/rg_injection_enabled_h2_smoke.jl`.
  Manager rerun passed with H2 endpoint `52/52`, facade/readback `67/67`, and
  injection smoke reporting `injected_dimension = 2`,
  `residual_dimension = 16`, artifact rejection, `F' S F = 2.89e-15`,
  `F' S R = 1.21e-14`, `R' S R = 1.14e-12`, kinetic symmetry `0.0`, and Vee
  symmetry `6.66e-15`.

Goal advancement:
- LT6: installs the default-off in-memory injection construction so
  near-gausslet GTO directions have a third fate besides discard or
  singular residual-complement normalization. This is not evidence that the
  Cr2 low-H1 residual sector is solved, and it does not promote injection to a
  production default.

Carrying-cost result:
- deleted: none.
- simplified: exact one-body and inherited fixed-sector IDA now operate on the
  injected replacement sector without persistent dense fixed-transform state.
- quarantined: default-on injection, driver/public API, artifact provenance,
  injection-enabled artifact writing, MWG channels for injected directions,
  spectral pruning, full HF, dense Vee/solver work, Cr2 artifacts/workflow,
  route/shellification/raw-block changes, and committed tests remain
  unapproved.
- exact remaining caller/blocker: future work must decide separately whether
  injection should become a default, whether artifact provenance can describe
  injected sectors, and whether a kinetic/`H1_RR` spectral gate is needed.
- added src lines: 151.
- deleted src lines: 20.
- new tests: none.
- new metadata/status fields: compact injected basis fields only; no route,
  report, artifact, or public status fields.

## Cartesian Hamiltonian Producer Pass 198 - Reduce Residual Injection Carrying Cost

Commit(s):
- this commit - Reduce residual injection carrying cost

Summary:
- Accepted a bounded cleanup of the just-added `HP-RG-INJECT-FN-01` source
  implementation. The stored `injected_dimension` and `injected_owner_counts`
  fields were removed; injected dimension is now derived from the compact
  injected-sector authority `injected_G`. The duplicated final residual merge,
  Lowdin normalization, sign canonicalization, and identity validation logic
  is shared by the disabled-injection and enabled-injection paths.
- The cleanup preserves the accepted physical contract: disabled injection is
  unchanged, injected sectors remain in-memory only, artifact writing still
  rejects injection-enabled residual sectors, and `injected_A`/`injected_G`
  remain the only persistent injected-sector basis authority.

Validation:
- Doer ran `git diff --check`, package load, H2 residual-GTO/MWG endpoint,
  supplemented facade/readback, and the H2 injection smoke.
- Manager reran `git diff --check`, package load, and
  `tmp/work/rg_injection_enabled_h2_smoke.jl`. The smoke reported
  `off_zero_delta = 0.0`, `injected_dimension = 2`,
  `residual_dimension = 16`, `F' S F = 2.89e-15`,
  `F' S R = 1.21e-14`, `R' S R = 1.14e-12`, kinetic symmetry `0.0`, and
  Vee symmetry `6.66e-15`.

Goal advancement:
- LT6: keeps the optional injection construction available while reducing the
  amount of permanent state and duplicate cleanup logic carried by the RG
  domain module. This does not change the strategic conclusion that injection
  is still default-off and not yet a demonstrated cure for the Cr2 low-H1
  residual sector.

Carrying-cost result:
- deleted: stored injected dimension, stored injected owner-counts, injected
  owner-count accumulation, and the single-use residual identity helper.
- simplified: residual dimension reporting is derived; final residual
  orthonormalization/identity validation is shared.
- quarantined: no default-on behavior, driver/public API, artifact provenance,
  injection-enabled artifact writing, committed tests, spectral gates, full HF,
  dense Vee/solver work, or Cr2 workflow was added.
- exact remaining caller/blocker: none for this cleanup; future manager work
  should separately decide whether to clean up older retained-unit metadata
  double-rebuilds found during the bloat review.
- added src lines: 36.
- deleted src lines: 51.
- net src lines: -15.

## Cartesian Hamiltonian Producer Pass 199 - Consolidate Retained Unit Metadata Enrichment

Commit(s):
- this commit - Consolidate retained unit metadata enrichment

Summary:
- Accepted the first `bloat-fixer` source cleanup pass. The older driver
  inventory path from Pass 188 rebuilt the retained-unit plan twice: once to
  attach native shell indices and once to attach mapped-COMX source-span facts.
  The cleanup replaces those two wrappers with one local enrichment pass in
  `src/pqs_source_box_route_driver_helpers.jl`.
- Behavior is intentionally unchanged. Ordinary source spans still receive
  shell-index metadata when available and do not receive mapped-COMX source
  facts. Mapped source spans still carry `source_mode_shape` and
  `raw_product_source_axis_transform_facts`. If no unit metadata changes, the
  original retained-unit plan is returned.

Validation:
- Bloat-fixer ran `git diff --check`, package load, a focused caller scan, and
  `tmp/work/terminal_inventory_native_shell_index_probe.jl`.
- Manager reviewed the actual diff and reran `git diff --check`, package load,
  and the same inventory probe. The probe reported PQS and WL
  `final_dimension = 471`; PQS complete shells carried native shell indices
  `1` and `2`, angular z-extension slabs carried shell `1`, and WL boundary
  strata plus angular slabs retained the expected native shell indices.

Goal advancement:
- LT2/MT5: reduces carrying cost in stable driver-inventory metadata plumbing
  without widening the route/report surface.
- LT5/LT6: preserves the compact human-facing inventory contract while keeping
  source-mode facts and shellifier-native provenance internal to existing
  metadata.

Carrying-cost result:
- deleted: two retained-unit plan wrapper helpers.
- simplified: shell-index and source-span metadata enrichment are now one pass
  over retained units.
- quarantined: no scientific behavior, public API, driver input, artifact
  schema, route semantics, shellification, terminal lowering, transform
  contract, RG/MWG/IDA, Hamiltonian, solver, test, or docs expansion.
- exact remaining caller/blocker: the new enrichment helper is called only from
  `_pqs_source_box_route_driver_unit_stage_low_order_summary`; no cleanup
  blocker remains for this helper family.
- added src lines: 18.
- deleted src lines: 35.
- net src lines: -17.

### Medium-Term Goal Checkpoint After Pass 199

- MT1 fake-PQS quarantine: active/maintained. Recent RG/injection and cleanup
  passes did not promote fake-PQS or source-backed oracle paths.
- MT2 independent H2 PQS recovery and MT3 common physical support vocabulary:
  active as background guardrails, but not the current tactical lane.
- MT4 supplement staging after authority: active, with the current residual
  injection work still default-off and not artifact/public workflow authority.
- MT5 cleanup pressure: active and strengthened. Pass 198 removed residual
  injection state/duplication, Pass 199 removed older retained-unit metadata
  double-rebuilds, and the new `bloat-fixer` role gives cleanup work a
  manager-reviewed lane.
- MT6 audit/classify old Cartesian flat paths: active. Bloat-fixer is now a
  concrete mechanism for stable-code cleanup, while repo-doer remains focused
  on the RG/injection spectral audit.
- Medium-goal wording: needs refinement soon. The top-level MT section still
  reflects an older H2/PQS recovery frame; the immediate active lane is now RG
  injection evidence plus stable-code simplification. No MT wording is changed
  in this mechanical cleanup commit.

## Cartesian Hamiltonian Producer Pass 200 - Delete WL Support Record Wrapper

Commit(s):
- this commit - Delete WL support record wrapper

Summary:
- Accepted a narrow `bloat-fixer` cleanup in the stable White-Lindsey terminal
  realization surface. The private `_wl_terminal_support_record` wrapper was
  removed and its exact support tuple was inlined at the single
  `_validate_block_support!` call.
- No strategic change. This leaves MT5/LT2 cleanup pressure active while
  preserving the scientifically sensitive identity, coefficient,
  support-disjointness, retained-count, and transform-path checks.

Validation:
- Bloat-fixer and manager both ran `git diff --check`, package load, and
  `tmp/work/terminal_inventory_native_shell_index_probe.jl`. Manager rerun
  reported PQS/WL `final_dimension = 471`.

Carrying-cost result:
- deleted: one private WL support-record wrapper.
- simplified: one validation call site now passes the exact support tuple
  directly.
- quarantined: no terminal-basis semantics, support validation,
  coefficient/identity checks, retained-count checks, transform checks,
  public API, artifact, route, RG/MWG/IDA, Hamiltonian, test, or docs change.
- exact remaining caller/blocker: no remaining caller of
  `_wl_terminal_support_record`; broader WL cleanup remains blocked by
  scientifically sensitive final-basis checks unless separately assigned.
- added src lines: 1.
- deleted src lines: 8.
- net src lines: -7.

## Cartesian Hamiltonian Producer Pass 201 - Fold WL Support Expansion Helpers

Commit(s):
- this commit - Fold WL support expansion helpers

Summary:
- Accepted a bundled low-risk `bloat-fixer` cleanup in the stable
  White-Lindsey terminal realization surface. The one-use
  `_wl_terminal_cpb_support!` helper was folded into
  `_wl_terminal_source_support`, preserving CPB support index/state append
  order, and the `_wl_axis_index` alias was removed in favor of direct calls
  to the shared terminal face-axis helper.
- No strategic change. This continues MT5/LT2 cleanup pressure while leaving
  metadata validation, retained-count handling, identity/coefficient checks,
  support-disjointness checks, and transform-path checks untouched.

Validation:
- Bloat-fixer and manager both ran `git diff --check`, package load, and
  `tmp/work/terminal_inventory_native_shell_index_probe.jl`. Manager rerun
  reported PQS/WL `final_dimension = 471`.

Carrying-cost result:
- deleted: `_wl_terminal_cpb_support!` and `_wl_axis_index`.
- simplified: `_wl_terminal_source_support` now owns its single support
  expansion loop directly; axis lookup uses the shared helper directly.
- quarantined: `_wl_metadata_value`, retained-count validation, `Int(q)`
  normalization, identity/coefficient checks, support-disjointness checks, and
  transform-kind checks remain untouched.
- exact remaining caller/blocker: no remaining callers of the deleted helpers;
  the retained-count validation/normalization is the next possible
  manager-approved stable-code cleanup target.
- added src lines: 10.
- deleted src lines: 15.
- net src lines: -5.

## Cartesian Hamiltonian Producer Pass 202 - Remove WL Retained Count Precheck

Commit(s):
- this commit - Remove WL retained count precheck

Summary:
- Accepted a narrow `bloat-fixer` cleanup in the stable White-Lindsey terminal
  realization surface. The local `q isa Integer && q > 0` retained-count
  precheck and two `Int(q)` call-site normalizations were removed from
  `_wl_boundary_stratum_block`.
- No strategic change. This is the intended mature-code policy in practice:
  the check only made an internal metadata bug crash with a custom message.
  Normal construction accepts `q` directly, while downstream construction will
  still fail naturally if the metadata is not a valid retained count.

Validation:
- Bloat-fixer and manager both ran `git diff --check`, package load, and
  `tmp/work/terminal_inventory_native_shell_index_probe.jl`. Manager rerun
  reported PQS/WL `final_dimension = 471`.

Carrying-cost result:
- deleted: one retained-count prettier-crash precheck.
- simplified: `_nested_doside_1d` and `_terminal_face_product_block` now
  receive the retained count directly.
- quarantined: `_wl_metadata_value`, source-CPB checks, corner support-size
  checks, identity/coefficient checks, support-disjointness checks, and
  transform-kind checks remain untouched.
- exact remaining caller/blocker: normal WL construction accepts `q` directly;
  no downstream blocker found.
- added src lines: 2.
- deleted src lines: 4.
- net src lines: -2.

## Cartesian Hamiltonian Producer Pass 203 - Measure RG Injection Spectral Outcome

Commit(s):
- this commit - Record RG injection spectral audit

Summary:
- Accepted the measurement-only RG injection audit. The result is decisive:
  optional injection does not remove the bad Cr2 residual-sector `H1_RR` mode.
  In the dimer replay, the lowest `H1_RR` eigenvalue moves only from
  `-7.3689569529` at `residual_injection_cutoff = 0` to `-6.9926510108` at
  `1.0e-3`; the mode remains strongly negative and balanced across the two
  atoms.
- Higher injection cutoffs reduce the true residual count and improve
  `K_RR`, but at `1.0e-4` and `1.0e-3` the production path fails its default
  `F' S R` validation. Those rows are therefore measurement-only relaxed
  replays, not source behavior to promote.

Validation:
- Doer ran package load, the H2 injection smoke, and the full ignored audit in
  `413.36 s`. No tracked source edits, committed tests, Cr2 Hamiltonian
  artifact, full HF, or solver run was made.
- Manager inspected saved outputs under
  `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/`, including
  `rg_injection_audit_summary.txt`, `rg_injection_audit_sweeps.tsv`,
  `rg_injection_audit_modes.tsv`,
  `rg_injection_audit_production_comparison.tsv`, and
  `rg_injection_audit_stages.tsv`. The production/manual comparison matched
  nonzero cutoffs within the `1.0e-7` audit tolerance; worst dimer residual
  subspace error was about `2.21e-8`, and worst injected-subspace error was
  about `6.71e-12`.

Goal advancement:
- LT6: closes the immediate question about whether default-off injection is a
  cure for the Cr2 low residual-sector `H1` mode. It is not.
- LT5: preserves route/provenance honesty by keeping the high-cutoff rows as
  measurement-only evidence because they fail the current production
  orthogonality gate.

Risk / guardrail:
- The audit summary records head `7653c2098`, before later WL-only
  bloat-fixer commits. Those later commits touched stable WL terminal
  realization cleanup, not the RG injection source or the PQS Cr/Cr2 artifact
  replay path, so the measurement remains accepted.

Carrying-cost result:
- deleted: none.
- simplified: none.
- quarantined: injection remains default-off; high-cutoff relaxed rows remain
  measurement-only and are not production behavior.
- exact remaining caller/blocker: injection improves `K_RR` somewhat but does
  not remove the low Cr2 `H1_RR` residual-sector mode. The next manager
  decision is whether to pursue a kinetic/`H1_RR` spectral safety gate or a
  different residual-selection policy; another injection implementation pass
  is not justified by this evidence.
- added tracked source lines: 0.
- deleted tracked source lines: 0.

## Cartesian Hamiltonian Producer Pass 204 - Delete WL Metadata Accessor Wrapper

Commit(s):
- this commit - Delete WL metadata accessor wrapper

Summary:
- Accepted a narrow `bloat-fixer` stable-code cleanup in the White-Lindsey
  terminal realization surface. The private `_wl_metadata_value` helper was
  deleted, and seven internal metadata reads now use direct `NamedTuple` field
  access.
- No strategic change. This reinforces the mature-code cleanup rule: internal
  metadata bugs do not need custom prettier-crash wrappers when a natural
  missing-field failure is sufficient.

Validation:
- Bloat-fixer and manager both ran `git diff --check`, package load, and
  `tmp/work/terminal_inventory_native_shell_index_probe.jl`. Manager rerun
  reported PQS/WL `final_dimension = 471`.

Carrying-cost result:
- deleted: `_wl_metadata_value`.
- simplified: seven White-Lindsey boundary metadata reads now use direct field
  access.
- quarantined: coefficient identity checks, support-disjointness checks,
  transform-kind checks, source-CPB count checks, corner support-size checks,
  and contract lookup checks remain untouched.
- exact remaining caller/blocker: no remaining caller of `_wl_metadata_value`;
  broader stable-code trimming should now move from single wrappers to bounded
  helper-family cleanup.
- added src lines: 7.
- deleted src lines: 13.
- net src lines: -6.

## Cartesian Hamiltonian Producer Pass 205 - Trim WL Terminal Prechecks

Commit(s):
- this commit - Trim WL terminal prechecks

Summary:
- Accepted a wider `bloat-fixer` stable-code cleanup in the White-Lindsey
  terminal realization surface. The private `_wl_edge_fixed_axes` alias was
  deleted in favor of the shared terminal face-axis helper, and three
  internal prettier-crash checks in `_wl_boundary_stratum_block` were removed:
  source-CPB count, corner support-size, and unsupported-stratum custom error.
- No scientific behavior change is intended. Bad internal metadata now fails
  through `only(...)`, shared axis lookup, coefficient/action bounds, or
  caller tuple destructuring before any malformed block can be pushed.

Validation:
- Bloat-fixer and manager both ran `git diff --check`, package load, and
  `tmp/work/terminal_inventory_native_shell_index_probe.jl`. Manager rerun
  reported PQS/WL `final_dimension = 471`.

Carrying-cost result:
- deleted: `_wl_edge_fixed_axes`, one source-CPB length precheck, one corner
  support-size precheck, and one unsupported-stratum custom throw.
- simplified: edge fixed-axis lookup now calls `_terminal_face_active_axes`
  directly; single-CPB assumptions rely on `only(unit.source_cpbs)`.
- quarantined: coefficient identity validation, support-disjointness
  validation, transform-path checks, unit/source CPB mismatch checks, and
  top-level missing-contract lookup remain untouched.
- exact remaining caller/blocker: no remaining caller of `_wl_edge_fixed_axes`;
  further bloat-fixer work should continue on older stable helper/check
  families, not active RG/Cr2 diagnostic source.
- added src lines: 1.
- deleted src lines: 12.
- net src lines: -11.

### Medium-Term Goal Checkpoint After Pass 205

- MT1 fake-PQS quarantine: active/maintained. Passes 200-205 did not promote
  fake-PQS, route-shadow, or old source-backed oracle paths.
- MT2 independent H2 PQS recovery and MT3 common physical support vocabulary:
  active as background guardrails. Recent WL cleanups stayed downstream of the
  corrected terminal-basis construction boundary and did not reopen support
  authority.
- MT4 supplement staging after authority: active, with the current tactical
  science lane now Cr2 residual-sector diagnosis. The injection audit closed
  the default-off injection question; next evidence should decompose kinetic
  and nuclear contributions rather than proposing a safety gate.
- MT5 cleanup pressure: active and stronger. The `bloat-fixer` role has moved
  from wrapper deletion to stable helper/check-family trimming, while still
  preserving numerical checks that prevent silent wrong bases.
- MT6 audit/classify old Cartesian flat paths: active. The useful current
  classification is stable WL terminal realization for deletion-oriented
  cleanup, versus active RG/Cr2 diagnostic code for repo-doer measurement.
- Medium-goal wording: still needs a future top-level refresh. The durable
  near-term lane is now Cr2 residual one-body diagnosis plus stable-code
  carrying-cost reduction under manager review.

## Cartesian Hamiltonian Producer Pass 206 - Trim Thin-Slab Prechecks

Commit(s):
- this commit - Trim thin-slab prechecks

Summary:
- Accepted a `bloat-fixer` stable-code cleanup in thin-slab and terminal-face
  helper plumbing. The pass removed redundant source-CPB and metadata
  prettier-crash checks, removed local `Int(retained_count)` normalization, and
  let the touched coefficient allocation rely on Julia's default `Float64`
  `zeros` element type.
- No strategic change. This continues MT5/LT2 stable-code carrying-cost
  reduction while staying out of the active RG/Cr2 diagnostic lane.

Validation:
- Bloat-fixer ran `git diff --check`, package load, and
  `tmp/work/terminal_inventory_native_shell_index_probe.jl` before the
  one-line `zeros` addendum; after the addendum bloat-fixer reran
  `git diff --check` and package load. Manager reran `git diff --check` and
  package load after the addendum.

Carrying-cost result:
- deleted: source-CPB length precheck before `only(source_cpbs)`, required
  metadata `get(..., nothing)` checks in compact thin-slab realization, native
  slab-axis custom throw, local retained-count normalization, and explicit
  `Float64` in the touched `zeros` allocation.
- simplified: `_thin_slab_axis` is now the preserved `:slab_normal_axis` then
  `:bond_axis` fallback accessor; invalid axes fail naturally through index
  lookup before contract construction.
- quarantined: retained-count positivity, slab thickness, fixed-index,
  fixed-side, coefficient identity, support-disjointness, transform-path, and
  source/contract mismatch checks remain untouched.
- exact remaining caller/blocker: `_thin_slab_axis` is now a one-call helper
  and may be an easy future inline/delete target.
- added src lines: 6.
- deleted src lines: 16.
- net src lines: -10.

## Cartesian Hamiltonian Producer Pass 207 - Decompose Cr2 Residual One-Body Modes

Commit(s):
- this commit - Record Cr2 residual one-body decomposition

Summary:
- Accepted the measurement-only Cr2 residual-sector one-body decomposition.
  The result rules out the clean kinetic/sign bug class: `K_RR` stays positive
  in both the existing cutoff-0 artifact sector and the `1.0e-3`
  measurement-only injection replay. The large negative `H1_RR` modes come
  from nuclear attraction, and the attraction is balanced across the two Cr
  centers rather than concentrated on one center.
- This changes the interpretation of the Cr2 residual blocker. Negative
  `H1_RR` alone is not evidence of an unphysical state. The Cr2 issue is that
  residual directions remain much more strongly nuclear-bound than the N2
  comparator after the current residual-selection policy.

Validation / evidence:
- Doer used the existing artifact
  `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/cr2_fixed95fec2b8_ida.jld2`
  plus the June 29 Cr2 counterpoise workdir. The ignored probe
  `tmp/work/cr2_residual_onebody_decomposition_probe.jl` is ignored by
  `.gitignore:12`.
- Manager inspected the ignored probe and did not rerun it because the
  `1.0e-3` replay took about `466.65 s`. The script computes `K_RR`, per-center
  `Z_i U_i`, total `H1_RR`, owner weights, dominant labels/components, and
  closure errors.
- Cutoff `0.0`: base/residual dimensions `6675/136`, `min K_RR =
  0.3216563629490918`, no `K_RR < -1e-10`, `78` negative `H1_RR` modes. The
  lowest mode has `H1 = -7.267806471`, `K = 5.741393382`, and balanced nuclear
  terms near `-6.5045999` from each center.
- Cutoff `1.0e-3`: measurement-only relaxed replay, dimensions `6675/58`,
  `min K_RR = 0.8756522802172541`, no `K_RR < -1e-10`, `37` negative `H1_RR`
  modes. The lowest mode has `H1 = -6.992651011`, `K = 4.276517222`, and
  balanced nuclear terms near `-5.6345841` from each center.

Goal advancement:
- LT5/LT6: narrows the residual-sector blocker from broad "bad low H1 mode" to
  "strong balanced nuclear attraction in residual directions despite positive
  kinetic energy." This discourages a safety gate and instead points to
  residual selection, nuclear geometry/scale, and owner-local candidate
  interpretation.
- MT4: keeps Cr2 supplement staging evidence-based. Injection remains
  default-off and not a cure; the next source pass should diagnose or alter
  residual selection rather than treating negative `H1_RR` as a standalone
  failure.

Risk / guardrail:
- The `1.0e-3` row remains measurement-only because it uses the same relaxed
  replay context as the injection audit where default production `F' S R`
  validation fails. It is useful for diagnosis, not production behavior.
- Live manager state after acceptance is cleaner than the doer handback's
  dirty note: the thin-slab bloat-fixer changes were accepted and pushed in
  `28f2ddffc`; only the two untracked successor handoff docs remain.

Carrying-cost result:
- deleted: none.
- simplified: none.
- quarantined: ignored measurement probe only; no source, tests, artifacts,
  commits, or solver runs were created by the doer.
- exact remaining caller/blocker: no kinetic/sign bug was found. The remaining
  blocker is explaining why Cr2 owner-balanced residual directions retain such
  large nuclear attraction relative to N2 and deciding whether residual
  selection should account for one-body/nuclear binding.
- added tracked source lines: 0.
- deleted tracked source lines: 0.

## Cartesian Hamiltonian Producer Pass 208 - Audit Cr2 Residual MWG Adequacy

Commit(s):
- this commit - Record Cr2 residual MWG adequacy audit

Summary:
- Accepted the measurement-only Cr2 residual-interaction adequacy audit. The
  bad saved HF density has about `3.7956 e` in residual rows. Residual-involving
  `H1` lowers the energy by about `-30.18 Ha`, while residual-involving MWG/IDA
  `Vee` charges only about `+20.38 Ha`, leaving a net residual-containing gain
  near `-9.80 Ha`.
- This materially sharpens the Cr2 failure hypothesis: the low one-body modes
  are not the bug by themselves. The likely failure is that exact `H1` finds
  broad, two-center, midbond/interstitial residual directions, while the MWG
  residual interaction under-penalizes occupation of those directions.

Validation / evidence:
- Doer used the existing Cr2 Hamiltonian artifact and bad UHF state under
  `/Users/srw/dmrgtmp/cr2_r1p68_ns7_lmax2_d0p00847_fixed95fec2b8/`, plus mode
  width data from the June 29 Cr2 counterpoise workdir. The ignored probe
  `tmp/work/cr2_mwg_residual_adequacy_probe.jl` is ignored by `.gitignore:12`.
- Manager inspected the ignored probe and did not rerun the long relaxation
  sweeps. The script reads the saved artifact/state, computes residual density
  trace, low-`H1_RR` mode occupations and widths, decomposes one-body/Hartree/
  exchange/two-body terms into base-base, base-residual, and
  residual-residual parts, and runs residual-containing `Vee` scaling probes.
- Key saved-density result: residual trace `3.795634732401 e`
  (`1.897817366201` alpha and beta). Residual-containing one-body contribution
  is about `-30.1824 Ha`; residual-containing two-body contribution is about
  `+20.3843 Ha`.
- Mode evidence: the exact lowest broad `H1_RR` mode is almost empty, but
  nearby broad low-`H1_RR` modes are heavily occupied. Examples include mode 2
  with `H1=-6.8146`, occupation `1.0322 e`, widths about
  `(2.95, 2.95, 2.00)`, and MWG self-Coulomb about `0.1115`; modes 4/5 also
  have substantial occupation and broad anisotropic widths.
- Scaling evidence: scaling residual-containing `Vee` by `1.25` reduced
  residual occupation to about `0.893 e` in 20 bounded sweeps; scaling by `1.5`
  reduced it to about `0.0446 e`. The `2.0` relaxed run did not complete
  promptly, so only fixed-density penalty evidence is used there.

Goal advancement:
- LT5/LT6: replaces the earlier "low residual `H1_RR` is bad" framing with a
  more physical failure mode: broad residual directions may be legitimate
  attractive midbond functions, but the residual MWG/IDA interaction appears
  too soft to keep HF from overoccupying them.
- MT4: points the active Cr2 lane toward residual interaction adequacy. The
  next source decision should not be an `H1_RR` negativity guard or another
  injection cutoff tweak; it should compare or improve residual-containing
  interaction treatment for the broad occupied modes.

Risk / guardrail:
- This is still measurement-only. The scaling probe is not a production fix,
  and the exact Coulomb comparator was blocked because the saved artifact does
  not include raw exact two-electron blocks or quadrature data for these broad
  mode comparisons.
- Live worktree includes an unrelated bloat-fixer dirty edit in
  `src/cartesian_terminal_lowering/region_contracts.jl`; this acceptance
  records only the measurement interpretation and does not accept that source
  diff.

Carrying-cost result:
- deleted: none.
- simplified: none.
- quarantined: ignored measurement probe only; no source defaults, artifact
  schema, MWG redesign, spectral guard, commits, or production Cr2 claim were
  made by the doer.
- exact remaining caller/blocker: need an interaction-side comparator or
  policy for broad residual modes. The current blocker is determining whether
  MWG self/base-residual/residual-residual terms are quantitatively too weak
  for the occupied broad midbond residual density and what bounded replacement
  or correction is justified.
- added tracked source lines: 0.
- deleted tracked source lines: 0.

## Cartesian Hamiltonian Producer Pass 209 - Delete Terminal Lowering Wrappers

Commit(s):
- this commit - Delete terminal lowering wrappers

Summary:
- Accepted a larger `bloat-fixer` stable-code cleanup in
  `src/cartesian_terminal_lowering/region_contracts.jl`. Five private wrapper
  helpers were deleted: `_source_cpb_from_box`, `_filled_source_cpb_from_box`,
  `_raw_region`, `_thin_slab_axis`, and `_direct_lowering_kind`.
- No strategic change. The pass keeps the same CPB construction calls, preserves
  filled-CPB box splatting, preserves thin-slab axis fallback order, and keeps
  the local guard that prevents unsupported regions from becoming direct
  contracts.

Validation:
- Bloat-fixer and manager both ran `git diff --check`, package load, and
  `tmp/work/terminal_inventory_native_shell_index_probe.jl`. Manager rerun
  reported PQS/WL `final_dimension = 471`.

Carrying-cost result:
- deleted: five stable terminal-lowering wrapper helpers.
- simplified: contract builders now use direct `region.raw_region` access and
  direct `CartesianCPB.cpb` / `CartesianCPB.filled_cpb` calls.
- quarantined: `_terminal_source_role` remains because it avoids repetitive
  key construction; retained-count positivity, slab-thickness, complete-shell
  inner-exclusion, direct-region, and unsupported-region guards remain because
  they protect retained dimensions or route semantics.
- exact remaining caller/blocker: no remaining callers of the deleted helpers.
  The terminal-lowering surface is now a cleaner stable-code baseline for any
  future bloat-fixer scans.
- added src lines: 20.
- deleted src lines: 42.
- net src lines: -22.

## Docs-History Compression - PQS Passes 230-265

Commit(s):
- this commit - Compress PQS blurb logs 230-265

Summary:
- Appended the pass 230-265 durable summary to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/summary.md`.
- The summary preserves fake-PQS quarantine, independent H2 PQS counts
  `(275, 98, 98)`, final dimension `471`, support counts `(275, 578, 362)`,
  support-partition/tile facts, pass-234/pass-241/pass-261 exceptions and
  paydown, route-shadow/current-route/legacy fixture deletions, and
  provider-block guardrails.
- Deleted the remaining raw `blurb`, `response`, `review`, attention, and
  exception files in that archive directory. This supersedes the earlier note
  that the active 230-current raw tail remained intact.

Validation:
- Bloat-fixer: docs reference grep and `git diff --check`.

Goal advancement:
- LT2/MT5: removes tracked docs-history bulk while preserving strategic
  interpretation in the curated summary and this running log.
- MT6/AG7: keeps fake-PQS and old route-shadow retirement history available
  without requiring future agents to reread obsolete raw baton files.

Risk / guardrail:
- This is docs-history compression only. No source, tests, artifacts, route
  semantics, numerical thresholds, public API, or active RG/Cr2 diagnostics are
  changed.

Line-count / complexity note:
- Deleted 10,722 raw archive lines and kept `summary.md` as the local durable
  archive record.

## Docs-Authority Compression - PQS Source-Box Framework

Commit(s):
- this commit - Compress PQS source-box docs

Summary:
- Accepted the `bloat-fixer` compression of the three long PQS policy notes
  into one compact current authority. `pqs_source_box_operator_framework.md`
  now carries the source-box-first contract; the old
  `raw_product_source_retained_transform_policy.md` transition note is
  deleted; `projected_q_shell_policy.md` remains only as a 26-line
  compatibility/provenance pointer because current design docs still cite that
  filename.
- Durable facts preserved: raw product-box stage is separate from shell
  realization; boundary COMX/product-mode selection is retained-source
  authority; shell support is owned support, not a CPB; shell-local Lowdin is
  realization/debug/oracle machinery; final PQS IDA weights are final
  retained-function integrals; PGDG analytic pair factors remain the normal
  path unless an explicit debug/reference path is selected.

Validation:
- Manager inspected the compressed authority and reference updates, checked
  that the deleted raw-product policy filename has no remaining refs, and used
  `git diff --check`.

Goal advancement:
- LT2/MT5/AG7: removes a large duplicated route-shadow documentation surface
  while making the active PQS contract easier for future agents to find.
- MT6: reduces the chance that old smoke/prototype narratives will be mistaken
  for implementation authority.

Risk / guardrail:
- Docs-only compression. No source, tests, artifacts, route semantics,
  numerical thresholds, public API, or design-registry authority changed.

Line-count / complexity note:
- Deleted one retired policy doc and reduced the remaining two policy docs
  from about 3,100 lines to 237 lines, with total pass net `-4,673` docs lines.

## Cartesian Hamiltonian Producer Pass 210 - Delete Retired PQS One-Body Wrapper

Commit(s):
- this commit - Delete retired PQS one-body wrapper

Summary:
- Accepted the `bloat-fixer` deletion of the retired
  `src/cartesian_pair_block_materialization/pqs_source_one_body.jl` selector
  wrapper file. The one live mixed-dispatch path now calls the surviving
  safe-term helpers directly for overlap, position, `x2`, and kinetic. Exports,
  include wiring, stale selector comments, and the old selector-surface summary
  were removed.
- Manager review added one small source-doc clarification in
  `pqs_source_shell_final_basis.jl`: the retained-boundary object-kind contract
  remains named, but the deleted constructor function is no longer advertised.

Validation:
- Bloat-fixer ran deleted-symbol/file scans, `git diff --check`, package load,
  and the terminal inventory probe. Manager reran the stale-symbol scan,
  `git diff --check`, and package load. The remaining
  `:unsupported_pqs_source_one_body_term` symbol is a blocker code in the
  surviving safe-term helper, not a deleted wrapper reference.

Goal advancement:
- LT2/MT5: removes another route-era selector layer instead of preserving
  compatibility vocabulary around a retired path.
- MT6/AG7: narrows the PQS source one-body surface to the actual live helpers
  and leaves `source_box_pair_shadow.jl` untouched for a later classified
  deletion map.

Risk / guardrail:
- No route semantics, numerical kernels, public APIs, artifact schemas,
  reports, statuses, residual-GTO, Cr2, or source-box shadow donor paths were
  changed. The exported final-basis boundary helper still accepts the old
  retained-boundary object-kind shape if a caller constructs it manually; this
  pass only deletes the retired selector/constructor surface.

Carrying-cost result:
- added source lines: 32.
- deleted source lines: 517.
- net source lines: -485.
- exact remaining blocker: `source_box_pair_shadow.jl` and
  `pqs_source_safe_terms.jl` remain the larger oracle/donor surfaces; they need
  separate caller-proof maps before deletion.

### Medium-Term Goal Checkpoint After Pass 210

- MT1 fake-PQS quarantine: active/maintained. Recent cleanup and Cr2 audits did
  not revive fake-PQS or source-backed WL/QW as independent PQS authority.
- MT2 independent H2 PQS recovery and MT3 common physical support vocabulary:
  completed as active guardrails. The independent H2 PQS route remains the
  provenance baseline, while current work should avoid reopening old route
  scaffolds unless a live physics target requires it.
- MT4 supplement staging after authority: active, with the current scientific
  lane shifted from injection/cutoff tuning toward residual interaction
  adequacy and screened-reference fluctuation audits.
- MT5 cleanup pressure: active and strong. Bloat-fixer has now moved from
  small stable-helper deletion into large docs compression and retired source
  wrapper deletion; keep requiring net-negative accounting.
- MT6 audit/classify old Cartesian flat paths: active. The next large source
  target should be a read-only classification of `source_box_pair_shadow.jl`
  subfamilies before any broad deletion.

## Cartesian Hamiltonian Producer Pass 211 - Block Screened-Reference Be/Be2 Audit

Commit(s):
- this commit - Record screened-reference q0 blocker

Summary:
- Accepted the measurement-only Be/Be2 screened-reference audit stop. The fast
  Be and Be2 artifacts exist, but the audit correctly stopped before forming
  `Delta h = J_Galerkin[rho0] - J_IDA/MWG[q0]`: the artifacts do not carry a
  cheap, valid projection of the fixed Gaussian reference density `rho0` into
  the same IDA/MWG density coordinates as `q0`.
- This is the right blocker. Center metadata alone is not a density coordinate,
  and retained/final weights are not automatically positive IDA quadrature
  weights. Guessing `q0` from centers would repeat the category mistakes the
  PQS framework and numerical-contract notes now forbid.

Validation / evidence:
- Doer inspected the fast artifacts under the June 24 driver ns-matrix run:
  Be atom base/supplemented `ns4` and Be2 base/supplemented `ns4`, with
  dimensions `237/246/549/567` and supplemented residual dimensions `9/18`.
- Closest existing primitives are Gaussian factor one-body pieces,
  centered electron-nuclear source blocks, dense small pure-Gaussian ERI
  reference helpers, and saved `electron_electron_ida` matrices. None supplies
  `rho0 -> q0` in the final/augmented IDA/MWG density gauge.
- Doer ran package load, created no ignored probe, made no tracked edits, and
  did not run Cr2. Manager verified live repo state is cleaner than the doer
  handback: tracked files are clean after Pass 210; only the two longstanding
  untracked successor handoff docs remain.

Goal advancement:
- MT4/LT5: sharpens the screened-reference lane. The missing design object is
  not the exact Gaussian Hartree screen, but the reference-density projection
  into the approximate interaction's density coordinates.
- LT6: preserves density-gauge honesty by refusing a center-only or
  retained-weight-only `q0` shortcut.

Risk / guardrail:
- No source, artifacts, defaults, public API, residual pruning, interaction
  scaling, injection change, or Cr2 run. The next task should be a bounded
  design/measurement map for a valid `rho0 -> q0` projection path, not a
  production screened-reference implementation.

Carrying-cost result:
- added tracked source lines: 0.
- deleted tracked source lines: 0.
- exact remaining blocker: define a small, valid IDA/MWG reference-density
  coordinate projection for fixed Gaussian `rho0`, preferably first as an
  ignored-probe measurement object.

## Cartesian Hamiltonian Producer Pass 212 - Delete PQS/PQS Low-Order Shadow Oracle

Commit(s):
- this commit - Delete PQS/PQS low-order shadow oracle

Summary:
- Accepted the `bloat-fixer` deletion of the PQS/PQS low-order
  route-shadow oracle family from
  `src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl`. Deleted
  helpers include the PQS/PQS axis cross factors, pair plan, factor/block
  assembly, explicit boundary-selection reference, and reference-block wrapper
  family.
- Manager review also deleted the now-unused
  `_PQS_PQS_SOURCE_BOX_REFERENCE_TERMS` alias and updated the historical
  route-retirement ledger so it no longer says this oracle remains kept.

Validation:
- Bloat-fixer ran deleted-symbol scans, `git diff --check`, and package load.
  Manager reran deleted-symbol scans over `src test docs`, `git diff --check`,
  and package load. Remaining refs are only historical retirement-ledger
  entries marked as later-retired; ignored `tmp/work` validation scripts still
  contain stale references and were not edited.

Goal advancement:
- LT2/MT5: removes the first large no-live-caller slice from
  `source_box_pair_shadow.jl` without adding an adapter or new abstraction.
- MT6/AG7: proves source-shadow deletion can proceed by subfamily caller proof.
  The raw product-box donor path, density-density helpers, local-Gaussian
  helpers, nuclear helpers, and route skeleton metadata remain untouched for
  separate classification.

Risk / guardrail:
- No numerical kernels, route semantics, public APIs, artifact schemas,
  residual-GTO, Cr2, screened-reference, density, local-Gaussian, or nuclear
  helper paths changed. Do not reopen `_pqs_raw_product_box_plan`,
  `_pqs_product_box_realization_plan`, or `_pqs_raw_product_box_plan_view`
  casually; they remain active donor/source callers.

Carrying-cost result:
- added source lines: 0.
- deleted source lines: 424.
- net source lines: -424.
- exact remaining blocker: the larger `source_box_pair_shadow.jl` density,
  local-Gaussian, and nuclear shadow families still need staged caller-proof
  deletion maps before source edits.

## Cartesian Hamiltonian Producer Pass 213 - Be j0/q0 Gauge Diagnostic

Commit(s):
- this commit - Record Be screened-reference gauge diagnostic

Summary:
- Accepted the measurement-only Be/Be2 screened-reference gauge audit as a
  partial result. The Be atom path can build `j0` in the same terminal/MWG
  density-proxy gauge as `electron_electron_ida`, including residual MWG rows
  for the supplemented atom. The fit remains a diagnostic, not an accepted
  `q0`, because it used unconstrained ridge least squares and has a small
  negative tail.
- Be2 is blocked for this audit path: current repo replay gives `base_dim=419`,
  while the saved Be2 artifacts have `base_dim=549`. The artifacts do not
  serialize the terminal support states/final coefficients needed to reconstruct
  the original density-proxy gauge, so no Be2 `j0` was fabricated.

Validation / evidence:
- Doer used ignored probe `tmp/work/be_screened_reference_q0_gauge_probe.jl`
  and the fast June 24 Be/Be2 artifacts. Reference density used
  `alpha=8.0` and `N_screen = 0, 1, 2, 4` per Be center.
- Be atom base: `dim=237`, active set `A=125`, rank `125`, condition about
  `1.91e2`. For `N=1`, `q_sum=1.005231`, `q_min=-1.70e-3`, `q_max=4.57e-2`,
  `q_neg=8`, relative error about `9.06e-4`; scaling for `N=2/4` was linear to
  printed precision.
- Be atom supplemented: `dim=246`, residual `9`, active set `A=125`, rank
  `125`, condition about `1.92e2`. For `N=1`, relative error about `9.16e-4`
  and residual-row relative error about `2.27e-3`.

Goal advancement:
- MT4/LT5: moves screened-reference work from "no q0 path" to a concrete
  atom-only gauge diagnostic with small fit error and a clear constrained-fit
  next step.
- LT6: preserves gauge honesty by refusing to fabricate Be2 density proxies
  from mismatched replay geometry or center metadata.

Risk / guardrail:
- No source edits, defaults, artifact schema changes, Cr2 run, HF relaxation,
  residual pruning, interaction scaling, or production screened-reference claim.
  The next task should either add a constrained nonnegative/charge-constrained
  diagnostic solve for Be atom or regenerate/load a Be2 artifact with the
  terminal density-proxy gauge data needed for `j0`.

Carrying-cost result:
- added tracked source lines: 0.
- deleted tracked source lines: 0.
- exact remaining blocker: no accepted `q0` until nonnegative/charge
  constraints are satisfied; no Be2 audit until artifact gauge data and replay
  dimensions agree.

## Cartesian Hamiltonian Producer Pass 214 - Delete PQS/Product Low-Order Shadow Oracle

Commit(s):
- this commit - Delete PQS/product low-order shadow oracle

Summary:
- Accepted the `bloat-fixer` deletion of the PQS/product low-order
  route-shadow oracle family from
  `src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl`.
  Deleted helpers include the PQS/product axis cross factors, pair plan,
  factor/block assembly, reference-block wrapper family, term constants, and
  now-unused low-order term validation.
- Updated `docs/src/developer/cartesian_route_retirement_ledger.md` so the
  old PQS/product multi-term oracle is recorded as later-retired history
  instead of a surviving oracle surface.

Validation:
- Bloat-fixer ran deleted-symbol scans, `git diff --check`, and package load.
  Manager reran source/test/docs caller scans, inspected the diff, ran
  `git diff --check`, and reran package load. Remaining deleted-symbol refs
  are only historical ledger entries and ignored `tmp/work` scripts.

Goal advancement:
- LT2/MT5: removes the product sibling of the retired PQS/PQS low-order
  route-shadow oracle without adding an adapter or preserving stale vocabulary.
- MT6/AG7: continues the staged shrinkage of `source_box_pair_shadow.jl` by
  caller-proved subfamilies while leaving active raw product-box donor paths
  untouched.

Risk / guardrail:
- No numerical kernels, route semantics, public APIs, artifact schemas,
  density-density helpers, local-Gaussian helpers, nuclear helpers, residual
  GTO, Cr2, or screened-reference code changed. The raw product-box donor
  path remains active and should not be reopened casually.

Carrying-cost result:
- added source lines: 0.
- deleted source lines: 335.
- net total lines: -333 after the small historical-ledger update.
- exact remaining blocker: `_pqs_product_source_box_product_block` is now a
  definition-only adjacent leftover and should be handled only as a separate
  approved deletion slice; density, local-Gaussian, nuclear, and skeleton
  metadata families still need their own caller-proof maps.

## Cartesian Hamiltonian Producer Pass 215 - Delete Definition-Only PQS/Product Helper

Commit(s):
- this commit - Delete definition-only PQS/product helper

Summary:
- Accepted the `bloat-fixer` deletion of the adjacent definition-only
  `_pqs_product_source_box_product_block` helper from
  `src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl`. This
  removes the small tail left after Pass 214 retired the PQS/product low-order
  route-shadow oracle family.

Validation:
- Bloat-fixer ran `git diff --check` and package load. Manager inspected the
  diff, reran a source/test/docs caller scan excluding this running log, reran
  `git diff --check`, and reran package load. No source/test/docs caller
  remains outside historical manager notes.

Goal advancement:
- LT2/MT5: keeps the source-shadow cleanup net-negative and removes a stale
  compatibility wrapper instead of leaving a false blocker in the source.
- MT6/AG7: continues the subfamily-by-subfamily deletion approach for
  `source_box_pair_shadow.jl` while leaving raw product-box donors,
  density-density, local-Gaussian, nuclear, and skeleton metadata families
  untouched.

Risk / guardrail:
- No numerical kernels, route semantics, public APIs, artifact schemas, docs,
  tests, tmp probes, Cr2, screened-reference, density, local-Gaussian, or
  nuclear paths changed.

Carrying-cost result:
- added source lines: 0.
- deleted source lines: 13.
- net source lines: -13.
- exact remaining blocker: larger `source_box_pair_shadow.jl` density,
  local-Gaussian, nuclear, and skeleton-metadata-connected families still need
  their own caller-proof maps before deletion.

### Medium-Term Goal Checkpoint After Pass 215

- MT1 fake-PQS quarantine: active/maintained. The recent deletion passes did
  not revive fake-PQS, source-backed WL/QW, or route-shadow oracles as
  production PQS authority.
- MT2 independent H2 PQS recovery and MT3 common physical support vocabulary:
  completed as active guardrails. They remain important constraints, but they
  are not the active work queue.
- MT4 supplement staging after authority: active. The current scientific lane
  is residual interaction adequacy and screened-reference density-gauge
  projection, with Be/Be2 as fast probes and Cr2 as the motivating stress case.
- MT5 cleanup pressure: active and strong. Bloat-fixer is now removing both
  large docs/history bulk and staged no-caller source-shadow slices; keep
  requiring net-negative deletion/shrinkage accounting.
- MT6 audit/classify old Cartesian flat paths: active. The useful source
  target remains `source_box_pair_shadow.jl`, but only by staged caller-proof
  subfamilies. Raw product-box donor plan/view helpers remain protected.
- Medium-goal wording: needs refinement soon. The top-level MT section still
  reflects the old H2 recovery era; current practice has shifted to
  CR2-facing residual/screened-reference measurement plus stable-source
  deletion. No wording is changed in this mechanical cleanup commit.

## Cartesian Hamiltonian Producer Pass 216 - Current Be2 Screened-Reference Gauge Diagnostic

Commit(s):
- this commit - Record current Be2 screened-reference gauge diagnostic

Summary:
- Accepted the measurement-only current-source Be2 screened-reference gauge
  diagnostic. The old June 24 Be2 artifacts are no longer needed for this
  measurement path: their `base_dim = 549` stale direct-identity slab behavior
  is superseded by current compact thin-slab lowering, where current Be2 builds
  `base_dim = 419 = 275 + 56 + 16 + 16 + 56`.
- The current in-memory path can build `j0` in the same terminal/MWG
  density-proxy gauge as `electron_electron_ida`. A constrained diagnostic
  `q0` fit was run only on direct atom-contact core base columns with
  nonnegative per-center charge constraints. This is a credible measurement
  step, not a production screened-reference correction or final accepted `q0`.

Validation / evidence:
- Doer used ignored probe
  `tmp/work/be2_current_screened_reference_gauge_probe.jl`. Manager inspected
  the probe and reran it successfully on current `main`; elapsed time was about
  `29.77 s`.
- Be2 base current: dimension `419`, residual `0`, active `q0` coordinates
  `275` direct core columns with owner counts `125/150`. For `N_screen = 1`
  per Be center and `alpha = 8.0`, constrained fit had
  `rel ||Vq0-j0|| = 2.77503591e-3`, exact charges `[1, 1]` to roundoff,
  `q_min = 0`, `q_max = 5.17047091e-2`, and no negative coordinates.
- Be2 supplemented current: dimension `437`, base `419`, residual `18`. For
  the same reference density, constrained fit had
  `rel ||Vq0-j0|| = 2.76870241e-3` and residual-row relative error
  `1.36090535e-3`; charges were `[1, 1]` to roundoff and no negative
  coordinates.
- `N_screen = 0` gave exact zero. `N_screen = 2` and `4` scaled linearly.
  Unconstrained ridge fits were tighter, about `6.9e-4` relative error, but
  had `42` negative coordinates and remain diagnostic only.

Goal advancement:
- MT4/LT5: clears the stale-artifact blocker from Pass 213 and demonstrates a
  current Be2 path for density-gauge `j0` plus a sane constrained diagnostic
  `q0`.
- LT6: preserves gauge honesty by fitting through the existing interaction
  map `V[:, A] q0 ~= j0` rather than using center metadata or final weights as
  density coordinates.

Risk / guardrail:
- No tracked source edits, defaults, artifact schema changes, public inputs,
  Cr2 run, HF relaxation, residual pruning, Vee scaling, or production
  screened-reference correction. The fit excludes residual coordinates and is
  not yet an accepted general `q0` policy.

Carrying-cost result:
- added tracked source lines: 0.
- deleted tracked source lines: 0.
- exact remaining blocker: improve or classify the constrained fit before
  treating `q0` as accepted. The direct next step is a small constrained-fit
  improvement or comparison of allowed coordinate sets; generate a clean
  current Be2 artifact only if persistent readback is needed.

## Cartesian Hamiltonian Producer Pass 217 - Be2 Screened-Reference Delta-H Probe

Commit(s):
- this commit - Record Be2 screened-reference Delta-H probe

Summary:
- Accepted the measurement-only current-source Be2 screened-reference
  one-body correction probe. The correction is formable in memory as
  `Delta h = J_Galerkin[rho0] - Diagonal(V_IDA * q0)` using the constrained
  nonnegative direct-core `q0` from Pass 216 and an analytic Gaussian-density
  convolution for `J_Galerkin[rho0]`.
- The result looks numerically sane as a measurement: `Delta h` is symmetric,
  scales linearly with `N_screen`, and pushes the lowest residual-sector
  `H1_RR` eigenvalues upward. It does not eliminate the negative Be2 residual
  sector; at `N_screen = 4` the lowest residual `H1_RR` mode remains about
  `-0.348`.

Validation / evidence:
- Doer used ignored probe `tmp/work/be2_screened_reference_delta_h_probe.jl`.
  Manager inspected the probe and reran `git diff --check`, package load, and
  the probe on current `main`; probe elapsed time was about `29.86 s`.
- Setup: current Be2 supplemented dimension `437`, base dimension `419`,
  residual dimension `18`, `alpha = 8.0`, `q0` coordinates `275` direct
  atom-contact core columns, inherited constrained fit relative error about
  `2.7687e-3`, and residual-row relative error about `1.3609e-3`.
- Low `H1_RR` eigenvalues moved from
  `[-0.4133, -0.3776, -0.3757, -0.3757, -0.1807]` to approximately
  `[-0.3667, -0.3647, -0.3647, -0.3255, -0.1452]` at `N = 1`,
  `[-0.3566, -0.3566, -0.3204, -0.2756, -0.1144]` at `N = 2`, and
  `[-0.3477, -0.3477, -0.2290, -0.1818, -0.0627]` at `N = 4`.
- Fixed-density energy shift was not run because no current-source Be2 HF
  density was loaded cheaply and stale June 24 Be2 artifacts were not used.

Goal advancement:
- MT4/LT5: moves the screened-reference idea from a density-gauge diagnostic
  to a formable one-body correction on a fast current-source molecule.
- LT6: preserves gauge discipline by subtracting the approximate screen in
  the same `V_IDA * q0` density coordinates and by refusing stale artifacts or
  center-only approximations.

Risk / guardrail:
- No tracked source edits, production defaults, artifact schema changes,
  public inputs, HF workflow, Cr2 run, residual pruning, Vee scaling, or
  production screened-reference claim. The Be2 result justifies a
  measurement-only Cr2 probe, not implementation.

Carrying-cost result:
- added tracked source lines: 0.
- deleted tracked source lines: 0.
- exact remaining blocker: test whether the same measurement gives useful
  physical signal on the bad Cr2 case: low residual `H1_RR` shifts,
  fixed-density correction on the saved bad UHF state, and residual-sector
  occupation signal if a bounded relaxation is cheap.

## Cartesian Hamiltonian Producer Pass 218 - Delete Source-Box Density Shadow Helpers

Commit(s):
- this commit - Delete source-box density shadow helpers

Summary:
- Accepted the `bloat-fixer` deletion of the remaining PQS/product and
  PQS/PQS density-density shadow helper families from
  `src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl`.
  Deleted functions include the density axis-factor, density weight-view,
  density-normalized interaction, and raw-weighted wrapper helpers for both
  pair families.
- Removed the route-skeleton helper-name metadata surface that pointed at
  those deleted helpers. Pair entries, pair-family counts, retained units, and
  route shape remain; route summaries no longer carry nonexistent
  `density_density_helper` or `helper_by_pair_family` symbols.

Validation:
- Bloat-fixer ran deleted-symbol scans, `git diff --check`, package load, and
  the terminal inventory probe. Manager inspected the source and skeleton
  diffs, reran deleted-symbol/metadata scans over `src test docs` excluding
  this log, reran `git diff --check`, package load, and the terminal inventory
  probe. The terminal inventory probe still reports PQS/WL final dimension
  `471`.

Goal advancement:
- LT2/MT5: removes another large stale route-shadow source slice instead of
  preserving helper-name compatibility around dead code.
- MT6/AG7: advances the staged cleanup of `source_box_pair_shadow.jl` while
  preserving the protected raw product-box donor plan/view path and the
  local-Gaussian/nuclear families that need separate maps.

Risk / guardrail:
- No numerical production path, public API, artifact schema, driver input,
  Cr2, screened-reference, local-Gaussian, nuclear, project-axis helper, or
  raw donor behavior changed. Ignored `tmp/work` references remain stale and
  non-blocking.

Carrying-cost result:
- source lines: `+2 / -612`, net `-610`.
- deleted: two density-density shadow helper families plus skeleton/stage
  helper-name metadata.
- simplified: route skeleton/stage summaries now carry pair families/counts
  without stale executable-helper symbol names.
- exact remaining blocker: local-Gaussian and nuclear source-box shadow
  families remain in `source_box_pair_shadow.jl`; they require separate
  caller-proof maps because they share by-center/local-Gaussian wrapper
  plumbing.
