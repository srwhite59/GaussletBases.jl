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
