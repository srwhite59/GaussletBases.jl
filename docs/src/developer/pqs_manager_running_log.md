# GaussletBases Cartesian/PQS Manager Running Log

This is the live manager decision ledger. It records current strategic state
and the most recent accepted passes; it does not replace doer reports,
canonical subsystem contracts, `authority.toml`, or `current.md`.

Read this live file before drafting a Cartesian/PQS blurb, accepting a pass, or
resuming manager work after compaction. Historical archives are task-gated
archaeology and are not normal startup reading.

## Archive Index

- [Initiation through Pass 379](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_through_pass_379.md)
  preserves the first `28,549` lines of the pre-rotation ledger verbatim
  (`SHA-256 7c8c72261786da0e09a3fc60bac3ea16b03b41ea4ac71b26a12a48e91f71af85`).
- [Passes 380 through 406](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_380_through_406.md)
  preserves the next `1,052` accepted ledger lines verbatim
  (`SHA-256 14ccf05eb960757cb3335351658cf54b707c6003e9bbb40916332b398e9a0767`).
- [Passes 407 through 429](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_407_through_429.md)
  preserves the next `927` accepted ledger lines verbatim
  (`SHA-256 83075d2d43762538fcad502a075cae85f03024f01a6a08e911b4b3f03ad42d12`).
- [Passes 430 through 450](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_430_through_450.md)
  preserves the next `972` accepted ledger lines verbatim
  (`SHA-256 1ae541c29af607853f637200a70fd1ba53938a1394ec8aac748def842c74bba3`).
- [Passes 451 through 474](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_451_through_474.md)
  preserves the next `952` accepted ledger lines verbatim
  (`SHA-256 e07d17fd739f8519511e79dca4ee994a7bf2a8bc61c7703ce147fca18d353169`).
- [Passes 475 through 495](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_475_through_495.md)
  preserves the next `926` accepted ledger lines verbatim
  (`SHA-256 a968d79f768462336d941309b6c0fae3262b14a3467c63a219290bf0e60beb7f`).
- [Passes 496 through 516](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_496_through_516.md)
  preserves the next `881` accepted ledger lines verbatim
  (`SHA-256 8ba3de7a83cc2e24191978ebd773e78e5763148df7427d604dbe85d16e4eba21`).
- [Passes 517 through 537](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_517_through_537.md)
  preserves the next `926` accepted ledger lines verbatim
  (`SHA-256 c0def005ed4292d5cbbc684d70298f25d9a9f03b3ea777bd2001cbfec6dc517a`).
- [Passes 538 through 566](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_538_through_566.md)
  preserves the next `1,082` accepted ledger lines verbatim
  (`SHA-256 fa5207b905c69acb7094599b5d0ec67c0303c5d8e5c3b341e0082023048edd83`).
- [Passes 567 through 591](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_567_through_591.md)
  preserves 753 accepted ledger lines verbatim
  (SHA-256 a54ecb97e95c369bec83d5d2c1a8fa6d993b4f67f2ac3b2f039d98f77936c08a).
- [Passes 592 through 624](designs/cartesian_hamiltonian_producer/history/manager_log/pqs_manager_running_log_passes_592_through_624.md)
  preserves 812 accepted ledger lines byte-for-byte
  (SHA-256 c8690eac982b3a2196a1eed15db23ba034b35fa9b1a6aa27869ba59e271f42b7).
- This live volume begins with Pass 625. Pass entries are preserved in accepted
  order; duplicate or nonmonotonic historical pass numbers are not rewritten.

## Current Strategic State

- The broad producer-documentation reorganization is complete. Schema-v3
  `authority.toml` is authoritative, generated registry/execution-whitelist
  views are checked one-way outputs, and authority CI is fail-closed.
- The historical first post-cutover static conformance audit covered `150`
  execution records: `107` matched, `11` documented gaps, `8` numerical gates,
  and `24` discrepancies. Pass 399 closed the atomic-packet fail-fast subset.
- Current screened-Hartree source and contracts use determinant orbitals for
  `P0/q0`, the density fit for `E0`, and the fitted potential as an approximate
  `J0` evaluator with reported consistency error. A determinant-exact `J0/E0`
  convention requires a separate scientific amendment.
- REQ-084 stopped correctly before molecular-full Cr2 interpretation because
  no repo owner constructs the complete Hartree field of a represented density
  containing terminal and supplement components. The direct internal producer
  design is deferred; Pass 634 grants only constructor validation repair.
  Pass 459 implementation preflight established that the
  missing operation is the contracted source pair-product action itself, not
  the existing target evaluators; source and full-space certification remain
  pending.
- The source-backed Cr2 composition/replay migration is closed. The matched
  full-parent/PQS/White-Lindsey H2+ and fixed-state H2 paper mechanics now use
  common aspect-aware shared-shell dimensions. Parent and PQS rows are
  unchanged; rebuilt WL rows await external same-density oracle
  interpretation before a method-accuracy claim.
- PRF mechanics retain semantic region discovery, consumer-targeted residual
  construction, additive/fixed-span composition, and the category-owned
  unscreened Hamiltonian as private diagnostic/provenance surfaces. Historical
  de-promotion removed six root exports; physical target selection,
  transition-density exchange, and PRF-to-GTO-residual interactions remain
  consumer questions. Released compatibility bindings restored in v0.2.1
  remain preservation-only, not renewed PRF promotion.
- Explicit charged electron sectors are implemented without changing bare
  basis or operator arrays. Specialized retained-GTO EGOI remains archived
  and deferred with no execution grant.
- The expert sliced hydrogen-chain producer is accepted with compact finite and
  periodic-template storage, analytic Galerkin H1, longitudinal
  IntegralDiagonal Vee, and checked long-range continuations. Solver and MPO
  policy remain consumer-owned.
- Immutable v0.2.0, RC1/RC2 and v0.2.1 transactions are closed. Final v0.2.1
  targets 2ad8441efe83; stable is pinned to its verified versioned docs. Preserve
  release-0.2.0 and all original folders. Registration, citation and later
  releases remain separate; normal main deployments update dev, not stable.
- A+B multipole arithmetic/full-span interaction repair is accepted in Pass 653;
  changed fixed-radial identities require owner review, not silent receipt edits.
  Pass 654 accepts C normalized harmonics through l=256 and one actual high-order
  table; quadrature D and consumer rebuilding remain separately unqualified.
- The source-layout moves, Lanczos consolidation, mapped-representation
  relocation, scheduled occupied-first coverage, and HFDMRG test-path repair
  are accepted. Generated execution authority is outside AGENTS.md. Public
  mechanical checks run in Supported floor; structural docs checks stay docs-only.
  Prose wording locks were removed, not rescheduled.
- September 4 review repairs have closed schema/onboarding, all 127 verified
  documented-reference gaps, quick onboarding examples, and algorithm pointers.
  The residual screening/defaults, duplicate-page, and example-guide packet is
  accepted at `719e19430`; the September 4 review is addressed. Absence from
  routine CI is not evidence of deadness.
- Production defaults, public workflows, corrected artifacts, and Cr2 endpoint
  claims remain unchanged unless separately authorized.

## Durable Long-Term Goals

1. Build reliable, provenance-complete Cartesian Hamiltonians and reusable
   consumer artifacts.
2. Keep route identity and scientific conventions explicit; never promote an
   oracle, compatibility path, or diagnostic fixture as production physics.
3. Prefer source-first, factorized, and one-dimensional construction over dense
   production paths while retaining bounded numerical oracles.
4. Maintain stable representation-transfer, residual-Gaussian, interaction,
   reference-density, and artifact boundaries for downstream consumers.
5. Reduce carrying cost and conceptual drift by deleting stale paths and tests
   when their live contract ends.
6. Use stratified validation: small contract tests during implementation,
   bounded physical endpoints at acceptance, and explicit terminal due
   diligence for every interpreted numerical result.

## Current Medium-Term Goals

**MT1 - Conformance remediation (active).** Resolve the bounded Pass 398
discrepancies and demonstrated exported-surface defects under explicit
authority. The ordinary-QW nested diatomic repair, PRF API de-promotion, and
four-name package/internal export repair are closed. The separately validated
centered/displaced arithmetic and radial exit/notification repairs are complete;
near-origin inverse-radius quadrature and finite-expansion limits remain.
Keep any new repair narrow and evidence-led.

**MT2 - Controlled Cr2 source migration (completed).** The source-backed
fixed-state and bounded replay reproduced the former consumer-local path. Any
new Cr2 endpoint, contraction, exchange, or solver interpretation is a
separate scientific choice.

**MT3 - Deferred proposals and blocked producer work.** Standard60 is an
unimplemented proposal deferred pending a named consumer and demonstrated
accuracy/cost benefit. Its exact identity and recovered evidence remain;
kernel arithmetic is repaired, but diffuse/long-range limits survive.
Represented molecular Hartree remains scaling-blocked on the nonmaterializing
contracted pair-product action and full-field certification. Its bounded
constructor-only grant is not evaluator authority. Fitting and Cr2 acceptance
remain separate. Completed H2+ controls, matched shells, charged sectors, and sliced
chain facilities are maintenance, not pending implementations. External
same-density oracle interpretation of corrected WL fingerprints and sliced-chain
HFDMRG adaptation remain consumer work. Specialized retained-GTO EGOI remains
archived/deferred with no execution grant.

**MT4 - Residual and protected-basis evidence (active).** Keep the residual
spectral audit measurement-only. Protected atoms, counterpoise, and any new
injection/localization policy remain separate future decisions. Parent residual
function mechanics, the onsite-calibrated Gaussian direct resource, and the
private diagnostic/provenance wrappers are implemented and source-backed.
Hooke owns the first Be `1s/2s` target study. Selection, transition-density
exchange, and PRF-to-GTO-residual interactions remain consumer or measurement
questions.

**MT5 - Documentation and authority maintenance (active).** v0.2.0/RC and
v0.2.1 lifecycles are closed; stable is pinned to verified v0.2.1 documentation.
Preserve all immutable identities/folders; normal main deployment updates dev.
Tag verification is repaired; RC-era counts/absent-stable requirements are
historical. Registration, citation and later releases remain separate.
Path-aware CI, the generated whole-file whitelist,
Supported-floor mechanical documentation checks, and prose-test cleanup are
implemented/maintenance. Full angular research remains outside per-push CI.
Separate PQS and screening surfaces and current group selections stay fixed.

**MT6 - Carrying-cost control (active).** Exact-order performance improvements,
shared Example 41 release execution, classified public documentation, source
layout, duplicate Lanczos removal, mapped-representation relocation, inert
sidecar retirement, and prose-test reduction are accepted. Cold reporting
barriers were regressions and remain closed without implementation authority.
The documented-reference gap is zero. Five undocumented exports remain an
exact reserved set for separately authorized next-minor namespace decisions;
the associated documented QWRG diagnostic type stays paired with its function.
No further source, compatibility, or numerical change follows from this review.
The September 4 reader review is addressed; further work needs a new target.

**MT7 - External Cartesian GTO interchange (completed/maintenance).** The
strict versioned reader, checkpoint-only PySCF exporter, frozen d-shell
fixture, and explicit determinant cleanup are implemented. The read-only C2
replay reproduced the accepted occupied subspace without a permutation ledger.
No solver, Hamiltonian payload, mandatory PySCF dependency, basis-only/live-mf
export, or release action belongs to this goal. Reader-facing documentation is
implemented and included in immutable RC2 and final v0.2.0.

**MT8 - Finite collinear PQS (construction and atomic-fit connection maintenance).** Pass 630 accepts small-chain
bases, complete operators and raw transfer with general ordered z/positive
charges, subject to mapping validity. Compact outer retention is explicitly
converged, not a long-chain or chemical-accuracy claim. Preserve MT7, all released
interfaces and the distinct sliced model; no legacy retirement is coupled.

Pass 633 closes the residual-GTO/MWG connection to maintenance. hchain-doer's
implementation prerequisite is cleared; consumer scope and supplemented
convergence remain separate. No bare scientific H-chain or H10/H20 campaign
is authorized by this repository closeout.

The occupation-one atomic-fit connection is accepted in Pass 636. Scientific-q,
cutoff, early-z and screening-memory repairs are accepted; Pass 652 closes stable
direct residual construction. None certifies H10 energy/virtual-space accuracy
or restarts hchain-doer. Molecular-field optimization is not a prerequisite.

## Manager Guardrails

- `authority.toml` grants execution; this ledger records interpretation only.
- Preserve exact basis, density, interaction, ordering, and energy-accounting
  conventions across producer and consumer boundaries.
- Fail closed on authority disagreement, malformed or unconverged scientific
  inputs, materially invalid metrics, nonfinite operators, and unsupported
  artifact conventions.
- Do not combine conformance repair with a scientific-policy change.
- Do not infer public, artifact, solver, or production authority from an
  implemented internal helper or successful ignored probe.
- Inspect terminal due diligence before accepting any endpoint interpretation.
- Work with concurrent changes; never absorb or revert unrelated WIP.

## Entry And Rotation Policy

- Append one compact strategic entry after each accepted substantive manager
  pass. Record commits, interpretation, validation actually run, goal movement,
  guardrails, remaining blockers, and carrying-cost impact where relevant.
- Every five accepted passes, append a medium-term checkpoint. Every 10-20
  passes or after a major correction, add a strategic compression entry.
- Do not duplicate the full doer handback or numerical tables already preserved
  in a report.
- Rotate the live volume after 25 additional accepted passes or before it
  exceeds roughly `2,000` lines. Move old entries verbatim to a task-gated
  archive, retain at least the latest 20 passes live, refresh this strategic
  preamble, and record the archive line count and SHA-256.
- `docs/check_manager_log.jl` enforces the `2,000`-line ceiling before every
  Documenter build; exceeding it is a CI failure, not an advisory warning.
- Never reorder, renumber, silently summarize, or delete accepted historical
  entries during rotation.

### Pass 625: Bounded v0.2.1 Candidate Authority

- Reviewed repo-manager's conversation-delivered readiness audit against
  57e61403262042923605f94a28e910dd2b543550 and released adfcaba32d.
  Project declarations and public examples are unchanged. Independently
  confirmed the seven-line export/type removal and the tag-namespace collision.
- Authorized HP-PQS-PUBLIC-V021-FN-01/TEST-01 for one candidate push:
  compatibility restoration, namespaced tag verification, honest unreleased
  version/changelog documentation and focused checks. No additional RC absent
  a concrete blocker. Superseded only the three-name absence maintenance
  condition; historical cleanup evidence remains intact.
- LT public-package reliability advances; MT7 remains completed/maintenance.
  Reuse accepted centered/displaced/radial evidence; the final candidate still
  requires one full matrix and isolated install/example. No numerical runs in
  this authority pass. Source restoration is seven lines, source documentation
  at most 36; bounded test/CI/documentation budgets live in the canonical heading.
- Deleted: obsolete absence requirement. Simplified: one candidate transaction.
  Quarantined: none. Not deleted because: released bindings need compatibility.
  Remaining blocker: implementation and independent candidate acceptance.
  Added/deleted src: 0/0; new tests/files/metadata: none in this pass.
- Validation: authority/self-test, generated views, package load, docs,
  Documenter and diff checks required before repo-manager proceeds. Tagging,
  publication and stable promotion remain separate; release-0.2.0 stays pinned.

### Pass 626: Candidate Context And Metadata Correction Authority

Reviewed 3d26994f8dc5b7d5116ea2bc113a3e4fdf58f7d7: exact 354-export
compatibility, restoration, archive hash and remote CI/Docs evidence agree.
Actual production tag-label evaluation still says unreleased; final acceptance
is withheld. The existing V021 pair is narrowed to five files for context-aware
labels and coherent changelog, landing-page and conditional installation text.
This supersedes initial candidate-only metadata instructions, not accepted
source or tag-verifier behavior. LT public reliability advances; MT7 stays
maintenance, with no new numerical or release goal. Root README/CHANGELOG
changes require a full matrix under the unchanged classifier, despite the
documentation-only substance; this authority commit uses docs-only checks.

Deleted/simplified: replace phrase-presence checking with actual label-context
evaluation; no parallel implementation. Quarantined: none. Retain released
bindings and prior numerical evidence, contingent on explicit unchanged
source/test comparison. Remaining blocker: corrected commit/tree/archive and
rendered tag-context verification, then independent candidate review.
Added/deleted src: 0/0; new tests/files/metadata: none in this authority pass.
Authority/self-test, generated views, docs/package/build/log/diff and remote
checks are required. No publication, tag or stable promotion is granted.

### Pass 627: Candidate Closeout And Conditional v0.2.1 Release

The user explicitly authorized this combined transition, overriding the usual
separate release-authority review boundary for this transaction only. Accepted
candidate 2ad8441efe8328b4f325bf42994925d8f2a491c9 has tree
4343ba110ada71fe7085fb9d2c504886a2e7a1e9; its 690-entry, 10,854,400-byte
archive hash is 8db862e58e083fca604b2e10ec297fede7ea69422f5543cd4f87060ff6eae198.
Independent archive/tree/exclusion, protected implementation and rendered-tag
review passed; docs 8/8 +162/162 +10/10 and exact-head CI 34429827030/Docs
34429827052 passed. Candidate FN/TEST grants are exhausted, not reopened.

The new temporary V021-RELEASE pair owns the exact annotated tag, frozen
1,677-byte body/final/latest/no-assets publication, then the frozen six-line
stable/label/selector substitution patch. Scratch substitution and actual
Documenter/label checks passed 11/11. Stable changes only after publication
and live versioned-doc verification; no copy through a symlink. One final
docs-only closeout removes the temporary grants and verifies pin persistence.
LT public-package reliability advances; MT7 remains maintenance. No new RC or
scientific goal: Standard60 and existing radial/molecular limitations remain.

Deleted: temporary candidate execution surfaces. Simplified: one ordered grant
without another approval cycle. Quarantined: none. Retained: all old release
objects and documentation snapshots. Exact blocker: successful ordered execution;
any failure preserves completed objects and stops. Added/deleted src: 0/0;
new tests/files/metadata: none in authority. Validation is authority/self-test,
generated views, package/docs/build/log/diff and docs-only remote checks.

### Pass 628: v0.2.1 Release And Stable Closeout

Accepted publication 386324250 (final/latest, zero assets) at candidate
2ad8441efe8328b4f325bf42994925d8f2a491c9, tree
4343ba110ada71fe7085fb9d2c504886a2e7a1e9; annotated object
c95f8868afc7a0b703f2fd97bb7b5c4af17416b8. The frozen 1,677-byte body
retains SHA-256 e13082c0bcff9e89c11157a94de31d1898ee02a50f8ed814f324bd691ae00c10.
Automatic tar/zip reconstruct the exact tree: hashes
afa1eb012dac8adb5a872296910b421c6a1a5ce16ce122dc5a0715d6b4a4294d and
fe03ec8fe3346d4b62b02f875c85d545774b52eedac314925f5eefeb11da82be.
Remote installation/load and Example 01 passed. Tag CI/Docs
34434227579/34434227566 passed. Promotion
1d440f8f78c052b41e86a93651d4268d5080388d matches the frozen +6/-6 patch;
CI/Docs 34485494842/34485494779 passed without numerical reruns.

Reviewed completed HTTP evidence: all 14 README destinations, canonical/footer
checks, reference source pins, stable/versioned byte equality, root, selectors,
dev and unchanged protected trees. Sandbox DNS and scratch-parser stops were
resolved by explicit read-only resumptions, never by changing published objects.
Evidence SHA-256: 70139f388a131790729230b3e06abb12b0d0fa820adbee89b96fe160fcfb4c2f.
LT public reliability advances; MT7 stays maintenance. Deleted: temporary
execution grants. Simplified: stable points directly to verified v0.2.1.
Quarantined: none; immutable release evidence is retained. Added/deleted source:
0/0; new tests/files/metadata: none. Remaining gate is this docs-only closeout's
normal-main deployment and stable-pin persistence check. Existing scientific
limitations remain; no successor task is opened.
Local closeout checks passed: authority/self-test and generated parity,
package load, docs 8/8 +162/162 +10/10, Documenter, log bound and diff checks.

### Pass 629: Finite Collinear PQS Implementation Authority

Independent review of both September 14 design/completion reports and their
passing scripts/logs supports one usable feature, not another general audit.
Existing all-nucleus angular calibration, shell/slab realization, GTO overlap,
kinetic, nuclear and IDA kernels suffice. General positive charges are retained
with unchanged mapping rejection. Complete enlarged H3 outputs (1265 squared)
were demonstrated at 30.33 MB retained payload and 2.36/0.485 s first/warm
operator build; these are fixture evidence, not end-to-end scaling promises.
No broad probes were repeated for authorization. Independent integer-set
inspection confirms the missing three-group termination fixture: 1029 group,
98 interior-gap and 98 exterior sites, totaling 1225 exactly once.

The new pair freezes two expert operations, a two-field private bare basis,
accumulated nuclear attraction and a three-field numerical return. User's
450-added-source-line exception is feature-specific; preferred 417 includes
relocated decomposition/docstrings. Existing public Cartesian tests have
170/220 preferred/hard additions; two existing reader pages have 40/55.
Separate px/py capture and projected H/IDA bounds retain the measured slab
asymmetry; they do not certify energy accuracy. Source-bearing implementation
must pass all three existing CI jobs; this authority uses docs-only checks.

MT8 opens for finite-chain construction; MT7 remains maintenance. Deleted:
duplicated outer decomposition must be removed when made callable. Simplified:
no odd/even recursion, dense global coefficient map or per-center dense list.
Quarantined: dense Q/S are bounded test oracles only. Not deleted: released
ordinary/sliced interfaces and private recursive paths still require separate
replacement evidence and external-owner clearance. Remaining blocker: bounded
implementation and frozen acceptance, including unmerged termination. Source
changes in this authority pass: 0/0; new files/tests/metadata: none. No release,
Standard60, solver or retirement authority follows from this grant.
Local authority/self-test, deterministic generated views, package load,
docs 8/8 +162/162 +10/10, Documenter, log bound and diff checks passed.
Repo-manager waits for this commit's docs-only CI and Docs before implementation.

### Pass 630: Finite Collinear PQS Closeout

Accepted cf608bd0e151f95ea2f88960959f463184325151 and normalization
3e32d7e50e9a232d8545fa313d2f06ac83565031, tree
6ad007c2d52bc9c095c23de350f1ea834199a256. Small finite chains now produce
complete H1/IDA matrices, nuclear repulsion and raw GTO transfer without a
global coefficient map or retained per-nucleus matrix list. Independent diff,
oracle and evidence review found no blocking issue; the focused owner passed
222/222 locally in 17.1 s. Exact-head CI 34897836745 and Docs 34897836789
passed. Reused accepted core, Cartesian, residual-GTO and matched-H2+ owners
without another angular suite or duplicate paper example.

Frozen dimensions 1265/619, inner coefficients, six old geometries and the
1029+98+98 three-group ownership partition passed. Reviewed equivalent due
diligence: enlarged parent 13x13x17, transverse bounds +/-6.43929, longitudinal
bounds +/-6.62135 bohr, requested padding 6/3, inner 617 and outer 648.
Measured warm construction/operators were .078/.557 s; operators allocated
210 MB with 30.3 MB retained. These are small-fixture costs, not scaling claims.

MT8 is maintenance; MT7 and release/stable identities remain unchanged.
Deleted: duplicate outer decomposition. Simplified: one shared decomposition
and direct nuclear accumulation. Quarantined: dense oracles stay test-only.
Not deleted: ordinary/sliced interfaces and private recursive paths; exact
remaining retirement blocker is matched replacement evidence and external-owner
clearance. Added/deleted source: 263/76 (265 cumulative additions across commits);
tests: 165 added lines (166 cumulative), no new owner; reader docs: 39 added.
New files/persistent metadata/status fields: none. Existing local shell records
reuse the established kernel contract. Frozen tolerances, diffuse parent loss,
px/py asymmetry and quadratic storage remain constraints, not energy accuracy
claims. No successor task or renewed feature-size exception is opened.
Local closeout checks passed: package load, docs 8/8 +162/162 +10/10,
authority/self-test, generated parity, Documenter, log bound and diff checks.
The closeout uses docs-only CI/Docs; numerical implementation evidence is reused.

### Pass 631: Finite Collinear Supplementation Authority

Independent inspection of the report, corrected GH/convolution oracle, logs,
carrier callers and numerical owners supports one four-file connection.
The actual restriction is old frontend/carrier composition, not new residual
mathematics. The two-center guard and old CartesianIDAHamiltonian behavior stay
intact. The new overload uses supplied potential/ownership nuclei independently
of the basis geometry, with exactly H1, IDA and nuclear repulsion in .hamiltonian.
Only the existing carrier gains a concrete Hamiltonian parameter; validation
recognizes two explicit representations. No unrelated consumer clearance follows.

Evidence H3/H-He-H final dimensions 243/235 and actual supplemented metrics
3.23e-12/7.07e-12 are distinct from raw-block 1e-15-level agreement and complete
transformed H1 2.08e-12/5.81e-12. Reviewed the corrected integral-normalized MWG
width convention; the initial squared-Gaussian oracle was wrong, not production.
Warm assembly .111/.071 s and about 2 MB retained justify bounded use, not
H10 scaling. Existing contracted cc-pVTZ checks remain; a shared compact
contracted-candidate oracle and different-potential case protect the new overload.
No numerical campaign was rerun for authority. Source/test/reader added-line
budgets are 100/130, 190/220 and 25/40 preferred/hard, respectively.

LT usable numerical construction advances; MT8 bare construction stays
maintenance while the supplementation connection is active. Deleted: no live
old facade; it owns released semantics. Simplified: reuse residual/MWG assembly
and stream nuclear raw blocks. Quarantined: scratch specialization is evidence,
not a consumer adapter. Exact remaining blocker: implementation and independent
acceptance before hchain-doer resumes. Added/deleted source in this authority:
0/0; new tests/files/metadata: none. No new kernels, approximation, solver,
screening, release, legacy retirement or scientific H10/H20 campaign is granted.
Local authority/self-test, generated parity, package load, docs
8/8 +162/162 +10/10, Documenter, log bound and diff checks passed.
Repo-manager waits for this docs-only authority commit's CI/Docs success.

### Pass 632: Supplement Keyword Load-Order Amendment

Accepted the reported declaration-order obstruction, not the draft implementation.
Independent inspection confirms cartesian_base_hamiltonian.jl loads at root
line 942, before the representation definition included at line 960. The
scratch study ran after package loading and could not expose this restriction.
Permit only a required unannotated supplement keyword followed by its exact
CartesianGaussianShellSupplementRepresentation3D assertion as the first body
statement. A separate minimal Julia probe passed 3/3: late type resolution,
wrong-type TypeError, and missing-keyword rejection. No numerical run was needed.

No strategic change: MT8's supplementation connection remains active and
hchain-doer paused. Existing atom/diatomic behavior, accepted representation,
four-file scope, budgets and numerical gates are unchanged. Deleted/simplified:
the invalid declaration-time annotation becomes an invocation-time assertion.
Quarantined: archived draft remains unvalidated. No source, tests, new files,
metadata or grant expansion in this amendment; no include reorder or type move.
Repo-manager waits for amendment CI/Docs, then resumes and validates the original
packet. No implementation acceptance or consumer authorization is implied.
Local package/docs 8/8 +162/162 +10/10, authority/self-test, generated parity,
Documenter, manager-log bound and diff checks passed. No numerical suite ran.

### Pass 633: Finite Collinear Supplementation Closeout

Accepted cb1b2c7829f4f10e3a653cd541c1998f1d65bf64, tree
287d9453d8cde0bcfa868f1f52ce01c577c67f01, under Passes 631-632.
Independent source/oracle review found no blocking issue; the focused
supplemented tests passed 77/77 locally in 35.6 s. Exact-head full CI
34934021767 and Docs 34934021540 passed. Reused accepted original owner
80/80, core/IDA, Cartesian 232/232 plus collinear 222/222, and old atomic
smoke evidence rather than repeat broad numerical owners.

The required first-body type assertion works; the existing carrier concretely
supports only the old Hamiltonian and new accumulated matrix representation.
Contracted candidates and different potential geometry are covered by shared
independent oracles. Actual supplemented metric, H1 congruence, MWG convolution,
unchanged base IDA and raw transfer pass frozen gates. Raw-block near-roundoff
agreement is not complete transformed-matrix precision or chemical accuracy.
Due diligence: H3 parent 9x9x13, padding 3/3, bounds +/-4.335003225 transverse
and +/-5.206143805 longitudinal; 1053 support sites, 231 base plus 12 residuals.
No construction warning; H-He-H final dimension 235. Warm supplementation
.102 s/111.9 MB allocated/2.10 MB retained is bounded evidence, not H10 scaling.

MT8 including supplementation is maintenance. The implementation prerequisite
for hchain-doer is cleared; scientific calculations still require supplemented
construction and separate consumer scope/convergence. Deleted: no live facade.
Simplified: existing residual/MWG composition with streamed nuclear raw blocks.
Quarantined: stopped draft is superseded evidence. Not deleted: old facades own
released artifact/sector/reweighting semantics. Remaining blocker: no repository
connection blocker; diffuse-space/observable convergence remains consumer work.
Added/deleted source 79/3; tests +220 at the hard cap; reader docs +34.
New files/types/fields/metadata: none; the existing carrier gains one concrete
type parameter. No numerical policy, release, stable or retirement change and
no successor task is opened.
Local closeout package/docs 8/8 +162/162 +10/10, authority/self-test,
generated parity, Documenter, log bound and diff checks passed. Remote closeout
uses docs-only CI/Docs; no broad numerical suite is repeated.

### Pass 634: Separate Residual Validity From Occupied-State Accuracy

Independently reviewed tmp/reviews/h10-gaussian-hartree-qualification-2026-09-16.md
at 5563bce8f5b917c44baf674ab770460f3c668e76. Source inspection confirms the
constructor wrongly applies its state tolerance to the whole residual metric;
the canonical contract already separates these checks. An independent bounded
probe passed 11/11 in 8.5 s including compilation: unused residual error
2.0e-8 satisfies its roughly 1.0e-7 owner bound, accurate occupied recovery is
exact, yet construction rejects. Occupying the perturbed direction fails the
unchanged state check. Invalid cross/identity cases remain independently tested.
Scratch: /private/tmp/pass634_constructor_review.jl; no H10 calculation ran.

Reused and restricted HP-REP-MIXDENS-HARTREE-FN-01/TEST-01 to one constructor
file and its existing test owner. Source 20/30 preferred/hard added lines;
tests 25/35. MT1 gains a contract-conformance repair; MT3 scaling and complete
field certification remain blocked. The report's source-contraction estimates
do not meet the consumer envelope and grant no evaluator or screening work.
hchain-doer's Hartree calculation stays paused. No new scientific policy.

Deleted (required implementation): conflated residual/state tolerance gate.
Simplified: reuse existing residual overlap and scale-aware owner rule.
Quarantined: broader contraction design is deferred, not executable authority.
Not deleted: bounded field oracle and strict spin recovery serve existing callers.
Exact remaining blocker: implementation/acceptance here; molecular field cost
and independent complete-field closure separately. Added/deleted source in this
grant 0/0; new tests/files/metadata none. Validation for implementation includes
bounded owners, unchanged valid outputs, constructor cost, full CI and Docs.
Repo-manager waits for the recorded grant and its checks; no successor is implied.
Local package/docs 8/8 +162/162 +10/10, authority/self-test, generated parity,
Documenter, log bound and diff checks passed. This authority uses docs-only
CI/Docs; implementation must use the unchanged full source-bearing matrix.

### Pass 635: Practical Atomic-Fit Collinear Connection

Accepted the bounded design in
tmp/reviews/h10-atomic-fit-screening-qualification-2026-09-16.md, not an H10
screening result. Independent code/contract review and 18/18 saved-input checks
(0.8 s including compilation) reproduced H normalization, analytic energy,
independent self/cross factors and signed consistency -1.03338551014e-5 Ha.
Both evidence hashes match. Frozen contractions differ from the legacy basis
entry: preserve exact arrays, not a basis-name substitution. Scratch review:
/private/tmp/pass635_atomic_fit_review.jl. No complete field or SCF was run.

The screening paper's translated spherical atomic fits are the practical route.
Density self/cross energies and potential-fit consistency remain independent;
do not force their nonzero difference to zero. Preserve the molecular 6Z(s,p)
reference for optional separate HF matching. The stopped molecular Gaussian
optimization evidence remains, but is not an H10 prerequisite. H needs a truthful
occupation-one input, not an RHF packet. Existing kernels and fit records suffice.

New paired authority owns exactly three source files (95/140 preferred/hard
additions, per-owner 65/25/50) and two existing test owners (65/105, caps 40/65).
Required deletion/simplification: one shared numerical fitter, no H copy; remove
unnecessary dense four-index cloud-self construction through existing pair terms.
Not deleted: RHF spec/packet/writer and all exact/fitted validation branches.
Quarantined: no molecular evaluator implementation. Added/deleted source here
0/0; no new type, file, metadata schema or numerical policy. MT8 gains the
bounded screening connection; MT3's represented-density scaling stays separate.

One field-only production acceptance may load frozen H10 artifacts under
600 s/16 GiB RSS/512 MiB additional-scratch limits, verify full fitted-field
consistency against independently predicted self/cross data, and record actual
cost. The forecast is not acceptance. Hchain-doer remains paused until independent
implementation acceptance and subsequent consumer authority. Pass 634 remains
separate; no new SCF, matching, basis rebuild, Gaussian reference or release work.
Local package/docs 8/8 +162/162 +10/10, authority/self-test, generated parity,
Documenter, log bound and diff checks passed. Repo-manager waits for this
docs-only grant's required CI/Docs; source implementation requires full CI.

### Pass 636: Atomic-Fit Connection Accepted

Accepted bf97c94b5652e6588f113122ff0a93662f718e45 after independent source/test
review, saved acceptance script/log/resource inspection and source-hash matching.
Exact-head CI 35115844724 and Docs 35115844994 passed. The complete H10 field
took 22.53 s (33.35 s monitored process), peak RSS 2.41 GiB, and 16,208 bytes
additional scratch. Signed consistency -1.03338551618e-5 Ha agrees with the
independently predicted self/cross result to 6.04e-14 Ha. This accepts assembly,
not exact fitted-field closure or H10 scientific accuracy. Density-fit,
potential-fit and finite-expansion errors remain separate. Frozen source hashes
bind the evidence to production; physical matrices and frozen inputs survived.

MT8's atomic-fit connection moves to maintenance; the one-shot H10 execution
grant is consumed, not renewable. Hchain-doer stays paused pending a separate
consumer assignment. Pass 634 constructor repair and MT3 scaling remain separate;
no SCF, HF matching, new reference, release work or successor task is granted.

Deleted: dense four-index density-cloud self-energy construction. Simplified:
one shared occupation-generic fitter and streamed raw-block accumulation.
Quarantined: none. Not deleted: RHF packet/writer and exact-field branches serve
live consumers. Exact remaining blocker: separate consumer authority, not atomic
field assembly. Added/deleted source 114/23; tests +89 in two existing owners;
new files/types/metadata fields none. The fixed betas/weights tuple is an
ephemeral existing-kernel argument, not a staged inventory. Independent manager
validation reuses the two small owners; no full H10 rerun or angular campaign.
Manager reruns passed 12+117 and 25+85 checks (28.7 s and 46.8 s). Local
package/docs 8/8 +162/162 +10/10, authority/self-test, generated parity,
Documenter, log bound and diff checks passed. Closeout uses docs-only CI/Docs.

### Pass 637: Private Recursive Chain Deletion Grant

Independent read-only audit at 5a770a571d9527838aaaf540fed515cd03380b97 accepts
the bounded repo-manager caller report. Main has no executable consumer outside
the three-carrier/seventeen-method closure. Repeated local HFDMRG source/test/
validation and detached snapshot, PQS/Hartree/angular paper, H-chain, multisliced,
grid-study and codex-tree scans found no matching caller. Archived high-order
callers remain in their own preserved source; the owner closure makes that lane
historical, not a future import target. This is not universal downstream proof.

The existing finite-collinear owner supplies the recorded H3/H4, contact,
mixed-charge, operator/transfer and three-group ownership evidence. MT8 does not
need H10 convergence to remove unused recursion. Freeze exactly 708 deletions,
zero additions in one file, with baseline/deleted/retained hashes. Preserve the
shared three-child helper and all square-lattice code. The remaining developer
reference is already labeled historical: no documentation cleanup is granted.

Deleted (required): private recursive chain closure. Simplified: no competing
odd/even chain construction. Quarantined: archives unchanged. Not deleted:
released product chains, prebuilt fixed-block consumers, sliced-chain capability
and square code. Exact remaining blocker: implementation and full acceptance;
stop on any newly found live caller. Added/deleted source here 0/0; planned
0/708; new tests/files/metadata none. Existing core/public Cartesian owners,
hash/caller checks and full source-bearing CI/Docs are required after deletion.
No H10, supplementation or screening rerun solely for orphanhood, no consumer
assignment, release or other retirement. Repo-manager waits for the recorded
grant and checks. Scratch static audit: /private/tmp/pass637_closure_audit.jl.
Local package/docs 8/8 +162/162 +10/10, authority/self-test, deterministic views,
Documenter, log bound and diff checks passed. No numerical owner ran for this
authorization; its remote checks use the docs-only route.

### Pass 638: Private Recursive Chain Retirement Closed

Accepted 4b6e1c9dab7696a4a6b9f8198b225450de8924a9: one file, source +0/-708,
three private carriers and seventeen methods deleted. Independent review matched
the frozen retained-file/deleted-range/helper hashes, reran the 12-check package
and x/y/z helper probe, and found only the preserved square-helper definition
and caller. Exact-head full CI 35136279056 and Docs 35136278883 passed; accepted
core, Cartesian 232/232 and collinear 222/222 evidence was not rerun.
No strategic change: MT8 construction/atomic-fit maintenance, separate consumer
authority and the constructor-repair boundary remain unchanged. Deleted: orphaned
recursive chain closure. Simplified: mixed owner 1488 to 780 lines. Quarantined:
none. Not deleted: square/shared helper, released product/sliced chains and
archives. Remaining blocker: none for this retirement; hchain-doer remains
paused for separate assignment. New tests/files/metadata/status fields: zero.
Both temporary grants are closed; no successor or numerical campaign is opened.
Local package/docs 8/8 +162/162 +10/10, authority/self-test, generated parity,
Documenter, log bound and diff checks passed. Closeout uses docs-only CI/Docs.

### Pass 639: General Scientific-q Chain Entrance

Steven's paper-first baseline settles the scientific direction; MT8 now needs
software translation, not an energy-tolerance decision or another fixture
family. Preserve the coarse anisotropic H10 runs as pilots. General q, not a
q5 branch, resolves hydrogen spacing 1.2/(q-1) and the established odd core
width. Fixed-parent comparisons remain explicitly distinct.

Independent production geometry/selector probes at 0386affc used H3/H10 R1.8,
q4/5/6, padding10, tail2.8 and scale1.4. H3 transverse orders drifted to q+1
on 4/7, 2/8 and 3/8 shells; H10 matched on every shell. No parent limitation,
but existing lower-band/upper-safe fallbacks occurred in all cases. The exact
longitudinal sequences and bounds are in the canonical amendment and ignored
standard-q-selector-2026-09-16 report. No terminal or complete operators ran.
The minimal correction explicitly sets transverse q while preserving selected L;
all prescribed dimensions fit the probed boxes. Actual realization still gates
implementation acceptance. First-call probe 1.94s; remaining cases .08-.56s,
not physical construction timings.

Authorize only two existing source owners (80 added lines), the existing public
test owner (80), and two reader pages (35). Outer count starts at q but remains
an independent convergence control. Deleted: none required. Simplified: one
entrance resolves the standard recipe; no duplicated constructor. Quarantined:
none. Not deleted: live explicit expert and released atom/diatomic paths.
Remaining blocker: implementation/qualification; no H10 operators, screening,
SCF or consumer assignment. Added/deleted source here 0/0; no new committed
tests/files/metadata. Tests will catch odd/even source-order drift that the
existing pilot coefficient checks cannot. Required implementation validation:
bounded H3 realization, existing owners and full CI/Docs. This authority pass
passed local package/docs 8/8 +162/162 +10/10, authority/self-test, generated
parity, Documenter, log bound and diff checks. Remote validation uses docs-only
CI/Docs; repo-manager waits for its commit and checks.

### Pass 640: Scientific-q Entrance Accepted

Accepted 2dd316f6aae255c8013833e73daa91579e5b4c42 after independent diff,
scope and exact-head remote review. General hydrogen q now resolves the scaled
parent and odd core recipe; explicit transverse (q,q) preserves selected L.
H3 q4/5/6 successfully realizes 797/1371/2569 columns, not just a rejection
guard. Manager rerun passed 126/126 full-overlap, ownership, positive-weight,
fixed-parent and input checks in 14.9s (13.12 GiB cumulative allocation,
compilation-dominated). Explicit expert blocks reproduced the frozen
419cf3be... serialized SHA-256 byte-for-byte. Production construction timings
remain the bounded implementation report's measurements, not H10 forecasts.

MT8 translation moves to maintenance; LT2's paper-first convention is now an
executable entrance. Scientific q remains distinct from angular calibration,
outer completion and a fixed-parent comparison. Preserve all pilot evidence.
Hchain-doer remains paused: supplemented, atomically screened H10 q5/R1.8,
padding10 needs a separate physical assignment; no energy-accuracy or diffuse
convergence claim follows. Constructor repair and MT3 remain separate.

Deleted: obsolete required-keyword signature/doc wording. Simplified: one
normalization path and a local transverse substitution. Quarantined: none.
Not deleted: live explicit expert and atom/diatomic routes. Remaining blocker:
physical qualification, not general-q translation. Source +44/-9, tests +52,
reader docs +26/-2; no new files, helpers, types, exports or metadata fields.
Mechanical diff and suspicious-addition review passed. Accepted existing owner
evidence plus full numerical CI 35154687776 and Docs 35154687653 are reused;
no full numerical rerun or H10 construction for closeout. Both FN/TEST records
return to maintenance. Local package/docs, authority/self-test, generated parity,
Documenter, log bound and diff checks plus docs-only CI/Docs gate this closeout.
No successor implementation or consumer grant is opened.

### Pass 641: Collinear Residual Cutoff Forwarding

Independent review of the H10 standard-q5 cutoff report, selector call chain,
probe source and hash-matched final log supports one missing facade keyword.
At 1e-8 the existing algorithm naturally retains 210 rather than 208 directions:
occupied Gram error 1.23e-14, atomic norm/charge errors 7.77e-15/1.95e-14,
residual metric max-entry 1.23e-8. All existing gates remain unchanged. The
full-metric row-sum is different and is not substituted for the max-entry gate.
Preserve both outcomes, intentional atomic overlap and frozen inputs. No H10
rerun; the report's 26.34s task/3.42 GiB peak RSS is selection evidence, not
complete-operator, field or HF acceptance. MT8 advances only the input connection.

Authorize the existing collinear supplemented overload (16 added source/docstring
lines), small public test owner (25) and manual (6). Default remains 1e-6;
finite nonnegative validation and direct forwarding suffice. A cutoff-sensitive
regression must catch silently ignored input; existing fixture/oracles suffice.
Deleted: none required. Simplified: expose an existing numerical choice without
a new selector. Quarantined: none. Not deleted: valid default and all other
interfaces. Remaining blocker: implementation acceptance and separate consumer
assignment; hchain-doer remains paused. Added/deleted source here 0/0; no new
files/helpers/types/metadata or physical grant. Local package/docs, authority/
self-test, generated parity, Documenter, log and diff checks plus docs-only
CI/Docs gate this authority; implementation requires normal full CI/Docs.

### Pass 642: Shared Residual Selection Default Amendment

Steven supersedes Pass 641's default-preservation rule: use 1e-8 as the standard
Gaussian-residual occupation cutoff, retain explicit tighter 1e-10 and deliberate
1e-6 comparisons. Forwarding already landed at 7f9dbfd1e; this amendment is not
implementation acceptance. Repo-manager is paused pending the amended grant's
checks. Audit found four default sites: collinear facade, terminal augmentation,
residual builder and protected-ladder recipe fallback. Atom/diatomic construction
inherits this numerical policy; signatures, guards and Hamiltonian meanings stay.
Explicit stored recipe cutoffs and numerical-complete 1e-10 routes stay unchanged.

Reuse the reviewed H10 208/210-direction evidence; no new scientific-default
campaign, operators or physical run. Historical Cr/Cr2 cutoff evidence remains
visible, not newly qualified. Merge, metric, representation, screening and solver
thresholds do not change. MT8 advances the standard selection connection, not
physical acceptance. Incremental hard budgets beyond 7f9dbfd1e: source/docstrings
16, two existing test owners 25, manual 6 added lines. Reuse small-fixture default
parity, explicit overrides and invalid-input checks; update only nested default
provenance assertions, not numerical expectations. Full implementation CI/Docs
and affected existing owners remain required; no H10 or broad angular rerun.

Deleted: obsolete current-default claims, not historical evidence. Simplified:
one shared policy across existing defaults. Quarantined: none. Not deleted:
explicit comparison and numerical-complete choices. Source changes here 0/0;
no new tests/files/helpers/metadata. Package/docs, authority/self-test, generated
parity, Documenter, log bound, diff and docs-only remote checks gate this grant.
Hchain-doer remains paused; after acceptance its separate continuation uses 1e-8
and the saved standard-q5 basis, preserving all 208-direction evidence.

### Pass 643: Residual Selection Connection Accepted

Accept forwarding 7f9dbfd1e48f8fdd080a1eeaf33420506643bab8 and default correction
b64499918bf7a43e806b37aef72dfa69ef49c2e4 after independent diff, scope and
exact-head remote review. The four shared defaults are 1e-8; stored choices,
explicit 1e-10/1e-6, first-body type assertion and all validity gates remain.
Default parity, overrides, rejection and raw-import coverage protect the actual
connection. Incremental source/docstrings +6/-6, tests +15/-12, manual +4/-4
fit the grant. No new file, helper, type, metadata/status field or algorithm.

Reuse reported supplementation 96/96 +80/80, Cartesian 232/232 +126/126 +222/222,
and nested 464/464 +64/64. Nested 487+18=505 dimensions, energies and self-Coulomb
expectations are unchanged; reviewed due-diligence retains axes9x9x15, padding4,
counts275/114/98 and existing warnings. Exact-head CI35185034012 ran all three
numerical jobs; Docs35185033953 passed. No closeout numerical rerun. Local
package/docs, authority/self-test, generated parity, Documenter, log and diff
checks plus docs-only remote gates validate this lifecycle change.

MT8's selection-connection prerequisite is complete, not H10 physical acceptance.
Both RG FN/TEST grants return to maintenance. Deleted: stale default literals
and pending-grant wording. Simplified: one shared selection policy. Quarantined:
none. Not deleted: deliberate overrides, stored recipes, validity checks and
historical Cr/Cr2 evidence. Remaining blocker: separate bounded Hchain assignment
at1e-8 using the saved standard-q5 basis; preserve208 evidence. Hchain-doer remains
paused. No successor implementation, physical run, release or stable change.

### Pass 644: Early Longitudinal Sum Integration Grant

Independent source and evidence review accepts the September17 q7 component
benchmark as integration evidence: nuclear GG948.57->94.64s, identical-H fitted
GG311.50->31.05s, not a measured complete-build speedup. Common transverse
factors permit exact distributivity before contraction; reassociation changes
raw entries near9e-15 and residual-transformed entries up to4.53e-10. Full signed
energy1e-8 Ha and relative action1e-10 gates passed; never substitute a large
cancelling GG-only denominator or claim bitwise equality. Snapshot/README hashes
and frozen q7 identities are in the canonical grant. No local numerical replay.

Authorize only two collinear call sites and the small existing placed-potential
GA/AA separation: three source owners, 40-55 preferred/65 hard added lines
including TimeG; two test owners, 20-35 preferred/45 hard. Existing weighted
nuclear and fitted post-transform oracles cover the live contract. Delete the
per-center collinear GG loops; do not retain a selectable duplicate. The general
placed wrapper remains for RHF/off-axis callers. No kernel, cache, public API,
fallback machinery, fit policy or approximation change. MT8 advances measured
operator construction cost, not new scientific accuracy or long-chain scaling.

One complete q7 supplemented operator+field acceptance is allowed through adviser
coordination on Mac Studio, using saved8509/8719 basis/system/field and unchanged
210-direction reference. Require input/runtime receipt and feasible memory plan;
45min/48GiB/2GiB hard resource caps. Separate construction stage timings from
validation/I/O. No basis rebuild, baseline replay or HF. Deleted: required old
GG loops in implementation. Simplified: sum before expensive contraction.
Quarantined: none. Not deleted: other placed-potential callers and atomic
self/cross accounting. No production edits here; no new metadata/status fields.
Remaining gate: grant checks, implementation, adviser-mediated acceptance and
independent closeout. Package/docs, authority/self-test, generated parity,
Documenter, log bound, diff and docs-only CI/Docs gate this authorization.

### Pass 645: Early Longitudinal Sum Accepted

Accept 6663d392258e5372e018fd8081cd19c3f9a6ed57 after independent diff/hash,
consumer-script/results and exact-head CI/Docs review. All five files match the
frozen candidate; kernel unchanged. MT8 advances complete construction cost,
not H10 convergence or scalable-chain accuracy. Saved q7/8509+210 inputs give
operators 1361.52->499.73s and fitted field/correction 358.78->75.90s: observed
2.99x subtotal, not a controlled cold/warm benchmark. Compilation/cache and
validation/I/O remain separately reported; nested TimeG is not additive.

Energy changes7.11e-15 Ha and maximum complete action9.92e-12 pass unchanged
1e-8/1e-10 gates. Residual-block H1/field differences6.05e-9/9.31e-10 are retained,
not hidden by bitwise claims. V, residual, transfer and fit arrays are exact;
nonzero fitted consistency remains. The interrupted scalar-audit typo attempt
is preserved, followed by the reviewed corrected replacement; cumulative951.93s,
peak15.83GiB and39.46MiB scratch fit original caps without reset.

Evidence: tmp/reviews/pass644-implementation-2026-09-17.md and consumer
h10_pass644_acceptance_20260917/replacement/README.md (SHA256
c85cb62e9afc5f577fa4c5f03fac044573b3be6d91b19d5dc42c8c1aed7bcf5e).
Reuse 1032 local numerical assertions and exact-head full CI35285875675 /
Docs35285875678; no closeout numerical replay. Local package/docs,
authority/self-test, deterministic views, Documenter/log/diff and docs-only
remote gates validate closeout. Deleted: per-center collinear GG loops.
Simplified: early z sum and one GA/AA owner. Quarantined: none. Not deleted:
general placed wrapper for live off-axis/RHF callers. Added/deleted src46/23;
tests21 added in existing owners; new files/metadata/status fields0.
Both grants return to maintenance; acceptance is consumed. Remaining blocker:
separate consumer assignment. Hchain-doer stays paused; no successor opened.

### Pass 646: Bounded Screening Memory Grant

Independent review at97db4dea3 accepts the two-function proposal, with the user's
explicit exception to net-negative refactor policy: expected source+17/-3,
hard25 added source lines; expected29/hard35 test lines in the existing owner.
No production edit here. Reproduced75 scratch and17 compact checks; combined
P0/q0 bit-identical, trace/trace-loss changes near1e-15, input/conversion errors
and ownership preserved. n768 additive allocation248.485->23.304MiB and dense
trace14,156,016->0B confirm the targeted waste; warmed single calls are not a
statistical timing or peak-RSS claim. Proposal and frozen patch hashes live in
Screening Memory Repair; reuse existing small-owner evidence, no heavy run.

MT8 advances economical screening assembly, not q9 scientific qualification.
Deleted: required per-atom dense densities and scalar-product temporaries.
Simplified: scalar per-block diagnostics and Float64 Frobenius reduction.
Quarantined: none. Not deleted: combined P0, _sym, validation/result copies and
copy-returning accessors with separate contracts. Exact remaining blocker:
q9 correction window~40.83GiB plus unreclaimed prior arrays may exceed46GiB;
keep6GiB reserve and separate adviser admission. No automatic follow-on cleanup.
New tests:17 compact checks, no new owner; new metadata/status fields0.
Existing correction/public screening owners and full source CI validate the
future implementation. Package/docs, authority/self-test, deterministic renders,
Documenter/log/diff and docs-only CI/Docs validate this grant. Hchain-doer remains
paused; no q9 loading/build, HF, checkpoints, release or numerical-policy change.

### Pass 647: Screening Memory Repair Accepted

Accept 8744b01c26faf29608677a7463554a21c837c376 after independent diff, frozen
patch/regression identity and exact-head remote review. Two source functions
only, +17/-3; existing owner +29 with17 checks. The user-approved net-positive
exception is consumed, not a standing policy change. MT8 reduces unnecessary
screening storage without changing combined P0/q0, thresholds or ownership.

Reviewed integrated n768 evidence: additive248.485->23.304MiB; dense trace
14,156,016->0B. Scalar differences near1e-12 are reassociation, not bitwise
equality; P0/q0 remain exact. Measurements are warmed single-call allocation
traffic, not q9 peak RSS or end-to-end timing. Evidence:
tmp/reviews/pass646-implementation-2026-09-17.md and linked allocation script/log.
Reuse nested132/132, public22/22 and integrated6/6; prior independent75/75
scratch checks apply to the exact patch. CI35311484055 ran all three numerical
gates; Docs35311483858 passed. No numerical rerun for lifecycle closeout.
Local package/docs, authority/self-test, generated parity, Documenter/log/diff
and docs-only remote gates validate this pass.

Deleted: per-atom full densities and trace copies/product. Simplified: scalar
Gram diagnostics and Frobenius reduction. Quarantined: none. Not deleted:
combined density, _sym, validation/result copies and public accessors.
Added/deleted src17/3; tests29 in existing owner; new files/metadata/status0.
Both grants return to maintenance. Exact remaining blocker: q9 correction
window~40.83GiB plus unreclaimed arrays is not a bound below46/48GiB; preserve
6GiB reserve and require separate adviser/application-owner admission.
Hchain-doer remains paused. No successor, heavy run or lifetime repair granted.

### Pass 648: Basic Gaussian Integral Repair Authorized

Task GB-BASIC-INTEGRAL-20260918 advances LT1/LT2: accurate shared overlap
inputs before interpreting residual accuracy. Steven's bounded manager-cycle
ARM authorizes this packet, not q9 admission. Add HP-GAUSSIAN-BASIC-ARITH-FN-01/
TEST-01 under the existing non-nuclear contract; the centered/displaced grants
and all residual policy stay unchanged. Baseline 68a05b6c1 and audit report
SHA-256 56cc10241f63bf31b872b3d20e24e2740a2757bceb4c54e27a05a5ea18f356c7
are frozen. Independent review reproduced the actual primitive's 1.07033e-11
error and relative-coordinate reduction to 2.85e-17. The final scratch moment
summary was unavailable: no claim it passed. Reuse completed audit evidence;
do not load q9 or replay it. The proposed repair removes absolute-center
cancellation while preserving signed accepted exponents, absolute moments and
kinetic semantics. Focused independent checks are required before acceptance.

Deleted/simplified: only the old damping/shift expressions in future source
work. Quarantined: none. Not deleted: polynomial/moment machinery and callers.
Budget: 30 added src lines, 40 added existing-core-test lines; new files,
helpers, metadata/status fields zero. Local package/docs 8/8 + 162/162 + 10/10,
authority/self-test, two equal external renders, Documenter/log/diff passed.
Remote docs-only CI/Docs must pass before handoff; focused caller gates and
full source CI are required for implementation.
Exact remaining blockers: amplified input errors, real base-metric defects and
residual merge arithmetic require separate physical requalification, as does
q9 resource admission. No q8 energy conclusion or consumer restart. MT8's
allocation work remains completed, not q9 admission. Hchain-doer stays paused;
this bounded cycle ends at kernel acceptance/closeout, no successor.

### Pass 649: Basic Gaussian Integral Repair Accepted

Accept 1dbc1856d7b699311dcf65ae6a68fb193f9c1009 after independent source/test
diff, hashes and exact-head CI/Docs review. Task GB-BASIC-INTEGRAL-20260918
ends here; both implementation grants are consumed and become maintenance.
LT1/LT2 advance through accurate overlap inputs, not a physical residual or
energy claim. Signed accepted exponents, absolute moments, prefactors, errors
and kinetic callers are preserved. No correction round was needed.

Evidence: tmp/reviews/pass648-basic-integral-implementation.md, SHA-256
5e03934068dfeaea0001e6c6be334c1f61fef1ff6a6db5d2eaace8ced9100b84.
Manager reran 81 committed checks plus 33 independent absolute-coordinate
256-bit checks in /Users/srw/dmrgtmp/pass648_design_review.jl (0.95 seconds),
reproducing primitive error 2.85e-17 versus old 1.07e-11. The previously lost
preflight summary remains unclaimed. Reuse doer core/IDA/Cartesian/residual,
atomic-packet/screening owners and full CI35370892362; Docs35370892353 passed.
Small warmed scalar/moment/kinetic timings overlap baseline ranges; allocation
counts are unchanged, not an end-to-end performance claim. Closeout local
package/docs 8/8 + 162/162 + 10/10, authority/self-test, two-render parity,
Documenter/log/diff passed. Remote docs-only CI/Docs must pass before the
completion notice; no full numerical replay.

Deleted: absolute quadratic and polynomial-center subtractions. Simplified:
relative-coordinate arithmetic. Quarantined: none. Not deleted: shared moment
machinery and callers. Added/deleted src10/10; tests35, 81 assertions in the
existing owner; new files/metadata/status0. Exact remaining blockers: accurate
input propagation, merge arithmetic and real base-metric effects need separate
physical residual requalification; q9 resource admission is also separate.
MT8 allocation work remains completed, not admission. No q8 energy conclusion,
q9 load, consumer restart or successor. Hchain-doer stays paused; send only the
authorized adviser completion notice after closeout checks pass.

### Pass 650: Residual Premerge Arithmetic Authorized

Task GB-H10-Q8-RESIDUAL-20260918, explicitly armed by adviser under Steven's
bounded cycle, is distinct from the closed integral repair. At aa2ac41c5,
fresh q5 passes and fresh q8 fails with identical270 ordinary s/p/d candidates;
both naturally retain270, not an assumed rank. Independently reviewed diagnosis
tmp/reviews/h10-q8-residual-diagnosis-2026-09-18.md (SHA-256
1e7e7495519d082222c89ecc676fa3cae915e703438644d50f46fcb2792e5233)
attributes the1.13e-7 accurate rounded-input defect to premerge cancellation,
not the inverse square root or a marginal owner cutoff. Physical spot checks
still fail; neither rounded-input success nor the old mixed-provenance replay
certifies the fresh basis. General D identity independently passed10 compact
256-bit checks in0.84s; no q8 rerun or candidate whitening by design-manager.

Authorize only general premerge reassociation under HP-RG-ORTHO-FN-01/TEST-01,
with unchanged final assertion. Deleted/simplified: one premerge four-term
evaluation, not the live overlap helper. Quarantined: none. Not deleted:
final independent validation, ordinary/injected callers and all safeguards.
Budget: eight added source lines,35 existing-misc-test lines; no new production
helpers/files/metadata/status. LT1/LT2 advance through accurate construction,
not acceptance by loosened checks. Diagnosis used87.029s/2.425GiB/140.668MiB;
remaining task-specific computation shares the20min/16GiB/2GiB envelope.
Require fresh q5/q8 unchanged gates and independent rounded/selected physical
checks, active owners/full source CI, plus docs-only checks for this grant.
Exact blocker on failure: this candidate is insufficient, not permission for
another policy or rank choice. Hchain-doer stays paused, no successor. This
does not reopen q9 resource admission or claim a q8 energy correction.
Authorization checks passed: package/docs8/8 +162/162 +10/10, authority and
adversarial self-test, two independent matching renders, Documenter/log/diff.
Require exact-head remote docs-only CI/Docs success before the implementation
handoff; no numerical owners were replayed for this authorization.

### Pass 650 Amendment: Occupation Snapshot Portability

Steven explicitly superseded the two occupation golden-value tolerance freezes
through AMEND-GB-H10-Q8-RESIDUAL-20260918-2. The original R3A failure remains
463/464: a 1.6896e-14 minimum-occupation difference exceeded atol1e-14, and
the later64-check section was not reached. Static upstream dataflow does not
prove its historical cause; no baseline campaign or causal claim is required.
Allow exactly two atol1e-12 substitutions and up to two portability-comment
lines in the existing R3A owner, preserving both values. Relative bounds near
2e-9/8.2e-11 remain strict at the two observed scales. No physical, occupation
cutoff, rank, metric, capture, action or energy gate changes.

Preserve the uncommitted source+3/-1 and misc+32 draft and the successful
rounded/selected-physical q5/q8 evidence in the implementation report (SHA-256
c89d539f171edb123679ee8f2c5987deab03bed1b0fcca6e702c6a615f79130f).
Deleted/simplified: only brittle snapshot precision; quarantined: none; no new
helper, file or metadata. LT1/LT2 and all consumer restrictions remain unchanged.
Require complete R3A including64 checks, full implementation CI/Docs and review;
no q5/q8 replay. Remaining blocker is implementation acceptance, not a waived
failure. Hchain-doer stays paused; no successor. This amendment changes only
documentation authority; the source/test draft is excluded from its commit.

### Pass 651: Residual Premerge Repair Accepted

Accept b974887980ef90bca6e753c049376d4548883dc0; close
GB-H10-Q8-RESIDUAL-20260918 and return both ORTHO grants to maintenance.
General premerge reassociation fixes demonstrated construction cancellation
without dropping injected D, changing rank/cutoffs or weakening final checks.
Fresh q5/q8 naturally retain270; accurate rounded maxima3.09e-9/4.45e-9 and
selected physical maxima5.24e-9/1.86e-8 meet the assigned bounds. Selected
physical checks are not full certification. Independent review inspected the
oracle code, signed contributions, log/report hashes and exact committed diff;
no q5/q8 replay. Report: tmp/reviews/pass650-q8-premerge-implementation.md,
SHA-256 b43f6c540e75ac31a71521c48a42e43a5a398278eeb48c36f43f4d4ca196b60e.

Original R3A463/464 failure remains evidence. Steven's narrow snapshot amendment
changed two bounds to1e-12, not their values or physical checks; amended R3A
464/464 and facade64/64 passed. Cause of the old final digits remains unclaimed.
One amended continuation sufficed. Manager reran misc67 and docs8+162+10;
reuse doer core/public owners, full exact-head CI35425393601 and
Docs35425393550. Closeout local package/docs, authority/self-test, two-render
parity, Documenter/log/diff passed; require docs-only remote checks before
adviser notification.

Deleted: one four-term premerge call. Simplified: premerge arithmetic.
Quarantined: none. Not deleted: final independent overlap assertion and injected
callers. Added/deleted src3/1; existing tests35/2, no new test files or
metadata/status. q8 premerge adds62.88MiB temporary allocation and4.77ms;
no material total construction regression observed or speedup claimed.
Bounded task used under211s, peak2.77GiB, scratch295MB. LT1/LT2 advance through
accurate construction; physical consumer admission remains separate. q9,
H10 operators/fields/HF, release work and successors remain excluded.
Hchain-doer stays paused; send the single final adviser notice after checks.

### Pass 652: Stable Terminal Residual Construction Authorized

Steven approves H10-SUPPLEMENT-IMPLEMENTATION-CYCLE-20260920 through Screening-
advisor, including the explicit240-added-source-line exception. Reuse all270
evidence rather than another design campaign: direct functions pass full
refined Gram/cross and selected analytic checks, while rounded final validation
still rejects them. This does not measure a full H10 energy error. Freeze a
four-grid order/tail schedule, residual-level stabilization and checked memory
admission before source work. Preserve symmetric Lowdin orientation with the
small QR-factor SVD; scratch QR gauge is not an IDA-equivalent replacement.

LT1/LT2 advance accurate usable supplementation, not another metric framework.
MT boundary is non-injected function-aware terminal construction and validation
together. Matrix-only/nonzero-D, injected/protected and operator definitions
remain unchanged. Canonical contract owns exact source30/210, test35/100 and
reader15/20 ceilings, selected operator gates and one cumulative60min/16GiB/8GiB
saved-H10 qualification. No HF, consumer restart or successor; at most two
in-scope correction rounds. Old evidence and both handoffs remain preserved.

Deleted/simplified target: cancellation-prone terminal validation and duplicated
entry routing, not selection. Quarantined: old formula only as diagnostic on
the repaired branch. Not deleted: matrix-only finalizer and injected callers.
Added source/tests/metadata in this authority pass: zero. Independent review
used source closure, scratch receipts and exact assignment/design hashes; no
numerical replay. Require local package/docs/authority/self-test/generated/build/
log/diff and exact-head remote docs-only checks before dispatch. Failure of
scope, budget, physical gates or admission stops the cycle for review.

Authorization local checks passed: package load, docs8+162+10, authority/check
self-test, two matching renders, Documenter with existing size warnings,
log1883/2000 and diff. Source/tests remain untouched. Remote exact-head
docs-only CI/Docs must pass before the single bounded implementation dispatch.

### Pass 652 Amendment: Reviewed Performance Resumption

Steven accepts the roughly0.6s easy-H2 overhead; preserve the original failed
0.1s gate and5.001s draft evidence. Local primitive specialization measured
0.613841s median with exact returned-field parity. Freeze0.80s warm median and
550MB cumulative-allocation bounds on that same fixture/settings; all physical,
grid/tail/GC and memory gates remain. Specialization must preserve supported
input types, not copy scratch Float64 assertions as a public restriction.

The reviewed saved-q9 construction and checks pass: all270 retained, full span
and occupied charge, independent selected kinetic2.46e-11Ha; rounded selected
kinetic/metric Ritz0.907microHa passes with limited margin. Not full-HF or IDA
accuracy. Charge390s/3600s, peak11.010GiB, scratch118,048,990bytes; reuse only
with final-code correspondence, no new campaign. Reports remain frozen under
tmp/reviews/h10-supplement-{primitive-specialization,saved-q9-qualification,
saved-q9-checks}-2026-09-20.md; current phase record is the separate resume report.

LT1/LT2 unchanged: finish usable, accurate supplementation; no consumer restart.
Remaining small operator/compatibility tests and independent acceptance remain
mandatory; zero of two correction rounds used. Deleted/simplified: only the
superseded task-local cost rule; quarantined:none. No source/test/metadata
additions in this amendment; preserve both source drafts outside its commit.
Original source240/test135/reader35 caps and whitelist remain; exact blocker
is remaining implementation validation. Require docs-only authority/build/CI
checks before dispatch; stop at accepted closeout, with no successor.

### Pass 652 Amendment: R3A Same-Function Contract

Steven authorizes the single R3A test-path reconciliation after460/464 failures;
two bitwise coefficient assertions and two cross-route operator comparisons
assumed the old shared finalizer. Their1.4461e-10/1.2551e-10 differences remain
failed evidence, not a proof of harmlessness or a kernel defect. Permit checked
full augmented-coordinate or independent same-function comparisons, retaining
1e-10 and tighter physical bounds. Preserve all snapshots/IDA/base/energy tests.

Reallocate135 added test lines: misc25/public80/R3A30. Source240/reader35 and
all numerical/memory/domain gates unchanged; required coverage cannot be dropped
to fit. Independent small operator tests and full R3A including64 facade checks
remain. Easy-cost0.61993s passes; q9 evidence reused only by correspondence.
LT1/LT2 unchanged. Actual occupied weights inform energy relevance, but do not
waive gates or certify virtual-space/correlation/relaxation accuracy. No HF/FCI.
Deleted/simplified target: obsolete identical-finalizer test pressure; no source
or test edits in this amendment, no new metadata. Exact remaining blocker is
test-contract implementation and independent regression acceptance. Zero of two
correction rounds used; no H10 replay, consumer restart or successor. Preserve
drafts and stop for new scope. Require checked docs-only grant before dispatch.

### Pass 652 Accepted: Stable Terminal Residual Functions

Accepted8b7536f72f6d8f32fc4c67f66b36192613907a6f after independent source,
contract, oracle and exact-head CI review. Direct residual formation retains
outside-parent Gaussian content; symmetric normalization and four-grid/tail
checks replace cancellation-sensitive terminal finalization. R3A compares the
same full functions without changing1e-10 limits, snapshots or IDA policy.
Reviewer public residual tests120+80 passed in72.68s including compilation;
inspected R3A465+64, compatibility/occupied/protected logs and full numerical
CI35551226941/Docs35551226944. Easy median0.61993s/452962584bytes passes the
user-amended gate. Prior performance/R3A failures and fixture mistakes remain.

LT1/LT2 advance to accepted construction, not consumer or energy certification.
Selected q9 Ritz0.907microHa has limited margin; tiny occupied-weight RR effects
exclude GR, nuclear/J/K, unselected virtuals, relaxation and false minima.
No H10 replay:390s/11.010GiB/118048990bytes remains frozen. Hchain-doer stays
paused. Implementation grants consumed; maintenance only, no successor.
Deleted:6 source lines. Simplified: shared selection/finalizer dispatch and
obsolete coefficient equality. Quarantined:none. Not deleted: matrix-only and
injected paths have live callers. Remaining blocker: separate consumer authority
and physical qualification, not this completed cycle. Added src170; tests96
added/9 deleted in existing owners; reader21; new files/metadata/status:none.
Mechanical diff gate and all per-owner budgets pass. No independent correction
round required after final handback; local docs-only closeout checks follow.

### Pass 653 Authorized: Radial Multipoles and Angular Cap

GB-ANGULAR-MULTIPOLE-20261001 reviews packet 01+02+improvements, not03/duplicate00.
Independent baseline erf-grid subsample reproduces L17 norm error9.17e-7 and
all-zero L18/L24; direct pair sums agree with selected256-bit checks below4e-16.
Freeze full-moment auto, explicit legacy cap, sorted points, old constructor
arities and truthful external unknowns. Caps:330 source/250 test/85 reader added
lines in existing owners only. LT1/LT2 accuracy and provenance advance; no
consumer rebuild, quadrature/high-l redesign, Cr2 shift claim or release.
Deleted/simplified target: global scaling/recovery and repeated radial sampling;
quarantined:none. Existing dense export remains; carrier compatibility stops
on a live unsupported binary consumer. No source/test change in this grant.
Required exact-candidate angular run once, remaining owners and full CI/Docs;
post-acceptance read-only cap inventory. Zero of two correction rounds used.
Grant checks: package load, docs8+162+10, authority/self-test, two matching
renders, Documenter and diff pass; log1977/2000. Source/tests remain unchanged.

Pass 653 continuation: Steven authorizes only the existing docs changelog check
to allow Unreleased while requiring first released v0.2.1. Six added/three
deleted docs-test lines fit the unchanged250 test cap (angular allocation114).
No strategic change: LT1/LT2 scientific gates and exclusions remain frozen;
eleven draft hashes and completed angular evidence preserved, no acceptance yet.

### Pass 653 Accepted: Multipole Arithmetic and Full Moment Span

Accepted35684b6ddeb0433f6f6e1d9ee324c9a6173d7cb9 after full diff/contract review; all per-owner budgets pass.
Independent367-point direct/256-bit oracle: L17/18/24 nonzero, norm errors1.45-1.81e-15;2.976s. Production changelog checks13/13.
Reused exact angular61834 (one optional HFDMRG skip), radial435/core2155/ida12058/public89/misc72; no angular replay.
Exact-head CI36922280281 ran all three numerical jobs; Docs36922280251 passed. Reviewer docs8+162+10, authority/self-test/renders and build pass.
LT1/LT2 advance arithmetic/provenance correctness, not continuum or molecular-energy certification; harmonics/quadrature remain separate.
Read-only inventory: Be15 stored manifest cap4 confirmed, actual-cap fields absent; remote Be ladder/Cr/Hooke caps remain partly source-derived or unknown.
Owners must review coherent rebuilds and changed identities/receipts; no mutation or consumer restart. Ne is closed-shell; its proxy is not a Cr2 shift.
Warm radial assembly2.713->1.264s, allocation3.708->1.722GB. Auto angular includes more physics; no application speedup claim.
Deleted:108 source lines, including two orphaned global-scale helpers. Simplified: local-ratio recurrence and reused radial sampling.
Quarantined:none. Not deleted: stored-cap and exact radial/Ylm routes have live callers. Remaining blockers: separate consumer review and C/D limitations.
Added src281/tests218 (3 deleted)/reader85; no new files/exports. One sample owner/field, five-field cap plan and four metadata keys were explicitly granted.
Mechanical diff gate passed. One approved placement continuation used; original failure preserved. Both task execution grants consumed, no successor.
Hchain-doer remains paused; immutable releases/stable, scientific tolerances and solver anchors unchanged. Closeout uses docs-only checks.

### Pass 654 Authorized: Normalized Real Spherical Harmonics

Steven explicitly activates C as GB-HIGH-L-HARMONICS-20261001 after closed A+B.
Patch03 is a candidate; old order-460 tables are not repaired integration evidence.
Independent scalar review 96/96, 0.860s: all m at 0/1/2/87/88/89/151/256,
generic/equatorial/pole/near-pole directions, worst scaled error 8.96e-14;
freeze 1e-12 accuracy with genuine zeros allowed. Tracked caller scan clears
the two obsolete Legendre helpers. Source hard 45/net decrease, existing angular
test hard 110, Unreleased reader hard 10. No source/test edit in this grant.
One order-460 profile/table qualification is mandatory; 45min/8GiB/512MiB cap,
no second long build. Full angular once on exact candidate. LT1/LT2 accuracy
and MT6 deletion advance; no continuum/energy certification or consumer restart.
Rotation archives 812 accepted lines verbatim, hash in index, retaining 625-653;
strategic preamble now reflects accepted release/residual decisions. No history
rewritten. Deleted/simplified target: unused factorial/polynomial chain, not A+B.
Quarantined: none. Exact remaining blocker: implementation and high-order table
qualification; D remains separate. No new source files, exports or metadata.
Require checked committed authority and docs-only grant CI/Docs before dispatch;
implementation requires normal full CI/Docs. Zero of two correction rounds used.
Local grant checks pass: docs 8/162/10, authority/self-test, two renders,
package load, Documenter without deployment, log/archive identity and diff.

### Pass 654 Accepted: High-L Harmonics And Actual Moments

Accepted 6a698b91609833e0b22259a1b42755211c93d47d after exact diff, hashes,
independent oracle and remote review. Reviewer 119/119 in 2.10s: all-m scalars,
real conventions, addition identity and selected saved L89/96 moments. Reuse
once-only angular 63541 (one optional HFDMRG skip), radial 435/core 2155/ida 12058/
public 89/misc 72 and docs 8+162+10. CI 36941084574 ran all three numerical jobs;
Docs 36941084454 passed. One order-460 qualification 1437.6s/1.081GiB/54.02MB
checks 9409 rows and 136 significant high-m rows, scaled error 2.993e-14;
independent first stop 96/lexpand 94. No profile or full-angular replay.
LT1/LT2 advance reliable coupling/table arithmetic, not continuum/energy
certification. Tail 1.991e-7 survives incremental stopping. Small coupling
roundoff 3.33e-16 permits identity changes; no downstream receipt repinning.
Warm scalar allocation 0; coupling 4160B unchanged, 66.875 versus 64.208us;
high-l failing baseline timing is not an application speedup claim.
Deleted: 35 source lines including two unused helpers. Simplified: normalization
inside one recurrence. Quarantined: none. Not deleted: live harmonic/coupling
and moment owners. Remaining blocker: quadrature D and owner-reviewed consumers,
not C. Added src 18/tests 64/reader 6; net source -17, no new files/metadata/status.
Mechanical gate passed; zero correction rounds. Both execution grants consumed.
Local closeout package/docs 8+162+10, authority/self-test/renders/build/log/diff
pass; remote docs-only checks required. Hchain paused; releases/stable preserved.

### Medium-Term Checkpoint After Pass 654

- MT1 active: demonstrated A+B+C arithmetic/default-cap repairs completed; remaining quadrature limits separate.
- MT2 completed: no new Cr2 scientific interpretation or work.
- MT3 deferred/blocked: standard60 awaits a named consumer; represented Hartree scaling remains blocked.
- MT4 active: residual/protected evidence and consumer qualification stay distinct.
- MT5 maintenance: fail-closed authority, immutable release identities and stable pin preserved.
- MT6 active: obsolete harmonic helpers deleted; no replacement framework or new cleanup grant.
- MT7 completed/maintenance: external transfer unchanged.
- MT8 maintenance: accepted collinear machinery does not certify Hchain energies or restart consumers.
