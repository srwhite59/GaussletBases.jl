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
- This live volume begins with Pass 538. Pass entries are preserved in accepted
  order; duplicate or nonmonotonic historical pass numbers are not rewritten.

## Current Strategic State

- The broad producer-documentation reorganization is complete. Schema-v3
  `authority.toml` is authoritative, generated registry/AGENTS views are
  checked one-way outputs, and authority CI is fail-closed.
- The first complete post-cutover static conformance audit covered all `150`
  execution records: `107` matched, `11` documented gaps, `8` numerical gates,
  and `24` discrepancies. Pass 399 closed the atomic-packet fail-fast subset.
- Current screened-Hartree source and contracts use determinant orbitals for
  `P0/q0`, the density fit for `E0`, and the fitted potential as an approximate
  `J0` evaluator with reported consistency error. A determinant-exact `J0/E0`
  convention requires a separate scientific amendment.
- REQ-084 stopped correctly before molecular-full Cr2 interpretation because
  no repo owner constructs the complete Hartree field of a represented density
  containing terminal and supplement components. A direct, typed internal
  producer is approved. Pass 459 implementation preflight established that the
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
  unscreened Hamiltonian as private diagnostic/provenance surfaces. Six
  unsupported root exports have been removed; physical target selection,
  transition-density exchange, and PRF-to-GTO-residual interactions remain
  consumer questions.
- Explicit charged electron sectors are implemented without changing bare
  basis or operator arrays. Specialized retained-GTO EGOI remains archived
  and deferred with no execution grant.
- The expert sliced hydrogen-chain producer is accepted with compact finite and
  periodic-template storage, analytic Galerkin H1, longitudinal
  IntegralDiagonal Vee, and checked long-range continuations. Solver and MPO
  policy remain consumer-owned.
- The immutable annotated `v0.2.0-rc1` tag is published at object
  `a4284f0bf448fb9d717de26ccbe1e9fc16db5ed2`, peeling to frozen target
  `1546c18d3058cce2b5051b50788cda3c12585e51`. Its canonical folder is live,
  `versions.js` lists RC1 and `dev`, `/stable/` remains absent, and both
  canonical folders are intact. Tag acceptance is closed. GitHub release
  `373460389` publishes the exact package-centered RC1 prerelease with no
  uploaded assets or latest-final status; archive, clean-install, and docs
  verification passed. PQS and reference-density screening remain distinct
  method and future citation surfaces. Release acceptance is closed.
- The exact `v0.2.0-rc2` candidate is prepared at commit `2b3c23970`. It adds
  the post-RC1 public residual-GTO/interchange surface, three reader-front-door
  links, and explicit RC2 selector retention without changing implementation,
  dependencies, workflows, or scientific policy. Its immutable annotated tag
  is accepted at object `7c8a21b99`, and package-centered GitHub prerelease
  `376503169` is published and validated. Direct final `v0.2.0` candidate
  `adfcaba32` is accepted without RC3 under a modest package-nearest-software
  claim. Its immutable annotated tag is accepted at object `722e8e875`, and
  final GitHub release `378216554` is published as the latest non-prerelease
  with zero uploaded assets. The fresh namespaced-ref recovery lane closed the
  frozen tag workflow's local-ref collision without rerunning numerical gates.
  Versioned and stable documentation, archives, and clean installation are
  accepted; registration and citation remain separate.
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
four-name package/internal export repair are closed without numerical changes.
Keep each subsequent repair narrow and evidence-led.

**MT2 - Controlled Cr2 source migration (completed).** The source-backed
fixed-state and bounded replay reproduced the former consumer-local path. Any
new Cr2 endpoint, contraction, exchange, or solver interpretation is a
separate scientific choice.

**MT3 - Approved pending producer facilities (active).** The private
matrix-free parent/PQS/WL H2+ controls, parent capture, supplemented H2
one-body preflight, fixed-state density-density gate, and matched PQS/WL
shared-shell construction and explicit charged electron sectors are
implemented. The external same-density oracle must use the corrected WL
fingerprints before a method-accuracy claim. RHF, scans, and ladders remain
separate. Standard60/canonical-driver exposure remains pending.
The specialized retained-GTO EGOI helper/test are deferred with no execution
grant; the archived WIP is not a mainline merge candidate.
The minimal sliced hydrogen-chain producer is implemented and in maintenance
as a separate G10-plus-transverse operator facility; it is not a
Cartesian/PQS route or solver lane. Its stable consumer fields are fixed in
the canonical contract; HFDMRG adaptation remains downstream work.
The represented mixed-density direct Hartree producer is approved for bounded
implementation with a streamed/tiled contracted PGDG/Gaussian pair action. Its
first pass is direct-only; fitted acceleration and Cr2 acceptance remain later,
separately reviewed gates.

**MT4 - Residual and protected-basis evidence (active).** Keep the residual
spectral audit measurement-only. Protected atoms, counterpoise, and any new
injection/localization policy remain separate future decisions. Parent residual
function mechanics, the onsite-calibrated Gaussian direct resource, and the
private diagnostic/provenance wrappers are implemented and source-backed.
Hooke owns the first Be `1s/2s` target study. Selection, transition-density
exchange, and PRF-to-GTO-residual interactions remain consumer or measurement
questions.

**MT5 - Documentation and authority maintenance (active).** The broad
reorganization and cutover are complete. The RC1 candidate, immutable tag,
versioned folder, explicit selector entry, and package-centered GitHub
prerelease are accepted. The RC2 candidate, immutable tag, versioned folder,
selector entry, and exact GitHub prerelease are also accepted. Path-aware CI is
implemented and accepted as a bounded release-efficiency refinement. The exact
final candidate, immutable annotated tag, final/latest GitHub release,
versioned/stable documentation, archives, and clean installation are accepted.
The tag-lane local-ref collision was recovered through exact namespaced manual
verification without a numerical rerun. Separate PQS and screening surfaces
remain fixed. The bounded Julia `1.10` Supported-floor extension to the live
`radial`, `misc`, and fixed-radial `angular_public` groups is implemented and
accepted. The complete `angular` research group remains outside per-push CI.
Future tag-lane repair, registration, citation metadata, and any later release
are separate decisions.

**MT6 - Carrying-cost control (active).** Remove stale helpers, compatibility
metadata, unsupported exports, and development-era tests as conformance work
identifies them. New scaffolding must advance a live physics or module contract
and account for what it replaces. The repaired ordinary-QW diatomic front doors
remain supported; the broader endpoint/correction audit remains a deferred
retirement map. Direct support-local PQS shell seeds, matched-H2+ parent-action
scratch reuse, terminal buffer reuse, and the exact-order four-element
Gaussian-sum reduction are implemented. The next performance investigation is
the cold reporting boundary; it has no implementation grant. Eight-lane batching,
further loop restructuring, compatibility cleanup, and release work remain
separate. The duplicate matched-H2+ release/example execution is removed by an
accepted test-only replacement that preserves Example 41 and all release
assertions. The fixed-radial angular public-owner extraction is implemented as
an `80`-line, `83`-assertion owner replacing a `163`-line block, with total
test delta `+88/-165` and net reduction of `77` lines; angular export
documentation and shell-key de-export remain separate. The bounded post-v0.2
export cleanup is implemented: the unused
timed nested wrapper is deleted, two undocumented diagnostic/QW names are
unexported, and the inexpensive qualified QW alias remains. The angular public
surface audit now authorizes documentation of seven durable experimental
profile/sequence bindings without treating the cache key as consumer API. One
bounded foundational documentation grant and four successor packets document
the supported basis/mapping,
function/stencil, partition/leaf, atomic-IDA reference, and final six
supported-public families without changing behavior or global documentation
policy. Their bounded sequence reduced the undocumented-export backlog from
`71` to `24`. The read-only follow-up retains
`bond_aligned_diatomic_geometry_payload` as expert/experimental and gives it a
bounded documentation grant. Three one-center nested names form one reserved
next-minor namespace transaction; `diagnose_qwrg_residual_space` and its
already-documented `QWRGResidualSpaceDiagnostics` type form another. No export
change is authorized. Accepted geometry documentation left `23` undocumented
exports. The angular packet must reduce that count exactly to `16` while
leaving `ShellLocalAngularProfileKey` exported, undocumented, and unreferenced
as the sole angular next-minor de-export candidate. That future namespace
change remains separate from documentation and does not reopen v0.2 release
work.

**MT7 - External Cartesian GTO interchange (completed/maintenance).** The
strict versioned reader, checkpoint-only PySCF exporter, frozen d-shell
fixture, and explicit determinant cleanup are implemented. The read-only C2
replay reproduced the accepted occupied subspace without a permutation ledger.
No solver, Hamiltonian payload, mandatory PySCF dependency, basis-only/live-mf
export, or release action belongs to this goal. Reader-facing documentation is
implemented and included in the immutable RC2 surface.

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

## Cartesian Hamiltonian Producer Pass 538 - Authorize Direct v0.2.0 Candidate

Commit(s):
- this docs-only authority amendment.

Summary:
- Added `HP-PQS-PUBLIC-V020-FN-01/TEST-01` for one direct final candidate from
  accepted post-RC2 `main`; no RC3 is required. The implementation may change
  only the package version, prepend a concise post-RC2 changelog section, move
  README installation/documentation links to the immutable final tag and
  `/stable/`, and update focused docs tests. RC2/RC1 history, source, exports,
  dependencies, examples, numerical policy, workflows, and `docs/make.jl`
  remain exact.
- The release claim is deliberately limited: v0.2.0 is the supported public
  package version closest to software used during the separate PQS and
  reference-density screening work, not an exact archive of either paper's
  computational history. Changelog scope is post-RC2 performance, fail-closed
  path-aware CI, and removal of duplicate release execution, with no manuscript
  or benchmark narrative.

Goal advancement / guardrail:
- MT5 advances to final-candidate preparation. This pass grants no tag or
  release. Candidate acceptance must freeze SHA, tree, archive, and an exact
  unpublished release body before a second docs-only pass may open one
  conditional tag-plus-publication transaction. The exact tag will use
  identity/install and independent Docs/Pages checks rather than rerunning the
  three numerical gates; a third docs-only pass closes the transaction.

Carrying-cost accounting:
- deleted: no code; simplified: no RC3 and one established three-pass release
  process; quarantined: registration, citation, paper metadata, and all new
  scientific/performance work; not deleted because: RC1/RC2 remain immutable
  release evidence.
- exact remaining caller/blocker: repo-manager must prepare and fully validate
  the candidate, then return exact freeze and release-body identities; added/
  deleted `src` lines `0/0`; new files none; new tests: focused existing docs
  owner only; new metadata/status fields none; validation: authority
  render/check/self-test, docs tests, Documenter, manager-log bound, scoped
  docs-only review, remote CI/Docs, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 539 - Accept Candidate And Authorize Final Transaction

Commit(s):
- `adfcaba32d4db06d9d796d947276433717bd2d89` - prepare the exact final
  candidate.
- this docs-only candidate closeout and conditional publication authority.

Summary:
- Accepted candidate tree `f64ba21e06ff57e2b5e78d91214398115afbe8de`
  with version `0.2.0` and only the four authorized file changes. RC2/RC1
  changelog history remains byte-identical. The clean archive has `677`
  entries, `10,137,600` bytes, and SHA-256
  `df09cc6fd7dc144daa168c9feb4a41be9b974ef450e1e81bf586787318ad1566`;
  it excludes a root manifest and both handoffs.
- Closed `HP-PQS-PUBLIC-V020-FN-01/TEST-01` to maintenance/completed and
  opened `HP-PQS-PUBLIC-V020-RELEASE-FN-01/TEST-01` for one ordered exact-hash
  transaction. The release body is frozen at `2,278` ASCII bytes including
  its final newline, SHA-256
  `e9ae9bcdad74b33bb66fb3e7e6a149d26285cb9bcc2f4c9555ac713be8bc90d2`.

Goal advancement / guardrail:
- MT5 advances from candidate preparation to final publication. The tag must
  target `adfcaba32` explicitly; remote annotated-object/install checks and
  independent Docs/Pages deployment must pass before the final/latest GitHub
  release is created. The numerical matrix is not rerun at the identical tag.
  Any partial tag or release is preserved and reported, never moved or edited.

Carrying-cost accounting:
- deleted: no code; simplified: candidate acceptance and publication opening
  share one docs-only review; quarantined: registration, citation, paper
  metadata, and all new source/performance work; not deleted because: RC1/RC2
  remain immutable evidence.
- exact remaining caller/blocker: repo-manager must execute and validate the
  ordered transaction; source delta `0/0`, new files/tests/metadata zero;
  candidate validation: Julia `1.10.12/1.12.6`, all three CI gates
  `33126022579`, Docs `33126022531`, examples 01/39/40/41, H2+ `18/18`,
  screening, export integrity, residual-GTO fixture, archive install,
  authority/docs/Documenter, and diff checks.

## Cartesian Hamiltonian Producer Pass 540 - Authorize Immutable-Tag Verification Recovery

Commit(s):
- this docs-only recovery authority amendment.

Summary:
- Accepted annotated tag object `722e8e8752a9d23f45e95d2f88e1749f9f3002e4`,
  message `GaussletBases v0.2.0`, peel `adfcaba32`, and tree `f64ba21e0` as
  exact. Fresh Julia `1.12.6` installation from the remote tag loaded version
  `0.2.0` with the frozen tree. Docs `33130411176` and Pages `33130489319`
  passed, with exact `/v0.2.0/`, real `/stable/`, intact prior folders, and the
  required selectors. The final GitHub release remains absent.
- Tag CI `33130411193` failed only because checkout created a local lightweight
  tag and the lane then fetched the annotated remote object to the same ref.
  Git correctly refused to clobber it. One fresh machine-local scratch lane may
  fetch the remote tag into a namespaced ref, prove annotated-object, message,
  peel, tree, version, API/remote, and clean-install identity, then recheck
  deployment before the already-frozen final release is published. No
  numerical rerun, tag mutation, workflow event, or workflow edit is part of
  this recovery.

Medium-term checkpoint:
- **MT1 active:** continue only bounded evidence-led conformance repairs.
- **MT2 completed:** controlled Cr2 source migration remains closed.
- **MT3 active:** represented-Hartree scaling, corrected-WL interpretation,
  and Standard60 remain separate pending work.
- **MT4 active:** residual/protected and consumer-owned PRF questions are
  unchanged.
- **MT5 active, refined:** the final candidate, tag, and documentation are
  accepted; publication is blocked only on the namespaced manual verification
  lane. Main-branch tag-lane repair is separate maintenance.
- **MT6 active:** accepted performance and test-efficiency work remains closed;
  this pass opens no optimization.
- **MT7 completed/maintenance:** external Cartesian-GTO interchange remains
  accepted in the final candidate.

Carrying-cost accounting:
- deleted: no code; simplified: one bounded manual lane replaces only the
  failed mechanical tag check; quarantined: future workflow repair,
  registration, citation, and all source/numerical work; not deleted because:
  immutable tag and release identity checks remain mandatory.
- exact remaining caller/blocker: repo-manager must pass the manual lane and
  publish/verify the exact release, or stop without mutation; source/test/
  workflow delta `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 541 - Close v0.2.0 Final Publication

Commit(s):
- immutable annotated tag object
  `722e8e8752a9d23f45e95d2f88e1749f9f3002e4`, peeling to candidate
  `adfcaba32d4db06d9d796d947276433717bd2d89` and tree `f64ba21e0`.
- this docs-only publication lifecycle closeout.

Summary:
- Accepted final GitHub release `378216554` with exact tag/title, final and
  latest status, zero uploaded assets, and the frozen `2,278`-byte body with
  SHA-256 `e9ae9bcdad74b33bb66fb3e7e6a149d26285cb9bcc2f4c9555ac713be8bc90d2`.
  Both automatic archives reconstruct the frozen tree, and a fresh Julia
  `1.12.6` installation from `rev = "v0.2.0"` loaded GaussletBases `0.2.0`.
- The manual namespaced-ref lane matched the tag object, message, peel, tree,
  version, `git ls-remote`, and GitHub API. `/v0.2.0/` and `/stable/` remain
  identical with the exact-version canonical URL; selectors retain stable,
  `v0.2`, RC2, RC1, and `dev`. No numerical gate was rerun.

Goal advancement / guardrail:
- MT5 closes final-v0.2 publication and returns both release records to
  completed/no-grant maintenance. The package claim remains modest and keeps
  PQS, screening, and interchange distinct. Do not mutate the final tag or
  release. Registration, citation, future tag-lane repair, and later releases
  require separate authority.

Carrying-cost accounting:
- deleted: the active publication grant; simplified: one completed release
  record now owns the durable evidence; quarantined: registration, citation,
  workflow repair, and later-release work; not deleted because: immutable
  release, archive, installation, and documentation evidence remains live.
- exact remaining caller/blocker: none for v0.2.0 publication; source/test/
  workflow delta `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 542 - Authorize Supported-Floor Radial And Misc Coverage

Commit(s):
- this docs-only authority amendment.

Summary:
- Reopened `HP-PUBLIC-PAPER-CI-FN-01/TEST-01` only to extend the existing Julia
  `1.10` Supported-floor selection from `core,ida,cartesian,examples` to
  `core,ida,cartesian,examples,radial,misc`. Clean evidence passed `radial`
  `322/322`, `misc` `59/59`, and combined `381/381` in about `70` seconds of
  testset time. Implementation owns one workflow-line replacement and one
  focused exact-selection policy assertion, with at most one angular-absence
  assertion.

Goal advancement / guardrail:
- MT5 improves compatibility-floor coverage without adding a CI row or changing
  any job name, Julia version, timeout, command, trigger, permission, path
  classifier, docs-only marker, PQS/Screening gate, or tag lane. `angular` is
  excluded: its first package-owned one-body fixture ran `13m49s` without an
  assertion and requires a separate fast-versus-acceptance audit. Weekly
  Cartesian, private occupied-first, blocked represented-Hartree, HFDMRG, docs,
  export, and release lanes remain separate.

Carrying-cost accounting:
- deleted: no code in this authority pass; simplified: existing live test
  groups join the established floor row; quarantined: angular and all unrelated
  owners; not deleted because: the three-row matrix and fail-closed routing are
  accepted release infrastructure.
- exact remaining caller/blocker: repo-manager must make only the two-file
  bounded replacement, run the expanded Julia `1.10` floor and unchanged remote
  matrix, then return for closeout; production-source delta `0/0`, workflow
  expected `+1/-1`, focused test additions preferred `1`, hard `2`, new files
  and metadata fields none.

## Cartesian Hamiltonian Producer Pass 543 - Authorize Post-v0.2 Nested/QW Export Cleanup

Commit(s):
- this docs-only authority amendment.

Summary:
- Reopened only `HP-PUBLIC-EXPORT-INTEGRITY-FN-01` for a two-file, seven-line
  source reduction. Implementation deletes the unused four-line
  `TimedNestedFixedBlockBuild` and removes three root exports. The internal
  `OneCenterAtomicNestedLayerStructure` type and one-line
  `QiuWhiteResidualGaussianOperators` alias remain defined for qualified use.
  The completed dynamic export-integrity test remains unchanged and supplies
  the acceptance owner.

Goal advancement / guardrail:
- MT6 removes unsupported public vocabulary without changing ordinary/QW/
  nested numerics. The deliberate alias retention avoids a needless qualified
  compatibility break. `ShellLocalAngularProfileKey` and the angular audit,
  all release objects, workflows, dependencies, and broader API classification
  remain separate.

Carrying-cost accounting:
- authorized for deletion: one orphaned wrapper and three export entries;
  simplified: the undocumented export backlog is expected to fall from `74`
  to `71`;
  quarantined: ignored historical probes and the later angular decision; not
  deleted because: the QW alias is one line and retains cheap qualified
  compatibility.
- exact remaining caller/blocker: repo-manager must confirm no committed
  caller, implement exact `+0/-7` scope, and return for lifecycle closeout;
  no test, file, helper, metadata, numerical, workflow, or release addition is
  authorized.

## Cartesian Hamiltonian Producer Pass 544 - Close Supported-Floor Radial And Misc Coverage

Commit(s):
- `15676153aec1569f5224ffa6ff5ed67b054c837f`
- this docs-only lifecycle closeout.

Summary:
- Accepted the exact Supported-floor extension from
  `core,ida,cartesian,examples` to
  `core,ida,cartesian,examples,radial,misc`. The implementation changed one
  workflow value and added one focused policy assertion, for `+2/-1` lines;
  source and all other workflow behavior remained unchanged.
- Julia `1.10` passed `radial` `322/322`, `misc` `59/59`, and their combined
  `381/381`. Remote CI run `33141930944` passed Supported floor, PQS, and
  Screening; Docs run `33141930932` passed. Authority/self-test, package load,
  docs `115/115` plus `10/10`, Documenter, and diff checks also passed.

Goal advancement / guardrail:
- MT5 returns `HP-PUBLIC-PAPER-CI-FN-01/TEST-01` to maintenance after improving
  compatibility-floor coverage. Preserve the three names/rows, path-aware
  routing, Julia versions, and tag lane. `angular` remains excluded pending its
  separate fast-versus-acceptance audit; no release or numerical policy follows.

Carrying-cost accounting:
- deleted: no owner or test; simplified: two live groups now run in the
  existing floor gate; quarantined: angular and unrelated internal suites; not
  deleted because: the three-gate public CI boundary remains live.
- exact remaining caller/blocker: none for Pass 542; production-source delta
  `0/0`, workflow/test delta `+2/-1`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 545 - Close Post-v0.2 Nested/QW Export Cleanup

Commit(s):
- `6e7bcbb7dae4e865dbdc0362b8f39ffd23f0a468`
- this docs-only lifecycle closeout.

Summary:
- Accepted exact `+0/-7` source reduction across the two authorized files.
  `TimedNestedFixedBlockBuild` is deleted. `OneCenterAtomicNestedLayerStructure`
  and `QiuWhiteResidualGaussianOperators` remain defined but unexported, with
  the one-line QW alias retained for qualified compatibility.
- The undocumented-export count fell from `74` to `71`. Complete core,
  export-integrity, and nested/QW regressions passed; docs passed `115/115`
  plus `10/10`. CI run `33143556046`, Docs run `33143556067`, authority
  check/self-test, package load, Documenter, and diff checks passed.

Goal advancement / guardrail:
- MT6 closes this evidence-led surface reduction without changing numerics or
  mutating v0.2.0. `ShellLocalAngularProfileKey` and the angular public-surface
  audit remain separate. No alias deletion, shim, warning, or replacement API
  follows from this pass.

Medium-term checkpoint:
- **MT1 active:** continue only narrow evidence-led conformance repairs.
- **MT2 completed:** controlled Cr2 source migration remains closed.
- **MT3 active:** represented-Hartree scaling, corrected-WL interpretation,
  and Standard60/canonical-driver exposure remain separate pending work.
- **MT4 active:** residual/protected and consumer-owned PRF questions are
  unchanged.
- **MT5 active/maintenance:** final v0.2.0 and expanded Supported-floor
  coverage are accepted; registration, citation, and future tag-lane repair
  remain separate.
- **MT6 active:** this export cleanup is closed; angular classification and any
  further carrying-cost work require their own evidence and authority.
- **MT7 completed/maintenance:** external Cartesian-GTO interchange remains
  accepted and unchanged.

Carrying-cost accounting:
- deleted: one orphaned wrapper and three export entries; simplified: the root
  surface and undocumented-export backlog; quarantined: angular and ignored
  historical probes; not deleted because: the one-line QW alias cheaply
  preserves qualified compatibility.
- exact remaining caller/blocker: none for Pass 543; source delta `+0/-7`, test
  delta `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 546 - Authorize Foundational Basis And Mapping Documentation

Commit(s):
- this docs-only authority amendment.

Summary:
- Added `HP-PUBLIC-FOUNDATION-DOC-FN-01/TEST-01` for ten existing exported
  basis/mapping generics. Implementation may add only declaration docstrings,
  one compact section in the existing bases/mappings reference, and focused
  checks in the existing docs owner. The map direction, derivatives, physical
  versus reference centers, and basis-integral meaning of `integral_weights`
  are explicit.

Goal advancement / guardrail:
- MT6 addresses a bounded reader-facing documentation deficit without adding
  API or preserving unsupported vocabulary. Acceptance requires the
  undocumented-export backlog to fall exactly from `71` to `61`, but forbids a
  global allowlist or `checkdocs` change. Function evaluation, stencils,
  partitions, atomic operators, angular work, and further API reduction remain
  separate.

Carrying-cost accounting:
- deleted: no code in this authority pass; simplified: ten existing public
  names receive one coherent reference family; quarantined: the remaining `61`
  undocumented exports and all unrelated API classification; not deleted
  because: these live bindings already have methods and supported semantics.
- exact remaining caller/blocker: repo-manager must remain within source-doc
  `70/100`, reference-page `20/35`, and docs-test `12/20` preferred/hard added
  line budgets, make no new file or behavioral change, run package/radial/docs/
  authority/Documenter/CI acceptance, and return for lifecycle closeout.

## Cartesian Hamiltonian Producer Pass 547 - Close Foundational Basis And Mapping Documentation

Commit(s):
- `03eb75c66e5cc8b24714d3769097ea81f0d74b15`
- this docs-only lifecycle closeout.

Summary:
- Accepted concise declaration docstrings and one compact reference family for
  the ten foundational mapping/basis generics. The exact implementation delta
  was `+58` source-docstring lines, `+19` reference lines, and `+8` docs-test
  lines, all within preferred limits. No definition, method, dispatch, export,
  return value, API, or numerical behavior changed.
- The undocumented exported-binding backlog fell exactly from `71` to `61`.
  Focused radial tests passed `322/322`; docs passed `118/118` plus `10/10`.
  Package load, authority check/self-test, Documenter, and diff checks passed;
  remote CI run `33181875262` and Docs run `33181875247` passed.

Goal advancement / guardrail:
- MT6 returns `HP-PUBLIC-FOUNDATION-DOC-FN-01/TEST-01` to maintenance after
  closing this reader-facing documentation deficit. Preserve the physical-x to
  reference-u convention, derivative meanings, physical/reference center
  distinction, and basis-integral meaning of `integral_weights`. The remaining
  `61` undocumented exports, angular work, API reduction, and global docs policy
  remain separate. No checkpoint is due.

Carrying-cost accounting:
- deleted: no implementation; simplified: ten existing public bindings now
  form one coherent documented family; quarantined: the remaining export
  backlog and unrelated classification; not deleted because: all ten bindings
  remain live supported interfaces.
- exact remaining caller/blocker: none for Pass 546; source-behavior delta
  `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 548 - Authorize Function Evaluation And Stencil Documentation

Commit(s):
- this docs-only authority amendment.

Summary:
- Added `HP-PUBLIC-FUNCTION-STENCIL-DOC-FN-01/TEST-01` for nine existing
  exported function-evaluation and stencil generics. Implementation may add
  only declaration docstrings, one compact section in the existing
  bases/mappings reference, and focused checks in the existing docs owner.
  Normal versus represented evaluation, derivative order, physical versus
  reference centers, per-function integration, and ordered stencil data are
  explicit.

Goal advancement / guardrail:
- MT6 addresses the next bounded reader-facing documentation deficit without
  adding API or preserving unsupported vocabulary. Acceptance requires the
  undocumented-export backlog to fall exactly from `61` to `52`, but forbids a
  global allowlist or `checkdocs` change. Basis metadata, partitions,
  operators, angular work, and further API reduction remain separate. No
  checkpoint is due.

Carrying-cost accounting:
- deleted: no code in this authority pass; simplified: nine existing public
  names receive one coherent function/stencil reference family; quarantined:
  the remaining `52` undocumented exports and unrelated API classification;
  not deleted because: these live bindings already have supported methods and
  semantics.
- exact remaining caller/blocker: repo-manager must remain within source-doc
  `60/90`, reference-page `18/30`, and docs-test `10/16` preferred/hard added
  line budgets, make no new file or behavioral change, run package/core/radial/
  docs/authority/Documenter/CI acceptance, and return for lifecycle closeout.

## Cartesian Hamiltonian Producer Pass 549 - Close Function Evaluation And Stencil Documentation

Commit(s):
- `624718972627b5ebf74beaa57ad73b8725aac587`
- this docs-only lifecycle closeout.

Summary:
- Accepted concise declaration docstrings and one compact reference family for
  the nine function-evaluation and stencil generics. The exact implementation
  delta was `+57` source-docstring lines, `+18` reference lines, and `+8`
  docs-test lines, all within preferred limits. No definition, method,
  dispatch, export, allocation policy, API, or numerical behavior changed.
- The undocumented exported-binding backlog fell exactly from `61` to `52`.
  Focused core/radial suites and radial `322/322` passed; docs passed `121/121`
  plus `10/10`. Package load, authority check/self-test, Documenter, and diff
  checks passed; remote CI run `33213490510` and Docs run `33213490543` passed.

Goal advancement / guardrail:
- MT6 returns `HP-PUBLIC-FUNCTION-STENCIL-DOC-FN-01/TEST-01` to maintenance
  after closing this reader-facing documentation deficit. Preserve normal
  evaluation versus represented-stencil evaluation, physical derivative and
  center meanings, per-function integration, and ordered stencil data. The
  remaining `52` undocumented exports, angular work, API reduction, and global
  docs policy remain separate. No checkpoint is due.

Carrying-cost accounting:
- deleted: no implementation; simplified: nine existing public bindings now
  form one coherent documented family; quarantined: the remaining export
  backlog and unrelated classification; not deleted because: all nine bindings
  remain live supported interfaces.
- exact remaining caller/blocker: none for Pass 548; source-behavior delta
  `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 550 - Authorize Partition And Leaf-Local Documentation

Commit(s):
- this docs-only authority amendment.

Summary:
- Added `HP-PUBLIC-PARTITION-LEAF-DOC-FN-01/TEST-01` for twelve existing
  exported partition-hierarchy and leaf-local accessors. Implementation may
  add only declaration docstrings, one compact section in the existing
  bases/mappings reference, and focused checks in the existing docs owner.
- The contract distinguishes copied `Float64` box blocks from stored ordered
  vectors and records that callers must treat as read-only. It also gathers
  the already-documented partition, hierarchy, leaf-PGDG, mapped-layer, and
  contraction context without authorizing new prose or behavior.

Goal advancement / guardrail:
- MT6 addresses the next bounded reader-facing documentation deficit without
  changing array ownership, indexing, API, or numerics. Acceptance requires
  the established undocumented-export backlog, excluding the module self-
  binding, to fall exactly from `52` to `40`; global allowlists and `checkdocs`
  changes remain forbidden.

Medium-term checkpoint:
- **MT1 active:** continue only narrow evidence-led conformance repairs.
- **MT2 completed:** controlled Cr2 source migration remains closed.
- **MT3 active:** represented-Hartree scaling, corrected-WL interpretation,
  and Standard60/canonical-driver exposure remain separate pending work.
- **MT4 active:** residual/protected and consumer-owned PRF questions are
  unchanged.
- **MT5 active/maintenance:** final v0.2.0 and expanded Supported-floor
  coverage remain accepted; registration, citation, and tag-lane repair are
  separate.
- **MT6 active, refined:** foundational and function/stencil documentation are
  closed; this bounded partition/leaf family is approved next. The remaining
  `40`-name backlog, angular classification, and API reduction remain separate.
- **MT7 completed/maintenance:** external Cartesian-GTO interchange remains
  accepted and unchanged.

Carrying-cost accounting:
- deleted: no code in this authority pass; simplified: twelve related accessors
  and their existing context receive one coherent reference family;
  quarantined: the remaining backlog and unrelated API classification; not
  deleted because: these bindings remain live public partition/local-layer
  interfaces.
- exact remaining caller/blocker: repo-manager must remain within source-doc
  `80/115`, reference-page `35/55`, and docs-test `10/18` preferred/hard added
  line budgets, make no new file or behavioral change, run package/core/docs/
  authority/Documenter/CI acceptance, and return for lifecycle closeout.

## Cartesian Hamiltonian Producer Pass 551 - Close Partition And Leaf-Local Documentation

Commit(s):
- `b24c0dbeeb0a6b5526fde9004cbc7319c1395ec4`
- this docs-only lifecycle closeout.

Summary:
- Accepted docstrings for all twelve partition-hierarchy and leaf-local
  accessors plus their compact reference section and family-scoped checks.
  Fresh `Float64` box-block copies remain distinct from stored vectors and
  records that callers treat as read-only. No definition, indexing, ownership,
  dispatch, export, or numerical behavior changed.
- Exact additions were `80` source-docstring, `35` reference, and `15`
  docs-test lines. These meet the preferred source/reference limits and the
  hard `18`-line test limit, with no new file. The undocumented-export backlog,
  excluding the module self-binding, fell exactly from `52` to `40`.

Validation / goal advancement:
- Focused core and docs `123/123` plus `10/10` passed. Package load, authority
  check/self-test, Documenter, and diff checks passed; remote CI run
  `33219963234` and Docs run `33219963250` passed.
- MT6 returns `HP-PUBLIC-PARTITION-LEAF-DOC-FN-01/TEST-01` to maintenance.
  Preserve ordering, hierarchy, ownership, and copy semantics. The remaining
  `40` undocumented exports, angular classification, API reduction, and global
  docs policy remain separate. No checkpoint is due.

Carrying-cost accounting:
- deleted: no implementation; simplified: twelve related public accessors now
  form one documented family; quarantined: the remaining export backlog and
  unrelated classification; not deleted because: all twelve bindings remain
  live supported interfaces.
- exact remaining caller/blocker: none for Pass 550; source-behavior delta
  `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 552 - Authorize Atomic IDA Inspection Documentation

Commit(s):
- this docs-only authority amendment.

Summary:
- Added `HP-PUBLIC-ATOMIC-IDA-DOC-FN-01/TEST-01` for ten existing exported
  atomic-IDA inspection and tiny two-electron reference functions. The grant
  permits only declaration docstrings, one compact section in the existing
  atomic/ordinary reference page, and focused checks in the existing docs
  owner. Existing carrier and builder docstrings provide the surrounding
  context without a new page or API.
- The contract distinguishes owner-specific stored ordering from a universal
  order, on-demand dense Gaunt/kernel reconstruction from retained dense
  storage, and the tiny dense one-up/one-down reference problem from a general
  solver. The small fully reorthogonalized Lanczos routine must document its
  controls and exact return fields without a production-eigensolver claim.

Goal advancement / guardrail:
- MT6 addresses one final coherent documentation family under the existing
  public surface. Acceptance requires the undocumented-export backlog to fall
  exactly from `40` to `30`, then stops for explicit classification rather
  than automatically treating the remainder as documentation debt. Source
  behavior, ordinary tests, angular research, storage/caching, generalized
  solvers, global docs policy, and API reduction remain separate. No checkpoint
  is due.

Carrying-cost accounting:
- deleted: no code in this authority pass; simplified: ten related bindings
  and six already-documented context names receive one bounded reference
  family; quarantined: the remaining `30` exports pending classification; not
  deleted because: the ten functions are live public inspection/reference
  interfaces.
- exact remaining caller/blocker: repo-manager must remain within source-doc
  `75/110`, reference-page `28/45`, and docs-test `10/18` preferred/hard added
  line budgets, add no file or behavioral change, run core/IDA/docs/authority/
  Documenter/CI acceptance, and return for lifecycle closeout.

## Cartesian Hamiltonian Producer Pass 553 - Close Atomic IDA Inspection Documentation

Commit(s):
- `8f84afb491ff7cc0f818a6ed982eca638a509b65`
- this docs-only lifecycle closeout.

Summary:
- Accepted docstrings for the ten atomic-IDA inspection and tiny two-electron
  reference functions, their compact atomic/ordinary reference section, and
  focused documentation checks. Stored records remain read-only and owner-
  ordered; dense angular tensors remain on-demand references; and the fully
  reorthogonalized Lanczos path remains a tiny one-up/one-down reference
  routine rather than a production solver.
- Exact additions were `70` source-docstring, `26` reference, and `15` docs-
  test lines, all within hard limits, with no new file or executable change.
  The undocumented-export count fell exactly from `40` to `30`.

Validation / goal advancement:
- The focused core group and docs `126/126` plus `10/10` passed. Package load,
  authority check/self-test, Documenter, manager-log bound, and diff checks
  passed; remote CI run `33230195333` and Docs run `33230195291` passed.
- MT6 returns `HP-PUBLIC-ATOMIC-IDA-DOC-FN-01/TEST-01` to maintenance and
  stops the automatic documentation sequence. The remaining names are now
  classified as six supported-public documentation candidates, nineteen
  expert/experimental bindings requiring accurate labels, and five future
  de-export-audit candidates; none is deletion-ready. No checkpoint is due.

Carrying-cost accounting:
- deleted: no implementation; simplified: ten related public functions now
  form one coherent documented reference family; quarantined: the remaining
  `30` names pending their separately approved classification work; not
  deleted because: every in-scope function remains a live inspection/reference
  interface.
- exact remaining caller/blocker: none for Pass 552; source-behavior delta
  `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 554 - Authorize Supported Public Surface Documentation

Commit(s):
- this docs-only authority amendment.

Summary:
- Added `HP-PUBLIC-SUPPORTED-SURFACE-DOC-FN-01/TEST-01` for the six
  supported-public names isolated by Pass 553. The grant moves the
  `BondAlignedDiatomicQWBasis3D` docstring from its unexported abstract parent,
  documents the finite callable Coulomb Gaussian representation and owner-
  specific metadata, keeps `cartesian_base_hamiltonian` limited to its stated
  H/H2 routes, and defines both external-GTO hashes as strict packet-integrity
  fingerprints rather than numerical-equivalence tests.
- Implementation is limited to five existing source docstring locations,
  three existing reference pages, and the existing docs-test family. Hard
  limits are `95` source-docstring, `35` reader-doc, and `18` docs-test added
  lines, with no new file or executable behavior.

Goal advancement / guardrail:
- MT6 authorizes the final supported-public documentation packet. Acceptance
  requires the undocumented-export count to fall exactly from `30` to `24`,
  then stop. The nineteen expert/experimental bindings and five future de-
  export candidates remain unchanged and require separate classification or
  audit authority. No checkpoint is due.

Carrying-cost accounting:
- deleted: no code in this authority pass; simplified: one misattached
  docstring must move rather than duplicate; quarantined: the remaining `24`
  classified names; not deleted because: all six in-scope names are live
  supported bindings.
- exact remaining caller/blocker: repo-manager must stay within the stated
  source/docs/test budgets, preserve every implementation and hash algorithm,
  run package/docs/authority/Documenter/CI acceptance, and return for
  lifecycle closeout.

## Cartesian Hamiltonian Producer Pass 555 - Close Supported Public Surface Documentation

Commit(s):
- `8161f131aa962fef979f8ef09c14d23231eb14e4`
- this docs-only lifecycle closeout.

Summary:
- Accepted accurate docstrings and curated placement for the final six
  classified supported-public bindings. The QW documentation moved from the
  unexported abstract parent to `BondAlignedDiatomicQWBasis3D` without
  duplication; the Coulomb expansion, owner-specific metadata, narrow H/H2
  facade, and strict non-invariant external-GTO fingerprints retain their
  authorized limits.
- Exact changes were `+52/-10` source-docstring lines, `+20` reader-doc lines,
  and `+18` focused docs-test lines, with no new file or executable behavior.
  All hard limits were met; the test addition used its hard ceiling exactly.
  The undocumented-export backlog fell exactly from `30` to `24`.

Validation / goal advancement:
- Docs passed `134/134` plus `10/10`; package load, authority check/self-test,
  Documenter, and diff checks passed. Remote CI run `33235720583` ran and
  passed Supported floor, PQS paper, and Screening paper; Docs run
  `33235720615` passed.
- MT6 returns `HP-PUBLIC-SUPPORTED-SURFACE-DOC-FN-01/TEST-01` to maintenance.
  The automatic supported-public documentation sequence is complete. The
  remaining nineteen expert/experimental bindings and five future de-export
  candidates remain classified but untouched; neither documentation nor
  removal follows without separate authority.

Medium-term checkpoint:
- **MT1 active:** continue only narrow evidence-led conformance repairs.
- **MT2 completed:** controlled Cr2 source migration remains closed.
- **MT3 active:** represented-Hartree scaling, corrected-WL interpretation,
  and Standard60/canonical-driver exposure remain separate pending work.
- **MT4 active:** residual/protected and consumer-owned PRF questions are
  unchanged.
- **MT5 active/maintenance:** final v0.2.0, documentation deployment, and the
  expanded Supported-floor gate remain accepted; registration, citation, and
  future tag-lane repair remain separate.
- **MT6 active, refined:** five bounded public-documentation families are
  closed at an undocumented-export count of `24`; expert labeling, de-export
  audits, angular classification, and further carrying-cost work each require
  separate evidence and authority.
- **MT7 completed/maintenance:** external Cartesian-GTO interchange remains
  accepted and unchanged.

Carrying-cost accounting:
- deleted: the misplaced abstract-type docstring; simplified: six live
  bindings now complete the supported-public reference set; quarantined: the
  remaining nineteen expert/experimental and five audit candidates; not
  deleted because: all six in-scope bindings remain supported public APIs.
- exact remaining caller/blocker: none for Pass 554; source-behavior delta
  `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 556 - Classify Remaining Surface And Authorize Expert Geometry Documentation

Commit(s):
- audit snapshot `0832bdce0f3ff9bdd07d8580bedcda0c47dadbe8`;
- this docs-only authority amendment.

Summary:
- Accepted the read-only de-export audit without granting an export change.
  `bond_aligned_diatomic_geometry_payload` remains an exported
  expert/experimental interface because it is the central constructor for the
  coherent bond-aligned geometry family with active high-order consumers.
- Reserved two inseparable next-minor namespace transactions: the three one-
  center nested diagnostic names, and `diagnose_qwrg_residual_space` together
  with its already-documented `QWRGResidualSpaceDiagnostics` type. Every
  implementation remains for qualified internal/research use. These
  reservations are namespace reductions, not deletion, performance, v0.2.x,
  release, or automatic v0.3 work.
- Added `HP-PUBLIC-EXPERT-GEOMETRY-DOC-FN-01/TEST-01` for one declaration
  docstring, one compact expert/experimental section in the existing bases and
  mappings reference, and focused existing-owner checks. Hard limits are `36`
  source-docstring, `28` reader-doc, and `14` docs-test added lines, with no
  new file or behavioral change.

Goal advancement / guardrail:
- MT6 refines the remaining public-surface classification instead of treating
  it as undifferentiated documentation debt. Acceptance must reduce the
  undocumented-export count exactly from `24` to `23`, leaving nineteen
  undocumented expert/experimental bindings and four original de-export
  candidates; the associated documented QWRG type remains coupled to the
  future residual transaction. No checkpoint is due.

Carrying-cost accounting:
- deleted: nothing; simplified: five formerly ambiguous candidates now have
  explicit retain/reserve dispositions; quarantined: both future next-minor
  namespace transactions; not deleted because: implementations remain active
  internal/research machinery and geometry remains a coherent expert family.
- exact remaining caller/blocker: repo-manager must remain within the stated
  docs-only budgets, preserve exports and behavior, run package/docs/authority/
  Documenter/CI acceptance, and return for lifecycle closeout.

## Cartesian Hamiltonian Producer Pass 557 - Close Expert Geometry Documentation

Commit(s):
- `312f69a4b757e8d559387e0ef14450f031691039`;
- this docs-only lifecycle closeout.

Summary:
- Accepted the expert/experimental documentation for
  `bond_aligned_diatomic_geometry_payload` and its seven already-documented
  companion records/operations. The text preserves the inspection-only,
  compressed-versus-source, same-geometry, read-only-data, and non-general-API
  boundaries without changing exports or behavior.
- Exact additions were `18` source-docstring, `21` reader-documentation, and
  `13` focused docs-test lines, all within preferred limits and with no new
  file. The undocumented-export count fell exactly from `24` to `23`.

Validation / goal advancement:
- Docs passed `138/138` plus `10/10`; radial passed `322/322`; package load,
  authority check/self-test, Documenter, and diff checks passed. Remote CI
  run `33326473189` passed all three gates, Docs run `33326473214` passed, and
  the live `/dev/` page was verified at the implementation commit.
- MT6 returns `HP-PUBLIC-EXPERT-GEOMETRY-DOC-FN-01/TEST-01` to maintenance.
  The remaining state is nineteen undocumented expert/experimental bindings
  plus four original next-minor de-export candidates. The two reserved
  namespace transactions remain non-executable and do not imply source
  deletion, v0.2.x compatibility work, or automatic v0.3 development. No
  checkpoint is due.

Carrying-cost accounting:
- deleted: nothing; simplified: one coherent expert geometry family now has
  one discoverable reader boundary; quarantined: both future namespace
  transactions; not deleted because: geometry has active high-order consumers
  and all reserved implementations remain useful for qualified research use.
- exact remaining caller/blocker: none for Pass 556; executable source delta
  `0/0`, new files and metadata fields none.

## Cartesian Hamiltonian Producer Pass 558 - Authorize Angular Optional-Consumer Control-Flow Repair

Commit(s):
- baseline `5f1708d989269562d9182ed630336104c258da2e`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-ANGULAR-TEST-CONTROL-FN-01/TEST-01` for one test-runner-only
  correctness repair. A missing optional HFDMRG checkout may visibly skip only
  its handshake; an existing checkout that fails to load/import must fail; the
  angular file may not return early; and explicit `angular` selection must run
  the complete owner after removal from the misleading `fast` set.
- MT1 advances a narrow fail-closed conformance repair. The implementation is
  limited to two existing runner files, about `+6/-4` lines with `12` added
  control-flow lines hard, plus at most `12` focused lines in the existing docs
  test. The complete roughly 16-minute angular owner is required once because
  it is the only acceptance proof that all tests after the former return run.

Guardrail / carrying cost:
- deleted: the early return, catch-all absence conflation, and `angular` fast
  membership; simplified: optional-consumer absence becomes local and visible;
  quarantined: fixed-radial extraction, seven angular exports, shell-key
  de-export, and CI ownership; not deleted because: all numerical angular tests
  and the explicit group remain live.
- exact remaining blocker: repo-manager implementation and acceptance evidence;
  source/API/workflow delta must remain `0/0`, with no new file or framework.

## Cartesian Hamiltonian Producer Pass 559 - Close Angular Optional-Consumer Control-Flow Repair

Commit(s):
- `22d25e8c5061abafa5a37ea35d4a51adfa4b9a72`;
- this docs-only lifecycle closeout.

Summary / validation:
- Accepted the two-file `+6/-9` repair. Missing HFDMRG now yields one visible
  handshake skip, while an existing checkout loads without a catch-all and
  therefore propagates failures. The file-level return is gone, explicit
  `angular` selection remains, and only `angular` membership in `fast` was
  removed. No numerical assertion, source, API, workflow, helper, file, or
  docs-test changed.
- The absent-checkout owner completed `61,907` passes and one skip in `929.68`
  seconds, proving later tests ran. Package load, docs `138/138` plus `10/10`,
  authority check/self-test, Documenter, and diff checks passed; remote CI
  `33475472072` passed all three numerical gates and Docs `33475471987` passed.

Goal / carrying cost:
- MT1 returns `HP-ANGULAR-TEST-CONTROL-FN-01` to maintenance and completes the
  optional TEST grant with no committed regression. The runner now fails closed
  without assigning the research suite to fast CI. No strategic change and no
  checkpoint is due.
- deleted: early return, catch-all, and false fast membership; simplified:
  absence is local; quarantined: all broader angular ownership/API work; exact
  remaining blocker: none for Pass 558.

## Cartesian Hamiltonian Producer Pass 560 - Authorize Fixed-Radial Angular Public CI Extraction

Commit(s):
- baseline `5d9befa246f38893449d28242ca242595a8f2695`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-ANGULAR-PUBLIC-CI-FN-01/TEST-01` for one separate, bounded
  extraction of the existing `[10,15,32]` fixed-radial angular sequence into a
  compact `angular_public` owner and the Julia `1.10` Supported-floor gate.
  The implementation must move and condense the old block, preserve its
  representative construction, profile, overlap, dense-payload, and one
  serialization identity contract, and reduce total tracked test LOC.
- `angular_public` becomes explicitly selectable and available through `all`,
  but both it and the complete roughly 16-minute `angular` research group stay
  out of `fast`. Only `angular_public` may be appended to the existing
  Supported-floor group list; every job identity, other selection, workflow
  behavior, and public gate remains fixed.
- MT1 advances public-contract coverage after the fail-closed angular runner
  repair. MT5 gains one bounded compatibility-floor owner without assigning
  the full angular suite to per-push CI. MT6 requires a net test reduction:
  preferred `60-80`, hard `95` owner lines, deletion of the whole superseded
  block, and at most `6` docs-policy lines.

Validation / guardrail:
- Implementation acceptance requires a Julia `1.10` focused run with runtime,
  one combined `angular,angular_public` run on one supported Julia version
  only, visible Supported-floor execution, all three unchanged CI gates, and
  package/docs/authority/diff checks. That expensive combined run must not be
  repeated for closeout. Angular export documentation,
  `ShellLocalAngularProfileKey` de-export, source, numerical policy, and release
  work remain separate.
- Rotated Passes 517-537 verbatim into the `926`-line archive with SHA-256
  `c0def005ed4292d5cbbc684d70298f25d9a9f03b3ea777bd2001cbfec6dc517a`;
  the live ledger now starts at Pass 538 and remains within its bound.

Carrying-cost accounting:
- deleted: no implementation in this authority pass; simplified: one focused
  endpoint will replace a repeated research-suite block; quarantined: full
  angular CI and angular API work; not deleted because: the remaining angular
  research owner stays directly selectable.
- exact remaining blocker: repo-manager implementation and the bounded
  acceptance evidence above; added/deleted source lines `0/0`, new metadata
  fields none, and exactly one planned test owner.

### Medium-Term Goal Checkpoint After Pass 560

- **MT1 - active:** fail-closed angular control flow is closed; the bounded
  fixed-radial public-owner extraction is approved and awaiting implementation.
- **MT2 - completed:** controlled Cr2 source migration remains closed.
- **MT3 - active:** pending producer facilities and corrected paper-oracle
  interpretation are unchanged.
- **MT4 - active:** residual/protected and consumer-owned PRF questions are
  unchanged.
- **MT5 - active/maintenance:** final v0.2 release maintenance and the three
  public CI identities remain fixed; only the bounded `angular_public` append
  is pending, while full `angular` remains excluded.
- **MT6 - active:** evidence-led test/API carrying-cost reduction continues;
  this extraction must be net-negative and adds no implementation framework.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange and its
  reader front door remain accepted and unchanged.

## Cartesian Hamiltonian Producer Pass 561 - Close Fixed-Radial Angular Public CI Extraction

Commit(s):
- `db156966a4a3a7bf2f685fa0f89312afca7b4280`;
- this docs-only lifecycle closeout.

Summary:
- Accepted the focused `[10,15,32]` fixed-radial sequence owner and its Julia
  `1.10` Supported-floor wiring. The `80`-line, `83`-assertion public owner
  replaced the complete `163`-line research-suite block; total tracked test
  delta was `+88/-165`, net `-77`. It remains explicitly selectable and
  included by `all`, while both `angular_public` and the complete `angular`
  research group remain outside `fast`; bare `angular` is not in per-push CI.

Validation / goal advancement:
- The focused Julia `1.10` owner passed `83/83` in `11.17` seconds. The one
  authorized combined Julia `1.12` `angular,angular_public` run passed in
  `930.34` seconds and was accepted without repetition during closeout. Remote
  CI run `33509253422` passed all three gates and visibly ran
  `angular_public` `83/83`; the initial Docs run `33509253439` failed only on
  the deliberately planned owner-path state.
- MT1 returns `HP-ANGULAR-PUBLIC-CI-FN-01` to implemented maintenance and its
  TEST record to completed maintenance. MT5 now includes the accepted bounded
  angular owner in compatibility-floor coverage. MT6 records the net test
  reduction without adding a framework. No checkpoint is due.

Carrying-cost accounting:
- deleted: the superseded `163`-line fixed-radial block; simplified: one
  table-driven owner and one existing CI selection; quarantined: the complete
  angular research suite remains direct-run only; not deleted because: its
  remaining research contracts are outside this public endpoint.
- exact remaining caller/blocker: none after the docs-only planned-to-existing
  reconciliation and passing Docs rerun; added/deleted production source lines
  `0/0`, new metadata fields none, and no numerical policy or release change.

## Cartesian Hamiltonian Producer Pass 562 - Authorize Experimental Angular Producer Documentation

Commit(s):
- baseline `628a46599dc7d1a85458a22a78fe095b7783b53a`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-PUBLIC-ANGULAR-PRODUCER-DOC-FN-01/TEST-01` for exactly seven
  durable shell-profile and fixed-radial-sequence bindings. The reader contract
  labels them experimental producer surfaces, documents exact-injected versus
  mixed-complement ordering, deterministic labels/gauge data, provenance-only
  identities, shell-independent source-target overlaps, adjacent versus direct
  upper-triangle sidecars, and one shared radial substrate across levels.
- MT6 advances the classified angular reader boundary without changing API or
  numerics. Final-basis orthonormality remains the working contract; sidecars
  are continuation/transfer data, not generalized-overlap authority.

Guardrail / carrying cost:
- Preferred/hard additions are `70/95` source-docstring, `30/45` combined
  reader/reference, and `12/18` docs-test lines, with no new file. Acceptance
  requires seven exported documented reference entries, exact undocumented-
  export reduction `23 -> 16`, and `angular_public` `83/83`; the full roughly
  `930`-second angular suite must not be repeated.
- deleted: nothing in this authority pass; simplified: seven related names get
  one coherent reader family; quarantined: restart orchestration, target lifts,
  Givens transforms, campaigns, and key de-export; not deleted because:
  `ShellLocalAngularProfileKey` remains live cache machinery.
- exact remaining caller/blocker: repo-manager implementation within the five
  listed files and budgets. The key must stay exported but undocumented and
  absent from the reference; any required behavioral/export/cache change is a
  stop-and-report outcome. No checkpoint is due.

## Cartesian Hamiltonian Producer Pass 563 - Close Experimental Angular Producer Documentation

Commit(s):
- `9b410ea14356cee6085ad2de51e73849695a0d97`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted documentation for all seven durable experimental shell-profile and
  fixed-radial-sequence producer bindings. The implementation added accurate
  docstrings, a curated export section, and a compact research-track reader
  explanation while preserving the producer-only boundary, exact/mixed mode
  ordering, deterministic labels and gauge data, provenance-only identities,
  shell-independent overlaps, adjacent/direct sidecars, and common radial
  substrate semantics.
- MT6 advances the classified angular reader surface to maintenance. The
  undocumented-export backlog fell exactly `23 -> 16` without an API or
  numerical change. `ShellLocalAngularProfileKey` remains exported,
  undocumented, and unlisted; its possible next-minor de-export stays separate.

Validation / carrying cost:
- Implementation commit `9b410ea14` passed `angular_public` `83/83` in `11.9`
  seconds, docs `143/143` and `10/10`, package load, authority check/self-test,
  Documenter, log bound, and diff checks. Remote CI `33520726731` passed all
  three gates; Docs `33520726763` passed and deployed. The full roughly
  `930`-second angular suite was intentionally not repeated.
- added: `58` source-docstring lines; changed: reader/reference `+42/-10` and
  docs tests `+18/-1`; deleted: no implementation; simplified: seven related
  exports now form one coherent documented family; quarantined: restart
  orchestration, target lifts, Givens transformations, campaigns, and key
  de-export; new files and metadata fields none.
- exact remaining blocker: none for this packet. Both IDs return to maintenance
  only; no workflow, cache, dependency, version, release, or strategic-goal
  change occurred. No checkpoint is due.

## Cartesian Hamiltonian Producer Pass 564 - Authorize Radial Paper-Parity Documentation

Commit(s):
- baseline `b845dea4cde4b58a9b4e8d8b4fc663c709402571`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-PUBLIC-RADIAL-PARITY-DOC-FN-01/TEST-01` for the four coherent
  `RadialBoundaryPrototype` paper-parity exports. The contract separates this
  frozen expert reproducibility route from the recommended `RadialBasisSpec`
  front door, requires read-only treatment of shared cached records, and
  distinguishes the exact prototype widths from rounded workflow presets.
- MT6 advances the classified expert reader surface without changing the
  radial implementation, artifact, cache, API, tests, or Supported-floor
  selection. Acceptance requires the exact undocumented-export reduction
  `16 -> 12`; every other expert family and next-minor de-export candidate stays
  outside this packet.

Guardrail / carrying cost:
- Preferred/hard additions are `50/70` source-docstring, `30/45` combined
  reader/reference, and `10/16` docs-test lines, with no new file. Existing
  radial tests and the unchanged Supported-floor gate remain the numerical
  acceptance owners; no new assertion or artifact test is authorized.
- deleted: nothing in this authority pass; simplified: four related exports
  receive one coherent reproducibility boundary; quarantined: all other expert
  surfaces and namespace reductions; not deleted because: the frozen prototype
  remains live, artifact-backed reproducibility infrastructure.
- exact remaining blocker: repo-manager implementation within the four listed
  files and budgets. Any required coefficient, artifact, cache, API, workflow,
  or numerical change is a stop-and-report outcome. No checkpoint is due.

## Cartesian Hamiltonian Producer Pass 565 - Close Radial Paper-Parity Documentation

Commit(s):
- `84ebdc5792c954b4c09f08b53ba86ef9332bf6b8`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted accurate documentation for all four radial paper-parity exports.
  The reader boundary keeps `RadialBasisSpec` as the ordinary front door,
  identifies the sole frozen prototype, distinguishes exact prototype widths
  from rounded workflow presets, and treats cached records and arrays as
  read-only inspection data rather than mutable construction state.
- MT6 advances the classified expert reader surface to maintenance. The
  undocumented-export backlog fell exactly `16 -> 12` without changing radial
  behavior, coefficients, artifacts, defaults, exports, or cache semantics.

Validation / carrying cost:
- Implementation commit `84ebdc579` passed the radial owner `322/322`, docs
  `148/148` and `10/10`, package load, authority check/self-test, Documenter,
  log bound, and diff checks. Remote CI `33527681191` passed all three gates;
  Docs `33527681116` passed, deployed, and the live reference was verified.
- added: `33` source-docstring and `32` reader-documentation lines; changed:
  docs tests `+14/-1`; deleted: no implementation; simplified: four related
  exports now have one explicit expert reproducibility contract; quarantined:
  all other expert families and namespace reductions; new files and metadata
  fields none.
- exact remaining blocker: none for this packet. The FN record returns to
  implemented maintenance and the TEST record to completed maintenance. No
  executable, workflow, dependency, version, release, or strategic-goal change
  occurred.

### Medium-Term Goal Checkpoint After Pass 565

- **MT1 - active/maintenance:** public endpoint ownership and fail-closed test
  selection remain accepted; no new numerical endpoint is opened here.
- **MT2 - completed:** controlled Cr2 source migration remains closed.
- **MT3 - active:** pending producer facilities and corrected paper-oracle
  interpretation are unchanged.
- **MT4 - active:** residual/protected and consumer-owned PRF questions remain
  unchanged.
- **MT5 - active/maintenance:** final v0.2 release state, three public CI
  identities, and bounded `angular_public` compatibility-floor coverage remain
  fixed.
- **MT6 - active:** classified reader and namespace carrying-cost work advances
  with the radial paper-parity family closed at twelve remaining undocumented
  exports; other expert families and next-minor de-export candidates still
  require separate decisions.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange and its
  reader front door remain accepted and unchanged.

## Cartesian Hamiltonian Producer Pass 566 - Authorize Experimental QW Geometry Diagnostics Documentation

Commit(s):
- baseline `6569fbb75107cdc24fa4213eb44dcd2a47578b7f`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-PUBLIC-QW-GEOMETRY-DOC-FN-01/TEST-01` for exactly three exported
  route-specific QW geometry diagnostics. The contract distinguishes exact
  nested-source reuse from the basis overload's normalized source construction
  and distinguishes lightweight existing-basis chain/square inspection from
  nested-diatomic construction.
- MT6 advances the classified expert reader surface without changing QW or
  nested implementations, API, numerics, workflows, or namespace. Acceptance
  requires the exact undocumented-export reduction `12 -> 9`; working-basis,
  sliced-chain, and next-minor de-export work remain separate.

Guardrail / carrying cost:
- Preferred/hard additions are `55/80` source-docstring, `25/40` combined
  reader/reference, and `10/16` docs-test lines, with no new file. Existing
  core tests and unchanged Supported-floor CI remain the numerical acceptance
  owners; no new numerical assertion or workflow change is authorized.
- deleted: nothing in this authority pass; simplified: three related route
  diagnostics receive one explicit experimental inspection boundary;
  quarantined: general geometry, serialization, Hamiltonians, solver inputs,
  arbitrary orientation, heteronuclear extensions, and namespace reductions;
  not deleted because: all three are exported live expert diagnostics.
- exact remaining blocker: repo-manager implementation within the five listed
  paths and hard budgets. If truthful documentation requires behavior, API, or
  excluded-surface work, implementation must stop and report. No checkpoint is
  due.

## Cartesian Hamiltonian Producer Pass 567 - Close Experimental QW Geometry Diagnostics Documentation

Commit(s):
- `c553c40c2a464f639c4385b6155d37963facdc3d`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted accurate expert-facing documentation for the nested diatomic,
  homonuclear chain, and axis-aligned square QW geometry diagnostics. The
  reader boundary distinguishes exact nested-source reuse from basis-overload
  construction and lightweight existing-basis chain/square inspection from
  nested construction while retaining read-only, route-specific semantics.
- MT6 advances the classified expert reader surface to maintenance. The
  undocumented-export backlog fell exactly `12 -> 9` without changing QW or
  nested behavior, defaults, API, numerics, workflow, or namespace.

Validation / carrying cost:
- Implementation commit `c553c40c2` passed relevant core geometry testsets
  `6/6`, `14/14`, `14/14`, and `35/35`, docs `151/151` and `10/10`, package
  load, authority check/self-test, Documenter, and diff checks. Remote CI
  `33545468306` passed all three gates; Docs `33545468193` passed, deployed,
  and the live reference was verified at the implementation commit.
- added: `31` source-docstring and `24` reader-documentation lines; changed:
  docs tests `+13/-1`; deleted: no implementation; simplified: three related
  diagnostics now share one explicit expert reference boundary; quarantined:
  general geometry, serialization, Hamiltonians, solvers, working-basis,
  sliced-chain, and namespace-reduction work; new files and metadata none.
- exact remaining blocker: none for this packet. The FN record returns to
  implemented maintenance and the TEST record to completed maintenance. No
  strategic-goal change occurred, and no checkpoint is due.
