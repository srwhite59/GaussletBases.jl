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
- This live volume begins with Pass 517. Pass entries are preserved in accepted
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
remain fixed. One bounded Julia `1.10` Supported-floor extension to the live
`radial` and `misc` groups is approved; `angular` remains outside CI pending a
separate runtime/ownership audit. Future tag-lane repair, registration,
citation metadata, and any later release are separate decisions.

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
assertions. One bounded post-v0.2 export cleanup is approved: delete the unused
timed nested wrapper, de-promote two undocumented diagnostic/QW names, and keep
the inexpensive qualified QW alias. The angular public surface remains a
separate audit.

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

## Cartesian Hamiltonian Producer Pass 517 - Authorize External GTO Reader Manual

Commit(s):
- this docs-only authority amendment.

Summary:
- Approved one reader-facing external Cartesian-GTO transfer page, its Manual
  and index links, one curated seven-binding API section, and focused checks in
  the existing docs owner. The page must show the checkpoint-only PySCF export,
  same-geometry working construction, overlap inspection, unchanged raw
  projection, and separately thresholded closest-determinant preparation.
- The bounds are `140/200` new-page lines, `24/36` combined existing-doc lines,
  and `20/30` docs-test lines, with exactly one new file and no source, example,
  fixture, dependency, workflow, binding, or release change.

Goal advancement / guardrail:
- MT7 remains completed/maintenance. This pass closes the authority gap for
  discoverability without reopening interchange numerics. The manual must say
  that capture loss is not repaired, no universal Gram threshold exists, and
  exporter invocation is an ordinary-HF Hamiltonian attestation rather than a
  fact provable from checkpoint fields. Basis-only/live-mean-field export and
  unsupported state classes remain outside v1. RC2 review is the next separate
  decision after reader acceptance.

Carrying-cost accounting:
- deleted: no code in this docs-only pass; simplified: one bounded page will
  replace developer-contract archaeology for ordinary use; quarantined: C2
  production evidence and paper identifiers remain outside reader docs; not
  deleted because: the developer contract remains the normative convention
  owner.
- exact remaining caller/blocker: implement and review the five authorized
  documentation surfaces; added/deleted `src` lines `0/0`; new tests none in
  this pass; new metadata/status fields: two authority records only;
  validation: public-surface/docstring/CLI reconciliation, authority rendering,
  docs tests, Documenter, manager-log bound, scoped diff, and remote Docs/CI.

## Cartesian Hamiltonian Producer Pass 518 - Accept External GTO Reader Manual

Commit(s):
- `a746038f0528c791314f273313e2f7142f3a03b0` - reader manual, navigation,
  curated reference, and focused documentation checks.
- this docs-only lifecycle closeout.

Summary:
- Accepted the checkpoint-to-bundle reader workflow, same-geometry working
  construction, final-by-source overlap inspection, unchanged raw projection,
  and explicitly thresholded closest-determinant preparation as reader-facing
  documentation. The page preserves capture-loss, occupation, Hamiltonian-
  attestation, unsupported-state, and no-universal-threshold boundaries.
- Actual additions were `126` manual lines, `18` existing-documentation lines,
  and `18` docs-test lines, all within preferred bounds. All seven curated
  bindings resolved; no source, API, fixture, workflow, or numerical behavior
  changed.

Goal advancement / guardrail:
- MT7 remains completed/maintenance, now with its discoverability prerequisite
  closed. RC2 review is the next separate decision; this pass grants no version,
  tag, release, registration, citation, basis-only/live-mean-field export, or
  solver-state work.

Carrying-cost accounting:
- deleted: stale approved-pending wording; simplified: ordinary users now have
  one bounded workflow page; quarantined: C2 production evidence and paper
  identifiers remain outside reader docs; not deleted because: the canonical
  developer contract remains the normative convention owner.
- exact remaining caller/blocker: none for reader documentation; RC2 requires
  separate authority; added/deleted `src` lines `0/0`; new tests: no new owner,
  only `18` focused existing-owner lines; new metadata/status fields none;
  validation: docs `94/94` and `10/10`, all seven bindings, authority
  check/self-test, package load, Documenter, live `/dev/` review, remote Docs
  `32748563028`, remote CI `32748562921`, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 519 - Authorize RC2 Candidate Preparation

Commit(s):
- this docs-only authority amendment.

Summary:
- Approved one bounded RC2 candidate-preparation pass because the public
  residual-GTO working-system facade and external Cartesian-GTO/PySCF transfer
  were added after immutable RC1. The pass may add three concise root-README
  links, bump only `0.2.0-rc1 -> 0.2.0-rc2`, prepend one post-RC1 changelog
  section, and retain exact RC2/RC1/dev documentation selectors with no
  `stable` alias.
- Candidate validation reuses the Julia 1.10 floor and Julia 1.12 paper gates,
  export-integrity regression, residual-GTO/interchange frozen fixture,
  examples 01/39/40/41, focused H2+ gate, and clean archive/install replay.

Goal advancement / guardrail:
- MT5 becomes active for RC2 preparation; MT7 remains completed/maintenance.
  The accepted C2 replay remains evidence and is not rerun without source
  change. Candidate preparation grants no tag, release, registration,
  citation, stable-doc, or final-v0.2 authority.

Carrying-cost accounting:
- deleted: no code; simplified: the root reader front door will name the three
  distinct public capabilities without duplicating tutorials; quarantined:
  tag/publication and final-release work remain separate; not deleted because:
  RC1 remains immutable release evidence.
- exact remaining caller/blocker: implement within existing files and stop on
  any source, numerical, test-owner, workflow, dependency, fixture-format, or
  manifest-policy change; added/deleted `src` lines `0/0`; new tests none in
  this pass; new metadata/status fields: two authority records. The authority
  checker gained only the exact existing-root-README docs-path case (`+2/-1`),
  with no schema field or broader path class. Validation: authority
  render/check/self-test, docs tests, Documenter, manager-log bound, scoped
  docs-only review, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 520 - Accept RC2 Candidate Preparation

Commit(s):
- `2b3c23970144aa030ae52b875a5cf01b32886b6e` - exact RC2 candidate.
- this docs-only lifecycle closeout.

Summary:
- Accepted version `0.2.0-rc2`, the concise post-RC1 changelog, three distinct
  README capability links, the exact RC2 selector, and focused docs checks.
  The five-file delta was `+48/-2`; the RC1 changelog section remained
  byte-identical and no source, API, numerical, dependency, example, workflow,
  fixture-format, or manifest policy changed.
- Julia `1.10.12` and `1.12.6` package loads passed. The residual-GTO/external
  transfer owner passed `80/80`, H2+ `18/18`, screening `22/22`, and local docs
  `99/99` plus `10/10`. Remote Docs `32780718822` and all three CI rows in
  `32780718923` passed. The clean `676`-entry archive excludes a root Manifest
  and both handoffs and has SHA-256
  `6728b80c1397f13b367c2d898fbdda3176c6cb87c39597817c590a0f41c1e2ac`.

Goal advancement / guardrail:
- MT5 moves to completed/maintenance for candidate preparation. RC1 remains
  immutable and MT7 remains completed/maintenance. RC2 tagging, GitHub release
  publication, registration, citation, `/stable/`, and final v0.2 require
  separate decisions; the accepted C2 replay remains canonical and was not
  repeated because implementation did not change.

Carrying-cost accounting:
- deleted: stale approved-pending lifecycle text; simplified: selector
  maintenance returns to the existing tag-deployment owner; quarantined: tag,
  publication, citation, registration, and final-release operations remain
  outside this candidate.
- not deleted because: RC1 remains immutable evidence and both candidate
  sections are reader history; exact remaining caller/blocker: none for RC2
  preparation, with tagging as the next optional decision; added/deleted `src`
  lines `0/0`; new tests: no new owner, only `20/1` focused docs-test lines;
  new metadata/status fields none; validation: exact diff and RC1 hash review,
  package/public gates, clean archive, remote Docs/CI, authority checks, docs
  tests, Documenter, manager-log bound, scoped diff, and `git diff --check`.

### Medium-Term Goal Checkpoint After Pass 520

- **MT1 - active:** continue narrow conformance repairs separately from release
  work.
- **MT2 - completed:** controlled Cr2 source migration remains closed.
- **MT3 - active:** pending producer facilities and corrected paper-oracle
  interpretation are unchanged.
- **MT4 - active:** residual/protected and consumer-owned PRF work is unchanged.
- **MT5 - completed/maintenance:** RC2 candidate preparation is accepted;
  tagging and publication are separate decisions.
- **MT6 - active:** evidence-led internal-test cleanup continues independently.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange and its
  reader front door are implemented, validated, and included in RC2.

## Cartesian Hamiltonian Producer Pass 521 - Authorize v0.2.0-rc2 Annotated Tag

Commit(s):
- this docs-only tag authority amendment.

Summary:
- Approved one annotated `v0.2.0-rc2` tag with message
  `GaussletBases v0.2.0-rc2`, frozen to candidate
  `2b3c23970144aa030ae52b875a5cf01b32886b6e`. Its tree is
  `7a4b51aec25f62436620f4ff938262d0f6b2fd62`; its accepted `676`-entry,
  `9,994,240`-byte archive has SHA-256
  `6728b80c1397f13b367c2d898fbdda3176c6cb87c39597817c590a0f41c1e2ac`.
- The local and remote tag names were absent during authorization. The
  candidate commit, tree, version, archive count, byte count, and hash matched
  exactly. The tag must target the candidate rather than current `main` or the
  Pass 520 lifecycle-closeout commit.

Goal advancement / guardrail:
- MT5 permits only immutable tag publication. Docs and all three public CI
  gates must pass afterward. A post-push failure preserves the tag and is
  reported; movement, deletion, recreation, or silent retry is forbidden.
  GitHub prerelease publication remains a separate decision.

Carrying-cost accounting:
- deleted: the blanket statement that RC2 tagging remained unauthorized;
  simplified: one exact target, tree, archive, tag name, message, and acceptance
  contract replace discretionary release procedure; quarantined: GitHub
  release, registration, citation, stable, final-v0.2, and tracked-file work;
  not deleted because: the existing Docs and public-CI workflows remain the
  validated tag-triggered invocation surfaces.
- exact remaining caller/blocker: repo-manager must create and push only the
  tag and return immutable-ref, workflow, and versioned-doc evidence for
  lifecycle closeout; added/deleted `src` lines `0/0`; new tests none; new
  metadata/status fields: two authority records only; validation: exact
  candidate and tag-absence preflight, authority render/check/self-test, docs
  tests, Documenter, manager-log bound, docs-only scope review, remote Docs/CI,
  and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 522 - Close v0.2.0-rc2 Tag Lifecycle

Commit(s):
- immutable tag object `7c8a21b998a838d245e0b5a7f4915910e2a091bc`;
- this docs-only lifecycle closeout.

Summary:
- Accepted the exact annotated `v0.2.0-rc2` tag and message. It peels to frozen
  candidate `2b3c23970144aa030ae52b875a5cf01b32886b6e` and tree
  `7a4b51aec25f62436620f4ff938262d0f6b2fd62`; the `676`-entry,
  `9,994,240`-byte archive retains SHA-256
  `6728b80c1397f13b367c2d898fbdda3176c6cb87c39597817c590a0f41c1e2ac`.
- Tag-triggered Docs `32798625043`, all three rows in CI `32798625045`, and
  Pages `32798719038` passed. RC2 and `/dev/` are live, selector order is
  RC2/RC1/dev, and `/stable/` remains absent. Only the tag ref was pushed.

Goal advancement / guardrail:
- MT5 returns to maintenance. The tag records are completed/no-grant. Preserve
  the immutable tag without movement, deletion, replacement, or recreation.
  GitHub prerelease publication remains a separate decision.

Carrying-cost accounting:
- deleted: the temporary tag execution grant and pending lifecycle language;
  simplified: immutable object, target, workflow, and deployed-doc evidence now
  replace creation instructions; quarantined: GitHub release, assets,
  registration, citation, stable, final-v0.2, and all tracked implementation;
  not deleted because: RC2/RC1 self-mappings remain live deployment policy.
- exact remaining caller/blocker: none for tag acceptance; added/deleted `src`
  lines `0/0`; new tests and metadata/status fields none; validation: local and
  remote tag identity, archive hash/size, workflow records, live canonical and
  selector checks, authority render/check/self-test, docs tests, Documenter,
  manager-log bound, docs-only scope review, remote Docs/CI, and diff checks.

## Cartesian Hamiltonian Producer Pass 523 - Authorize RC2 GitHub Prerelease

Commit(s):
- this docs-only publication authority amendment.

Summary:
- Added `HP-PQS-PUBLIC-RC2-RELEASE-FN-01/TEST-01` for one GitHub prerelease
  attached to the existing immutable RC2 tag. Preflight found no release for
  that tag; remote object `7c8a21b99` still peels to candidate `2b3c23970` and
  tree `7a4b51aec`.
- Froze the exact `2,200`-byte ASCII body, including its final newline, with
  SHA-256 `a2cbaaa2a349857e897d6d58fb728c6ccd9d731c7371bb39997bfc0360f3653a`.
  The body adds the external Cartesian-GTO transfer surface while retaining
  separate PQS and screening sections and package/solver boundaries.
- Publication requires `--verify-tag --prerelease --latest=false`, no uploaded
  assets, and no tracked mutation. Failed preflight makes no release; a
  post-publication discrepancy preserves and reports both release and tag.

Goal advancement / guardrail:
- MT5 permits exactly one RC2 GitHub prerelease publication. Final `v0.2.0`,
  registration, citation metadata, stable documentation, repository metadata,
  tag mutation, and broader release machinery remain unauthorized.

Carrying-cost accounting:
- deleted: the blanket RC2-release prohibition; simplified: one exact tag,
  title, body, command, and acceptance contract; quarantined: custom assets,
  RC1 edits, final release, registration, citation, and stable policy; not
  deleted because: candidate, tag, changelog, and versioned docs remain the
  immutable release identity and reader sources.
- exact remaining caller/blocker: repo-manager must publish once and return
  GitHub API/page, source-archive, clean-install, docs-selector, stable-absence,
  and unchanged-branch evidence; added/deleted `src` lines `0/0`; new tests
  none; new metadata/status fields: two authority records only; validation:
  synchronized-state/release-absence preflight, exact tag identity, body
  byte/hash check, authority render/check/self-test, docs tests, Documenter,
  manager-log bound, docs-only scope review, remote Docs/CI, and diff checks.

## Cartesian Hamiltonian Producer Pass 524 - Close RC2 GitHub Prerelease Lifecycle

Commit(s):
- GitHub release `376503169` at immutable tag `v0.2.0-rc2`;
- this docs-only lifecycle closeout.

Summary:
- Accepted the package-centered RC2 prerelease with exact title and tag,
  `draft=false`, `prerelease=true`, no latest final release, and zero uploaded
  assets. Its `2,200`-byte body is byte-identical to the canonical notes and has
  SHA-256 `a2cbaaa2a349857e897d6d58fb728c6ccd9d731c7371bb39997bfc0360f3653a`.
- Tag object `7c8a21b99` still peels to candidate `2b3c23970` and tree
  `7a4b51aec`. The `2,564,018`-byte tarball and `2,957,570`-byte zipball have
  recorded SHA-256 identities and both reconstruct that tree. A fresh isolated
  Julia `1.12.6` installation loaded GaussletBases `0.2.0-rc2` with the same
  manifest tree. RC2, RC1, and `dev` remain live; `/stable/` remains absent.

Goal advancement / guardrail:
- MT5 returns to completed/maintenance. Both RC2 release records are
  completed/no-grant. Preserve the immutable release and tag without editing,
  deletion, recreation, retargeting, asset upload, or narrative change. Final
  `v0.2.0`, registration, citation metadata, and stable documentation remain
  separate decisions.

Carrying-cost accounting:
- deleted: temporary publication permission and pending status; simplified:
  immutable release, archive, install, and documentation evidence now replace
  execution instructions; quarantined: final release, registration, citation,
  stable, and repository-metadata work; not deleted because: the exact body is
  durable reader and release evidence.
- exact remaining caller/blocker: none for RC2 prerelease publication;
  added/deleted `src` lines `0/0`; new tests and metadata/status fields none;
  validation: GitHub API/page state, tag peel/tree, archive reconstruction,
  clean install, live documentation and selector checks, authority
  render/check/self-test, docs tests, Documenter, manager-log bound, docs-only
  scope review, remote Docs/CI, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 525 - Authorize Path-Aware Release CI

Commit(s):
- this docs-only authority amendment.

Summary:
- Reopened `HP-PUBLIC-PAPER-CI-FN-01/TEST-01` only to implement one fail-closed
  classifier in the existing CI workflow and focused policy checks in the
  existing docs test. Pull requests and all three numerical job names remain
  unchanged. Main pushes default to the full matrix; only an exact fetched
  `github.event.before` object, successful nonempty `--no-renames` diff, and
  the four-path allowlist may select lightweight package/docs validation.
- The separate tag lane must fetch and prove the annotated
  `refs/tags/<tag>^{tag}` object before checking its commit peel, tree/version,
  remote installation, and package load. Deployed documentation remains an
  independent Docs/Pages verification by repo-manager before any separately
  authorized publication.

Goal advancement / guardrail:
- MT5 remains active but is refined: future candidate acceptance may combine
  exact candidate closeout with one version-specific conditional
  tag/publication grant, followed by one final closeout. This pass grants no
  tag, release, final-v0.2, stable, registration, citation, workflow chaining,
  or release framework. The parallel performance audit and possible duplicate
  H2+ example execution change remain separate.
- Transition validation is explicit: this authority commit incurs the legacy
  full matrix; the workflow implementation changes `ci.yml` and must run the
  full matrix; the docs-only closeout must show lightweight checks plus three
  visibly skipped-success numerical jobs and independent Docs success; the tag
  commands are rehearsed read-only against immutable RC2.

Carrying-cost accounting:
- deleted: future duplicate tag/release authorization ceremony from the
  accepted process model, but no current record or implementation; simplified:
  one existing workflow and authority pair own the routing repair;
  quarantined: performance/H2+ rewiring and all release operations remain
  separate; not deleted because: the three scientific gates remain the exact
  candidate/code acceptance boundary.
- exact remaining caller/blocker: repo-manager must implement within
  `.github/workflows/ci.yml` and `test/docs/runtests.jl`, then return staged
  transition evidence for lifecycle closeout; added/deleted `src` lines `0/0`;
  new tests: no file or numerical assertion, focused docs-policy additions
  only; new metadata/status fields none; validation: authority
  render/check/self-test, docs tests, Documenter, manager-log bound, scoped
  docs-only review, remote legacy CI/Docs, and `git diff --check`.

### Medium-Term Goal Checkpoint After Pass 525

- **MT1 - active:** continue bounded conformance repairs independently of
  release mechanics.
- **MT2 - completed:** controlled Cr2 source migration remains closed.
- **MT3 - active:** represented-Hartree scaling and external same-density
  interpretation remain separate producer/scientific work.
- **MT4 - active:** residual/protected and consumer-owned PRF evidence is
  unchanged.
- **MT5 - active, refined:** RC2 is fully published; implement and close the
  path-aware CI repair before a separately authorized final-v0.2 candidate.
- **MT6 - active:** carrying-cost and performance audits remain independent;
  no duplicate H2+ test-wiring decision is made here.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange remains
  included in immutable RC2 with no numerical or release-scope change.

## Cartesian Hamiltonian Producer Pass 526 - Authorize Direct Local PQS Shell Seeds

Commit(s):
- this docs-only authority amendment.

Summary:
- Reopened `HP-FN-00` and `HP-MCOMX-TERM-FN-01` only to replace the terminal
  PQS shell seed's global parent-row/source-column materialization with direct
  construction of the final support-row by retained-boundary-mode matrix. The
  accepted R=2 H2+ audit found eight shell seeds allocating roughly
  `969--1,028 MB`; the direct-local prototype allocated `13.75 MB` and produced
  byte-identical coefficients. The replacement must delete both global-full
  helpers, retain exact ordinary and carried-axis-fact arithmetic/order, and
  keep all Lowdin, metric, topology, provenance, due-diligence, and public
  physics behavior unchanged.
- Source authority is exactly one existing file, with `35/50` preferred/hard
  added-line limits and required net source nonincrease. No test edit, new ID,
  helper owner, fallback, workspace, cache, API, type, metadata, dependency, or
  validation-policy change is granted. Four clean-baseline/candidate block
  comparisons and the matched H2+ endpoint are mandatory; any one-bit
  difference or need for a second path stops without an implementation commit.

Goal advancement / guardrail:
- MT6 advances a measured carrying-cost/performance target without changing the
  PQS construction. Path-aware CI implementation remains a disjoint release
  lane. The compile-specialization breakdown, terminal Gaussian-sum nuclear
  loop, buffer presizing, workspace reuse, Gram policy, compatibility planning,
  and parent caching remain separately unauthorized.

Carrying-cost accounting:
- deleted: approved implementation must remove the two global-full seed
  helpers; simplified: ordinary and carried source matrices feed one direct
  local loop; quarantined: compile and nuclear optimization targets remain
  separate; not deleted because: `_nested_product_coefficients` remains a live
  general donor outside this shell-seed use.
- exact remaining caller/blocker: repo-manager must implement the one-file
  bitwise replacement and return four-route parity, endpoint, due-diligence,
  order-controlled performance, and CI evidence; added/deleted `src` lines in
  this authority pass `0/0`; new tests none; new metadata/status fields none;
  validation: source/authority review, render/check/self-test, docs tests,
  Documenter, manager-log bound, scoped diff, remote Docs/CI, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 527 - Path-Aware CI Transition Probe Rejected

Commit(s):
- `3cce96d40f9e4f06f23a190f782834271fca884b` - implement fail-closed
  path-aware main routing and the annotated-tag identity/install lane.
- this docs-only lifecycle closeout.

Summary:
- Reviewed the bounded Pass 525 implementation and attempted to return
  `HP-PUBLIC-PAPER-CI-FN-01` and `HP-PUBLIC-PAPER-CI-TEST-01` to maintenance.
  Pull requests retain the three scientific jobs. Main pushes default to the
  full matrix and may
  take the lightweight package/docs lane only after the exact fetched base and
  no-renames diff prove exclusive membership in the four-path documentation
  allowlist. Unknown or failed classification remains full. Version tags use a
  separate annotated-object, peel/tree/version, remote-install, and package-load
  lane; deployed-documentation acceptance remains in the independent Docs/Pages
  workflow.
- The implementation commit correctly classified itself full. CI run
  `32895536531` passed `Supported floor` in `5m58s`, `PQS paper` in `17m52s`,
  and `Screening paper` in `45s`; Docs run `32895536562` passed. Read-only RC2
  rehearsal verified tag object `7c8a21b99`, frozen commit `2b3c23970`, tree
  `7a4b51aec`, and a fresh remote-tag install. This closeout push is the required
  first live `docs_only` transition gate. Classifier and independent Docs run
  `32899962286` passed, but CI run `32899962285` rejected lifecycle closure:
  GitHub collapsed the unexpanded false matrix to one skipped `matrix.gate`,
  and the lightweight docs test failed `105/106` because its isolated
  Documenter fixture lacked an instantiated docs environment.

Goal advancement / guardrail:
- MT5 advances the accepted three-pass release model without weakening exact
  candidate/code validation. MT6 performance work and duplicate H2+ execution
  remain separate. No tag, publication, stable alias, source, numerical owner,
  or release authority is added.

Carrying-cost accounting:
- deleted: none because transition acceptance failed; simplified: the
  classifier and tag lane remain valid while the job graph needs one bounded
  repair; quarantined: publication, performance, and H2+ wiring remain separate;
  not deleted because: all three scientific jobs remain mandatory for PRs and
  every non-proven main push.
- exact remaining caller/blocker: expand three named no-op matrix rows on
  explicit docs-only pushes and instantiate the docs environment in the
  lightweight lane;
  added/deleted `src` lines `0/0`; new tests: no numerical tests, `22` focused
  documentation-policy lines in the implementation; new metadata/status fields
  none; validation: implementation diff and budget review, full implementation
  CI/Docs, local authority render/check/self-test, docs tests, Documenter,
  manager-log bound, scoped docs-only review, transition CI/Docs, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 528 - Authorize Path-Aware CI Transition Repair

Commit(s):
- this docs-only repair authority amendment.

Summary:
- Reopened `HP-PUBLIC-PAPER-CI-FN-01` and
  `HP-PUBLIC-PAPER-CI-TEST-01` only for the two transition failures observed in
  run `32899962285`. The existing three-row matrix must expand for every
  non-tag event. A proven documentation-only push runs one visible marker in
  each named row and skips all numerical steps; every unknown, failed, PR, or
  non-documentation main state retains the full numerical path. The lightweight
  lane must instantiate the existing docs environment before running the docs
  group. The independent Docs workflow already passed as run `32899962286`.
- Workflow and focused test repairs are bounded at `18/28` and `12/20` added
  lines respectively. No new workflow, job, row, action, group, helper, test
  file, numerical assertion, source, release operation, or tag change is
  permitted. The immutable RC2 rehearsal and tag lane remain accepted and
  untouched.

Goal advancement / guardrail:
- MT5 remains active until one full repair implementation run and one genuine
  docs-only closeout run prove the intended job graph. The failed transition is
  treated as evidence, not waived. MT6 performance and duplicate H2+ execution
  remain separate.

Carrying-cost accounting:
- deleted: none in this authority pass; simplified: the approved repair uses
  the existing matrix and Docs environment operation; quarantined: all release,
  numerical, performance, and test-wiring work remains separate; not deleted
  because: the three gate identities remain stable required checks.
- exact remaining caller/blocker: repo-manager must implement the two bounded
  workflow repairs and focused policy assertions, then return full-matrix and
  docs-only transition evidence; added/deleted `src` lines `0/0`; new tests: no
  file or numerical test, focused docs-policy edits only; new metadata/status
  fields none; validation: authority render/check/self-test, docs tests,
  Documenter, manager-log bound, scoped diff, YAML inspection, full repair
  CI/Docs, subsequent docs-only CI/Docs, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 529 - Close Corrected Path-Aware CI

Commit(s):
- `9ddc689c1bc806c7ec899cac7a39d77cb7fad3bf` - preserve the three named
  documentation-only gates and instantiate the docs environment.
- this docs-only lifecycle closeout.

Summary:
- Accepted the exact Pass 528 repair and returned
  `HP-PUBLIC-PAPER-CI-FN-01` to implemented maintenance and
  `HP-PUBLIC-PAPER-CI-TEST-01` to completed maintenance. Every non-tag event now
  expands the existing matrix. A proven documentation-only main push runs one
  visible marker per named gate and skips all numerical steps; every other
  non-tag state remains full. The lightweight job now instantiates both the
  root and docs environments before running the existing docs group.
- The repair used `+13/-1` workflow lines and `+9` focused policy-test lines,
  within both preferred budgets. Because it changed the CI workflow, run
  `32901325992` correctly took the full path and passed Screening in `50s`, the
  Supported floor in `6m`, and PQS in `17m01s`; Docs run `32901326008` passed.
  This closeout push remains the final live transition proof: its lightweight
  lane and Docs must pass, and all three named jobs must succeed through their
  marker while numerical steps remain skipped.

Goal advancement / guardrail:
- MT5 closes the path-aware routing repair while preserving exact candidate and
  code validation. Final-v0.2 authority remains separate. MT6 performance and
  duplicate H2+ execution remain untouched.

Carrying-cost accounting:
- deleted: one pre-expansion matrix-suppression condition; simplified: the same
  three jobs now cover full and marker-only paths; quarantined: release,
  performance, and H2+ rewiring remain separate; not deleted because: all three
  scientific gate identities remain stable checks.
- exact remaining caller/blocker: only the live docs-only closeout evidence;
  added/deleted `src` lines `0/0`; new tests: no file or numerical assertion,
  `9` focused policy-test lines; new metadata/status fields none; validation:
  implementation diff/budget review, full CI/Docs, authority
  render/check/self-test, docs `110/110` and `10/10`, package load, YAML,
  Documenter, manager-log bound, scoped diff, corrected transition CI/Docs, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 526 Closeout - Direct Local PQS Shell Seeds

Commit(s):
- `e1d5ca2ddb3a39134fddb476d00029ec590c431f` - replace global sparse
  shell-seed materialization with direct support-local assembly.
- this docs-only lifecycle closeout.

Summary:
- Accepted the one-file `+20/-47` implementation and returned `HP-FN-00` and
  `HP-MCOMX-TERM-FN-01` to implemented maintenance. All four ordinary/mapped
  and atomic/diatomic block sets are byte-identical; fingerprints, topology,
  dimensions, captures, energies, warnings, and column accounting are
  unchanged. Eight seeds improved from `134--216 ms` and `969--1,028 MB` to
  `2.56--2.61 ms` and `13.26--13.51 MB`; warm PQS improved from
  `0.622--0.729 s` and `1.756--1.764 GiB` to `0.280--0.384 s` and `0.874 GiB`.
  Fresh and post-WL allocations also fell substantially.

Goal advancement / guardrail:
- MT6 advances a measured production simplification without changing physics.
  The limited fresh-process time reduction leaves cold specialization,
  terminal nuclear/product buffers, Gram-policy changes, and compatibility
  cleanup as independent, unauthorized questions.

Carrying-cost accounting:
- deleted: two obsolete global-full helpers; simplified: one direct local loop
  now serves ordinary and carried axis matrices; quarantined: the four later
  optimization questions above; not deleted because: Lowdin and independent
  post-Lowdin validation remain required numerical checks.
- exact remaining caller/blocker: none for direct local seed assembly; added/
  deleted `src` lines `20/47`; new tests none; new metadata/status fields none;
  validation: exact four-route comparisons, public Cartesian `232/232`, H2+
  `18/18`, mapped/source-q `528/528`, bounded groups, authority/self-test,
  docs `110/110` and `10/10`, Documenter, diff checks, CI run `32909315394`,
  and Docs run `32909315469`.

## Cartesian Hamiltonian Producer Pass 530 - Authorize Matched-H2+ Parent Workspace

Commit(s):
- this docs-only authority amendment.

Summary:
- Reopened only `HP-PQS-PUBLIC-MATCHED-FN-01` for call-local reusable output
  and scratch in the private matrix-free parent product action and its
  high135 nuclear accumulation. The accepted audit measured `34,106` product
  calls at `10.652 s / 21.914 GiB`; `raw_nuclear` accounted inclusively for
  `21.538 GiB`, with `21.512 GiB` from its nested actions. Each comparison must
  own disjoint lexical storage; global/task-local caches, shared pools, locks,
  escaping scratch, public types, and callback-owner changes are forbidden.
- Source authority is exactly `src/pqs_matched_h2plus.jl`, with `40/60`
  preferred/hard added-line limits, at most `30` net added lines, and no test
  edit. Deterministic action parity and every dimension, fingerprint, topology,
  energy, capture, residual, symmetry, warning, and release tolerance remain
  fixed. Paired baseline/candidate measurements must remove at least `18 GiB`
  from the complete comparison without a greater than `10%` warmed-time
  regression.

Goal advancement / guardrail:
- MT6 advances a measured allocation bottleneck after the direct-local seed
  replacement. Terminal Gaussian-sum assembly, compilation/provenance cleanup,
  Gram policy, compatibility planning, route construction, release work, and
  duplicate-example policy remain independent and unauthorized.

Medium-term checkpoint:
- **MT1 - active:** continue only bounded evidence-led conformance repairs.
- **MT2 - completed:** controlled Cr2 source migration remains closed.
- **MT3 - active:** represented-Hartree scaling, corrected-WL interpretation,
  and Standard60 remain separate pending work.
- **MT4 - active:** residual/protected and consumer-owned PRF questions remain
  unchanged.
- **MT5 - active, refined:** RC2 and path-aware CI are accepted; final-v0.2,
  registration, and citation decisions remain separate.
- **MT6 - active:** direct-local shell seeds are closed; parent-action scratch
  reuse is the sole newly opened performance lane.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange and its
  reader documentation remain accepted in RC2.

Carrying-cost accounting:
- deleted: implementation must remove the allocating private hot path rather
  than preserve parallel actions; simplified: one per-comparison workspace
  serves product and nuclear accumulation; quarantined: all excluded
  performance/release lanes remain separate; not deleted because: the
  matrix-free parent oracle and existing release owner remain live.
- exact remaining caller/blocker: repo-manager must return exact parity,
  independent-workspace concurrency, isolated warmed apply, paired complete-
  comparison, due-diligence, owner, package, CI, and Docs evidence; added/
  deleted `src` lines in this authority pass `0/0`; new tests none; new
  metadata/status fields none; validation: authority render/check/self-test,
  docs tests, Documenter, manager-log bound, scoped diff, remote Docs/CI, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 531 Closeout - Matched-H2+ Parent Workspace

Commit(s):
- `fbef95e60ca3aadafe2082b4a63f8522a724e7be` - replace per-action
  temporaries with one lexical workspace per comparison.
- this docs-only lifecycle closeout.

Summary:
- Accepted the one-file `+55/-32` implementation and returned
  `HP-PQS-PUBLIC-MATCHED-FN-01` to implemented maintenance. The old allocating
  product helpers were removed. Product and complete parent actions matched
  exactly (`max_abs = 0.0`), independent tasks used disjoint workspaces, and
  the public rows remain exactly `12789/1285/1285` with unchanged energies,
  captures, topology, residuals, symmetry, and warnings.
- Paired Julia `1.12.6` comparisons reduced cumulative allocation by
  `21.606 GiB` fresh and `21.963 GiB` warm. Warm time improved from `33.074`
  to `32.795 s`; the isolated warmed action allocated `1088` bytes per call
  with no parent-sized allocation.

Goal advancement / guardrail:
- MT6 closes this allocation target. Terminal Gaussian-sum assembly,
  cold-compilation/provenance cleanup, Gram policy, compatibility planning,
  route construction, release work, and duplicate-example policy remain
  separate and unauthorized.

Carrying-cost accounting:
- deleted: allocating `_pqs_h2plus_mode_product` and
  `_pqs_h2plus_product_apply`; simplified: parent and capture actions now use
  explicit mutating storage; quarantined: all excluded optimization/release
  lanes remain separate; not deleted because: the matrix-free parent oracle
  and exact release comparison remain live.
- exact remaining caller/blocker: none for parent-action scratch reuse; added/
  deleted `src` lines `55/32`; new tests none; new metadata/status fields none;
  validation: exact action parity, independent-task probe, H2+ `18/18`, due
  diligence, package load, authority/self-test, docs `110/110` and `10/10`,
  Documenter, diff checks, CI run `32916092379`, and Docs run `32916092431`.

## Cartesian Hamiltonian Producer Pass 532 - Authorize Call-Local Terminal Buffers

Commit(s):
- this docs-only authority amendment.

Summary:
- Reopened `HP-FN-03` and `HP-DRV-STAGE-FN-01` only for pre-sizing and reuse of
  the existing call-local action, tile, and block buffers in terminal one-body
  assembly and its base low-order caller. One lexical set may serve the three
  kinetic terms and a separate lexical set may serve all nuclear centers.
  Ownership must remain disjoint and nonescaping; the 64 MiB tile bound,
  kinetic/block/center order, `mul!` and accumulation order, and exact 135-term
  scalar loop are unchanged.
- The controlled prototype preserved every matrix bit and reduced construction
  by `0.599 s / 3.942 GiB` combined: `0.397 s / 1.550 GiB` in unit-nuclear and
  `0.202 s / 2.392 GiB` in kinetic work. The implementation has a combined
  `25/35` preferred/hard added-source budget across exactly two existing files,
  with no test or file addition. Acceptance requires both phase allocations to
  fall, at least `3.0 GiB` combined reduction, and no greater than `10%`
  repeated-warm complete-comparison time regression.

Goal advancement / guardrail:
- MT6 reopens only this measured allocation lane after closing parent-action
  scratch reuse. Exact kinetic and by-center nuclear matrices, matched-H2+
  dimensions/fingerprints/energies/captures/residuals/topology/warnings, and
  release tolerances are immutable. Scalar-loop optimization, cold compilation
  and provenance cleanup, Gram policy, compatibility deletion, route work,
  duplicate-example policy, and release work remain separate.

Carrying-cost accounting:
- deleted: implementation should eliminate repeated buffer growth rather than
  retain parallel allocating call paths; simplified: existing lexical buffers
  are sized once and reused through existing internal workspace forms;
  quarantined: all excluded optimization and release lanes; not deleted because:
  the blockwise terminal kernels, 64 MiB tiling, and exact arithmetic remain
  live numerical owners.
- exact remaining caller/blocker: repo-manager must return bitwise PQS/WL
  kinetic and per-center nuclear parity, fresh/warm phase and complete-comparison
  measurements, due diligence, package/owner/CI/Docs evidence, and the scoped
  two-file diff; added/deleted `src` lines in this authority pass `0/0`; new
  tests none; new metadata/status fields none; validation: authority
  render/check/self-test, docs tests, package load, Documenter, manager-log
  bound, scoped docs-only diff, transition CI/Docs, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 533 Closeout - Call-Local Terminal Buffers

Commit(s):
- `3419da6132810d8c4454f5b013c6302ef7842cb3` - pre-size and reuse lexical
  terminal operator buffers.
- this docs-only lifecycle closeout.

Summary:
- Accepted the exact two-file `+25/-10` implementation and returned
  `HP-FN-03` and `HP-DRV-STAGE-FN-01` to implemented maintenance. PQS and
  White-Lindsey kinetic and both by-center nuclear matrices remained bitwise
  identical; the exact 135-term scalar loop was untouched. Parent axes remain
  `21 x 21 x 29`, terminal dimensions remain `1285/1285`, and column accounting
  remains `275 + 960 + 50` across eight semantic shells and two slabs.
- PQS kinetic allocation fell from `2.741 GB` to `0.140 GB`, unit-nuclear
  allocation from `1.926 GB` to `0.267 GB`, and warmed complete-comparison
  allocation from `6.816 GB` to `2.536 GB`. The combined stage reduction was
  `3.968 GiB`; warmed time improved from `32.774` to `32.028 s`.

Goal advancement / guardrail:
- MT6 closes the measured terminal-buffer lane. Scalar-loop optimization,
  cold-compilation cleanup, Gram policy, compatibility deletion, and release
  work remain separate and unauthorized. No checkpoint is due.

Carrying-cost accounting:
- deleted: repeated zero-sized production-buffer initialization; simplified:
  one pre-sized lexical buffer set serves three kinetic terms and another set
  serves physical centers; quarantined: all excluded performance/release lanes;
  not deleted because: blockwise kernels, 64 MiB tiling, and exact arithmetic
  remain live numerical owners.
- exact remaining caller/blocker: none for call-local terminal-buffer reuse;
  added/deleted `src` lines `25/10`; new tests none; new metadata/status fields
  none; validation: matrix parity, matched H2+ `18/18`, public Cartesian
  `232/232`, residual-GTO `80/80`, bounded core/IDA/examples, due diligence,
  authority/self-test, docs `110/110` and `10/10`, package load, Documenter,
  diff checks, CI run `32944547034`, and Docs run `32944547012`.

## Cartesian Hamiltonian Producer Pass 534 - Authorize Exact-Order Four-Element Gaussian Sum

Commit(s):
- this docs-only authority amendment.

Summary:
- Reopened only `HP-FN-03` for a one-file replacement of the 135-term scalar
  reduction in `_fill_terminal_gaussian_sum_action!`. Complete groups contain
  exactly four independent output elements with four explicit accumulators;
  each element retains term order and the existing left-associated arithmetic.
  A scalar remainder covers at most three elements. One call-local preweighted
  `coefficients * fx` table may be reused across existing block pairs.
- Scratch report SHA-256
  `6c33a5d264c92a45cd2d66e419ceacca901259463a6515dc6cc0a7c8fc8b703a`
  records bitwise parity and isolated improvements from `4.9875` to `2.7444 s`
  for PQS and `4.1457` to `2.1969 s` for White-Lindsey, with `491,600` bytes of
  call-local preweighting and about `0.15 s` additional fresh compilation.

Goal advancement / guardrail:
- MT6 opens only this measured scalar-loop target. Acceptance requires bitwise
  PQS/WL matrices for both centers, every release and due-diligence invariant,
  paired fresh/warm isolated and complete benchmarks, and bounded compile and
  allocation cost. The `45/50` source budget replaces the old loop and permits
  no helper, fallback, test, tuple, cache, API, type, file, or other lane.

Carrying-cost accounting:
- deleted: implementation must remove the old scalar inner loop; simplified:
  one lexical preweight and four independent accumulators replace repeated
  coefficient multiplication; quarantined: compilation, Gram, compatibility,
  route, duplicate-example, and release work; not deleted because: blockwise
  tiling, buffers, and exact termwise arithmetic remain live.
- exact remaining caller/blocker: repo-manager must stop without commit if
  parity, ownership, budget, isolated speed, complete speed, compilation, or
  allocation gates fail; added/deleted `src` lines in this authority pass
  `0/0`; new tests none; new metadata/status fields none.

## Cartesian Hamiltonian Producer Pass 535 Closeout - Exact-Order Four-Element Gaussian Sum

Commit(s):
- `94ec277d954b5435a04b0ad68ae352c95b0434c7` - batch exactly four
  independent terminal Gaussian-sum outputs while preserving each element's
  135-term arithmetic order.
- this docs-only lifecycle closeout.

Summary:
- Accepted the one-file `+50/-14` implementation at, but not beyond, the hard
  added-line limit and returned `HP-FN-03` to implemented maintenance. The old
  scalar loop was deleted; no fallback, helper, tuple machinery, API, test,
  cache, metadata, or file was added. Both nuclei's PQS and White-Lindsey
  matrices matched their frozen hashes bitwise, and all endpoint facts stayed
  unchanged.
- Warm isolated time improved `46.50%` for PQS (`5.021 -> 2.687 s`) and
  `47.09%` for White-Lindsey (`4.125 -> 2.182 s`). Warm complete time improved
  `25.35--25.60%`; fresh complete time improved `12.85--14.58%`. Compilation
  increased at most `0.134 s`; scalar accumulation remained allocation-free,
  and call-local preweighting used `491,600` bytes.

Medium-term checkpoint:
- **MT1 active:** continue only bounded evidence-led conformance repairs.
- **MT2 completed:** controlled Cr2 source migration remains closed.
- **MT3 active:** represented-Hartree scaling, corrected-WL interpretation,
  and Standard60 remain separate pending work.
- **MT4 active:** residual/protected and consumer-owned PRF questions are
  unchanged.
- **MT5 active:** RC2 and path-aware CI remain accepted; final-v0.2,
  registration, and citation decisions remain separate.
- **MT6 active, refined:** direct shell seeds, parent scratch, terminal buffers,
  and four-element reduction are closed. The next investigation is cold
  reporting, not eight-lane batching; neither has source authority here.
- **MT7 completed/maintenance:** external Cartesian-GTO interchange remains
  accepted in RC2.

Carrying-cost accounting:
- deleted: old scalar inner loop; simplified: four explicit accumulators share
  one lexical preweight; quarantined: eight-lane and all unrelated optimization
  or release lanes; not deleted because: exact-order remainder handling,
  blockwise tiling, and current buffers remain live.
- exact remaining caller/blocker: none for Pass 534; added/deleted `src` lines
  `50/14`; new tests none; new metadata/status fields none; validation:
  augmented `464/464`, release `18/18`, public Cartesian `232/232`, residual-
  GTO/MWG `80/80`, docs `110/110 + 10/10`, package load, authority/self-test,
  Documenter, diff checks, CI `33023091521`, and Docs `33023091519`.

## Cartesian Hamiltonian Producer Pass 536 - Authorize One Matched-H2+ Release Execution

Commit(s):
- this docs-only authority amendment.

Summary:
- Reopened `HP-PQS-PUBLIC-MATCHED-TEST-01` and
  `HP-PUBLIC-PAPER-CI-TEST-01` only to remove the duplicate complete H2+
  comparison. The release owner and Example 41 currently execute the same
  public comparison separately, costing `110.56 s / 19.342 GiB` together.
  Example 41 may now return its already-computed comparison as its final value;
  the release owner must apply all unchanged `18` assertions to that value, and
  the runner may delete only the duplicate subprocess smoke.
- Allowed implementation is exactly the existing Example 41, release owner,
  and root test runner. Example output remains the same three-row,
  eight-column TSV plus concise summary. Required-check identity, result types,
  tolerances, accounting, public API coverage, and nonfinite rejection remain
  exact. Fresh validation should approach `55 s / 9.7 GiB` and prove one full
  comparison rather than merely deleting coverage.

Goal advancement / guardrail:
- MT6 advances test/runtime carrying-cost reduction without reopening numerical
  performance work. No source, workflow, docs, API, release, fixture, cache,
  helper, or framework is authorized. Cold-reporting specialization remains
  closed after regressions; eight-lane batching and its measured at-most `2%`
  opportunity remain out of scope.

Carrying-cost accounting:
- deleted: implementation must remove the duplicate Example 41 subprocess;
  simplified: one public execution supplies both reader output and the release
  assertions; quarantined: screening, workflow, numerical, cold-reporting, and
  release lanes; not deleted because: standalone Example 41 and all `18`
  release checks remain live contracts.
- exact remaining caller/blocker: repo-manager must stop if the preferred
  include/return arrangement changes output, drops an assertion, executes the
  comparison more than once, fails to reduce fresh cost materially, or needs
  more than eight added lines; source delta must remain `0/0`, test/example
  delta net negative, new tests/files/metadata zero; validation: direct Example
  41, focused release owner, complete `pqs_release`, all three unchanged remote
  gates, authority/self-test, docs tests, package load, Documenter, manager-log
  bound, and diff checks.

## Cartesian Hamiltonian Producer Pass 537 - Close One-Execution H2+ Release Gate

Commit(s):
- `b0dbd9ea37317590334a24883ef0667bdb0195a5` - return Example 41's
  comparison to the release owner and remove its duplicate subprocess.
- this docs-only lifecycle closeout and bounded-ledger rotation.

Summary:
- Accepted the exact three-file `+4/-5` implementation and returned
  `HP-PQS-PUBLIC-MATCHED-TEST-01` and `HP-PUBLIC-PAPER-CI-TEST-01` to
  completed maintenance. Example 41 remains independently executable with the
  same concise output and three-row/eight-column TSV; its final value is the
  comparison consumed by all unchanged `18` release assertions. The ordinary
  successful gate now constructs the complete comparison exactly once.
- Fresh cost fell from `110.56 s / 19.342 GiB` to
  `54.91 s / 9.716 GiB`, saving about `55.65 s / 9.63 GiB`. Dimensions remain
  `12789/1285/1285`, frozen capture/energy tolerances pass, and terminal due
  diligence remains `21 x 21 x 29`, eight shells, two slabs, and exact
  `275 + 960 + 50 = 1285` accounting. TSV SHA-256 is
  `9a540cffff3b9ed87076063ed75626c55c2fda90457f4ac5a01e303ab3289d83`.

Goal advancement / guardrail:
- MT6 closes the duplicate test/runtime cost without numerical, source,
  workflow, API, fixture, or release change. Cold-reporting specialization and
  eight-lane batching remain closed; neither is a follow-on implementation
  target from this pass.

Carrying-cost accounting:
- deleted: duplicate Example 41 subprocess; simplified: one public execution
  supplies output and release assertions; quarantined: screening, workflow,
  cold-reporting, numerical, and release work; not deleted because: standalone
  Example 41 and all `18` assertions remain live.
- exact remaining caller/blocker: none for the single-execution gate; source
  delta `0/0`, test/example delta `+4/-5`, new tests/files/metadata zero;
  validation: local `18/18`, output/hash and due-diligence review, package load,
  authority/self-test, docs `110/110 + 10/10`, Documenter, diff checks, CI
  `33079914157`, and Docs `33079914164`.
- ledger rotation: Passes 496-516 moved verbatim to an `881`-line archive with
  SHA-256 `8ba3de7a83cc2e24191978ebd773e78e5763148df7427d604dbe85d16e4eba21`;
  the live ledger retains Passes 517-537.

### Strategic Compression After Pass 537

- RC2, external Cartesian-GTO interchange, reader documentation, and separate
  PQS/screening release gates are accepted maintenance surfaces.
- Path-aware CI is fail-closed and preserved; paper gates remain distinct.
- Four bounded PQS performance repairs and the duplicate release execution are
  closed with exact numerical parity and lower measured cost.
- Final v0.2, citation/registration, represented-Hartree scaling, corrected-WL
  interpretation, and any new cold-reporting work remain separate decisions.

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
