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
- This live volume begins with Pass 567. Pass entries are preserved in accepted
  order; duplicate or nonmonotonic historical pass numbers are not rewritten.

## Current Strategic State

- The broad producer-documentation reorganization is complete. Schema-v3
  `authority.toml` is authoritative, generated registry/execution-whitelist
  views are checked one-way outputs, and authority CI is fail-closed.
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
are separate decisions. The execution whitelist now lives in its generated
whole-file view, and `docs_fast` protects mechanical public-surface integrity
on source pushes. The next bounded maintenance step removes only prose wording
locks from the full docs owner while retaining all mechanical policy checks;
that bounded cleanup is now accepted and in maintenance.

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
work. The audited source-layout sequence, shared private Lanczos extraction,
and mapped-representation ownership inversion are complete. The accepted
documentation-test reduction removes 220 net lines without weakening the
mechanical public-surface boundary. Further carrying-cost work requires a new
audit; no source, API, workflow, or reader-prose change follows automatically.

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

## Cartesian Hamiltonian Producer Pass 568 - Authorize Experimental Sliced Hydrogen-Chain Operator Documentation

Commit(s):
- baseline `4dcb648d32d5368f19d0104b55425e582070cb3c`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-PUBLIC-SLICED-HCHAIN-DOC-FN-01/TEST-01` for exactly the three
  undocumented sliced-chain operator functions `sliced_h1`, `sliced_vee`, and
  `sliced_row!`, curated with the three already-documented companion exports.
  The contract preserves lazy read-only H1 and two-index density-density views,
  private concrete view types, represented H1 bandwidth semantics, caller-owned
  row buffers, and bounded-only dense materialization.
- MT6 advances the classified experimental reader surface without reopening the
  accepted compact producer or numerical subsystem. Acceptance requires the
  exact undocumented-export reduction `9 -> 6`; working-basis documentation and
  every next-minor namespace transaction remain separate.

Guardrail / carrying cost:
- Preferred/hard additions are `40/55` source-docstring, `35/55` combined
  reader/reference, and `10/16` docs-test lines, with no new file. Existing IDA
  validation and unchanged Supported-floor CI remain the numerical owners; no
  numerical assertion or workflow edit is authorized.
- deleted: nothing in this authority pass; simplified: six related exports gain
  one coherent expert operator boundary; quarantined: general Hamiltonians,
  solvers, MPOs, electron sectors, four-index ERIs, dense production storage,
  private view types, internal bands/work arrays, working-basis work, and
  namespace reductions; not deleted because: all six names form the accepted
  live expert producer interface.
- exact remaining blocker: repo-manager implementation within the three listed
  documentation/source paths, one existing docs-test owner, and hard budgets.
  If truthful documentation requires behavior, storage, allocation, API, or
  excluded-surface changes, implementation must stop and report. No checkpoint
  is due.

## Cartesian Hamiltonian Producer Pass 569 - Close Experimental Sliced Hydrogen-Chain Operator Documentation

Commit(s):
- `ec3b8fa12a68118c922ae87ec3cacae821633d5b`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted one coherent expert reference for all six sliced hydrogen-chain
  exports. The three operator functions now document lazy read-only H1 and
  two-index density-density views, separate nuclear repulsion, represented-band
  semantics, private concrete view types, caller-owned row overwrite/return,
  and zero-allocation compatible-buffer behavior.
- MT6 advances the classified experimental reader surface to maintenance. The
  undocumented-export backlog fell exactly `9 -> 6` without reopening compact
  storage, numerical algorithms, allocation policy, API, or workflow behavior.

Validation / carrying cost:
- Implementation commit `ec3b8fa12` passed the minimal sliced hydrogen-chain
  owner `72/72`, docs `154/154` and `10/10`, package load, authority
  check/self-test, Documenter, and diff checks. Remote CI `33565530332` passed
  all three gates; Docs `33565530379` passed, deployed, and the live reference
  was verified at the implementation commit.
- added: `36` source-docstring and `32` reader-reference lines; changed: docs
  tests `+12/-1`; deleted: no implementation; simplified: all six exports now
  share one explicit expert operator boundary; quarantined: general molecular
  Hamiltonians, solvers, MPOs, electron sectors, four-index ERIs, dense
  production storage, private view types and internal work arrays, working-basis
  work, and namespace reductions; new files and metadata none.
- exact remaining blocker: none for this packet. The unused optional ordinary-
  branch path is removed from ownership, the FN record returns to implemented
  maintenance, and the TEST record to completed maintenance. No strategic-goal
  change occurred, and no checkpoint is due.

## Cartesian Hamiltonian Producer Pass 570 - Authorize Expert Staged Cartesian Construction Documentation

Commit(s):
- baseline `83db3f50d8f3e44a0d3d31a6742648db99e2c575`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-PUBLIC-BASE-WORKING-DOC-FN-01/TEST-01` for the retained exported
  `cartesian_base_working_basis` constructor. The contract documents its
  resolved-input, Coulomb, parent, terminal, provenance, inventory, and due-
  diligence stage while keeping the complete Hamiltonian facade, artifacts,
  solvers, corrections, and PRF wrappers separate.
- MT6 advances the final retained undocumented expert export without making
  its staged return a compatibility schema. The current six-name runtime audit
  is exact: this constructor plus five reserved namespace candidates.

Guardrail / carrying cost:
- Preferred/hard additions are `24/36` source-docstring, `20/32` reader/
  reference, and `10/16` docs-test lines, with no new file. The count-only test
  must become exact equality to the five reserved names, preventing a newly
  undocumented export from silently substituting for an intended candidate.
- deleted: nothing in this authority pass; simplified: the last retained expert
  export receives one explicit unstable staged-construction boundary;
  quarantined: return-field compatibility, PRF promotion, complete operators,
  artifacts, solvers, paper workflows, all namespace transactions, and v0.3;
  not deleted because: active drivers, source consumers, and validation use the
  staged constructor.
- exact remaining blocker: repo-manager implementation within the three listed
  source/docs/test paths and hard budgets. If truthful documentation requires
  behavior or API changes, implementation must stop without a commit and report
  the mismatch.

### Medium-Term Goal Checkpoint After Pass 570

- **MT1 - active/maintenance:** public endpoint ownership and fail-closed test
  selection remain accepted; no numerical endpoint or test owner changes here.
- **MT2 - completed:** controlled Cr2 source migration remains closed.
- **MT3 - active:** pending producer facilities and corrected paper-oracle
  interpretation are unchanged.
- **MT4 - active:** residual/protected machinery remains internal where
  classified; PRF-specific wrappers stay unexported.
- **MT5 - active/maintenance:** final v0.2 release state, path-aware CI, and the
  three public gate identities remain fixed.
- **MT6 - active:** reader classification reaches the last retained expert
  export. After implementation, only the exact five reserved namespace names
  may remain undocumented; their next-minor transactions require separate
  authority and do not independently begin v0.3.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange and its
  reader front door remain accepted and unchanged.

### Strategic Compression After Pass 570

- The broad undocumented-export count is no longer an actionable backlog.
  Supported and retained expert surfaces have bounded reader contracts; the
  remaining names are explicit namespace-policy reservations.
- Public Hamiltonian construction stays centered on
  `cartesian_base_hamiltonian`. The working-basis constructor is retained only
  as expert staged inspection/composition state, with no stable field schema.
- PRF de-promotion remains durable: internal provenance and diagnostics survive,
  but no later expert constructor silently restores their public status.
- Final v0.2 tags/releases, public numerical gates, and path-aware docs-only CI
  are stable maintenance infrastructure rather than active release work.
- The next lane after implementation closeout is an explicit next-minor
  namespace decision, not automatic documentation, source deletion,
  performance work, or v0.3 development.

## Cartesian Hamiltonian Producer Pass 571 - Close Expert Staged Cartesian Construction Documentation

Commit(s):
- `5c4f5f363cb5c955a7046e9f287c0dbb67f19abb`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted the final retained expert-export documentation. The staged
  constructor is explicitly expert/unstable, directs ordinary users to
  `cartesian_base_hamiltonian`, and preserves its current return fields as
  non-schema inspection/composition state. PRF definitions remain private.
- MT6 moves this documentation/classification lane to completed maintenance.
  The remaining undocumented exports are exactly the five separately reserved
  next-minor namespace candidates; this closeout grants no de-export or v0.3
  work.

Validation / carrying cost:
- The implementation added `21` source-docstring and `18` reader-reference
  lines and changed docs tests by `+13/-1`, within all hard budgets. The former
  count-only assertion was deleted in favor of exact five-name set equality;
  no API, behavior, numerical code, schema, helper, cache, or file was added.
- Package load, docs `157/157` and `10/10`, public Cartesian `232/232` and
  `80/80`, authority check/self-test, Documenter, and diff checks passed.
  Remote CI `33579964076` passed all three gates; Docs `33579964078` passed and
  deployed the commit-pinned live reference.
- deleted: count-only namespace coverage; simplified: the last retained expert
  export now has one explicit reader boundary; quarantined: stable staged
  schemas, PRF promotion, namespace transactions, and v0.3; exact remaining
  blocker: none for this packet. FN is implemented/maintenance and TEST is
  completed/maintenance.

## Cartesian Hamiltonian Producer Pass 572 - Authorize Projected-Q-Shell Staged Descriptor Retirement

Commit(s):
- baseline `1e8e31377efb5778df8e48c9d06c4e2237fd8a81`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-RETIRE-PQS-STAGED-DESCRIPTOR-FN-01/TEST-01` for one bounded
  deletion of the inert projected-q-shell staged descriptor, prototype
  contractions, duplicate stored construction facts, and false descriptor
  consumption/status fields. The actual boundary-mode projection, symmetric
  Lowdin cleanup, coefficient matrix, supports, metric packet, layer,
  numerical diagnostics, and consumed provenance remain authoritative.
- The live tracked closure is confined to `src/cartesian_nested_faces.jl` and
  two propagation lines in `src/cartesian_nested_diatomic.jl`. No tracked
  test, driver, tool, example, or current contract calls it. A completed
  REQ-094 worker log is historical evidence; its retained executable replay
  does not call the accessor.

Guardrail / carrying cost:
- The source pass is expected to delete approximately `597` lines and may add
  at most `12` direct-control-flow lines, with net source decrease required.
  No replacement object, tuple, alias, compatibility surface, copied
  diagnostic, test edit, helper, cache, metadata vocabulary, or new file is
  permitted. The generic by-center sidecar and active product/factorized
  machinery remain out of scope.
- Existing owners plus transient baseline/candidate probes must establish
  exact coefficient, packet, and fingerprint parity; unchanged ordinary H2
  and matched H2+ endpoints; and before/after allocation/time evidence. Any
  live executable caller, parity loss, material performance regression, test
  change, or budget overrun stops implementation without a commit.
- deleted: nothing in this authority pass; simplified: one future-sidecar
  experiment is designated for complete removal; quarantined: all adjacent
  sidecars, kernels, high-order work, API, numerical policy, workflow, and
  release work; exact remaining blocker: repo-manager implementation and
  review under these two records. No medium-term goal changes and no
  checkpoint is due.

## Cartesian Hamiltonian Producer Pass 573 - Close Projected-Q-Shell Staged Descriptor Retirement

Commit(s):
- `8fcda0086f73dbd1348aa6b261c7862f1cf64bb3`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted complete removal of the inert staged-descriptor type, prototype
  contractions, duplicate construction facts, false availability/consumption
  flags, and diatomic propagation. The active boundary-mode projection,
  symmetric Lowdin cleanup, coefficient matrix, supports, metric packet,
  numerical diagnostics, and consumed provenance remain authoritative.
- LT5 advances through a net `591`-line source reduction with no replacement
  vocabulary. No medium-term goal changes: this closes one classified carrying-
  cost item without reopening PQS physics, API, release, or high-order work.

Validation / carrying cost:
- The implementation changed only the two authorized source files by
  `+2/-593`. Coefficients, metric packets, ordinary-QW H2 matrices, and packet/
  basis fingerprints were byte-identical. Matched H2+ dimensions, topology,
  energies, captures, residuals, warnings, fingerprints, and accounting were
  unchanged.
- Core, public Cartesian `232/232`, residual-GTO `80/80`, nested augmented H2
  `464/464` plus facade `64/64`, matched release `18/18`, docs `157/157` and
  `10/10`, package load, authority check/self-test, Documenter, log bound, and
  diff checks passed. Remote CI `33644285365` passed all three numerical gates;
  Docs `33644285334` passed.
- Warm layer cost changed `2.197 ms / 9.082 MB -> 2.154 ms / 8.986 MB`; warm
  matched comparison changed `22.573 s / 2.535 GB -> 22.570 s / 2.534 GB`, with
  no material cold or warm regression. Deleted: the complete inert closure;
  simplified: coefficient assembly now returns the active result directly;
  quarantined: generic by-center and active product/factorized machinery;
  exact remaining blocker: none. FN is retired/no-grant and TEST is
  completed/no-grant. The separately audited by-center retirement requires its
  own authority, and no checkpoint is due.

## Cartesian Hamiltonian Producer Pass 574 - Authorize Generic By-Center Sidecar Retirement

Commit(s):
- baseline `01ad223d3c822c54289d5d8843b3d71adbf4c520`;
- this docs-only authority amendment.

Summary / goal advancement:
- Added `HP-RETIRE-BYCENTER-SIDECAR-FN-01/TEST-01` for complete removal of the
  uninstalled generic staged by-center carrier, its builder/attach closure,
  selector/path branches, generic nuclear overload, and two stale current-doc
  claims. Tracked source has no installer or direct caller; active polymorphic
  consumption resolves through the separately preserved product-staged owner.
- LT5 advances through another bounded carrying-cost deletion without changing
  the ordinary factorized-final path, product-staged endcap-panel path, dense
  fallback, fixed-block cache/storage, or any numerical contract. No medium-
  term goal changes, and the next checkpoint remains Pass 575.

Guardrail / carrying cost:
- The expected reduction is approximately `214` source lines. Added source is
  capped at zero, no test edit or new file is permitted, and current-doc repair
  is limited to the two exact stale claims. Ignored old probes and a completed
  paper log remain historical evidence, not compatibility obligations.
- Existing owners plus a transient external probe must establish exact
  baseline/candidate by-center, summed, and complete H1 parity for ordinary
  factorized-final and product-staged endcap-panel routes, while retaining the
  existing product-staged dense-oracle tolerance. Any live constructor, sibling
  dependency, parity loss, source addition, replacement carrier, or test change
  stops implementation without a commit.
- deleted: nothing in this authority pass; simplified: one inert sidecar route
  is designated for deletion; quarantined: high-order/deferred milestones,
  performance, namespace, workflow, release, and all replacement machinery;
  exact remaining blocker: repo-manager implementation under these two grants.

## Cartesian Hamiltonian Producer Pass 575 - Defer Generic By-Center Sidecar Retirement

Commit(s):
- authorization `9cd74a86c893b83ddcdd49e2adbefdeda28d577e`;
- sibling consumer `ed43ff241b16d1e95ea258843017c6638166a940`;
- this docs-only corrective decision.

Summary / goal advancement:
- Repo-manager correctly stopped before editing. The Pass 574 main-only caller
  proof omitted the active `high-order/manager-lane` worktree, where four
  integrated-recipe paths install the generic sidecar and the diatomic owner
  validates its carrier, route, parity, and malformed inputs.
- The high-order contract identifies generic staged units as part of intended
  general-q mainline import. Classifying that lane as historical would be
  false. `HP-RETIRE-BYCENTER-SIDECAR-FN-01` is therefore deferred/no-grant and
  its TEST record is completed/no-grant. Main retains the approximately
  `214`-line closure and current documentation unchanged.

Guardrail / next step:
- LT5 does not advance through this candidate. A separate high-order owner must
  either accept the generic carrier as maintained import machinery or authorize
  complete migration of its installers, call paths, and tests. Only then may a
  fresh cross-worktree caller audit reopen retirement. No source, test,
  high-order worktree, API, numerical, workflow, namespace, or release change
  is granted here.
- deleted: no implementation; simplified: the false caller-free claim is
  withdrawn; quarantined: all migration/replacement design; exact blocker: the
  active high-order non-product consumer. The original failure rule remains
  effective rather than being weakened after discovery.

### Medium-Term Goal Checkpoint After Pass 575

- **MT1 - active/maintenance:** public endpoint ownership and fail-closed test
  selection remain unchanged.
- **MT2 - completed:** controlled Cr2 source migration remains closed.
- **MT3 - active:** pending producer facilities and paper-oracle interpretation
  are unchanged.
- **MT4 - active:** residual/protected and high-order research machinery remain
  internal where classified; the newly confirmed sidecar dependency requires
  its own owner decision.
- **MT5 - active/maintenance:** final v0.2 release state, path-aware CI, and the
  three public gate identities remain fixed.
- **MT6 - completed/maintenance:** reader classification remains closed with
  exactly five reserved namespace candidates.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange remains
  accepted and unchanged.

## Cartesian Hamiltonian Producer Pass 576 - Reauthorize Generic By-Center Sidecar Retirement

Commit(s):
- current-main audit baseline `b4e632f1417f15e7b03b4b9a165ebf519302ce22`;
- archived high-order snapshot `ed43ff241b16d1e95ea258843017c6638166a940`;
- this docs-only authority amendment.

Summary / goal advancement:
- Accepted the high-order owner closure in
  `chatarchive/reports/software_reviews/high_order_lane_closure_and_generic_sidecar_release_2026-09-02.md`.
  The old lane is archived, is not a future merge target, and releases its
  claim on the generic sidecar. Its worktree remains physically preserved only
  because `523 MB` of ignored scratch, including July AddNest evidence, has a
  separate archive gate.
- A fresh caller audit found no current-main installer and no active paper,
  external-work, test, driver, tool, or non-archived worktree consumer. The one
  current CR2 script using the general path reporter exercises preserved
  factorized/product-staged behavior. Old detached/release copies and stale
  refs are evidence snapshots, not compatibility owners.
- `HP-RETIRE-BYCENTER-SIDECAR-FN-01/TEST-01` are reactivated for the original
  narrow deletion. LT5 may advance by approximately `214` source lines without
  changing current PQS or high-order scientific policy.

Guardrail / next step:
- Source remains confined to the two existing owners, with zero added source
  lines; two stale current-doc claims may be corrected, and committed tests may
  not change. Product-staged support, factorized-final and dense routes, exact
  by-center/complete H1 parity, and all public endpoints remain mandatory.
- deleted: nothing in this authority pass; simplified: archived worktree
  retention is separated from source compatibility; quarantined: all high-order
  and AddNest evidence plus replacement design; exact blocker: repo-manager
  implementation under the two temporary grants. Any newly found non-archived
  caller, source addition, replacement need, or parity loss stops the pass.

## Cartesian Hamiltonian Producer Pass 577 - Close Generic By-Center Sidecar Retirement

Commit(s):
- implementation `f13e946dfffd47ca9bc06e003f76f54c25316415`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted the exact `+0/-210` source deletion of the orphaned generic staged
  by-center carrier, its builder/attach closure and route branches, and its
  matching nuclear overload. Only the two authorized current-document claims
  changed; no replacement, source addition, test, API, cache, dependency,
  workflow, version, or release change was introduced.
- Baseline probes established bitwise identity for factorized-final dimension
  `433` and product-staged dimension `397`: every per-center matrix, summed
  nuclear matrix, kinetic matrix, coefficients, and complete H1 matched. The
  product-staged dense-oracle error remained `8.88e-16`, and representative
  timing/allocation measurements showed no material regression.
- Core, public Cartesian `232/232`, nested PQS `464/464`, supplemented facade
  `64/64`, matched H2+ `18/18`, and docs `157/157` plus `10/10` passed. CI run
  `33664855002` and Docs run `33664855006` were green. LT5 advances through a
  net conceptual and source reduction with current PQS numerics unchanged.

Guardrail / next step:
- `HP-RETIRE-BYCENTER-SIDECAR-FN-01` is retired/no-grant and its TEST record is
  completed/no-grant. Product-staged, factorized-final, dense fallback,
  fixed-block cache storage, and shared contraction kernels remain active.
- deleted: the complete generic sidecar closure; simplified: by-center route
  vocabulary now names only live paths; quarantined: deferred historical
  milestones and all high-order/AddNest evidence; not deleted because: the
  product-staged and factorized paths serve current endpoints; exact remaining
  blocker: none for this retirement. The archived high-order worktree remains
  physically untouched until its separate evidence-preservation gate closes.

## Cartesian Hamiltonian Producer Pass 578 - Authorize Package-Root Path Indirection

Commit(s):
- audit/authority baseline `5b844e619b317633d6747974af54d8893df653f7`;
- this docs-only authority amendment.

Summary / goal advancement:
- Read the complete 2026-09-02 external follow-up review and its 87-row metrics
  ledger as audit evidence, not authority. Independent source inspection
  corrects Step 0 to five containing-file-relative data paths plus one nested
  include, all confined to the four assigned source files. No additional path
  construction belongs in this transaction.
- Runtime inspection resolved all five data paths to their tracked package
  targets, loaded `122` full angular orders and curated orders `[15,32,51]`,
  and confirmed the radial high-precision family table is loaded. SHA-256
  baselines are fixed in the canonical contract. A separate Julia scope probe
  confirmed that explicit `parentmodule(@__MODULE__)._PACKAGE_ROOT` access is
  valid; the nested module may not assume inheritance of a private parent name.
- `HP-PACKAGE-ROOT-PATH-FN-01/TEST-01` authorize only one private root constant,
  one private data helper, six substitutions, and unchanged validation owners.
  This reduces relocation coupling without authorizing any relocation or a
  broader source-layout program; no new medium-term goal is opened.

Guardrail / next step:
- Raw source additions are capped at `10` because Git counts the six replacement
  lines; only two new declaration lines are permitted, and net growth is capped
  at `5`. Tests, files, include order, APIs, data, dependencies, workflows,
  numerics, releases, and Steps 1-4 remain unchanged.
- deleted: nothing in this authority pass; simplified: one future source-root
  dependency boundary; quarantined: every other review recommendation and all
  file moves; not deleted because: existing path functions and nested module
  remain live; exact blocker: implementation must prove identical paths and
  bytes. Any mismatch, test edit, elaborate scope workaround, or need for a
  move stops the pass without an implementation commit.

## Cartesian Hamiltonian Producer Pass 579 - Close Package-Root Path Indirection

Commit(s):
- authority `d6b19f4d597d594bd9c7ea1eb3a814eeb4611d51`;
- implementation `cadf02f7c6785fdf8fa792358837c5e83010376b`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted the exact `+8/-8` implementation across the four authorized source
  files. It adds the two approved private declarations, replaces five repeated
  data constructions, and gives the radial nested module an explicit
  parent-qualified package root. No API, test, data, file, include order,
  dependency, cache, workflow, numerical behavior, or release state changed.
- Manager rerun of the external probe confirmed all six resolved targets and
  SHA-256 values, all `122` full angular orders, curated orders `[15,32,51]`,
  and the radial prototype. Unchanged core, radial `322/322`, angular-public
  `83/83`, misc `59/59`, and docs `157/157` plus `10/10` passed. CI run
  `33670968148` and Docs run `33670968155` were green.
- `HP-PACKAGE-ROOT-PATH-FN-01` moves to implemented/maintenance and its TEST
  record to completed/no-grant. Step 0 is complete without opening a new
  medium-term goal or granting any source relocation.

Guardrail / next step:
- Preserve the private root/helper, exact six targets and hashes, and
  parent-qualified nested-module access. Any file move, path framework, test
  edit, or later source-layout step requires a fresh transaction.
- deleted: five ad hoc data-root constructions and one source-local include
  assumption; simplified: package data location now has one owner; quarantined:
  every later proposal in the external review; not deleted because: all data
  loaders and the nested high-precision module remain live; exact remaining
  blocker: none for Step 0.

## Cartesian Hamiltonian Producer Pass 580 - Authorize Radial And Atomic Source Relocation

Commit(s):
- baseline `5643e3376bd9f8a1a17ded24f05c4d2c363f89cc`;
- this docs-only authority amendment.

Summary / goal advancement:
- Accepted Step 1 of the external follow-up review only after independent
  inspection. The exact surface is `14` tracked files and `5,469` lines, not
  the review's stale `5,471`; Pass 578 left no containing-file-relative data
  dependency. The sole nested include is already parent-qualified, and no test
  directly includes an old path.
- `HP-SOURCE-LAYOUT-RADIAL-ATOMIC-FN-01/TEST-01` authorize byte-identical
  `git mv` relocation into `src/radial/` and `src/atomic/`, plus only `14`
  in-place root include substitutions and exact current-pointer reconciliation.
  The review's one-record governance estimate is corrected to two maintained
  path owners, and the current documentation closure is `12` pointers across
  seven documents.

Guardrail / carrying cost:
- Production content remains unchanged: `14` complete renames and root
  `+14/-14` are required. The implementation may mechanically reconcile the
  exact authority path entries in the same commit so Docs can pass, but may not
  change lifecycle, grant, behavior, tests, modules, APIs, include order,
  dependencies, workflows, numerics, releases, empty directories, or any later
  layout step. Historical evidence remains historical.
- deleted: no code; simplified: radial/atomic file discovery; quarantined:
  angular, foundation, ordinary, Cartesian, submodule, and all later layout
  proposals; not deleted because: every moved owner remains live; exact
  blocker: implementation must prove all hashes and `14` 100% renames, then
  pass unchanged owners and normal CI/Docs.

### Medium-Term Goal Checkpoint After Pass 580

- **MT1 - active/maintenance:** no conformance or public-surface policy changes.
- **MT2 - completed:** controlled Cr2 source migration remains closed.
- **MT3 - active:** pending producer facilities and oracle interpretation are
  unchanged.
- **MT4 - active:** residual/protected research boundaries are unchanged.
- **MT5 - active/maintenance:** final-v0.2 release, path-aware CI, authority,
  and documentation contracts remain fixed; this relocation must preserve
  their current pointers.
- **MT6 - active/maintenance:** the bounded move improves source discovery at
  zero production-content growth, but does not create a broader source-layout
  program or authorize review Steps 2-4.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange remains
  accepted and unchanged.

## Cartesian Hamiltonian Producer Pass 581 - Close Radial And Atomic Source Relocation

Commit(s):
- authority `53e292ee7b20796f983482e6208ea53180d09197`;
- implementation `812a44d4476c598448bfa3ac733054f009c55b89`;
- this docs-only lifecycle closeout.

Summary / goal advancement:
- Accepted the exact path-only Step 1 implementation: all `14` files are
  `100%` renames with the frozen SHA-256 values, the root include list changed
  exactly `+14/-14` in place, and the implementation tree is
  `a4f50fcfaadfdb5a19ecaf384878477e4fb923a0`. The two maintained radial path
  owners and all `12` current documentation pointers were reconciled. No source
  body, test, module, API, dispatch, dependency, workflow, numerical behavior,
  or release state changed.
- Manager inspection independently confirmed the rename objects, hashes,
  include diff, path closure, and exact remote heads. The unchanged owners,
  package/resource loads, docs `157/157` plus `10/10`, authority checks,
  Documenter, and diff checks passed; CI run `33679026378` and Docs run
  `33679026359` were green.

Guardrail / next step:
- `HP-SOURCE-LAYOUT-RADIAL-ATOMIC-FN-01` is implemented/maintenance and its
  TEST record is completed/no-grant. This is a discovery/layout improvement
  under MT6 with no strategic change and no authorization for review Steps
  2-4 or another source-layout program.
- deleted: no code; simplified: radial and atomic source discovery;
  quarantined: all later review moves and refactors; not deleted because: every
  relocated owner remains live; exact remaining blocker: none for Step 1.

## Cartesian Hamiltonian Producer Pass 582 - Clarify Lifecycle Transition Economy

Commit(s):
- baseline `dbfd701a773997b6400803642a6d71d0a1b176e4`;
- this policy-only clarification.

Summary / guardrail:
- Independent review of G1/G2 and the 55-commit ledger confirmed substantial
  transition overhead, while current atomic authority maintenance already
  permits multiple lifecycle changes in one commit; Pass 556 supplies the
  working precedent. `AGENTS.md` now states the bounded operating rule without
  changing deny-by-default authority, independent review, budgets, validation,
  the two-role division, schema, checkers, generated views, or canonical
  contracts.
- A fully evidenced closeout may share a commit with a separately explicit and
  independently reviewed next authorization, but a ready closeout is never
  delayed for that opportunity. Failures, deviations, blockers, tags, releases,
  public API, compatibility, and numerical policy retain separate decisions.
  Mechanical closeouts stay factual and compact; new or rewritten scopes
  normally stay near 60 words without a hard enforcement mechanism.
- MT6 remains active with no strategic or source-layout grant. The next angular
  relocation may use a combined transition only if its own review is complete
  when a prior pass is ready to close; otherwise that pass closes immediately
  and angular relocation receives separate authority. Review Step 2 remains
  unauthorized here.

## Cartesian Hamiltonian Producer Pass 583 - Authorize Angular Source Relocation

Commit(s):
- baseline `82b8d97951fa08bd2b295f8d4a9cf6316d6d3b6c`;
- this docs-only authority amendment.

Summary / goal advancement:
- Independently verified Step 2 of the external review: the five flat angular
  owners total `4,750` lines and match all supplied SHA-256 values. Their root
  includes occupy five fixed positions; none has a containing-directory path.
  `_lanczos_ground_state_apply` remains a deliberate cross-track dependency
  from matched H2+ and must move only with its unchanged containing file.
- `HP-SOURCE-LAYOUT-ANGULAR-FN-01/TEST-01` authorize five byte-identical
  `git mv` operations, five in-place include substitutions, and exact current
  path reconciliation. The corrected governance closure is three authority
  entries, three milestone pointers, and two canonical pointers.

Guardrail / next step:
- MT6 remains active. The full angular owner runs once during implementation;
  docs-only closeout must not repeat it. No test edit, helper relocation, file
  body, API, numerical, workflow, release, namespace, foundation, ordinary,
  Cartesian, or later-layout work is granted. Any non-100% rename, helper
  adaptation, test edit, or parity failure stops the implementation.
- deleted: no code; simplified: angular source discovery only; quarantined:
  every later review step; not deleted because: all five owners and the Lanczos
  helper remain live; exact blocker: implementation and unchanged validation.

## Cartesian Hamiltonian Producer Pass 584 - Close Angular And Authorize Foundation Relocation

Commit(s):
- angular implementation `f8421494e19b3a4553a7f89e20049520c38e3b2c`;
- accepted tree `61aa2190735cbb31267195722a1197ff93ddaebf`;
- this combined docs-only transition.

Summary / goal advancement:
- Accepted Step 2 exactly: five byte-identical `100%` renames, five in-place
  include substitutions, and only current authority/generated/documentation
  path updates (`+31/-31`). The Lanczos helper and both call sites remain
  unchanged. Package/resources, angular-public `83/83`, the sole complete
  angular run `61,812/61,812`, isolated HFDMRG adapter `66/66`, matched H2+
  `18/18`, docs `157/157` plus `10/10`, authority, Documenter, and diff checks
  passed. CI `33711417479` and Docs `33711417480` passed after push.
- Independent Step 3 preflight remains valid at the accepted head and corrects
  the external proposal: the architectural destination is `src/foundation/`,
  with `16` files and `6,001` lines. Thirteen must remain byte-identical; three
  may change only four live path comments. `HP-SOURCE-LAYOUT-FOUNDATION-FN-01`
  and `TEST-01` separately authorize that exact relocation, fifteen root
  include paths, one parent-qualified radial include, two direct test paths,
  and current-pointer reconciliation. This uses the Pass 582 combined-transition
  rule without delaying the completed angular closeout.

Guardrail / next step:
- MT6 remains active. Repo-manager must wait for this commit and its checks
  before Step 3. No executable source body, symbol, include order, assertion,
  workflow, numerical behavior, API, release state, Step 4, or full-angular
  rerun is authorized. Any hash, diff, package-target, or test deviation stops
  implementation without a commit.
- deleted: no code; simplified: angular discovery is complete and foundation
  discovery is the next bounded move; quarantined: cross-track dependency
  repair and all later layout work; not deleted because: all foundation owners
  remain live; exact blocker: Step 3 must satisfy the frozen map and checks.

## Cartesian Hamiltonian Producer Pass 585 - Close Foundation And Authorize Docs-Fast

Accepted foundation implementation `f2e13ff3c788551c43eb2de269d4723de880ae8e` and tree `20faffb58324a8533acd0041119cdc39235db2e2`. All `16` owners moved within the authorized
`+20/-20` source-path and `+2/-2` test-path boundary: `13` are byte-identical and three contain only the four approved comments. Package/resources, core/radial `2,174/2,174`, docs
`157/157` plus `10/10`, authority, Documenter, CI `33713199108`, and Docs `33713199099` passed. No strategic change: MT6 retains the path-only result and Step 4 remains excluded.

The same transition reopens `HP-PUBLIC-PAPER-CI-FN-01/TEST-01` only for G4: one shared mechanical public-surface owner, a `docs_fast` group appended to the existing Julia `1.10`
Supported-floor selection, and a test-local 1.10 Docs semantic fallback cross-checked against Julia `1.12.6`. Existing prose checks, `fast`, classifiers, three jobs/rows, paper/tag
lanes, and numerics remain fixed. Repo-manager must wait for this commit and its checks; equivalence, alias, exact-five-set, once-only execution, or group-selection failure stops the pass.

## Cartesian Hamiltonian Producer Pass 586 - Close Docs-Fast And Authorize Whitelist Relocation

Accepted `891bf3e84653197b90b930d7d5787be17d8fb998`: the shared owner is `98` lines, tracked tests are net `+31`, Julia 1.10
passed `6/6`, Julia 1.12 passed `8/8`, and CI `33716146165` ran all three unchanged gates with `docs_fast` visible. Docs
`33716146139` failed only on the expected planned-path guard now reconciled. No strategic change: G4 enters maintenance.
The same transition authorizes only the independently reviewed G3 whole-file whitelist relocation under
`HP-AUTHORITY-EXECUTION-WHITELIST-FN-01/TEST-01`; MT5/MT6 retain deny-by-default authority, and Step 4 remains excluded.

## Cartesian Hamiltonian Producer Pass 587 - Close Whitelist And Authorize Ordinary Relocation

Accepted G3 implementation `24512b0beda75174572fc2b728010b432331c05c`
and tree `702c005b58c324e01d44478b1da808ce0cfab7bb`. Independent review
confirmed all `226` execution IDs unchanged, sorted, and unique; whole-file
registry/whitelist rendering; stable fail-closed startup instructions; and
preserved missing, stale, modified, nongenerated, and stale-registry failures.
Authority/self-test, two independent renders, docs `8/8`, `149/149`, and
`10/10`, package load, Documenter, CI `33721567583`, and Docs `33721567704`
passed. G3 is implemented/completed maintenance under MT5 with net `-31` lines
and no source, API, numerical, workflow, or authority-model change.

The independently repeated Step 4 preflight remains exact at this head. The
Ordinary half is cleanly separable from Cartesian/PQS: `15` frozen files,
`14,430` lines, unchanged blobs/hashes, no location-sensitive include/data
dependency, and no direct source/test consumer. `HP-SOURCE-LAYOUT-ORDINARY-FN-01`
and `TEST-01` authorize only byte-identical moves into `src/ordinary/`, fifteen
in-place root include substitutions, and classified current-pointer updates.
Repo-manager must wait for this combined transition and its checks. Any hash,
rename score, body, include-order, test, pointer-count, or validation mismatch
stops without a commit; Cartesian/PQS Step 4B remains unauthorized.

No strategic change: MT6 advances source discoverability only. The bounded
ledger rotated Passes 538-566 verbatim into a `1,082`-line archive with SHA-256
`fa5207b905c69acb7094599b5d0ec67c0303c5d8e5c3b341e0082023048edd83`;
the live volume retains Passes 567-587.

## Cartesian Hamiltonian Producer Pass 588 - Close Ordinary And Authorize Cartesian/PQS Relocation

Accepted Step 4A implementation `96802b00cc380707306106b529ff22a09ccb415e`
and tree `c09786c2fd73ddd7a4271c59860a95b6359a20cc`: all `15` Ordinary
owners are byte-identical `R100` moves, the root loader has exactly `15`
in-place path substitutions, and classified current pointers plus both generated
views were reconciled. Package load, unchanged Supported-floor owners, docs,
authority/self-test, two renders, Documenter, CI `33767227348`, and Docs
`33767227357` passed. The machine-local Julia 1.10 precompile fault is not a
package failure because the clean remote Julia 1.10 gate completed. The FN
record enters maintenance and TEST closes with no grant.

Independent Step 4B preflight remains exact at this head. It authorizes `76`
files (`40,736` lines) under `src/cartesian/`, with one prerequisite path-neutral
docstring edit followed by `76` `R100` moves in a second local commit and one
push. Root include order, nested includes/imports, executable bodies, tests,
APIs, numerics, workflows, and release state remain fixed. No strategic change:
MT6 advances source discovery only. Deleted: no implementation. Simplified:
Ordinary discovery is complete. Quarantined: decomposition, Lanczos movement,
dependency inversion, orphan/conditioning work, and later layout. Exact blocker:
any hash, rename, dependency, pointer, or unchanged-owner mismatch stops Step 4B.

## Cartesian Hamiltonian Producer Pass 589 - Correct Step 4B Root-Include Count

Independent comparison of the frozen map with the root loader established that
the `76` moved files comprise `39` direct `src/GaussletBases.jl` includes and
`37` byte-identical nested sibling includes. `CartesianParentAxisFactors.jl` is
the sole mapped flat owner in the latter set and is included by
`CartesianParentGaussletBases.jl`. The authority and canonical count are
corrected from `40` to `39`; the map, hashes, scope, prerequisite commit
`381743d13269f5d75f58fe2ffc9b483248451bc2`, and all failure rules are unchanged.
No strategic change: MT6 remains active, and Step 4B may resume only after this
correction and its checks pass.

## Cartesian Hamiltonian Producer Pass 590 - Correct Step 4B Derived Counts

Independent checks against the staged relocation established seven authority
scope substitutions, not ten, and recomputed the post-Ordinary normalized
88-entry loader digest as
`8ab555c688804b36f384f0516d75407e618b9fe78b875060a685e0aa454f5f7e`.
The canonical contract now records those exact facts. The frozen 76-file map,
39 root substitutions, 37 nested includes, 340 authority paths, 264 current
documentation pointers, file blob/content hashes, scope, validation, and
failure rules are otherwise unchanged. No strategic change: MT6 and the paused
Step 4B implementation remain active.

## Cartesian Hamiltonian Producer Pass 591 - Authorize Current Evidence-Path Reconciliation

The staged relocation exposed three current `repo_path` evidence values that
the preflight counted outside `records.paths` and scope text. They belong to
`HP-RHO0-FAPP-FN-01` and `HP-RHO0-JANCHOR-FN-01` and identify retained,
source-backed dormant implementations, so they must continue resolving to
tracked files rather than remain as historical spelling. Step 4B may substitute
only those three values for their mapped `src/cartesian/` targets. The 76-file
map, 39 loader edits, 340 authority paths, seven scope literals, 264 current
documentation pointers, tests, behavior, and failure rules remain unchanged.
No strategic change: MT6 remains the active source-discovery goal.

## Cartesian Hamiltonian Producer Pass 592 - Close Cartesian Relocation And Authorize Occupied-First Maintenance

Accepted Step 4B implementation `64d8d65b6891d5e56a0b09f63ff934270b31a9fc`
and tree `739aa1d036b1b9b6f0bc83f7e5588271b5c7d15a`: all `76` owners are
`R100` moves under `src/cartesian/`, the `39` root include substitutions retain
the normalized `8ab555c688804b36f384f0516d75407e618b9fe78b875060a685e0aa454f5f7e`
identity, and path closure is exact. Package load, PQS `18/18`, Screening
`23/23`, four maintenance owners, docs `167/167`, authority, Documenter, CI
`33798470616`, and Docs `33798470735` passed. The external review's Step 0-4
source-layout sequence is complete; no runtime, API, numerical, workflow, or
release behavior changed.

The orphan-owner audit remains deliberately split. Occupied-first passed
`64/64` in `27.0` seconds at the relocated head and uniquely protects physical
Be/Ne packet-to-PQS capture and selection beyond the `misc` synthetic checks.
The existing maintenance FN/TEST pair is reopened only to add it after atomic
packet as a fifth separate process and add a focused workflow-policy assertion;
the `15`-minute ceiling and every other workflow boundary remain fixed.
Represented molecular Hartree passed its bounded standalone oracle, but
`HP-REP-MIXDENS-HARTREE-*` remains approved/implementation with the planned
contraction owner absent, so it stays unwired pending that owner decision.

## Cartesian Hamiltonian Producer Pass 593 - Close Occupied-First And Authorize Shared Lanczos

Accepted commit `9367afd5ac47b81afe4ad8e2170f274fd36f362f` and tree
`8a3458d63c1625369e4ddfcb1d57ffcf5c45573d`. The unchanged occupied-first
owner now follows atomic packet in the scheduled lane; local `64/64`, manual
run `33801669502` (all five suites in `10m17s`), CI `33801664774`, and Docs
`33801664806` passed. The `+2` workflow and `+10` policy-test delta changed no
owner, numerical code, timeout, trigger, permission, or public CI. Both
maintenance records are closed; represented Hartree remains unwired.

Independent preflight also found two structurally identical 81-line reference
Lanczos implementations. `HP-FOUNDATION-LANCZOS-FN-01/TEST-01` authorize one
exact private-helper move to `src/foundation/lanczos.jl`, one root include,
and a thin behavior-exact atomic wrapper plus at most `20` focused existing-
owner test lines. The helper's three callers remain unchanged. Default,
supplied-`v0`, maximum-step, breakdown, malformed-problem, return, error,
allocation, compilation, angular, and PQS behavior must retain exact parity;
the source hard budget is `+110/-162` and must remain net negative. MT6 advances
only if duplicate solver code is actually removed. No full angular rerun,
public solver abstraction, policy, dependency, cache, workflow, represented-
Hartree work, or release change is authorized.

## Cartesian Hamiltonian Producer Pass 594 - Close Shared Lanczos And Authorize Mapped Representation Ownership

Accepted commit `99c4890261c5fae8bdc9ed1d1c66e9270d2f66a2` and tree
`0c846fb25f6eab04ac536bb516c9c81142bfca77`. The exact 81-line private
Lanczos kernel moved to Foundation with its frozen SHA, the atomic public method
became a thin wrapper, and the duplicate recurrence was deleted. Production
source is `+91/-149`, net `-58`; the existing IDA owner gained `19` focused
lines. Independent manager replay reproduced the exact old/new snapshot and
public errors; warm allocation fell from `30,733,984` to `30,313,696` bytes,
with neutral time and no compilation regression. IDA, angular small ED, matched
H2+, and CI `33825308690` passed. Docs `33825308655` failed only on the expected
planned-to-existing path transition. MT6 advances through real duplicate-code
removal; no solver or policy surface changed.

Independent preflight remains valid at the accepted head for the last known
Foundation-to-Ordinary inversion. `HP-MAPPED-ORDINARY-REP-FN-01/TEST-01`
authorize only the frozen 25-line mapped representation block, SHA-256
`20bd94d4ce9f5a5a250b37ae1cb383a185838f5892d651f1c5074f81d1ecfb78`,
to move byte-identically from `foundation/primitive_sets.jl` immediately after
the Ordinary proxy-localization owner. The sole Cartesian caller and load order
stay unchanged; source is exactly `+25/-25`. A `10/16`-line extension of the
existing small core fixture must exercise the specialized public
representation. Any caller adaptation, parity difference, broader machinery,
or budget failure stops without an implementation commit. Complete angular,
represented-Hartree, diagnostics, optimization, and release work remain out of
scope.

## Cartesian Hamiltonian Producer Pass 595 - Close Mapped Ownership And Authorize Docs-Test Reduction

Accepted commit `3f173f63f190a4869739f070c9fc100195413bc9` and tree
`e5567e5020bbd4184d21eb15a9eb7897a22bc18e`. The exact 25-line mapped
working-basis representation block moved byte-identically from Foundation to
Ordinary with SHA-256
`20bd94d4ce9f5a5a250b37ae1cb383a185838f5892d651f1c5074f81d1ecfb78`;
the sole Cartesian caller and package load order remain unchanged. Source is
exactly `+25/-25`; the existing core fixture gained 14 focused lines. Manager
inspection confirmed the block identity, caller closure, and load order, then
replayed the complete core owner. Serialized representation parity, core,
Cartesian `232/232`, residual-GTO `80/80`, CI `33833428267`, and Docs
`33833428264` passed. Both mapped-ownership records are maintenance-only.

Independent G4 classification authorizes only
`HP-PUBLIC-DOCS-PROSE-CLEANUP-FN-01/TEST-01`. The implementation may change
`test/docs/runtests.jl` from 717 to 497 lines (`+19/-239`, net `-220`) by
retaining all mechanical policy, placement, authority, release, workflow, and
path checks, reducing four durable reader boundaries to headings, and deleting
descriptive phrase lists. The shared public-surface owner stays byte-identical.
Any documentation/source/workflow edit, lost mechanical check, net deletion
below 200 lines, or numerical-matrix execution fails closed. MT5 and MT6
advance through narrower test ownership without weakening documentation review.

### Medium-Term Goal Checkpoint After Pass 595

- **MT1 - active:** continue only evidence-led conformance repairs.
- **MT2 - completed:** no Cr2 source-migration work is reopened.
- **MT3 - active:** represented mixed-density Hartree remains separately blocked.
- **MT4 - active:** residual and protected-basis boundaries remain unchanged.
- **MT5 - active/maintenance:** generated authority views and mechanical public-surface checks are stable; bounded prose-test cleanup is next.
- **MT6 - active:** source-layout and ownership inversions are closed; the authorized net-negative docs-test cleanup is the current carrying-cost target.
- **MT7 - completed/maintenance:** external Cartesian-GTO interchange remains unchanged.

## Cartesian Hamiltonian Producer Pass 596 - Close Documentation Prose-Test Cleanup

Accepted commit `e65075e789d66e4169da4b6e88a68a30cac36e43`. The full docs
owner is exactly 497 lines after the authorized `+19/-239` change, net `-220`
lines and 31 fewer test sites. Descriptive phrase locks and their unused helper
and reads are gone; all classified mechanical checks and four durable heading
boundaries remain. The shared public-surface owner stayed byte-identical at Git
blob `b3dfc94a9362997aec84a31acab835ca7fef919c`.

Manager replay passed `docs_fast` `8/8`, full docs `122/122`, policy `10/10`,
and package load. Authority/self-test, manager-log bound, Documenter, and diff
checks passed; CI `33836202784` used only the docs lane and three visible
markers, and Docs `33836202788` passed. No numerical matrix executed and no
reader prose, source, API, workflow, helper, or file changed. MT5/MT6 record a
completed net-negative ownership cleanup; no successor transaction is opened.

## Cartesian Hamiltonian Producer Pass 597 - Authorize Portable HFDMRG Test Path

Independent inspection found the remaining portability defect confined to
`_local_hfdmrg_module()` in `test/runtests.jl`: its only machine-specific input
is a hard-coded source directory, while the accepted cache, visible absence
skip, propagated import failures, later-test continuation, and angular group
ownership are already correct. A clean detached HFDMRG source at the accepted
Pass 583 revision loaded successfully when its `src` directory was added using
the existing `LOAD_PATH`/`using` sequence.

The maintenance amendment therefore authorizes only
`GAUSSLETBASES_HFDMRG_SRC`, with blank-as-absent and a clear non-directory
failure. Preferred/hard implementation budgets are `+4/-2` and `+7/-2` in the
single existing test runner. Five fresh-process cases must distinguish absence,
bad configuration, propagated import failure, and a valid module; the valid
case reuses the isolated 66-check adapter owner. No full angular rerun,
production source, workflow, dependency, discovery, fallback, cache change, or
new test is authorized. MT1 gains a bounded conformance repair; MT5/MT6 and all
release state remain unchanged. The exact blocker is now repo-manager's
single-function implementation and validation.

## Cartesian Hamiltonian Producer Pass 598 - Close Portable HFDMRG Test Path

Accepted implementation `9004fd5087c0a24bfe253f0b166a8cf91fae813f`.
`_local_hfdmrg_module()` now reads only `GAUSSLETBASES_HFDMRG_SRC`: blank is
optional-consumer absence, a non-directory raises a clear `ArgumentError`, and
configured-directory load/import failures propagate through the unchanged
cache and handshake. The exact test-runner delta is `+4/-2`; no production
source, dependency, workflow, helper, fallback search, or numerical assertion
changed.

Repo-manager's five fresh-process cases passed, and clean HFDMRG revision
`19e11fc2b12a142138dc3324f585d8aa92d46098` passed the isolated adapter owner
`66/66` in `14m07.8s`. Manager replay independently reproduced the four cheap
control-flow cases and resolved the valid module to that clean checkout. The
angular owner blob is unchanged; its full 16-minute suite was correctly not
repeated. Package/docs/authority/Documenter checks, full source-bearing CI
`33876664840`, and Docs `33876664823` passed at the implementation SHA. MT1's
portability defect is closed; MT5/MT6 and release state are unchanged. The
remaining guardrail is explicit configuration only: no personal fallback or
filesystem discovery may return.
