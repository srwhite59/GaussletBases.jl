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
- This live volume begins with Pass 380. Pass entries are preserved in accepted
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
- The current consumer target remains a controlled Cr2 screened-Hartree off/on
  fixed-density comparison in one high-accuracy numerical-complete basis, with
  the same imported occupied state and consumer-owned solver interpretation.
- Fixed-parent shell-q/PRF measurements now also support an internal
  parent-residual basis facility and direct-only parent-backed Gaussian
  interaction resource. Physical target selection, transition-density
  exchange, and PRF-to-GTO-residual interactions remain consumer questions.
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
discrepancies under existing authority. The immediate fail-fast sequence is
closed; continue with remaining consumer/readback correctness, due-diligence
warning shape, stale-path removal, and missing validation gates. Keep each
repair lane narrow.

**MT2 - Controlled Cr2 measurement (active).** Complete the authorized
numerical-complete screened-Hartree fixed-density comparison and the separate
fixed-parent PRF contraction study without changing each comparison's parent,
Coulomb policy, imported state, solver interpretation, or reference convention
mid-comparison.

**MT3 - Approved pending producer facilities (active).** Implement Standard60
and canonical-driver Coulomb exposure separately from the controlled Cr2 run.
The retained-GTO EGOI helper remains pending and must not absorb the unrelated
`hamiltonian_corrections.jl` WIP without review.

**MT4 - Residual and protected-basis evidence (active).** Keep the residual
spectral audit measurement-only. Protected atoms, counterpoise, and any new
injection/localization policy remain separate future decisions. Parent residual
function mechanics and the onsite-calibrated Gaussian direct resource are
approved pending implementation after comparison with MWG, full parent IDA,
and bounded continuum direct oracles. Transition-density exchange and
PRF-to-GTO-residual interactions remain measurement-only questions.

**MT5 - Documentation and authority maintenance (maintenance).** The broad
reorganization and cutover are complete. Update machine authority atomically,
keep startup documents compact, and attach documentation edits to concrete
contract or source findings rather than beginning another migration campaign.

**MT6 - Carrying-cost control (active).** Remove stale helpers, compatibility
metadata, and development-era tests as conformance work identifies them. New
scaffolding must advance a live physics or module contract and account for what
it replaces. Keep the ordinary-QW endpoint/correction audit as a deferred
retirement map rather than an immediate deletion program. Current capability
work should determine whether its unmatched correction, branch/fragment,
counterpoise, chain/square, and oracle concepts are reusable, superseded, or
obsolete.

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

## Cartesian Hamiltonian Producer Pass 380 - Compress External-GTO Representation Authority

Commit(s):
- `0452ab581` - implement the general external-GTO importer;
- `702aa4a62` - implement protected composition and the standalone sidecar;
- this commit - reconcile and compress their canonical documentation.

Summary:
- Reduced `external_gto_orbital_import.md` from `387` to `298` lines while
  preserving exact packet/result fields, `S_FG*C_G` and `S_LG*C_G`,
  validation-only `S_GG`, metric-aware full-span capture, direct
  unorthonormalized spin imports, and the complete native-order v1 sidecar key
  families. Removed the proposed packet reader and optional ladder-accessor
  story: current source consumes in-memory packets and member fields directly.
- Split two combined registry headings into four individually addressable
  lifecycle records and reduced AGENTS to the cross-overlap, orientation,
  no-transform, and no-generalized-metric guardrails.

Validation / evidence:
- Compared committed structs, exports, helpers, exact GTO handoff, member
  dependency, sidecar writer/reader/reimport, key inventory, focused test, and
  Passes 337, 363, and 365. ID-heading counts, source/test/path and stale-plan
  scans, `git diff --check`, and full local Documenter passed. No numerical test
  was needed because source, packet, and sidecar behavior are unchanged.

Goal advancement / guardrail:
- MT5: removes about `97` net documentation lines before this entry and gives
  every external-representation ID one metadata record. MT4 remains the active
  controlled numerical-complete Cr2 measurement. External import remains
  representation infrastructure, not solver or physics authority. Coulomb,
  driver, shell-resolution, protected-additive, and EGOI contracts are
  untouched.

### Medium-Term Goal Checkpoint After Pass 380

- **MT1 fake-PQS quarantine - completed/maintained.** Representation transfer
  adds no route or alternate-Hamiltonian path.
- **MT2 independent PQS recovery - completed for supported cells.** This pass
  changes no producer construction or recovery behavior.
- **MT3 common physical support vocabulary - active refinement.** Current
  shell-resolution questions remain separate; external AO identity and native
  `L` ordering do not decide shell policy.
- **MT4 supplement/reference staging - active.** The numerical-complete Cr2
  screened/unscreened measurement remains the physics target. External import
  supplies representation and capture only; consumers own orthonormalization,
  solver state, and interpretation.
- **MT5 documentation and carrying-cost cleanup - active.** Canonical ownership
  and individually addressable registry IDs continue to replace combined
  implementation-era prose.
- **MT6 old flat-path classification - maintenance.** Protected-additive,
  retained-GTO EGOI, Coulomb, driver, and shell-resolution lanes remain
  independent and were neither activated nor reinterpreted here.

## Cartesian Hamiltonian Producer Pass 381 - Canonicalize Neutral Gaussian Raw Blocks

Commit(s):
- `5da4c8a6e`, `47d9b2a3e`, `82b3f697f`, and `7c9afa8bd` - accepted nuclear
  extraction, reuse, and stable arithmetic;
- `00d052e29`, `9fa0cc16d`, `806b37e32`, and `71a89433c` - accepted
  non-nuclear extraction and family reuse;
- this commit - reconcile and compress their canonical documentation.

Summary:
- Replaced implementation plans and profiling chronology with source-backed
  nuclear and non-nuclear contracts: exact return fields, uncharged attractive
  nuclear sign, consumer-applied charges, deterministic supplement ordering,
  term/family reuse, and Residual-Gaussian/Qiu-White caller boundaries.
- Reconciled all ten IDs into individual registry records. `HP-CGAI-FN-01`
  remains unused optional authority because its proposed in-place helper never
  landed. Retained Qiu-White atomic, factor, sidecar, probe, and provider
  helpers remain live and are not deletion targets under the neutral owner.

Validation / evidence:
- Compared committed helpers, exact returned fields and dimensions, charge/sign
  assembly, operator families, direct and retained callers, implementation
  commits, focused tests, and manager Passes 078-086B. Heading/path/lifecycle
  scans, `git diff --check`, and the full local Documenter build passed. No
  Julia numerical test was needed because source and numerical behavior are
  unchanged.

Goal advancement / guardrail:
- MT5: the two canonical pages now own the numerical contract and the registry
  is metadata-oriented, removing `237` net documentation lines including
  this entry. MT4 remains parallel. Mixed-Hartree, terminal `G-G`, residual
  transforms, Qiu-White provider semantics, routes, artifacts, and solvers are
  unchanged and remain separate owners.

## Cartesian Hamiltonian Producer Pass 382 - Canonicalize Residual Core And Cutoff Policy

Commit(s):
- `93ac83ab5`, `7c0651be5`, and `7e45d5da1` - implemented the RG domain
  object, exact transforms, and MWG owner;
- `76453cc39`, `47e56593c`, `f0b662dca`, and `1f7f04e56` - implemented the
  robustness and successive default policies;
- this commit - reconcile and compress current RG authority.

Summary:
- Made the 324-line domain page the sole core-algorithm contract and the
  171-line robustness page the sole numerical-precedence contract. They now
  record the exact live object fields, owner-local metric, one final merge,
  hard negative/near-singular failures, exact operator formulas, final-basis
  MWG descriptors, weight-aware `V_GM`, and live compatibility callers.
- Converted all sixteen registry IDs to metadata records. Core and ORTHO
  remain implemented maintenance surfaces; CUTOFF-02 owns production `1e-6`.
  IDTOL and CUTOFF-01 remain addressable historical evidence but no longer
  grant source/test work. Explicit numerical-complete `1e-10` remains opt-in.

Validation / evidence:
- Reconciled source defaults, all 24 object fields, functions, callers,
  implementation commits, H2 assertions, and manager Passes 065-185. ID
  heading/whitelist, stale-language, path, cutoff-precedence, and staged-scope
  scans passed, as did `git diff --check` and the full local Documenter build.
  No numerical Julia test was needed because source and defaults are unchanged.

Goal advancement / guardrail:
- MT5 removes `545` net documentation lines including this entry. MT4
  remains parallel. No cutoff, tolerance, orientation, MWG, numerical-complete,
  injection/protected, artifact, solver, or Cr2 behavior changed; eigenvalue
  flooring and global residual selection remain forbidden.

## Cartesian Hamiltonian Producer Pass 383 - Reconcile R3 Exact-Operator Optimizations

Commit(s):
- `fb9b0414a`, `5cd9e15a6`, `79e5cd474`, `b9ad881df`, and `720912ca4` -
  implemented terminal product reuse, unit-nuclear assembly, trusted base-block
  reuse, and canonical-driver wiring;
- this commit - canonicalize their implemented contracts and lifecycle.

Summary:
- Rewrote the terminal `G-G` and unit-nuclear `U_GG` pages as compact
  source-backed contracts and added one narrow canonical home for trusted
  same-construction base kinetic/unit-nuclear reuse. Exact recomputation remains
  live, while function-local workspaces remain distinct from persistent caches.
- Converted all eight optimization IDs to individual metadata records and
  compressed startup/navigation repetition. The adjacent `HP-R3REM-AUDIT-01`
  record is now completed historical measurement with no source permission;
  its old prohibition no longer contradicts the later implementation.

Validation / evidence:
- Reconciled current helpers, fallback and trusted call paths, facade, driver,
  protected-ladder consumer, focused H2 test, implementation commits, and
  manager Passes 087-110. Heading, source/caller, stale-lifecycle, link, and
  scope scans passed, as did `git diff --check` and the full local Documenter
  build. No numerical test was run because source and behavior are unchanged.

Goal advancement / guardrail:
- MT5: removes `149` net documentation lines including this entry and makes
  exact-operator ownership machine-addressable. MT4 remains parallel. No
  kernel, allocation, caller,
  driver, residual/MWG, artifact, public API, solver, Coulomb, or Cr2 behavior
  changed. Same-construction trust is local call-graph responsibility, not a
  persisted proof or permission to reuse merely dimension-compatible matrices.

## Cartesian Hamiltonian Producer Pass 384 - Batch Canonicalize Producer Authority

Commit(s):
- this commit - canonicalize driver, mapping, White-Lindsey/shell, route-stage,
  packet, additive-reference, and protected-artifact documentation.

Summary:
- Four disjoint canonical-contract batches replaced implementation-era plans
  and copied registry prose with source-backed contracts. A new
  `route_stage_metadata_contract.md` now owns the current inventory, plan,
  carrier, ordering, and stage-handoff semantics; five old cleanup pages are
  retained only as short historical pointers.
- Registry authority is now structurally normalized to `230` unique IDs and
  `230` individual headings, with explicit lifecycle for every ID and no
  combined headings. All `181` active whitelist IDs resolve. Independent
  review restored exact per-ID source ceilings, separated dependency owners
  from executable permission, corrected validation ownership, and added the
  two previously inline rho0 candidate headings.

Validation / evidence:
- Five parallel read-only topology audits, four disjoint editing workers, and
  two independent final reviews were reconciled against committed source.
  Lifecycle, source/test path, whitelist, heading, formula, driver/facade,
  due-diligence-field, and stale-wording checks passed, as did local
  Documenter and `git diff --check`. No numerical Julia test was run because
  this pass changes documentation only.

Goal advancement / guardrail:
- MT5 advances substantially: active documentation falls by roughly `2.7k`
  lines before this entry while preserving equations, failure behavior, exact
  source ceilings, and scientific distinctions. MT4 remains unchanged. No
  source, test, API, artifact, shell policy, numerical default, EGOI, solver,
  or Cr2 authority changed; the open longitudinal shell-resolution question
  remains a separate design lane.
- Next: normalize explicit permission and canonical-owner fields for the
  remaining `50` mature registry records, then introduce non-authoritative
  shadow metadata only after exact parity. Pass 385 is the next medium-term
  checkpoint.

## Cartesian Hamiltonian Producer Pass 385 - Normalize Registry Metadata

Commit(s):
- this commit - make every registry record self-contained for later shadow
  metadata generation.

Summary:
- Added explicit permission to the remaining `50` records and explicit local
  canonical/history ownership where records still depended on parent-section
  context. All `230` unique IDs now have individual headings, lifecycle,
  permission, owner/canonical fields, and local canonical/history document
  links; all `181` whitelist IDs resolve.
- Preserved lifecycle and execution boundaries. Independent review caught and
  corrected nine evidence-only test IDs that must remain `Permission: none`
  rather than becoming validation-maintenance lanes. Closed audits,
  superseded helpers, retired workflows, and unapproved candidates remain
  non-executable. Standard60, NUMCOMP, direct-G compatibility, SPECTRAL,
  XPAIR, MIXH/FEXACT, and pending EGOI retain their prior exact boundaries.

Validation / evidence:
- Five family audits, three follow-up ownership audits, and two independent
  final reviews were reconciled. Structural checks found no missing or
  duplicate headings, lifecycle, permission, owner, or canonical links;
  whitelist parity, `git diff --check`, and local Documenter passed. No
  numerical Julia tests were run because source and behavior are unchanged.

Goal advancement / guardrail:
- MT5 completes prose-registry normalization and moves to shadow metadata.
  The pass adds about `371` net registry/status lines before this log entry;
  that deliberate carrying cost makes each authority record independently
  parseable and removes hidden section context. No source, test, API, artifact,
  numerical policy, EGOI WIP, solver, or Cr2 authority changed.
- Next: define non-authoritative `authority.toml` metadata and exact parity
  checks. Do not generate or cut over `registry.md`/`current.md` authority in
  the same pass.

### Medium-Term Goal Checkpoint After Pass 385

- **MT1 fake-PQS quarantine - completed/maintained.** Metadata normalization
  restores no retired route or wrapper path.
- **MT2 independent PQS recovery - completed for supported cells.** No
  construction or endpoint behavior changed.
- **MT3 common physical support vocabulary - active refinement.** The open
  shared-shell longitudinal COMX-resolution rule remains a separate design
  question and was not inferred from documentation metadata.
- **MT4 supplement/reference staging - active.** The controlled
  numerical-complete Cr2 screened/unscreened measurement remains consumer
  work; registry normalization changes no physics permission.
- **MT5 documentation and carrying-cost cleanup - active, next phase.**
  Canonical contracts and the prose registry are normalized. The next bounded
  step is non-authoritative shadow metadata plus parity tooling.
- **MT6 old flat-path classification - maintenance.** Retired and superseded
  records now carry explicit `Permission: none` and history ownership without
  compatibility restoration.

## Cartesian Hamiltonian Producer Pass 386 - Add Non-Authoritative Authority Shadow

Commit(s):
- this commit - add generated registry/whitelist shadow metadata and a
  fail-closed docs-build parity gate.

Summary:
- Added `registry_whitelist_shadow.toml` as an explicitly generated,
  non-authoritative, authorization-incomplete mirror of all `230` normalized
  registry records and raw membership in the `181`-ID `AGENTS.md` whitelist.
  The shadow preserves normalized lifecycle/status, permission, and
  ownership/link text plus record/source hashes; it does not infer effective
  authorization, callers, tests, dependencies, or semantic enums.
- Added a deterministic Julia generator/checker and made local/CI Documenter
  builds reject missing records, ambiguous fields, broken canonical links,
  whitelist drift, content drift, or noncanonical serialization. Markdown plus
  `AGENTS.md` remain the authority. Also repaired the missing local canonical
  link on `HP-WLDIAT-PARITY-TEST-01` found by the parser audit.

Validation / evidence:
- Standalone generation/check, byte-determinism, deliberate field/record drift,
  duplicate-state and missing-permission negative checks, full local Documenter,
  and `git diff --check` passed. Independent tooling and authority-safety
  reviews found no unaddressed authority cutover or numerical/source change.

Goal advancement / guardrail:
- MT5 advances from normalized prose to checked shadow metadata. The deliberate
  carrying cost is one roughly `400`-line checker and a generated `156` KiB
  shadow; it replaces no human contract yet. No source, test, API, artifact,
  physics, default, EGOI WIP, solver, or Cr2 authority changed.
- Next: decide separately whether a normalized lifecycle/permission taxonomy is
  useful. Do not generate authority views or cut over execution authority in
  that decision pass.

## Cartesian Hamiltonian Producer Pass 387 - Approve Semantic Shell Source-q Refinement

Commit(s):
- this commit - approve a private semantic complete-shell source-q override and
  bounded validation lane.

Summary:
- Added `HP-PQS-SHELLQ-OVERRIDE-FN-01` and
  `HP-PQS-SHELLQ-OVERRIDE-TEST-01`. The input identifies an atom-local or
  shared complete shell by actual terminal role plus positive native
  `shell_index`, requires symmetric `owner = :all`, and accepts only a
  shell-local integer `source_q > route_q`. Route metadata `q` remains fixed;
  shared selector retention, `nside`, and `selected_q` use `source_q`.
- The override acts after common shellification and before lowering,
  retained/support, and transform records freeze. Atom-local pairs receive
  cubic source dimensions; shared shells rerun the existing angular-band
  selector for `L`. Parent axes, route `q`, shell boxes/support, direct cores,
  slabs, and unmatched regions remain unchanged. The first lane is ordinary-
  span, in-memory, numerical-complete additive composition only.
- Consumer comparison imports the same external occupied determinant into each
  variant independently. No dense final-final overlap, artifact round-trip,
  mapped-COMX behavior, public/default input, or Cr2-specific source branch is
  approved.

Validation / evidence:
- Read-only source-seam and authority audits confirmed the existing
  post-shellification rewrite point and authoritative shape propagation through
  retained units, support, transforms, due diligence, and realization. Focused
  ID/link/scope scans, shadow parity, local Documenter, and `git diff --check`
  passed. No numerical test ran because this pass changes docs authority only.

Goal advancement / guardrail:
- MT3 gains a bounded diagnostic for testing local source resolution without
  changing common shell geometry or the ordinary aspect policy. MT4 may use it
  only after H2 and padded Be2 source gates pass, then compare CR2 shell 7 and
  shell 8 independently before any pair or HF run.
- Next: implement only the three approved private source surfaces and existing
  nested test. Only the numerical-complete additive call site may forward the
  recipe value; sibling/empty-reference builders reject it. Stop rather than
  edit multilayer consumers or add provenance, result fields, compatibility
  metadata, or transfer APIs.

## Cartesian Hamiltonian Producer Pass 388 - Implement Semantic Shell Source-q Refinement

Commit(s):
- `adc98163f` - implement private semantic per-shell PQS source-q overrides and
  bounded validation;
- this commit - reconcile lifecycle and record acceptance.

Summary:
- The numerical-complete additive-reference composer can now refine symmetric
  atom-local or shared complete shells by semantic role and native shell index.
  Parent geometry, support, route `q`, cores, slabs, ordering, and unmatched
  regions remain unchanged; shared shells reuse the existing angular selector
  to derive `L`. Omitted and empty requests retain exact matrix parity.
- Manager review added fail-closed collection validation, real sibling/empty
  placement/mapped-COMX composer rejection, and bounded due-diligence summary
  comparison. No public input, artifact field, mapped-COMX path, HF behavior,
  or Cr2-specific source branch was added.

Validation / evidence:
- Package load passed. The focused H2 gate passed `146/146` and the existing
  supplemented facade gate passed `69/69`; shared shell 1 changed from
  `(5,5,6)` to `(6,6,7)`, while both atom-local shell-1 owners changed from
  `(5,5,5)` to `(6,6,6)`. The doer's physically padded Be2 gate also passed
  packet recovery, metric, correction, no-reference parity, symmetry, and
  low-mode checks. `git diff --check` passed.

Goal advancement / guardrail:
- MT3 now has an implemented bounded source-resolution diagnostic; MT4 may
  proceed to the authorized CR2 shell-7 versus shell-8 fixed-density
  comparison. Source changed `+174/-13`; focused tests added `144` lines. The
  override remains private, in-memory, ordinary-span, and numerical-complete
  additive only. Production defaults and the common aspect policy are
  unchanged.

## Cartesian Hamiltonian Producer Pass 389 - Reconcile Execution And Validation Authority

Commit(s):
- this commit - reconcile whitelist grants and exact validation ownership.

Summary:
- Completed a lifecycle/grant/path audit of all `232` producer IDs, then
  reconciled prose authority without changing physics or source. Removed `33`
  IDs from the execution whitelist: ten contradictory no-grant records and 23
  validation IDs backed only by completed evidence, ignored probes, reports,
  or docs-only record maintenance. The whitelist now contains `150` exact
  executable IDs.
- Retained test authority only where ownership is concrete: mapping `s_factor`
  owns its public-facade assertions, WL terminal validation inherits one exact
  one-center smoke, and mixed-Hartree GAAA inherits placed `A-A` kernel checks.
  Pending retained-GTO EGOI now names one exact planned focused test path.
  Probe/report records remain documented but non-whitelisted.
- Retired the never-landed `HP-CGAI-FN-01` proposal rather than preserving an
  overlapping dormant grant. The regenerated shadow remains checked but
  explicitly non-authoritative.

Validation / evidence:
- Six family audits and three focused follow-ups reviewed registry text,
  canonical pages, committed Git objects, callers, tests, probes, and manager
  evidence. Exact-path, whitelist/grant, lifecycle, canonical-link, shadow,
  Documenter, staged-scope, and `git diff --check` checks passed. No numerical
  Julia test was needed for this docs-only pass.

Goal advancement / guardrail:
- MT5 removes ambiguous execution permission before candidate-v2 metadata. No
  source, test, API, artifact, numerical policy, solver, EGOI WIP, shell-q
  behavior, or Cr2 workflow changed. Next: encode the reviewed
  lifecycle/grant/surface/path decisions in a non-authoritative semantic
  candidate; authority cutover remains a later atomic decision.

## Cartesian Hamiltonian Producer Pass 390 - Add Reviewed Semantic Authority Candidate

Commit(s):
- this commit - add explicit candidate-v2 authority metadata and semantic
  validation without changing the authority source.

Summary:
- Added one explicit record for each of the `232` producer IDs with separate
  lifecycle, grant, and surface enums; raw `AGENTS.md` membership for all `150`
  executable IDs; exact canonical document headings; present source/test paths;
  planned tests; dependency IDs; scope; and raw registry-record hashes.
- Four read-only family audits reviewed every record. Manager reconciliation
  resolved mixed implemented/pending Coulomb and manifest grants,
  preservation-only compatibility, dormant-but-present rho0 helpers,
  measurement/probe records, and the pending retained-GTO EGOI test.
- The candidate is hand-reviewed rather than generated. It remains
  `authoritative = false` and `authorization_complete = false`; Markdown plus
  `AGENTS.md` continue to grant work.

Validation / evidence:
- Candidate positive checks and negative mutation checks cover enums,
  whitelist/grant compatibility, closed lifecycles, canonical headings,
  present/planned paths, dependencies, sorting, hashes, and deterministic TOML.
  Independent reviews found and closed path traversal/symlink escape and
  semantic-binding gaps; canonical owner mappings, canonical document hashes,
  and one reviewed semantic digest now fail closed. Raw-shadow parity, local
  Documenter, staged scope, and `git diff --check` also passed. No numerical
  Julia test was needed because source behavior is unchanged.

Goal advancement / guardrail:
- MT5 reaches reviewed dual representation but not authority cutover. The
  temporary carrying cost is one roughly 4,000-line candidate plus a focused
  checker; it is justified only while proving parity and should replace, not
  permanently duplicate, hand-maintained authority after a separate decision.
- Next: decide separately whether to perform one atomic authority cutover. Do
  not generate registry/current/whitelist views or flip authority flags under
  this pass.

### Medium-Term Goal Checkpoint After Pass 390

- **MT1 fake-PQS quarantine - completed/maintained.** No retired route or
  materialization authority was restored.
- **MT2 independent PQS recovery - completed for supported cells.** This pass
  changes no numerical construction or endpoint behavior.
- **MT3 common physical support vocabulary - active refinement.** Semantic
  shell-q overrides are implemented; the shared-shell longitudinal physical
  resolution rule remains a separate design question.
- **MT4 supplement/reference staging - active.** The controlled
  numerical-complete Cr2 screened/unscreened comparison remains consumer work.
- **MT5 documentation and carrying-cost cleanup - active, cutover decision
  next.** Canonical contracts, normalized prose records, raw shadow, and the
  independently reviewed semantic candidate agree under fail-closed checks.
  Authority remains Markdown.
- **MT6 old flat-path classification - maintenance.** Superseded, retired,
  preservation-only, evidence-only, and dormant-present cases are explicit and
  remain fail-closed.

## Cartesian Hamiltonian Producer Pass 391 - Harden The Authority Transition Rehearsal

Commit(s):
- this commit - separate candidate validation from legacy parity and harden the
  non-authoritative rehearsal path.

Summary:
- Migrated the reviewed `232`-record candidate to self-contained schema v3 with
  typed document roles, owned path kinds/states, evidence references, lifecycle,
  grant, surfaces, dependencies, and scope. It derives the same `150` execution
  IDs and remains `authoritative = false` and
  `authorization_complete = false`.
- Kept the prose registry and marked `AGENTS.md` block authoritative. A separate
  transition checker/snapshot now binds complete candidate bytes to current
  registry records and the exact marked whitelist block. Deterministic previews
  can be written only outside the repository and carry explicit rehearsal
  warnings.
- Closed review gaps for the active residual spectral probe, rejected
  `HP-CHANGE-01` qualification, and `HP-MCOMX-OBJ-01` ownership. Hardened
  traversal/symlink containment, raw-block parity, single-snapshot rendering,
  Markdown structure, path typing, atomic writes, and concurrent-change checks.

Validation / evidence:
- All shadow, candidate, and transition `--check` and negative `--self-test`
  gates passed. Two independent `/private/tmp` rehearsals were byte-identical,
  with the current checker hash recorded in each manifest. Local Documenter and
  `git diff --check` passed. The focused docs group remains `55/58`; its three
  failures are pre-existing stale prose assertions over unchanged pages, while
  every new authority assertion passed.

Goal advancement / guardrail:
- MT5 now has a fail-closed, CI-checked dual representation suitable for a
  second independent read-only rehearsal. No authority cutover, generated live
  registry/whitelist, physics/source/API/artifact change, or solver work is
  approved. This transition pass is `+4,007/-2,368` lines (net `+1,639`), mostly
  explicit candidate metadata and fail-closed tooling. That is temporary
  carrying cost: a later atomic cutover must retire the raw shadow, transition
  snapshot and checker, transition CI job, and duplicate Documenter hooks
  rather than preserve both authority systems.

## Cartesian Hamiltonian Producer Pass 392 - Second Authority Transition Rehearsal

Commit(s):
- this commit - record the completed no-go rehearsal; no candidate or authority
  semantics changed.

Summary:
- Generated two byte-identical external rehearsals from exact HEAD
  `321194193`, candidate digest `0028eadd...`, `232` records, `44` document
  entries, and `150` derived execution IDs. Three independent reviewers then
  covered disjoint ranges totaling `232/232` records; a fourth reviewed tooling
  and cutover readiness.
- Execution membership matched exactly, but semantic review found 19 candidate
  records with missing or over-broad owned paths, missing dependencies,
  incorrect canonical/history roles, or wrong evidence. Tooling review found
  full-document whitelist-context, manifest-binding, generated-link, preview
  warning, and cutover-transaction gaps that structural checks do not cover.
- Added one compact historical review record and updated live status; this
  docs-only pass is `+165/-7` lines. Prose authority remains unchanged.

Validation / evidence:
- Shadow, candidate, and transition checks passed before review; rehearsal
  outputs and manifests matched byte-for-byte. Reviewers covered `80 + 80 + 72`
  records with no unreviewed semantic subset. Manager triage confirmed the
  no-go categories and separated non-owned dependencies from edit authority.
  No source, test, candidate, registry, whitelist, artifact, API, or numerical
  behavior changed. Local Documenter passed. The focused docs group remains
  `55/58`; the same three pre-existing stale prose assertions failed and no new
  failure appeared.

Goal advancement / guardrail:
- MT5 advances from structural dual representation to a complete semantic
  defect inventory, but cutover is blocked. Next is one bounded metadata and
  tooling reconciliation, followed by another independent rehearsal. Do not
  treat transition-snapshot regeneration as review, convert dependency paths
  into owned surfaces, or reuse the external preview writer as an atomic
  authority installer.

## Cartesian Hamiltonian Producer Pass 393 - Reconcile And Bind The Authority Candidate

Commit(s):
- this commit - reconcile the reviewed candidate metadata and make rehearsals
  bind one captured transition state.

Summary:
- Corrected all 19 Pass 392 candidate discrepancies and one follow-up live
  dependency. Owned paths now exclude consumer dependencies, missing exact
  source/test surfaces and dependency IDs are restored, historical/canonical
  document roles and evidence are corrected, and supersession remains explicit
  scope rather than a false dependency edge. Candidate digest is
  `5af669e1517ccbb3a8cc35589541660320f7b6c74878dcaba64c758711bf86fd`.
- Hardened path/link and full-document whitelist parsing, warned standalone
  previews, and moved CI rehearsal generation behind the transition checker.
  The writer now captures candidate, transition snapshot, registry, `AGENTS.md`,
  checkers, and Git HEAD before validation; derives parity from that exact byte
  set; rejects candidate/snapshot mismatch; and rechecks every input after
  rendering.
- The candidate and transition snapshot remain non-authoritative and
  authorization-incomplete. Their semantic status remains
  `pending_independent_rehearsal`; the execution set is unchanged at 150 IDs.

Validation / evidence:
- Candidate, shadow, and transition `--check` and negative `--self-test` gates
  passed. Two transition-bound `/private/tmp` rehearsals were byte-identical.
  Independent reviewers confirmed the metadata reconciliation and stress-tested
  candidate mismatch plus post-capture registry mutation. Local Documenter and
  `git diff --check` passed. The focused docs group remains `55/58`; only the
  same three pre-existing stale prose assertions fail.

Goal advancement / guardrail:
- MT5 now has a reconciled, structurally bound rehearsal candidate ready for a
  new independent semantic review. It does not authorize cutover. The temporary
  net growth is concentrated in fail-closed transition tooling; do not expand
  this dual system further before the next rehearsal decides whether an atomic
  cutover can retire the shadow and transition machinery.

## Cartesian Hamiltonian Producer Pass 394 - Independent Authority Rehearsal And Reconciliation

Commit(s):
- this commit - complete the semantic rehearsal and its confined two-record
  reconciliation without promoting machine authority.

Summary:
- Generated two byte-identical transition-bound rehearsals from exact HEAD
  `5226ad711`, candidate digest `5af669e1...`, 232 records, 44 hashed
  documents, and 150 execution IDs. Three fresh reviewers covered disjoint
  `80 + 80 + 72` ranges; a fourth reviewed tooling and cutover readiness.
- Tooling and records 81-232 passed. Records 1-80 initially exposed two
  semantic omissions: `HP-COMP-FACEPROD-FN-01` lacked both terminal-realizer
  consumers, and `HP-COMP-THINSLAB-FN-01` lacked those realizers plus its
  conditionally authorized native metadata/route-summary support surfaces.
  All other 230 records and the full execution set agreed.
- At user direction, corrected those two records in the same pass. Parsed diff
  confinement showed that only their paths and the thin-slab conditional scope
  changed. The original finder and a fresh reviewer both accepted the fix,
  exact 150-ID parity, and deterministic corrected rehearsals. Candidate digest
  is now `30cf4ed840b00c09da39ba4e15b3cb6c3d2c1376263877a2e7778f5a15bef716`.
- Added a durable review report and marked manual semantic review complete. The
  candidate, whitelist, and prose permissions remain non-authoritative/unchanged;
  no source, test, numerical, artifact, or workflow behavior changed.

Validation / evidence:
- Candidate, shadow, and transition checks/self-tests passed before review and
  after final status recording; manifests and output hashes matched and each
  independent render pair was deterministic. Review covered every record plus
  tooling. Local Documenter and `git diff --check` passed. The aggregate docs
  group remains `55/58` only because of the same three stale prose assertions.

Goal advancement / guardrail:
- MT5 now has reviewed semantic parity for all 232 records and clean transition
  tooling, but still no cutover authority. Next is a separate atomic
  cutover-design decision covering machine-source promotion, generated live
  views, dual-system archival/removal, and whole-commit rollback/fail-closed
  behavior. Do not execute cutover in that design pass. This pass is
  `+225/-20` lines: 134 lines preserve review evidence, while the candidate
  growth is confined to the missing path inventory and one scope sentence; no
  new checker or schema layer was added.

## Cartesian Hamiltonian Producer Pass 395 - Review Atomic Authority Cutover

Commit(s):
- this commit - add the reviewed one-commit cutover design; no cutover executed.

Summary:
- Added `authority_atomic_cutover_plan.md` for the future transition from the
  reviewed schema-v3 candidate to one machine-owned record authority. The plan
  pins base `9b283e16c`, candidate digest `30cf4ed8...`, 232 records, 44 hashed
  documents, and 150 execution IDs. Exact promotion preserves every record and
  changes only artifact kind plus the two authority flags; expected promoted
  digest is `6057ef50...`.
- The design permits one direction only: `authority.toml` generates the full
  registry and marked `AGENTS.md` block. It removes active per-ID prose
  duplication, all shadow/transition files and reverse parsing, and replaces
  three temporary checkers with one permanent fail-closed checker.
- Activation is one Git commit/ref update prepared and validated in clean
  machine-local worktrees. The plan defines exact files, CI dependency,
  deterministic rendering, file-link/heading validation, expected docs-test
  cleanup, line-reduction gate, source freeze, and whole-commit rollback. It
  grants no execution authority.

Validation / evidence:
- Three independent read-only audits covered transaction safety, generated
  views/prose retirement, and permanent tooling/CI carrying cost. Written-plan
  review found and closed parent-pin, frozen-document, negative-test,
  external-render, AGENTS-range, file-allowlist, rollback, and terminology
  gaps. All three reviewers returned go for the plan only. Cutover execution
  remains no-go; no source or live authority changed.

Goal advancement / guardrail:
- MT5 moves from semantic parity to a concrete atomic transaction design. The
  next action is explicit user approval pinning this Pass 395 commit, not an
  inferred cutover.
  The implementation should remove roughly 5,200 or more net lines; a result
  below that gate requires review for retained dual machinery.

### Medium-Term Goal Checkpoint After Pass 395

- **MT1 fake-PQS quarantine - completed/maintained.** No retired materializer
  or compatibility authority is restored.
- **MT2 independent PQS recovery - completed for supported cells.** This pass
  changes no basis, operator, or endpoint behavior.
- **MT3 common physical support vocabulary - active refinement.** Semantic
  shell-q overrides are implemented; shared-shell longitudinal resolution
  remains a separate physics/design question.
- **MT4 supplement/reference staging - active.** Controlled numerical-complete
  Cr2 consumer measurements remain separate from documentation authority.
- **MT5 documentation and carrying-cost cleanup - active, cutover design
  reviewed.** Record parity and the transaction plan are reviewed; execution
  remains explicitly unapproved pending a new user/design-manager handoff.
- **MT6 old flat-path classification - maintenance.** Closed and historical
  records remain explicit and fail-closed.

## Cartesian Hamiltonian Producer Pass 396 - Activate Machine Authority

Commit(s):
- this commit - atomically promote the reviewed authority and remove the dual
  prose/shadow system.

Summary:
- Promoted the exact schema-v3 candidate to authoritative `authority.toml`
  with digest `6057ef50...`, preserving all 232 records, 44 hashed documents,
  and 150 execution IDs. The generated registry and marked `AGENTS.md` block
  are now checked one-way views; manually authored policy and contracts may
  restrict but never broaden machine authority.
- Frozen subsystem wording that the registry "owns" an ID is now explicitly
  interpreted as the exact generated human view, never a second prose source;
  machine authority controls any discrepancy.
- Replaced three temporary checkers with one permanent fail-closed checker,
  removed the raw shadow and transition snapshot, and made the docs workflow
  depend on the read-only `cartesian-authority` job. The completed cutover
  record preserves activation and whole-commit rollback rules.

Validation / evidence:
- The exact authorized parent `a8076f2fa` passed all six legacy gates and the
  reviewed `30cf4ed8...` candidate, `232/44/150` counts, whitelist digest, and
  two AGENTS range digests matched. Permanent check/self-test, repeated fresh
  rendering, focused docs tests (`53/53` plus `10/10`), local Documenter, link
  and legacy-absence scans, staged scope, detached-checkout replay, and
  whole-commit revert/tree proof form the activation gate.

Goal advancement / guardrail:
- MT5 machine-authority cutover is completed. Future authority amendments must
  edit `authority.toml`, regenerate both views externally, and land them in one
  checked commit; reverse parsing or prose fallback is forbidden.
- The transaction changes no source, numerical contract, artifact, driver, or
  hashed canonical document and removes more than 5,600 net lines. Any failed
  post-activation authority/docs CI requires a whole-commit revert before
  producer work resumes.

## Cartesian Hamiltonian Producer Pass 397 - Compress Startup Guidance

Commit(s):
- this commit - remove duplicated live status from producer startup pages and
  record the post-cutover maintenance procedure.

Summary:
- Replaced the long subsystem-by-subsystem `README.md` orientation with a
  stable six-stage pipeline and task-directed contract map. Compressed
  `current.md` to live facilities, active work, the Cr2 target, and current
  blockers; architecture guardrails remain in `invariants.md`.
- Added the compact one-way authority maintenance sequence: edit the machine
  record and contracts, update explicit document hashes, render externally,
  replace only the generated views, and validate before one atomic commit.
  No authority record, generated view, canonical contract, source, test,
  artifact, or numerical behavior changed.

Validation / evidence:
- Permanent authority `--check` and `--self-test`, focused documentation tests
  (`53/53` and `10/10`), local Documenter, link/path review, staged docs-only
  scope review, and `git diff --check` passed.

Goal advancement / guardrail:
- MT5 documentation cleanup remains active after the authority cutover. Startup
  reading is materially shorter without weakening fail-closed authority or
  endpoint due-diligence rules. Further AGENTS or hashed-contract compression
  remains a separately reviewed lane. This pass changes `+86/-220` lines, a
  net reduction of 134 lines.

## Cartesian Hamiltonian Producer Pass 398 - Audit Execution Conformance

Commit(s):
- this commit - preserve the first complete post-cutover conformance audit and
  harden authority dependency semantics.

Summary:
- Eight disjoint static reviews covered all `150` execution records at exact
  baseline `c8c1a4911`. Results are `107` matched, `11` accurately documented
  gaps, `8` requiring numerical validation, and `24` source/test/authority or
  wording discrepancies. The durable review separates fail-fast correctness,
  stale compatibility code, missing endpoint gates, and wording-only fixes so
  later passes remain bounded.
- The audit confirmed that current screened-Hartree source and contracts agree:
  determinant orbitals define `P0/q0`, the density fit defines `E0`, and the
  fitted potential approximates `J0`. Determinant-exact `J0/E0` would be a new
  scientific policy, not a conformance repair.
- Fixed the sole authority dependency cycle by retaining only the z-axis
  extension's dependency on the generic R3U facade. The checker now rejects
  dependency cycles and execution dependencies on closed/no-grant records.

Validation / evidence:
- Permanent authority check/self-test, focused docs tests, deterministic
  external rendering, generated-view parity, local Documenter, report-matrix
  completeness, staged scope review, and `git diff --check` form the gate. No
  producer source or numerical endpoint changed in this pass.

Goal advancement / guardrail:
- MT5 broad documentation reorganization is complete and moves to maintenance.
  Next work is conformance remediation under existing IDs, beginning with
  packet/one-body/injection fail-fast checks and misleading completed-test
  claims. Do not combine that work with a screened-reference policy change.
  Carrying cost is concentrated in one 224-line evidence report and a 79-line
  checker hardening; generated-view changes are three digest/dependency rows
  and no producer source is added.

## Cartesian Hamiltonian Producer Pass 399 - Fail Fast At Atomic Packet Boundaries

Commit(s):
- `d1bf11051` - validate atomic packet writes and direct consumers before using
  unconverged or malformed reference data.

Summary:
- Direct `P0/q0`, density-fit Hartree, and potential-fit Hartree consumers now
  reject packet data unless RHF convergence is explicit. Packet writes validate
  before creating an output file, so invalid in-memory packets leave no partial
  artifact.
- One ordinary-potential-fit validator now checks the complete row consumed by
  readback: compact scaffold count, coefficient/exponent shape and finiteness,
  positive exponents, fixed/trimmed term algebra, SVD-rank bounds, radial and
  tail diagnostics, charge/self-energy agreement, and determinant consistency
  algebra. The finite fitted-potential consistency error remains diagnostic;
  no `1e-8 Ha` magnitude gate or packet-physics change was introduced.
- Manager review caught and fixed a narrower first draft that bounded retained
  rank only by the 45-term source scaffold. Validation now uses the actual
  pre/post-trim fitted dimension and packet validation shares the same guard.

Validation / evidence:
- Package load passed in `7.0 s`; atomic packet tests passed `117/117` in
  `28.6 s`; screened-Hartree correction tests passed `93/93` in `64.3 s`;
  staged `git diff --check` passed. The source/test change is `+186/-12` lines.

Goal advancement / guardrail:
- This closes the Pass 398 atomic-packet fail-fast discrepancy under the
  existing packet IDs and advances MT4 reliability without changing the
  determinant/density-fit/potential-fit role split. Remaining immediate
  conformance work is the nonfinite one-body coefficient boundary, direct
  injection Gram/identity guards, and protected one-body diagnostics.

## Cartesian Hamiltonian Producer Pass 400 - Rotate The Manager Ledger

Commit(s):
- this commit - archive the pre-Pass-380 ledger verbatim and establish a compact
  live manager volume.

Summary:
- Moved the first `28,549` lines of the prior ledger unchanged into the
  task-gated manager-log history. The preserved-prefix hash is recorded in the
  live index, while accepted Passes 380-399 remain verbatim in this file.
- Replaced the stale initiation-era current-goal preamble with the post-cutover
  strategic state, active conformance and Cr2 goals, durable guardrails, and a
  bounded rotation policy. No authority, subsystem contract, source, test,
  artifact, workflow, or numerical behavior changed.

Validation / evidence:
- Verified the archive and retained-tail hashes against the exact pre-rotation
  split, reconstructed the original ledger byte-for-byte from those two pieces,
  checked links and pass boundaries, ran permanent authority check/self-test,
  local Documenter, and `git diff --check`.

Goal advancement / guardrail:
- MT5 startup carrying cost falls from `29,344` lines to roughly one thousand
  live lines while the full historical ledger remains available. Future
  rotations are mechanical evidence moves, not permission to rewrite accepted
  decisions or duplicate authority.

### Medium-Term Goal Checkpoint After Pass 400

- **MT1 conformance remediation - active.** Atomic packet boundaries are
  closed; one-body/injection/protected diagnostics are next.
- **MT2 controlled Cr2 measurement - active.** This docs-only pass changes no
  construction or comparison input.
- **MT3 pending producer facilities - active.** Standard60 and retained-GTO
  EGOI remain separate approved work.
- **MT4 residual/protected evidence - active.** No policy promotion occurred.
- **MT5 documentation/authority maintenance - maintenance.** Broad migration is
  complete; the live ledger is now bounded and historical evidence task-gated.
- **MT6 carrying-cost control - active.** Mandatory manager reading is reduced
  without deleting evidence or adding a tooling layer.

## Cartesian Hamiltonian Producer Pass 401 - Reject Nonfinite Terminal Coefficients

Commit(s):
- `af44e43b2` - fail before terminal Gaussian-sum one-body assembly can consume
  a nonfinite coefficient.

Summary:
- The shared terminal Gaussian-sum accumulator now rejects `NaN`, positive
  infinity, and negative infinity after coefficient conversion and before
  factor preparation or destination mutation. This closes the Pass 398
  `HP-FN-03` discrepancy for base unit-nuclear, residual-GTO, and mixed-Hartree
  callers without duplicating checks at each caller.
- Finite signed coefficients remain valid. No positivity, normalization,
  magnitude, workspace, factor, Coulomb-policy, or matrix convention changed.
  Terminal IDA retains its separate existing guard at its distinct entry
  boundary.

Validation / evidence:
- Package load passed in `0.45 s`; the public Cartesian base test passed
  `134/134` in `39.5 s`; malformed inputs left a prefilled destination exactly
  unchanged. Omitted and explicit compact matrices remained exactly equal, and
  existing compact/high H/H2 endpoints, artifacts, finiteness, and symmetry
  checks passed. Compact/high due diligence retained bounds
  `[-2.46794,2.46794]^3`, axes `7x7x7`, dimension `79`, direct-core plus
  complete-shell topology, and no warning flags. `git diff --check` passed.

Goal advancement / guardrail:
- MT1 closes its second fail-fast item with only `+2` source and `+15` test
  lines. The next immediate conformance target is direct-injection
  negative-Gram/final-identity validation, followed by protected one-body
  diagnostics; neither should be combined with a residual or screening policy
  change.

## Cartesian Hamiltonian Producer Pass 402 - Enforce The Live Ledger Bound

Commit(s):
- this commit - fail documentation builds when the live manager ledger exceeds
  `2,000` lines.

Summary:
- Added one evidence-policy checker, separate from Cartesian execution
  authority, and invoked it before Documenter setup. The live log still rotates
  by manager judgment, but ignored rotation can no longer pass docs CI.
- Focused fixtures prove that exactly `2,000` lines pass and `2,001` fail. No
  authority record, generated view, subsystem contract, producer source,
  artifact, workflow semantics, or numerical behavior changed.

Validation / evidence:
- The standalone live check, focused docs tests, authority check/self-test,
  local Documenter, and `git diff --check` passed. The actual live ledger
  remains well below its ceiling.

Goal advancement / guardrail:
- MT5 gains automatic carrying-cost enforcement without coupling historical
  evidence to `authority.toml` or adding a manager-log generator. Archive
  timing and content selection remain reviewed manager decisions; CI enforces
  only the hard cognitive-load ceiling.

## Cartesian Hamiltonian Producer Pass 403 - Validate Direct Injection Geometry

Commit(s):
- `8d7c8332f` - reject indefinite injected-mode geometry and validate the
  complete implicit replacement sector.

Summary:
- The default-off direct-`G` compatibility path now checks the global
  injected-mode Gram spectrum before rank cleanup. Materially negative
  directions can no longer disappear when the positive subspace is retained.
- After global cleanup, the implementation validates the three blocks of
  `F' S F`: `Y' S_AA Y`, `B' Q_perp`, and `Q_perp' Q_perp`. It uses the
  existing scale-aware `identity_atol`, handles an empty injected sector, and
  creates no persistent dense transform, diagnostic field, or new caller.
- Manager review moved this validation immediately after fixed-sector
  construction; the first draft checked only after residual construction had
  already consumed the replacement geometry.

Validation / evidence:
- Package load passed in `7.0 s`; the augmented H2 test passed `155/155` in
  `47.0 s`, and the supplemented facade passed `69/69` in `8.0 s`. A synthetic
  Gram with eigenvalues `[-1,3]` rejects; a valid `nG=3` case has one injected
  and one residual mode, zero blockwise `F' S F` and `F' S R` errors at printed
  precision, and `R' S R` error `4.44e-16`. Existing dimensions `487/18/505`,
  MWG self-Coulomb `0.4574161883692301`, readback parity, and due-diligence
  warnings remain unchanged. `git diff --check` passed. The pass is
  `+44/-1` lines.

Goal advancement / guardrail:
- MT1 closes the direct-injection fail-fast discrepancy without promoting this
  path or changing residual policy. The next immediate conformance target is
  the missing protected one-body finite, dimension, symmetry, and geometry
  diagnostics; protected physics and interaction semantics remain unchanged.

## Cartesian Hamiltonian Producer Pass 404 - Validate Protected One-Body Transforms

Commit(s):
- `3e0df456a` - validate protected fixed-sector operator inputs and return
  truthful pre-symmetrization diagnostics.

Summary:
- Exact protected fixed-sector kinetic, per-center unit-nuclear, and assembled
  `H1_F` construction now rejects malformed dimensions, nonfinite raw blocks,
  nonfinite charges, and nonfinite transformed matrices at its owned boundary.
  The existing matrix-returning operator helper and charge-weighted one-body
  convention remain unchanged.
- One compact diagnostics record reports dimensions, raw and transformed
  symmetry errors, traces, and the four existing same-geometry acceptance
  facts. Symmetry is measured before roundoff cleanup can hide it; no new
  numerical threshold, low-spectrum production solve, caller, or persistent
  artifact shape was introduced.
- Manager review added a post-symmetrization finite check because averaging two
  finite extreme values can overflow even when the pre-cleanup matrix is
  finite.

Validation / evidence:
- A transient synthetic smoke passed `25/25`, including malformed GG/GA/AA
  dimensions, nonfinite blocks and charges, finite-input symmetrization
  overflow, exact legacy matrix parity, and diagnostic values. The augmented
  H2 test passed `155/155` in `46.6 s`; the supplemented facade passed `69/69`
  in `7.9 s` with dimensions `487/18/505`, self-Coulomb and artifact readback
  parity unchanged. Due diligence retained axes `9x9x15`, padding `4`, rows
  `[275,114,98]`, and its existing advisory warnings. `git diff --check`
  passed. Source carrying cost is `+65/-9`; no committed test was authorized or
  added.

Goal advancement / guardrail:
- MT1 closes the final immediate fail-fast item identified in Pass 398. The
  next conformance work should move to consumer-data correctness, beginning
  with protected ladder readback/trace-loss or due-diligence warning shapes;
  protected physics, interaction semantics, and Cr2 measurement inputs remain
  unchanged.

## Cartesian Hamiltonian Producer Pass 405 - Reconcile Protected Ladder Readback

Commit(s):
- `c46b9e9ff` - return the complete written v1 manifest facts and compute
  restart trace loss from the actual source density trace.

Summary:
- Protected ladder readback now exposes the geometry, electron-count, and
  basis facts already written by the manifest: symbols, charges, locations,
  spin counts, `core_spacing`, `s_factor`, basis name, and `lmax`. Original v1
  bundles without `s_factor` retain their historical `1.0` default, while
  missing Coulomb provenance remains unavailable rather than inferred.
- Alpha and beta transfer loss are now each `source_trace - target_trace`.
  This remains signed, so a transfer-induced norm increase is visible instead
  of being hidden by occupied-column count. The obsolete private `nup/ndn`
  aliases were removed. No writer key, schema, transfer matrix, restart matrix,
  fixed-density formula, or sidecar convention changed.

Validation / evidence:
- The `/tmp` readback/trace smoke passed `37/37` in `1.6 s`, including exact
  current-manifest roundtrip with nonunit `s_factor`, two legacy manifests,
  absent Coulomb provenance, nonorthonormal source orbitals, and a negative
  beta trace loss. Package load and `git diff --check` passed. No Hamiltonian
  was built, so terminal due diligence was not applicable. Source carrying
  cost is `+18/-7`, with two stale private fields deleted and no committed test
  added.

Goal advancement / guardrail:
- MT1 closes the first consumer/readback discrepancy from Pass 398. The next
  bounded correctness target is due-diligence warning-shape reconciliation;
  ladder construction, transfer physics, solver behavior, and the controlled
  Cr2 inputs remain unchanged.

### Medium-Term Goal Checkpoint After Pass 405

- **MT1 conformance remediation - active.** The four immediate fail-fast items
  and protected ladder readback/trace-loss are closed. Continue with warning
  shape, stale compatibility paths, and truthful missing validation claims.
- **MT2 controlled Cr2 measurement - active.** Passes 401-405 changed no basis,
  imported state, screening convention, or solver interpretation.
- **MT3 pending producer facilities - active.** Standard60, canonical-driver
  Coulomb exposure, and retained-GTO EGOI remain separate authorized lanes.
- **MT4 residual/protected evidence - active.** No injection, localization,
  cutoff, or protected-basis policy was promoted by conformance work.
- **MT5 documentation/authority maintenance - maintenance.** The generated
  authority views and bounded manager ledger remain current; no new migration
  campaign is needed.
- **MT6 carrying-cost control - active.** This pass deletes obsolete aliases;
  subsequent reconciliation should continue removing stale paths rather than
  adding adapters or compatibility vocabulary.

## Cartesian Hamiltonian Producer Pass 406 - Reconcile Due-Diligence Warning Shapes

Commit(s):
- `9faece815` - make row-level warning fields match the canonical reporting
  contract.

Summary:
- Terminal due-diligence rows now carry `warning_flags` as the ordered tuple of
  active stable symbols, using `(:none,)` when clean, and `warning_summary` as
  concise comma-joined text. Boolean warning predicates remain local to row
  construction. The canonical driver consumes the text directly and deletes
  its redundant tuple formatter.
- The report-level `due.warnings` boolean `NamedTuple` remains a separate
  unchanged contract. No warning predicate, vocabulary, threshold, row order,
  shellification fact, retained count, source shape, matrix, artifact, or
  numerical construction behavior changed.

Validation / evidence:
- The doer package load passed. A `/tmp` row-contract smoke exercised one clean
  row and two warned rows in `18.34 s`; all flags were symbols and every summary
  exactly matched its ordered flags. The canonical H2 driver completed in
  `22.11 s`, wrote and read its artifact, and printed concise warning text.
  Due diligence retained bounds `x/y = +/-4.87075`, `z = +/-6.10974`, axes
  `9x9x15`, dimension `487`, rows `[275,114,98]`, direct atom-contact core plus
  two complete shells, and unchanged report warnings. Manager package load and
  `git diff --check` passed. The pass is `+8/-6`, net `+2` lines.

Goal advancement / guardrail:
- MT1 closes the Pass 398 due-diligence consumer-shape discrepancy. The fix
  removes one duplicate formatter and adds no compatibility representation.
  Future COMX nuclear-resolution diagnostics, warning vocabulary, and all
  basis or screening policy remain separate authority questions.

## Cartesian Hamiltonian Producer Pass 407 - Reconcile Compact WL And Thin-Slab Semantics

Commit(s):
- this commit - remove stale identity metadata from compact White-Lindsey and
  thin-slab records.

Summary:
- White-Lindsey boundary contracts, retained units, and transform contracts
  now consistently use the existing
  `:white_lindsey_boundary_stratum_product` vocabulary. Identity embedding is
  reserved for true direct/core units. Obsolete direct-identity slab kinds,
  compatibility branches, and counters were removed, together with one
  caller-free lowering-count helper.
- This is a semantic and carrying-cost correction, not a basis change. The
  active terminal realizer already built compact products; the repaired
  records now tell the same story as its coefficients.

Validation / evidence:
- The doer public base gate passed `134/134`. Manager reran the source-backed
  H2/WL parity check: dimension `471`, one direct block, `52` WL boundary
  blocks, two compact slabs, exact pre/post `H1` and `Vee`, identical ordered
  supports and column ranges, and maximum compact Gram error `4.77e-15`.
  Bounds were `x/y = +/-4.87075`, `z = +/-6.77742` on axes `9x9x15`; native
  inventory confirmed the direct core, two shells, and two angular-extension
  slabs. Focused stale-symbol scans and `git diff --check` passed.

Goal advancement / guardrail:
- MT1 closes the WL/thin-slab stale-path discrepancy and MT6 gains `+5/-75`
  source lines, net `-70`, without compatibility glue. WL row-level
  due-diligence remains empty even though native terminal inventory is
  complete; treat that as a separate reporting gap, not permission to reopen
  realization, shellification, or numerical policy.

## Cartesian Hamiltonian Producer Pass 408 - Authorize QW/High-Order Cluster Retirement

Commit(s):
- this commit - approve deletion of the obsolete four-file QW receipt,
  chain/square consumer, and high-order doside/IDA cluster.

Summary:
- Caller review confirmed that the eight top-level experimental exports have
  no committed source or test consumers, and that the 24-name carried-space
  receipt submodule is consumed only by the retiring chain/square file. The
  current `nesting=:wl` producer uses preserved adjacent kernels.
- Added `HP-RETIRE-QW-DONOR-FN-01` and `HP-RETIRE-QW-DONOR-TEST-01` with
  retirement-only authority for the four `6,008`-line files, their
  `GaussletBases.jl` includes/exports/generics, and bounded existing validation
  files. No compatibility stubs, aliases, moved implementations, or new tests
  are authorized.
- Compressed the active receipt and high-order plans into historical evidence.
  The doside/tensor-shell, FSB/FBU, He controls, distorted-parent uncertainty,
  Cr capture, and migrated endcap/panel lessons remain recorded. Atomic chains
  remain a future scientific target requiring a new design.

Validation / evidence:
- Focused committed-caller and preserved-owner scans, authority render/check
  and self-test, generated registry/AGENTS parity, local Documenter,
  manager-log bound, docs-only staged scope, and `git diff --check` form the
  gate. No source or test file changes in this pass.

Goal advancement / guardrail:
- MT6 now has a deletion-ready, fail-closed boundary for approximately `6,000`
  obsolete source lines. Pass 409 must delete the complete cluster together,
  preserve chain/square basis constructors and current WL behavior, inspect
  WL terminal due diligence, and stop if any live caller appears.

## Cartesian Hamiltonian Producer Pass 409 - Retire QW/High-Order Experimental Cluster

Commit(s):
- this commit - delete the obsolete QW carried-space, chain/square operator,
  high-order doside, and experimental IDA implementation cluster.

Summary:
- Deleted the four authorized files in full and removed their four includes,
  eight exports, and six empty generic declarations from `GaussletBases.jl`.
  The source change is `+0/-6,026`; no compatibility aliases, replacement
  implementation, or tests were added.
- Focused caller scans found no committed source or test consumer of the
  retired surface. Surviving chain/square basis constructors and geometry
  diagnostics remain covered by core tests, while the current
  `nesting=:wl` producer continues to use adjacent preserved kernels.
- Closed the machine-authority lifecycle in the same accepted batch so deleted
  paths are not recorded as existing: the function ID is now `retired`, the
  test ID is `completed`, and both have no grant, surface, path ownership, or
  execution-whitelist entry. This is closure, not new authority.

Validation / evidence:
- The bloat-fixer package load passed, the core group passed `432` checks, and
  the public Cartesian base gate passed `134/134`. Deleted-symbol and
  preserved-owner scans plus `git diff --check` passed. Manager review
  confirmed the exact deletion boundary and clean staging separation; package
  load, authority check/self-test, docs tests `56/56` and `10/10`, the
  manager-log bound, and local Documenter also passed.
- WL due diligence retained the one-center `ns=5`, derived-`q=3` endpoint:
  axes `7x7x7`, dimension `223`, one direct core plus one complete shell, and
  native retained counts `125 + 54 + 36 + 8`. The formal WL row table remains
  empty; its complete 27-row native inventory remains the equivalent audit.

Goal advancement / guardrail:
- MT6 closes a coherent `6,008`-line donor retirement and removes another 18
  lines of public/package surface. Atomic chains remain a future scientific
  target, not a preserved experimental API. Active WL, adjacent geometry,
  ordinary-QW kernels, and `cartesian_carried_spaces.jl` remain separately
  owned; the empty WL due-diligence row table is still a reporting gap.

## Cartesian Hamiltonian Producer Pass 410 - Authorize Carried-Space Adapter Retirement

Commit(s):
- `839395f63` - completed the Pass 409 QW/high-order cluster retirement.
- this commit - approve deletion of the orphaned internal carried-space
  adapter and its sole root include.

Summary:
- Caller and history review found that `src/cartesian_carried_spaces.jl` has no
  committed source or test consumer. It was introduced by `e0ca22c2d`, gained
  its sole production consumer in `231331ff8`, lost standalone tests in
  `bc425ce67`, and became orphaned when Pass 409 removed that consumer. The only
  remaining caller is an ignored stale Dropbox conflicted test copy.
- Added `HP-RETIRE-CARRIED-SPACE-FN-01` and
  `HP-RETIRE-CARRIED-SPACE-TEST-01`. The follow-on pass must delete the complete
  `266`-line adapter and one include, add no compatibility surface or tests,
  and preserve current representation, overlap, transfer, parent, hybrid,
  chain/square, and WL/PQS owners.

Validation / evidence:
- Focused caller/history/path scans, authority render/check and self-test,
  generated registry/AGENTS parity, docs tests, local Documenter, manager-log
  bound, docs-only staged review, and `git diff --check` form this pass's gate.
  No source or test file changes are authorized here.

Goal advancement / guardrail:
- MT6 remains active for one final `267`-line deletion. Qualified access to the
  internal unadvertised submodule is not a compatibility obligation. Current
  physics, WL/PQS behavior, ordinary-QW capability review, and EGOI work are
  unaffected.

### Medium-Term Goal Checkpoint After Pass 410

- **MT1 conformance remediation - active.** This retirement does not alter
  correctness findings, warning contracts, or endpoint interpretation.
- **MT2 controlled Cr2 measurement - active.** Current Hamiltonians, imported
  states, screening convention, and physics targets are unchanged.
- **MT3 pending producer facilities - active.** Standard60, canonical-driver
  Coulomb exposure, and retained-GTO EGOI remain independent lanes.
- **MT4 residual/protected evidence - active.** No residual, injection,
  protected-basis, or interaction policy changes.
- **MT5 documentation/authority maintenance - maintenance.** The amendment is
  bounded to exact retirement authority and generated views.
- **MT6 carrying-cost control - active.** Pass 409 is complete; delete the
  final orphaned adapter next, then close this retirement sublane. Ordinary-QW
  capability review remains separate and must begin from surviving callers.

## Cartesian Hamiltonian Producer Pass 411 - Retire Cartesian Carried-Space Adapter

Commit(s):
- this commit - delete the orphaned internal carried-space adapter and its sole
  root include.

Summary:
- Deleted `src/cartesian_carried_spaces.jl` in full and removed only its
  `GaussletBases.jl` include: `+0/-267` source lines. No alias, stub,
  deprecation, compatibility module, replacement adapter, or test was added.
- Pre/post scans confirmed that the retired module, type, constructor, and five
  accessors had no committed source or test consumer. The one ignored reference
  remains a stale Dropbox conflicted test copy and is not a compatibility
  obligation.
- Closed the machine records in the same reviewed batch: the function ID is
  `retired`, the validation ID is `completed`, and both now have no grant,
  surface, path ownership, dependency, or whitelist entry.

Validation / evidence:
- The bloat-fixer package load passed and the unchanged core group passed
  `432/432`; deleted-symbol scans and `git diff --check` passed. Manager review
  confirmed the exact two-file source boundary and preserved owner files.
  Authority check/self-test, generated-view parity, docs tests, manager-log
  bound, and local Documenter form the closure gate.

Goal advancement / guardrail:
- MT6 completes the carried-space/QW donor retirement sublane with no behavior
  change. Representation, overlap, transfer, parent, hybrid, chain/square,
  ordinary-QW, WL, and PQS ownership remains unchanged. The next carrying-cost
  decision is the separate ordinary-QW endpoint/correction capability audit;
  it must begin from live shared-kernel callers rather than file-level labels.

## Cartesian Hamiltonian Producer Pass 412 - Remove Dead QW Midpoint Blocks

No strategic change. Deleted the caller-free midpoint cross-block route and
its four now-orphaned sampling/support helpers from
`ordinary_qw_raw_blocks.jl`, reducing source by `132` lines while preserving
all live `_qwrg_*` kernels and current PQS/WL callers. Package load, the core
group (`432/432`), deleted-symbol scans, and `git diff --check` passed; MT6 now
returns to the separate public endpoint and correction-capability decisions.

### Strategic Direction After Pass 412 - Defer Broader Ordinary-QW Retirement

The ordinary-QW capability audit is retained as an ownership and carrying-cost
map, not as an instruction to begin broad deletion. Near-term work should
enhance the capabilities wanted in the current PQS/WL line. That evidence will
show which old correction, hydrogenic-core/ESOI, branch/fragment,
counterpoise, chain/square, residual-oracle, and public-endpoint ideas should be
reused, independently reimplemented, historically preserved, or retired.

Do not keep old code merely because a concept may matter later, but do not
force a replacement or retirement decision before current physics work makes
the ownership clear. The completed dead midpoint cleanup stands; shared
`_qwrg_*` kernels remain protected, and no broader ordinary-QW source or API
retirement should start without a new user-directed decision.

### Strategic Clarification After Pass 412 - Parent-Supported Residual Pieces

Shell-local parent-completion functions should not be treated as new terminal
IDA sites. They are intended to split broad residual content into more local
pieces that lie in omitted parent-shell spans. Their one-body matrices can use
the exact parent transform, while their near-zero linear IDA weights are normal
for residual-like functions. The proposed `abs2` parent-coefficient charge
model is therefore a candidate replacement for MWG only on interaction rows
involving these additions; it must be compared with current MWG and a bounded
exact Coulomb oracle. Existing base-base interaction entries remain unchanged.

This is a partial Cr2 basis/interaction improvement, not a substitute for
screened Hartree or later correction physics. A central efficiency target is
to lower selected shell-local source order from the current `q=7` level toward
`q=5` or `q=6` on the same adequate parent lattice, then add fewer localized
parent-completion modes than the ordinary contraction columns removed. Use
`source_q` or retained-order language for that shell-local experiment; public
global `ns` and parent resolution are separate unless explicitly varied. No
producer or interaction authority follows from this clarification; fixed-
density, dimension, locality, MWG-versus-parent-charge, and exact-oracle
evidence comes first.

## Cartesian Hamiltonian Producer Pass 413 - Authorize Shell-q Coarsening

Commit(s):
- this commit - amend the implemented semantic per-shell PQS source-order
  authority to permit bounded coarsening as well as refinement.

Summary:
- Source and history review confirmed that `source_q > route_q` was the narrow
  original refinement question, not a numerical requirement. The existing
  post-shellification path already rewrites one authoritative shape, preserves
  symmetric atom ownership, and reruns the shared-shell longitudinal selector.
- `HP-PQS-SHELLQ-OVERRIDE-FN-01` now permits non-Boolean integer
  `source_q >= 3` with `source_q != route_q`. Values below route q coarsen only
  the selected shell contraction; global `ns`, route q, parent axes, support,
  ownership, cores, slabs, and route metadata remain unchanged. Equal values
  remain redundant errors.
- The amended test gate requires route-q `7` to source-q `6` and `5` retained-
  count reduction, orthonormal columns, full finite/symmetric construction,
  exact omitted/empty parity, and strict malformed/below-3/equal/asymmetric/
  unmatched rejection. No lowering-plan accessor or parent-completion mode is
  approved.

Validation / evidence:
- Canonical-contract and source-signature review, authority render/check and
  self-test, generated registry/AGENTS parity, docs tests, local Documenter,
  manager-log bound, docs-only staged review, and `git diff --check` form the
  gate. No source or test file changes occur in this pass.

Goal advancement / guardrail:
- MT4 gains the bounded contraction-coarsening diagnostic needed before any
  parent-completion design. MT2 may consume reviewed static CR2 measurements,
  but parent completion, interaction changes, HF claims, public controls, and
  durable accessors remain separate authority.

## Cartesian Hamiltonian Producer Pass 414 - Implement Shell-q Coarsening

Commit(s):
- this commit - allow fixed-parent semantic shell source order below route q.

Summary:
- The existing override normalizer now accepts non-Boolean integer
  `source_q >= 3` with `source_q != route_q`; the complete post-shellification
  shape-rewrite path required no accessor or downstream repair. Route-q `7`
  H2 shared-shell dimensions changed `855 -> 789 -> 735` for source q
  `7 -> 6 -> 5`, while the atom-local pair changed
  `1759 -> 1627 -> 1519`.
- A padded Be2 route-q `7`, source-q `6` gate changed two shell rows from
  `218` to `152`, saving `132` base functions and producing dimensions
  `4335 + 42 = 4377`. Parent geometry, support, ownership, unmatched regions,
  packet traces/capture, screened-Hartree anchors, and due-diligence warning
  classes remained valid.
- Manager review replaced two brittle fixed `1e-10` residual-identity test
  assertions with the established scale-aware `1e-7` ceiling. The measured H2
  errors remained near `1e-11`; the Be2 error was `1.63e-9`.

Validation / evidence:
- The doer package load and padded Be2 gate passed. The focused H2 test passed
  `288/288` and the facade regression passed `69/69`; manager reran both after
  the test-bound correction. `git diff --check` and manager-log bounds passed.

Goal advancement / guardrail:
- MT4 now has the fixed-parent coarsening control needed for CR2 occupancy/H1
  and later parent-completion measurements. Source change is only `+2/-2`;
  focused tests add `137` net lines. This pass does not authorize completion
  modes, parent-charge interactions, HF, public controls, or Cr2 claims.

## Cartesian Hamiltonian Producer Pass 415 - Authorize Parent Residual Functions And Gaussian Direct Blocks

Commit(s):
- this commit - approve internal PRF mechanics and direct-only parent-backed
  Gaussian interaction blocks.

Summary:
- Added one canonical contract and four execution IDs. Consumers retain shell,
  target, orientation, mode-count, and physical-state policy; repo source may
  project and validate supplied parent targets, build exact PRF one-body blocks,
  and evaluate tiled PRF-PRF/PRF-G direct blocks without changing existing G-G
  operators or interactions.
- The Gaussian resource uses mapped parent centers and same-expansion positive
  parent-IDA onsite values. It is explicitly direct-only. Transition-density
  exchange, PRF-to-GTO-residual interactions, Hamiltonian integration,
  artifacts, public controls, and Cr2-specific behavior remain unapproved.
- Reconciled semantic shell-q coarsening as implemented and completed its test
  lifecycle after Pass 414.

Validation / evidence:
- Reviewed the live terminal realization, factorized one-body/IDA owners, the
  accepted shell-q implementation, and the CR2 occupancy, q-ladder, full-parent
  IDA, continuum, Gaussian-distance, and complete 16-by-6783 PRF-G reports.
  Authority render/check and self-test, generated registry/AGENTS parity, docs
  tests, local Documenter, manager-log bound, staged docs-only scope, and
  `git diff --check` form the amendment gate. No producer source or numerical
  test changes occur in this pass.

Goal advancement / guardrail:
- MT4 advances from measurement-only PRFs to a bounded implementation seam
  because direct evidence is now sufficient. The ongoing orbital-contracted
  exchange audit may quantify a later need but cannot broaden this direct
  authority. The two compact internal objects are the only new persistent
  source shapes approved; no file/module or metadata cloud is authorized.

### Medium-Term Goal Checkpoint After Pass 415

- **MT1 conformance remediation - active.** PRF authority does not reinterpret
  or close any remaining Pass 398 discrepancy.
- **MT2 controlled Cr2 measurement - active.** Screened-Hartree and fixed-parent
  PRF studies remain separate controlled comparisons with consumer-owned state
  selection and interpretation.
- **MT3 pending producer facilities - active.** PRF mechanics/Gaussian direct,
  Standard60/driver exposure, and retained-GTO EGOI are independent pending
  source lanes.
- **MT4 residual/protected evidence - active with refinement.** PRF exact
  one-body and direct blocks are approved; orbital-contracted exchange and
  PRF-to-GTO-residual interactions remain measurement-only.
- **MT5 documentation/authority maintenance - maintenance.** This is a bounded
  atomic machine-authority amendment, not a documentation migration.
- **MT6 carrying-cost control - active.** Implementation must reuse existing
  terminal/Coulomb owners, add no file or module, and report any helper or test
  pressure it cannot avoid.

## Cartesian Hamiltonian Producer Pass 416 - Implement Parent Residual Functions And Direct Blocks

Commit(s):
- `5b46ae073` - add internal PRF mechanics, exact one-body blocks, and
  parent-backed Gaussian/direct-oracle infrastructure.

Summary:
- Consumer-supplied support-local parent targets now produce explicit PRF
  blocks without dropping columns or changing the terminal `G` basis. Exact
  kinetic, per-center unit-nuclear, assembled `H1`, position, and second-moment
  `G-R`/`R-R` blocks support multiple owner/shell PRF blocks.
- Added the in-memory onsite-calibrated Gaussian resource and tiled PRF-G/
  PRF-PRF direct blocks, plus the bounded full-parent IDA comparator. Manager
  review added collection-wide terminal/metric validation, cross-block
  one-body terms, exact Coulomb-expansion identity checking, and direct-block
  symmetry guards.

Validation / evidence:
- Package load passed. Core tests passed `440/440`; the final focused H2 gate
  passed `358/358` and the supplemented facade passed `69/69`. Exact one-body
  oracle errors remained at roughly `1e-15`, existing `H1`/`Vee` stayed
  unchanged, authority check and `git diff --check` passed, and terminal due
  diligence was inspected. No Cr2 fixture, HF, or endpoint claim was added.

Goal advancement / guardrail:
- MT4 now has source-backed PRF mechanics and direct-only interaction
  resources. Selection policy, transition exchange, PRF-GTO interactions,
  Hamiltonian integration, artifacts, and public controls remain outside the
  lane. The change is `+976/-2` across seven source/test files; no new file or
  module was introduced. Machine lifecycle/current-status reconciliation is
  the remaining documentation closure step.

## Cartesian Hamiltonian Producer Pass 417 - Correct Angular-Style Injection Interpretation

Commit(s):
- this commit - reconcile the repo injection descriptions with the primary
  angular-injection construction.

Summary:
- Clarified that angular-style injection is not complete after constructing
  the replacement span `F = Y + (G intersect Y_perp)`. Every old localized
  `G` seed is projected into `F`, and the projected seeds are symmetrically
  Lowdin-orthogonalized to obtain the final injected localized basis.
- Distinguished span-only occupied/direct-G helpers from the complete
  relocalization already used by protected-localized construction. No source
  lifecycle, numerical threshold, or execution grant changes in this pass.

Validation / evidence:
- Compared the formulas with the local primary `Injection.tex`, including its
  clean geometric definition and practical project-then-Lowdin construction.
  The coefficient-space nullspace form implements the stated
  `G intersect Y_perp` exactly and preserves the full-target-rank stop rule.
- Authority document hashes, generated-view parity, local Documenter, docs
  tests, manager-log bound, and `git diff --check` form the acceptance gate.

Goal advancement / guardrail:
- MT4 gains an explicit construction invariant before parent-backed injected
  composition is authorized. This is a documentation correction only: it does
  not wire occupied-first injection, add a new injection policy, or authorize
  Hamiltonian, artifact, solver, or Cr2-specific behavior. Net carrying cost
  is limited to the operational refresher and focused canonical clarifications.

## Cartesian Hamiltonian Producer Pass 418 - Authorize Parent-Backed Injected Composition

Commit(s):
- this commit - approve private fixed-span injection, mixed interaction
  assembly, and screened-Hartree delegation.

Summary:
- Added independent implementation/test authority for composing
  consumer-selected parent targets into angular-style injected terminal rows,
  retaining the exact parent-backed complement, and residualizing the explicit
  GTO supplement against that complete span. The old localized terminal seeds
  must be projected before symmetric Lowdin; dimension and span are invariant.
- The interaction contract keeps categories explicit: rebuilt terminal IDA,
  parent-Gaussian direct blocks, terminal/residual MWG, and moment-derived
  parent-residual/external-residual MWG. Screening is rebuilt only after the
  complete native `B`, `H1`, and `Vee` are fixed.
- Reconciled the four PRF/direct IDs as implemented/completed under
  `5b46ae073` and Pass 416.

Validation / evidence:
- Reviewed committed PRF/direct APIs and bounded tests plus the July 15 q7 Cr2
  report. That report supports the mechanics but remains consumer measurement:
  q7 fixed-span interaction error improved while q6 remained unsuitable, so no
  shell, source-q, cutoff, PRF count, localization, or Cr2 default is promoted.
- Authority render/check/self-test, generated-view parity, docs tests, local
  Documenter, manager-log bounds, staged docs-only review, and
  `git diff --check` form the amendment gate.

Goal advancement / guardrail:
- MT4 advances to a reusable composition seam. Selection policy, exact
  PRF-GTO interactions, transition exchange, public workflow, artifacts,
  solvers, and production endpoints remain outside authority. No source or
  numerical test changes occur in this pass; one compact internal result object
  is the maximum newly approved persistent shape.

## Cartesian Hamiltonian Producer Pass 419 - Implement Parent-Backed Injected Basis Composition

Commit(s):
- `cdd2c27af` - compose the private fixed-span parent-backed injected basis.

Summary:
- Implemented native `B = [Ginj,Rnew,RGexternal]` composition without changing
  the parent-backed span or dimension. Every old terminal seed is projected
  into the injected span before symmetric Lowdin, the exact parent complement
  is retained, and the explicit Gaussian supplement is residualized against
  the complete parent-backed basis through the existing numerical-complete
  builder.
- Exact parent-backed kinetic, nuclear, physical one-body, position, and
  second-moment blocks now feed the existing augmented-operator path. Manager
  review removed an independently supplied one-body bypass, tied augmentation
  overlaps to their raw blocks, rejected center/block truncation, and required
  all consumers to revalidate actual terminal and PRF geometry.

Validation / evidence:
- Manager reran the focused H2 gate (`417/417`) and supplemented facade
  (`69/69`). Target recovery, old/new span singular values, metric identities,
  exact one-body oracle parity, and omitted-path parity passed. The doer also
  completed the padded Be2 gate at dimensions `1729 + 2 + 42 = 1773`; packet
  recovery and native ordering passed, with only the established shell/slab
  due-diligence warnings. Machine authority and `git diff --check` passed.

Goal advancement / guardrail:
- MT4 now has the reusable fixed-span basis composition required by the
  successful consumer recipe. The source/test change is `+791/-22`; it adds no
  file, module, export, artifact, solver, interaction, screening, or Cr2
  policy. The next bounded lane is separated native interaction assembly and
  screened-Hartree delegation under the existing interaction IDs.

## Cartesian Hamiltonian Producer Pass 420 - Implement Parent-Backed Injected Interaction

Commit(s):
- `006432e9d` - assemble the private category-owned interaction and delegate
  additive screened Hartree.

Summary:
- The native `[Ginj,Rnew,RGexternal]` basis now receives a freshly rebuilt
  interaction: terminal IDA for `Ginj-Ginj`, the onsite-calibrated parent
  Gaussian model for blocks involving `Rnew`, and the established moment-
  matched Gaussian convention for blocks involving external residuals. No old
  interaction is copied or rotated.
- Packet-backed construction represents each atomic occupied block in the
  complete basis, evaluates the ordinary fitted atomic fields there, and calls
  the existing additive screened-Hartree API only after `H` and `V` are fixed.
  The correction remains separate and leaves both matrices unchanged.
- Manager review bound all six moment matrices to the actual residual/raw-block
  construction and removed retained category matrices that would have
  duplicated roughly one full terminal interaction at Cr2 scale.

Validation / evidence:
- Manager reran the focused H2 gate (`434/434`) and supplemented facade
  (`69/69`). Every native interaction slice matched its lower-level owner,
  stale external moments failed, fresh terminal rebuilding differed from the
  old matrix as expected, and legacy facade matrices/readback remained exact.
- The doer reran the padded Be2 gate at dimensions `1729 + 2 + 42 = 1773`.
  Packet traces/capture, fitted-potential diagnostics, additive decomposition,
  derivative anchor, correction recomputation, and unchanged `H/V` passed.
  Terminal due diligence showed only the established shell/slab warnings.
  Machine authority and `git diff --check` passed.

Goal advancement / guardrail:
- MT4 now has the complete private source mechanics needed for controlled CR2
  consumer validation. The change is `+366/-9`; it adds no type, export, file,
  artifact, solver, exchange correction, or physical selection policy. PRF
  targets, shell/source order, RDM, cutoffs, orientations, and interpretation
  remain consumer-owned.

### Medium-Term Goal Checkpoint After Pass 420

- **MT1 conformance remediation - active.** This physics facility does not
  close or reinterpret remaining Pass 398 discrepancies.
- **MT2 controlled Cr2 measurement - active.** The next step is a consumer-side
  fixed-state validation using the source-backed basis, interaction, and newly
  rebuilt correction before any relaxed HF interpretation.
- **MT3 pending producer facilities - active.** Standard60/driver exposure and
  retained-GTO EGOI remain independent work.
- **MT4 residual/protected evidence - source milestone complete, validation
  active.** PRF construction, fixed-span injection, exact one-body assembly,
  category-owned density interaction, and correction delegation are now
  implemented. Automatic selection, transition-density exchange, exact
  PRF-GTO interactions, and production defaults remain out of scope.
- **MT5 documentation/authority maintenance - maintenance.** Reconcile the four
  composition/interaction IDs to their source-backed lifecycle in one bounded
  docs pass after this implementation commit.
- **MT6 carrying-cost control - active.** Review removed duplicate Cr-scale
  component matrices and added no compatibility carrier. The net source/test
  increase is justified by the current controlled physics target; broader API
  or artifact growth remains forbidden.

## Cartesian Hamiltonian Producer Pass 421 - Reconcile Parent-Backed Composition Lifecycle

Commit(s):
- `cdd2c27af` - implemented fixed-span parent-backed basis composition.
- `006432e9d` - implemented separated interaction and screened-Hartree
  delegation.
- this commit - reconcile machine authority and current status.

Summary:
- Marked both source IDs implemented with maintenance permission and both test
  IDs completed with maintenance permission. Maintenance ownership is narrowed
  to the six source files actually carrying the facility and the single
  committed H2 test; unused prospective owners were removed.
- The canonical contract now records native `[Ginj,Rnew,RGexternal]` order,
  projected old seeds, fixed parent span/dimension, numerical-complete external
  residuals, exact one-body reconstruction, category-owned interaction blocks,
  represented packet fields, and correction delegation after `H`/`V` are fixed.

Validation / evidence:
- Passes 419-420 and the implementation diffs support every lifecycle claim.
  Manager rerun evidence is H2 `434/434` plus facade `69/69`. Padded Be2 at
  `1729 + 2 + 42 = 1773` passed packet/correction accounting and unchanged
  `H`/`V`, but remains external evidence rather than a committed slow fixture.
- Authority check/self-test, generated-view parity, docs tests, local
  Documenter, manager-log bound, staged docs-only review, and
  `git diff --check` form this pass's gate.

Goal advancement / guardrail:
- MT4 source work is complete; MT2 now owns consumer-side CR2 fixed-state
  validation. Selection, transition exchange, exact PRF-GTO interaction,
  public workflow, artifacts, solvers, and endpoint acceptance remain outside
  authority. This pass changes no source, test, numerical, or workflow behavior.

## Cartesian Hamiltonian Producer Pass 422 - Close Source-Backed CR2 Validation

Commit(s):
- `da607b855` - source-backed lifecycle baseline.
- this commit - record fixed-state and bounded replay evidence.

Summary:
- The repo construction reproduced the former CR2-local native basis and
  Hamiltonian path at dimension `6915 + 16 + 138 = 7069`. Fixed screened error
  remained `+1.576416 mHa` versus matched PySCF.
- Eight source-backed sweeps reached `-2086.524053675786 Ha`, differing from
  the prior result by `-1.09e-11 Ha`; the maximum sweep-energy difference was
  `6.22e-9 Ha`. No collapse or spin-basin change occurred. Final PRF/external
  occupations were `2.377826e-6 / 1.841130e-3 e`. Strict convergence was not
  declared, so this closes a bounded replay rather than a production endpoint.

Validation / evidence:
- Reviewed both July 16 reports and independently verified restart SHA-256
  `5e85af3caffe129a00593221233d70cf414bce9a09e3ac931dff0625d1b989f3`.
  The largest raw external `H1` delta was `3.28e-7 Ha` on a roughly `57336 Ha`
  diagonal, `5.72e-12` relative, with fixed-state expectation below
  `2.4e-13 Ha`. Full trajectory parity accepts this replay only; it creates no
  generic or state-dependent tolerance policy.
- Authority checks, docs tests, Documenter, manager-log bound, and
  `git diff --check` form the docs-only gate.

Goal advancement / guardrail:
- MT2 source-migration validation is complete; old CR2-local mechanics are
  historical. Selection and interpretation remain consumer-owned. The
  approximately 17-minute moment-only localization rebuild friction is
  deferred because the accepted restart removes it from the immediate path.
  Pause for the next scientific choice; do not automatically open q6, helper,
  solver, artifact, exchange, or production-endpoint work.
