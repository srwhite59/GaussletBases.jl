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
- This live volume begins with Pass 407. Pass entries are preserved in accepted
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
- The source-backed Cr2 composition/replay migration is closed. The matched
  full-parent/PQS/White-Lindsey H2+ paper point at `R=2` is box-converged under
  the accepted padding-10/20 control; padding `10` remains the campaign value.
  The finer `tail_spacing=2.0` control is complete and favors retaining `2.8`;
  the fixed parent-state capture diagnostic is complete. The next bounded step
  is a neutral-metadata H2 one-body preflight with the frozen contracted
  H/cc-pVTZ `s,p` supplement and production residual cutoff.
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

**MT2 - Controlled Cr2 source migration (completed).** The source-backed
fixed-state and bounded replay reproduced the former consumer-local path. Any
new Cr2 endpoint, contraction, exchange, or solver interpretation is a
separate scientific choice.

**MT3 - Approved pending producer facilities (active).** The private matched
matrix-free parent/PQS/WL H2+ point and padding-10/20 control are complete;
padding `10` is retained. The approved `tail_spacing=2.8/2.0` control is also
complete and favors `2.8`; the parent-ground-state capture diagnostic at the
accepted point is complete. One frozen neutral-H2 supplemented one-body
preflight is approved next, using contracted H/cc-pVTZ `s,p`, no width filter,
and the production `1e-6` residual cutoff. H2 interaction/RHF, supplement and
cutoff scans, and fixed-parent ladders remain separate. Standard60 and
canonical-driver Coulomb exposure remain separate. The retained-GTO EGOI
helper remains pending and must not absorb the unrelated
`hamiltonian_corrections.jl` WIP without review.

**MT4 - Residual and protected-basis evidence (active).** Keep the residual
spectral audit measurement-only. Protected atoms, counterpoise, and any new
injection/localization policy remain separate future decisions. Parent residual
function mechanics and the onsite-calibrated Gaussian direct resource are
implemented and source-backed. Selection, transition-density exchange, and
PRF-to-GTO-residual interactions remain consumer or measurement questions.

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

## Cartesian Hamiltonian Producer Pass 423 - Authorize Private PQS/WL Paper Driver

Commit(s):
- this commit - approve one private matched H2+/H2 paper-validation command.

Summary:
- Added `HP-PQS-PAPER-H2-DRV-FN-01` for exactly
  `bin/pqs_paper_h2_driver.jl`. The script must use the current one-center
  staged working-basis route, genuine PQS and White-Lindsey terminal
  realizations, arbitrary-center nuclear attraction, and the frozen
  high-accuracy paper setup. It does not replace the canonical driver.
- The first gate is only the `R=2` H2+ lowest terminal `H1` state for both
  methods, with no `Vee` or RHF. H2 fixed-state/closed-shell work remains
  ordered after that acceptance.

Validation / evidence:
- Read the July 24 route-readiness audit and reconciled the live staged
  constructors, terminal realizers, arbitrary-center nuclear helper,
  high-accuracy IDA owner, matrix contractions, canonical-driver boundary, and
  machine authority. Authority checks, generated-view parity, focused conflict
  scans, docs links, manager-log bound, and `git diff --check` form the gate.

Goal advancement / guardrail:
- MT3 gains one bounded paper-validation driver. Preferred/hard size limits are
  `125/150` added `bin` lines, with zero `src` lines and no committed test.
  Scratch parent/density oracles, generalized-overlap machinery, AddNest,
  public APIs, schemas, broad driver infrastructure, and physics-default
  changes remain forbidden. The single scientific command is the acceptance
  test; failure within the one-file budget means no implementation commit.

Carrying-cost accounting:
- deleted: none; this pass grants authority only.
- simplified: none; the canonical driver remains unchanged.
- quarantined: scratch parent/density oracles and AddNest reporting remain
  ignored evidence, not implementation inputs.
- not deleted because: the new driver does not exist yet.
- exact remaining caller/blocker: no stable repo command joins the existing
  staged kernels for the paper endpoint.
- added/deleted `src` lines: `0/0`.
- new tests: none.
- new metadata/status fields: none; one machine authority record and exact
  planned driver path were added.
- validation: authority check/self-test, docs tests, Documenter, focused
  path/conflict scans, manager-log bound, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 424 - Correct Private PQS/WL Paper Route

Commit(s):
- this commit - replace the false one-center paper-driver premise with the
  existing two-center staged H2+ construction.

Summary:
- `HP-PQS-PAPER-H2-DRV-FN-01` remains approved for exactly
  `bin/pqs_paper_h2_driver.jl`, but the public neutral facade is no longer part
  of its contract. The driver must carry the actual `nup=1`, `ndn=0`,
  two-center system through `cartesian_system` to `cartesian_transforms`.
- Both methods use one `ns=5`-derived positive combined-inverse-sqrt parent
  with `core_spacing=0.30`, `s_factor=1`, `tail_spacing=2.8`, at least `10`
  bohr padding at `R=2`, and high135 Coulomb data. PQS/WL route-local orders
  are `q=5/3`; dimensions and topology are live outputs, not fixed assertions.
- The first gate remains H2+ `H1` only. Its readable report must carry the
  complete existing structured due-diligence object for both methods. Missing
  White-Lindsey shell/slab rows are a source-backed blocker, not permission
  for a driver-local reconstruction.

Validation / evidence:
- Reconciled the neutral facade rejection, multicenter combined-inverse-sqrt
  parent owner, staged PQS/WL realizers, arbitrary-center nuclear owner,
  terminal inventory, and due-diligence implementation against committed
  source. Authority check/self-test, generated-view parity, docs tests,
  Documenter, manager-log bound, scoped diff review, and
  `git diff --check` form the gate.

Goal advancement / guardrail:
- MT3 remains active with corrected implementation authority. Historical
  fixed dimensions, `core_spacing=0.15`, radius `6`, midpoint-parent language,
  old energies, and a required PQS/WL energy ordering are removed. The
  `125/150` line budget, zero-source/test rule, and no-framework/public/schema
  boundaries remain unchanged.

Carrying-cost accounting:
- deleted: false numerical assumptions and historical acceptance values from
  the live driver contract.
- simplified: one physical parent and one padding convention now govern both
  terminal methods.
- quarantined: the current untracked driver draft is neither read nor accepted
  as implementation evidence.
- not deleted because: the approved driver path remains planned until the
  corrected endpoint and full due diligence fit the hard budget.
- exact remaining caller/blocker: implementation must prove that the existing
  White-Lindsey route returns complete shell/slab due-diligence rows.
- added/deleted `src` lines: `0/0`.
- new tests: none.
- new metadata/status fields: none.
- validation: authority check/self-test, generated-view parity, docs tests,
  Documenter, manager-log bound, scoped docs-only review, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 425 - Complete White-Lindsey Due Diligence

Commit(s):
- this commit - repair the existing terminal-report join under
  `HP-DRV-SHELLDD-FN-01`.

Summary:
- White-Lindsey reporting now follows its native retained-unit order instead
  of the PQS-only support-rule records. The corrected two-center `R=2` route
  exposes `211` realized rows in `11` physical regions; each of its eight
  complete shells sums to `98` columns and all rows close the `1109`-column
  terminal dimension.
- PQS keeps its support-record join unchanged. No basis, coefficient, H1,
  interaction, artifact, or public-driver behavior changed.

Validation / evidence:
- Package load and the public Cartesian base regression passed. The corrected
  private H2+ command completed both high135 H1 solves with shared parent
  fingerprint, finite symmetric matrices, and eigen-residuals below `7e-14`.
  The source patch is `60/25` added/deleted lines and passes
  `git diff --check`.

Carrying-cost accounting:
- deleted: WL dependence on unavailable PQS support-rule records.
- simplified: native retained-unit order now owns WL report/block alignment.
- quarantined: key probe and paper-gate outputs remain under `/tmp`.
- not deleted because: the PQS support-record path remains its live owner.
- exact remaining caller/blocker: implement and lifecycle-close the separately
  authorized private paper driver.
- added/deleted `src` lines: `60/25`.
- new tests: none.
- new metadata/status fields: none.
- validation: package load, focused public-base regression, corrected H2+
  endpoint, scoped diff review, and `git diff --check`.

### Medium-Term Goal Checkpoint After Pass 425

- **MT1 conformance remediation - active.** The WL report gap is closed; other
  Pass 398 discrepancies remain separate.
- **MT2 controlled Cr2 source migration - completed.** Unchanged.
- **MT3 approved pending facilities - active.** The paper-driver source blocker
  is cleared; its own implementation/lifecycle commit is next.
- **MT4 residual/protected evidence - active.** Parent-backed mechanics are
  implemented; selection and transition/exchange questions remain scientific.
- **MT5 authority maintenance - maintenance.** No new authority was needed.
- **MT6 carrying-cost control - active.** This pass adds only the bounded
  report join and no test, payload, metadata, or compatibility layer.

## Cartesian Hamiltonian Producer Pass 426 - Implement Private H2+ Paper Driver

Commit(s):
- `3267d0b80` - complete native White-Lindsey due diligence.
- this commit - add the private H2+ command and move its path to existing.

Summary:
- Added the `150`-line `bin/pqs_paper_h2_driver.jl`. It constructs actual
  two-center H2+ independently through source-box-first PQS and current
  White-Lindsey stages, using one `ns=5`, `core_spacing=0.30`, high135 parent.
  It builds only terminal H1, solves the lowest state, and writes TSV/text
  evidence; H2, Vee, RHF, and artifacts remain disabled.
- At `R=2`, PQS/WL dimensions are `1285/1109` and electronic energies are
  `-1.1019722712680262/-1.0990733343694545 Ha`. Both use the same
  `21x21x29`, dimension-`12789` parent and fingerprint. Live dimensions are
  evidence, not an equality or ordering policy.

Validation / evidence:
- The sole acceptance command passed with high135, H1 symmetry errors
  `7.11e-15 Ha`, eigen-residuals below `7e-14`, and diagnostic combined-process
  peak RSS about `2.21 GB`.
  Full due diligence records bounds near `+/-11.10` transverse and
  `+/-11.96` longitudinal, more than `10.95` bohr beyond each nucleus,
  direct-core/eight-shell/two-slab topology, and all warning rows.
- No test edit was accepted: the existing public test still masks empty due rows
  through its inventory fallback. Durable WL-row assertions require a separate
  docs-only test-authority amendment.
- The anti-bloat scan flags only fixed two-nucleus and fixed three-axis
  `Tuple` conversions; neither is a basis-size or runtime route inventory.

Carrying-cost accounting:
- deleted: no numerical owner; the driver joins existing stages only.
- simplified: no scratch parent, density oracle, or reporting framework copied.
- quarantined: H2/Vee/RHF and optional artifact behavior remain unimplemented.
- not deleted because: the private command is the reproducible paper endpoint.
- exact remaining caller/blocker: paper-manager review before the separate H2
  endpoint pass.
- added/deleted `src` lines: `0/0`; added `bin` lines: `150`.
- new tests: none.
- new metadata/status fields: none.
- validation: package load, authority checks, public-base regression, sole
  H2+ command, TSV/report readback, scope/line review, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 427 - Authorize Full-Parent H2+ Row

Commit(s):
- this commit - extend the existing private-driver authority to the mandatory
  matrix-free full-parent reference row.

Summary:
- `HP-PQS-PAPER-H2-DRV-FN-01` remains the sole driver ID and
  `bin/pqs_paper_h2_driver.jl` remains its sole path. For `R=2`,
  `method=:both` must now report full parent, PQS `q=5`, and White-Lindsey
  `q=3` from the shared `21x21x29` PGDG parent.
- The parent calculation may use live axis overlap/kinetic data, analytic
  high135 factors, symmetric axis orthogonalization, small file-local
  mode-product applies, and the existing apply-based Lanczos routine. A dense
  `12789^2` matrix, source helper, generalized eigensolver, quadrature, and
  reusable oracle framework remain forbidden.
- The report adds independent-reference, parent-resolution, and contraction
  errors against `-0.6026342144949465 Ha`; the parent residual and
  `T/U_left/U_right` decomposition must be independently recomputed.

Validation / evidence:
- Reconciled the 150-line accepted driver, live parent dimensions, PGDG axis
  data, and existing apply-based solver. Authority check/self-test, generated
  views, docs tests, Documenter, manager-log bound, scoped diff review, and
  `git diff --check` form the docs gate.
- Rotated Passes 380-406 verbatim into a 1,052-line task-gated archive with
  SHA-256 `14ccf05eb960757cb3335351658cf54b707c6003e9bbb40916332b398e9a0767`;
  Passes 407-426 remain live.

Goal advancement / guardrail:
- MT3 advances to the three-row H2+ comparison. H2, Vee/IDA, RHF/SCF, scans,
  artifacts, public surfaces, dependencies, source/test/tool changes, PRFs,
  screening, and EGOI remain outside this pass. Preferred/hard final driver
  size is `225/250` lines.

Carrying-cost accounting:
- deleted: stale live-ledger startup pressure through verbatim rotation.
- simplified: one driver and one report own all three H2+ rows.
- quarantined: historical ledger entries remain in the linked archive.
- not deleted because: the accepted 150-line terminal driver is the extension
  base.
- exact remaining caller/blocker: implement the parent row readably within the
  250-line hard limit.
- added/deleted `src` lines: `0/0`; new tests: none.
- new metadata/status fields: none; only private TSV/report columns are
  authorized later.
- validation: authority/docs checks, manager-log bound, scoped review, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 428 - Close Full-Parent H2+ Gate

Commit(s):
- `171e2f368` - add the matrix-free full-parent row to the private paper
  driver.

Summary:
- The `248`-line driver now reports the shared `12789`-function parent, PQS
  `q=5` (`1285` functions), and White-Lindsey `q=3` (`1109` functions).
  Their total energies are `-0.602240572951`, `-0.601972271268`, and
  `-0.599073334369 Ha`.
- Against the frozen reference, parent resolution contributes
  `+0.393642 mHa`; PQS adds `+0.268302 mHa` contraction error and
  White-Lindsey adds `+3.167239 mHa`. The archived PQS/WL dimensions,
  decompositions, and energies reproduce to numerical precision.

Validation / evidence:
- A detached clean worktree at `171e2f368` emitted three `git_dirty=false`
  rows with one parent fingerprint. Parent residual, matrix-free symmetry, and
  overlap-identity errors were `9.57e-11`, `5.69e-16`, and `2.78e-14`;
  left/right nuclear expectations agreed to `4.4e-16 Ha`.
- Package load, authority check/self-test, manager-log bound, an independent
  rectangular tensor-index oracle, archived-row comparison, due-diligence and
  TSV/report readback, anti-bloat scan, staged scope review, and
  `git diff --check` passed.

Goal advancement / guardrail:
- MT3 closes the first parent-versus-contraction `R=2` H2+ gate. Padding,
  tail-spacing, geometry, H2, Vee/IDA, RHF/SCF, artifacts, and publication
  timing remain separate later decisions.

Carrying-cost accounting:
- deleted: the review rejected a runtime-sized operator-validation tuple.
- simplified: overlap whitening is shared around the summed parent action.
- quarantined: the parent oracle remains private to the removable paper
  driver.
- not deleted because: the terminal rows and complete due diligence remain
  the active paper comparison.
- exact remaining caller/blocker: paper-manager review before separately
  authorized convergence controls or H2 work.
- added/deleted `src` lines: `0/0`; added/deleted `bin` lines: `110/12`.
- new tests: none.
- new metadata/status fields: none; only approved private report columns.
- validation: the clean scientific gate and mechanical checks listed above.

## Cartesian Hamiltonian Producer Pass 429 - Authorize H2+ Padding Control

Commit(s):
- this commit - amend the private paper-driver authority for one padding-only
  convergence control.

Summary:
- The accepted padding-10 parent/PQS/White-Lindsey outputs were preserved with
  matching hashes. A clean padding-20 attempt stopped before solving because
  the driver still required the padding-10 `21x21x29` axes.
- `HP-PQS-PAPER-H2-DRV-FN-01` now authorizes only the `R=2`,
  `padding=10.0/20.0` pair. The follow-on replaces that one fixed-shape
  assertion with live positive odd axis counts, dimension-product consistency,
  finite positive weights, and PQS/WL parent-object and fingerprint parity.

Goal advancement / guardrail:
- MT3 advances to a bounded box-convergence check without changing the parent
  mapping, high135 interaction, route-local `q=5/q=3`, matrix-free solve,
  decomposition, report columns, or endpoint physics. Each of the three energy
  shifts has a provisional `0.01 mHa` tolerance. This is not a general padding
  or tail-spacing scan and does not authorize H2, Vee/IDA, RHF/SCF, artifacts,
  tests, helpers, or source outside the existing driver.

Carrying-cost accounting:
- deleted: none in this docs-only authority pass.
- simplified: fixed-axis policy becomes one live parent-validity contract.
- quarantined: wider padding, tail-spacing, and geometry scans remain outside
  authority.
- not deleted because: the accepted 248-line private driver owns the paper
  comparison.
- exact remaining caller/blocker: replace one driver assertion and run the
  clean padding-20 command.
- added/deleted `src` lines: `0/0`; new tests: none.
- new metadata/status fields: none; output columns remain unchanged.
- validation: authority/docs checks, scoped review, manager-log bound, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 430 - Close H2+ Padding Control

Commit(s):
- `67351aca6` - replace the padding-10 axis equality with live parent-dimension
  validation.

Summary:
- A tracked-clean worktree at the implementation commit completed the exact
  `R=2`, padding-20 parent/PQS/White-Lindsey command. Padding-20 axes are
  `31x31x39`; parent/PQS/WL dimensions are `37479/1969/1599`, compared with
  padding-10 dimensions `12789/1285/1109`.
- Padding-20 total-energy shifts are `-1.55e-9`, `-5.25e-7`, and
  `-7.50e-6 mHa` for parent, PQS, and White-Lindsey. All are far below the
  provisional `0.01 mHa` threshold, so padding `10` remains the campaign
  value. The `0.393642 mHa` parent-resolution conclusion is not a finite-box
  effect at this scale.

Validation / evidence:
- Independent TSV/full-report readback verified three clean rows, one
  per-padding parent fingerprint, exact contraction-error identities,
  contiguous final-column accounting, finite symmetric operators, and
  residuals `8.91e-11`, `9.77e-14`, and `1.14e-13`.
- Both paddings retain direct-core/complete-shell/two-z-slab topology. The
  larger box adds five outer shells (`8 -> 13`) with no new warning class and
  no weight, padding, or large-identity warning. Diagnostic peak RSS rises
  from `2.32` to `3.59 GB`; row times rise to `54.7/116.1/108.4 s`.

Goal advancement / guardrail:
- LT3 and MT3 gain a bounded factorized parent convergence result without a
  dense parent matrix. One finer tail-spacing control is next, but requires
  separate authority; H2, Vee/IDA, RHF/SCF, geometry curves, and broad scans
  remain closed.

Carrying-cost accounting:
- deleted: the fixed `21x21x29` padding-10 assumption.
- simplified: one live odd-axis/product check now covers both authorized boxes.
- quarantined: comparison code and padding-20 evidence remain under `/tmp`.
- not deleted because: the private driver remains the reproducible paper route.
- exact remaining caller/blocker: docs-only authority for one finer
  tail-spacing control.
- added/deleted `src` lines: `0/0`; added/deleted `bin` lines: `1/1`.
- new tests: none.
- new metadata/status fields: none.
- validation: package load/parse, authority check/self-test, clean endpoint,
  due-diligence comparison, anti-bloat/scope review, and `git diff --check`.

### Medium-Term Goal Checkpoint After Pass 430

- **MT1 conformance remediation - active.** Unchanged; this scientific control
  does not consume another Pass 398 discrepancy.
- **MT2 controlled Cr2 source migration - completed.** Unchanged.
- **MT3 approved pending facilities - active.** The H2+ box control is complete
  and retains padding `10`; authorize only one finer tail-spacing control next.
- **MT4 residual/protected evidence - active.** Unchanged; no residual,
  protected, injection, or interaction policy entered this pass.
- **MT5 authority maintenance - maintenance.** Pass 429 kept the control
  record-local; lifecycle closure should remain compact.
- **MT6 carrying-cost control - active.** The implementation is one replacement
  line, adds no test or schema, and removes one stale fixed-size assumption.

## Cartesian Hamiltonian Producer Pass 431 - Authorize H2+ Tail Control

Commit(s):
- this commit - amend the existing private-driver authority for one finer
  tail-spacing control.

Summary:
- Padding-20 raw evidence was preserved with hashes matching the accepted
  temporary outputs; repository files were unchanged. The box control remains
  evidence at `tail_spacing=2.8`, and padding `10` remains the paper setting.
- `HP-PQS-PAPER-H2-DRV-FN-01` now permits exactly three `R=2` combinations:
  `(10.0,2.8)`, `(20.0,2.8)`, and `(10.0,2.0)` for
  `(padding,tail_spacing)`. The follow-on adds one private input, defaults it to
  `2.8`, and rejects `padding=20.0`, `tail_spacing=2.0`.

Goal advancement / guardrail:
- MT3 advances to one finer-tail parent/PQS/White-Lindsey comparison while
  preserving `ns=5`, `core_spacing=0.30`, `s_factor=1`, high135, ordinary
  source span, route-local `q=5/q=3`, matrix-free parent application, complete
  due diligence, and existing report columns. The three energy shifts are
  reported without an ordering or convergence threshold.

Carrying-cost accounting:
- deleted: none in this docs-only authority pass.
- simplified: the existing recorded tail-spacing value becomes one bounded
  resolved input rather than a new result surface.
- quarantined: other padding/tail combinations, capture, `q` ladders, geometry
  curves, H2/Vee/IDA/RHF, and enrichment remain outside authority.
- not deleted because: the accepted 248-line private driver owns the paper
  comparison.
- exact remaining caller/blocker: implement the bounded input/combination gate
  within the 250-line limit and run the clean `tail_spacing=2.0` command.
- added/deleted `src` lines: `0/0`; new tests: none.
- new metadata/status fields: none; the existing TSV/report column is reused.
- validation: authority/docs checks, scoped review, manager-log bound, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 432 - Close H2+ Tail Control

Commit(s):
- `33357da4f` - add the exact private tail-spacing input and combination gate.

Summary:
- A tracked-clean worktree completed the `R=2`, padding-10,
  `tail_spacing=2.0` parent/PQS/White-Lindsey comparison. Axes changed from
  `21x21x29` to `23x23x31`; parent/PQS/WL dimensions changed from
  `12789/1285/1109` to `16399/1431/1207`.
- Relative to `2.8`, the full-parent total improves by `0.007154 mHa`, but PQS
  and White-Lindsey totals worsen by `0.044124` and `0.341989 mHa`.
  Contraction errors increase by `0.051278` and `0.349143 mHa`. The evidence
  therefore favors retaining `tail_spacing=2.8`; paper-manager acceptance is
  still required before calling it frozen.

Validation / evidence:
- Independent TSV/full-report readback verified three clean rows, one live
  parent fingerprint, exact contraction identities, contiguous terminal
  columns, high135, and residuals `9.81e-11`, `8.37e-14`, and `8.08e-14`.
- Both values retain direct-core/complete-shell/two-z-slab topology. The finer
  tail adds one shell (`8 -> 9`), with no weight, padding, or large-identity
  warning. Diagnostic row times are `15.3/40.0/19.1 s`, and peak RSS is
  `2.45 GB`.

Goal advancement / guardrail:
- MT3 closes the bounded tail control without H2, Vee, SCF, capture, or a
  general scan. The nonmonotonic terminal result is retained as scientific
  evidence, not converted into a new default or source policy.

Carrying-cost accounting:
- deleted: the hard-coded `tail_spacing=2.8` route value and broad positive-R
  padding gate.
- simplified: one exact three-combination input gate owns the accepted controls.
- quarantined: comparison code and tail-2.0 outputs remain under `/tmp`.
- not deleted because: the private driver remains the reproducible paper route.
- exact remaining caller/blocker: paper-manager tail-spacing acceptance, then
  separate authority for physical-target capture.
- added/deleted `src` lines: `0/0`; added/deleted `bin` lines: `6/6`.
- new tests: none.
- new metadata/status fields: none; the existing tail-spacing column is reused.
- validation: package load/parse, forbidden-combination rejection, authority
  check/self-test, clean endpoint, due-diligence comparison, staged scope, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 433 - Authorize Parent-State Capture

Commit(s):
- this commit - amend the existing private H2+ driver authority for one
  full-parent-state capture diagnostic.

Summary:
- The accepted tail-2.0 outputs were preserved under
  `/Users/srw/Dropbox/Papers/PQS/provenance/h2plus_R2_tail2_clean_2026-07-26/`
  with matching hashes and no repository edit. The bounded tail decision
  retains `tail_spacing=2.8`.
- `HP-PQS-PAPER-H2-DRV-FN-01` now permits capture only at `R=2`,
  `padding=10.0`, `tail_spacing=2.8`, and `method=:both`. The driver projects
  the one normalized matrix-free parent ground state into the frozen PQS
  `q=5` and White-Lindsey `q=3` terminal spans using existing axis overlap
  factors and support-local terminal blocks.

Goal advancement / guardrail:
- MT3 advances from box/tail controls to physical-target retention. Exact
  support partition, native due-row/block order, terminal orthogonality,
  global norm closure, and physical shell/slab regional closure are mandatory.
  No capture ordering, terminal reconstruction change, enrichment, H2, Vee,
  SCF, or broader scan is authorized.

Carrying-cost accounting:
- deleted: none in this docs-only authority pass.
- simplified: one blockwise factorized projection reuses the existing parent
  state and terminal realization; no dense parent-by-terminal object is added.
- quarantined: variable regional details stay in the readable report; capture
  at every other campaign point remains outside authority.
- not deleted because: the private driver remains the reproducible paper route.
- exact remaining caller/blocker: implement the diagnostic within the
  preferred `300`/hard `315` line limit and run the clean baseline command.
- added/deleted `src` lines: `0/0`; new tests: none.
- new metadata/status fields: none; exactly four private TSV/report columns are
  approved.
- validation: authority check/self-test, generated-view parity, docs tests,
  Documenter, manager-log bound, staged scope, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 434 - Implement H2+ Parent-State Capture

Commit(s):
- `1e373d0c1` - add the bounded blockwise parent-state capture diagnostic.

Summary:
- A clean detached replay at the frozen `R=2`, padding-10, tail-2.8 point
  retained `0.9999709356338488` of the normalized parent ground state in PQS
  and `0.9991007020512097` in White-Lindsey. Lost norms are
  `2.9064366e-5` and `8.9929795e-4`; WL therefore loses `30.94` times more.
- PQS loss is dominated by shell 1, while WL loss is dominated by shells 2
  and 3. Direct-core loss is below `5e-28` for both. Regional support and
  column totals close exactly at `12789` and `1285/1109`.
- Review caught and removed an x-fast local index reconstruction. The
  matrix-free parent apply now uses canonical z-fast ordering and consumes
  each terminal block's actual support indices. Energies changed by at most
  `5.1e-14 Ha`; accepted parent/PQS/WL totals remain
  `-0.6022405729510050/-0.6019722712680262/-0.5990733343694545 Ha`.

Goal advancement / guardrail:
- MT3 now has physical-target retention evidence after box and tail controls.
  Capture is diagnostic, with no ordering gate, enrichment decision, H2,
  interaction, or SCF authority.

Carrying-cost accounting:
- deleted: the redundant and incorrectly ordered local flat-index helper.
- simplified: one canonical parent ordering lets projection use live support
  indices directly.
- quarantined: variable regional rows remain readable-report-only.
- not deleted because: the 299-line private driver remains the paper route.
- exact remaining caller/blocker: paper-manager interpretation before any
  fixed-parent dimension ladder; H2 remains deferred.
- added/deleted `src` lines: `0/0`; added/deleted `bin` lines: `56/5`.
- new tests: none.
- new metadata/status fields: none; four approved private report columns.
- validation: package load, authority check/self-test, endpoint gates,
  independent TSV/regional closure, clean-worktree replay, scoped diff, and
  `git diff --check`.

## Cartesian Hamiltonian Producer Pass 435 - Authorize Supplemented H2 One-Body Preflight

Commit(s):
- this commit - amend the existing private paper-driver authority for one
  frozen neutral-H2 supplemented one-body endpoint.

Summary:
- Pass 434 closed the parent-state diagnostic at the accepted H2+ point:
  PQS/WL capture fractions are `0.9999709356338488` and
  `0.9991007020512097`, with exact regional/support closure and no enrichment
  decision. The next gate therefore adds the established Gaussian supplement
  while preserving one-body isolation.
- `HP-PQS-PAPER-H2-DRV-FN-01` now permits exactly neutral H2 metadata at
  `R=2`, padding `10`, tail spacing `2.8`, matched PQS `q=5`/WL `q=3`, and the
  bundled contracted H/cc-pVTZ `s,p` supplement with no width filter.
  Residual selection remains owner-local at the production `1e-6` cutoff.

Goal advancement / guardrail:
- LT3 and MT3 advance from terminal capture to exact augmented one-body
  behavior. Parent/terminal supplement capture, owner residual spectra, exact
  augmented kinetic/by-center nuclear matrices, overlap identity, and lowest
  one-body states are required. The full parent remains capture-only; no
  generalized parent augmentation, `Vee`, IDA/MWG interpretation, RHF/SCF,
  artifact, public input, or supplement/cutoff scan is authorized.

Carrying-cost accounting:
- deleted: none in this docs-only authority pass.
- simplified: the private driver reuses the existing supplemented facade,
  residual builder, and exact augmented operators rather than adding a source
  helper or second composition path.
- quarantined: variable supplement/capture/rank detail remains readable-report
  only.
- not deleted because: the 299-line driver remains the bounded reproducible
  paper route.
- exact remaining caller/blocker: implement the five-row preflight within the
  preferred `400`/hard `425` line limit and run its clean `R=2` gate.
- added/deleted `src` lines: `0/0`; new tests: none.
- new metadata/status fields: none; ten private TSV/report fields are approved.
- validation: authority check/self-test, generated-view parity, docs tests,
  Documenter, manager-log bound, staged scope, and `git diff --check`.

### Medium-Term Goal Checkpoint After Pass 435

- **MT1 conformance remediation - active.** Unchanged; this paper endpoint does
  not close a Pass 398 discrepancy.
- **MT2 controlled Cr2 source migration - completed.** Unchanged.
- **MT3 approved pending facilities - active.** Box, tail, and capture controls
  are complete. The next gate is only the frozen supplemented H2 one-body
  preflight; interaction and orbital relaxation remain closed.
- **MT4 residual/protected evidence - active.** The endpoint consumes the
  production residual policy without changing cutoff, injection, protection,
  or interpretation.
- **MT5 authority maintenance - maintenance.** The new scope stays on one
  existing ID, canonical page, and driver path.
- **MT6 carrying-cost control - active.** Zero source/test/helper files are
  approved; implementation must reuse current numerical owners and remain
  within the private-driver line ceiling.

## Cartesian Hamiltonian Producer Pass 436 - Implement Supplemented H2 One-Body Preflight

Commit(s):
- `119bc17c7` - implement the bounded five-row neutral-metadata supplemented
  one-body endpoint.
- this commit - reconcile the existing ID to implemented maintenance.

Summary:
- A detached clean replay at `R=2`, padding `10`, and tail spacing `2.8`
  reproduced the accepted parent/PQS/WL dimensions and energies exactly, then
  appended PQS/WL supplemented dimensions `1303/1127`. Their total one-body
  energies are `-0.6024777336542515` and `-0.6024638359415411 Ha`,
  improvements of `0.505462` and `3.390502 mHa` over the matching bare rows.
- Both routes used the same 18-function contracted H/cc-pVTZ supplement and
  parent fingerprint. Parent capture has minimum singular value
  `0.9991415015`; terminal minima are `0.9781753530/0.7793743583`. Each owner
  retained all nine residual directions at the production `1e-6` cutoff.
- Full augmented-overlap errors are `3.64e-11/2.05e-12`; H1 symmetry errors
  are zero and lowest-state residuals remain below `9.4e-14`. Complete PQS/WL
  due diligence retained the `21x21x29` parent, `11.10/10.96` bohr transverse/
  longitudinal padding, exact column accounting, and only existing advisory
  shell/source-shape warnings.
- The final driver is `416` lines, nine below the hard cap. Its `+128/-11`
  delta exceeds the initial expected addition band by 18 lines, but adds no
  source, test, helper, payload, artifact, or persistent metadata surface.

Goal advancement / guardrail:
- LT3/MT3 now have a source-backed supplemented one-body preflight. This is
  repository acceptance, not independent paper acceptance. Vee, IDA/MWG
  interpretation, same-state interaction, RHF/SCF, artifacts, and scans remain
  outside authority.

Carrying-cost accounting:
- deleted: the old H2 rejection and H2+-only output literals; total driver
  deletion was 11 lines.
- simplified: one supplement load and identity record serve both routes, while
  existing residual and exact augmented-operator owners remain authoritative.
- quarantined: variable spectra, primitive identity, and operator fingerprints
  remain readable-report-only.
- not deleted because: the accepted bare and parent-state controls remain the
  comparison baseline in the same private driver.
- exact remaining caller/blocker: independent clean scientific replay; any H2
  interaction endpoint requires a separate authority decision.
- added/deleted `src` lines: `0/0`; added/deleted `bin` lines: `128/11`.
- new tests: none.
- new metadata/status fields: none; exactly ten approved private report fields.
- validation: package load, authority check/self-test, exact clean-worktree
  endpoint, independent TSV/report and archived-bare-row comparison, staged
  scope/anti-bloat review, docs tests `56/56` and `10/10`, manager-log bound,
  local Documenter, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 437 - Authorize Fixed-State H2 Interaction Gate

Commit(s):
- this commit - amend the existing private paper-driver authority for one
  five-row fixed-state density-density interaction comparison.

Summary:
- The independent supplemented one-body replay accepted dimensions
  `12789/1285/1109/1303/1127` and PQS/WL augmented one-body errors
  `0.156/0.170 mHa`. The next gate now holds the accepted full-parent H2+
  orbital fixed after separate projection into the parent, bare, and
  supplemented spaces.
- Existing numerical owners suffice: the parent uses a matrix-free high135
  ordinary-IDA contraction, bare rows use terminal/site IDA, and supplemented
  rows retain terminal IDA in `G-G` plus existing MWG in `G-R/R-R`. No
  interaction rotation or new source seam is authorized.
- Repo output owns projection, category, symmetry, and energy-accounting
  checks. An external paper-validation script must reconstruct the exact same
  densities and converge the direct Coulomb oracle before any production
  interaction error is interpreted.

Goal advancement / guardrail:
- LT3/MT3 advance from one-body readiness to the bounded same-state
  interaction question without allowing relaxation. RHF/SCF, continuum-exact
  exchange, scans, artifacts, public controls, and an in-repo numerical oracle
  remain closed.

Carrying-cost accounting:
- deleted: none in this docs-only authority pass.
- simplified: the existing live five-row construction is reused rather than
  adding a companion driver or script-to-script state carrier.
- quarantined: variable block, factor-resource, and matrix diagnostics remain
  readable-report-only; the independent oracle remains external scratch work.
- not deleted because: the 416-line driver is the accepted reproducible state
  owner and already retains every object required by this endpoint.
- exact remaining caller/blocker: implement within preferred `500`/hard `525`
  lines, preserve all existing one-body fields, pass the clean repo gate, then
  complete the external same-density oracle replay.
- added/deleted `src` lines: `0/0`; new tests: none.
- new metadata/status fields: none; thirteen private TSV/report fields are
  approved.
- validation: source/authority seam audit, external contract reconciliation,
  authority check/self-test, generated-view parity, docs tests, Documenter,
  manager-log bound, docs-only staged scope, and `git diff --check`.

## Cartesian Hamiltonian Producer Pass 438 - Implement Fixed-State H2 Interaction Gate

Commit(s):
- `7d2b6dc61` - implement the bounded five-row fixed-state interaction
  measurement.
- this commit - reconcile the existing ID to implemented maintenance.

Summary:
- The clean `R=2` replay retained dimensions
  `12789/1285/1109/1303/1127` and projected the same full-parent H2+ orbital
  into the parent, bare PQS/WL, and supplemented PQS/WL spaces. Captures are
  `1.0/0.9999709356/0.9991007021/0.9999972239/0.9999933045`.
- Matrix-free parent high135 IDA, fresh terminal IDA, and augmented
  terminal-IDA/residual-MWG produced `J =
  0.6613194326/0.6614182833/0.6632012296/0.6614042003/0.6625621285 Ha`.
  Native coordinate identity and augmented G-G parity are exact; maximum
  state-norm and interaction-symmetry errors are `1.04e-14` and `4.89e-15`.
- Independent review caught a missing identity gate between the separately
  staged bare and supplemented working terminal coordinates. The accepted
  patch now compares block order, support indices, column ranges, and
  coefficient matrices before combining projected states, H1, and Vee.

Goal advancement / guardrail:
- LT3/MT3 now have the repo-owned same-state density-density measurement.
  These values establish production construction and accounting only. Method
  accuracy remains uninterpreted until the external replay matches captures,
  state fingerprints, H1, and an independently converged direct-Coulomb
  oracle; RHF/SCF remains closed.

Carrying-cost accounting:
- deleted: seven superseded driver lines while extending the existing route.
- simplified: one common target projection and the existing category owners
  serve all five rows; no companion path or interaction adapter was added.
- quarantined: variable matrix/resource diagnostics remain report-only and
  the direct-Coulomb oracle remains external.
- not deleted because: the private driver remains the reproducible owner of
  the accepted H2+ controls and five-row H2 construction.
- exact remaining caller/blocker: external clean same-density oracle replay
  before any PQS/White-Lindsey interaction-accuracy interpretation.
- added/deleted `src` lines: `0/0`; added/deleted `bin` lines: `109/7`.
- new tests: none.
- new metadata/status fields: none; exactly thirteen approved private report
  columns.
- validation: package load, doer endpoint replay, independent code review,
  exact clean-worktree endpoint at `git_dirty=false`, deterministic
  pre-existing/new-field parity, due-diligence inspection, staged scope,
  authority check/self-test, docs tests, Documenter, manager-log bound, and
  `git diff --check`.
