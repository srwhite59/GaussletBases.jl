# Residual Gaussian Injection Hybrid Memo

Status: design memo plus approved measurement-only audit authority under
`HP-RG-INJECT-AUDIT-01`, historical default-off `G`-injection authority under
`HP-RG-INJECT-FN-01`, and the current protected-original compact-main design
authority under `HP-RG-PROTECT-INJECT-DESIGN-01`. The protected-original design
is the current direction for compact-first RG/injection work. This document
does not approve source edits for that new design, a production default,
artifact schema changes, driver inputs, public API, full HF, or Cr2 workflow.

## Motivation

The current Residual Gaussian path classifies supplement directions by the
owner-local residual metric

```text
M = S_AA - X'X
```

after projection against the orthonormal terminal gausslet space. Retaining a
direction with small residual occupation `lambda` then forms a normalized
residual roughly like

```text
r = (I - P_G) a / sqrt(lambda)
```

This is exactly the operation that becomes fragile when a GTO direction is
already almost represented by the gausslets. Tightening
`residual_occupation_cutoff` can discard marginal directions, but evidence
after `HP-RG-CUTOFF-FN-02` still showed a low two-owner residual-sector mode
in the Cr2 residual-only audit. The likely failure mode is not merely "the
cutoff is too loose"; it is that near-gausslet GTO directions have only two
current fates:

- discard them, losing exact GTO-span information;
- keep and normalize their tiny residual complements, risking ghost residual
  functions.

The proposed third fate is injection: represent near-gausslet directions
exactly in the one-body basis by replacing the corresponding gausslet subspace,
without creating residual-Gaussian/MWG channels for them.

This is conceptually analogous to the injection construction in the Angular
Gausslet manuscript: an approximate subspace of a localized Gaussian span is
replaced by an exact target subspace, and the final localized basis is obtained
by orthonormalizing inside the injected span. Here the injected target is not
low-`l` spherical harmonics; it is the set of supplement modes whose residual
norm against the terminal gausslet span is small.

## Optional Switch

The proposed user/internal policy knob is:

```text
residual_injection_cutoff = lambda_inj
```

with the simple off rule:

```text
lambda_inj <= 0    injection disabled; current RG behavior
lambda_inj > 0     near-gausslet modes with lambda <= lambda_inj are injected
```

For the first audit policy, `lambda_inj` should be at least the active
`residual_occupation_cutoff`. Otherwise there is an ambiguous band of modes
that are neither injected nor retained as true residuals. This is not a
permanent mathematical requirement: a later design could intentionally approve
a discard buffer between injection and true-RG retention, but that policy must
be explicit rather than accidental.

The first practical sweep values should be treated as audit choices, not
defaults. A plausible starting range is:

```text
lambda_inj = 0          off/current behavior
lambda_inj = 1.0e-6    near current residual cutoff
lambda_inj = 1.0e-5
lambda_inj = 1.0e-4
```

## Candidate GTO Orthonormalization

The first step must stabilize the raw owner-local supplement span before
classifying injection or residual content. For each physical owner atom:

```text
S_AA = owner-local candidate overlap
S_AA = U diag(s) U'
keep candidate metric modes with s above the candidate-overlap threshold
A_tilde = A U_keep diag(s_keep)^(-1/2)
```

The candidate-overlap threshold removes linearly dependent GTO candidate
combinations. It is not a residual occupation cutoff and has no direct
physical meaning as residual content. The threshold should be recorded
separately, for example as an absolute/relative rule:

```text
keep s_i > max(candidate_overlap_atol,
               candidate_overlap_rtol * maximum(s))
```

The simple audit policy can begin with a relative scale near `1.0e-8` when
the candidates are individually normalized, but a future source lane must name
the actual threshold and failure rule explicitly.

## Local Classification, Global Injection

The tricky part is that injected functions must be orthonormal as one global
set, while true residual Gaussians should remain as owner-local as possible
for MWG descriptors.

The proposed split is:

```text
For each owner:
  1. build stable orthonormal candidate modes A_tilde;
  2. compute C = G' S A_tilde;
  3. diagonalize M = I - C'C;
  4. mark modes with lambda <= lambda_inj as provisional injected modes;
  5. keep modes with lambda > lambda_inj as provisional residual candidates.

Across all owners:
  6. concatenate provisional injected modes;
  7. globally orthonormalize/merge them in the S metric;
  8. drop duplicate injected directions by a separate injected-overlap rank
     threshold;
  9. obtain one global orthonormal injected subspace Y_inj.

Then:
 10. build the injected gausslet sector F = Y_inj + (G cap Y_inj^perp);
 11. for each owner separately, residualize its remaining candidate modes
     against F;
 12. apply owner-local residual occupation selection;
 13. perform the final inter-owner residual merge.
```

Classification is applied to the owner-local orthonormal candidate principal
modes, not to raw contracted GTO columns. If

```text
M v_i = lambda_i v_i
```

then the classified supplement-space mode is

```text
y_i = A_tilde v_i
```

These `y_i` modes are the provisional injected modes or provisional residual
candidates.

This deliberately does not diagonalize one global residual metric over all
atoms to decide residual Gaussians. Global injection is acceptable because
injected functions do not become MWG residual channels. Global residual
selection would reintroduce nonlocal residual rotations and would be the wrong
owner model for MWG.

## Injected Gausslet Sector

If `Y_inj` is the global orthonormal injected subspace, the injected gausslet
sector is:

```text
F = Y_inj op (G cap Y_inj^perp)
```

It has dimension `nG`, not `nG + dim(Y_inj)`. The injected functions replace
the corresponding approximate directions in the gausslet sector. They are not
added on top of the original gausslets.

A concrete construction should use the projection of the global injected
subspace into the original gausslet coefficient space. Let `Y` be the global
orthonormal injected modes and

```text
B = G' S Y
```

Then build an orthonormal complement `Q_perp` inside the original gausslet
coefficient space:

```text
Q_perp' Q_perp = I
B' Q_perp = 0
```

The injected gausslet sector can then be represented as

```text
F = [Y, G Q_perp]
```

This makes the `nG`-dimensional replacement explicit and gives a direct rank
and conditioning diagnostic for whether the injected modes lie stably in the
gausslet span.

Required guards for any future source lane:

- `dim(Y_inj) < nG`;
- the projection `B = G' S Y_inj` has full rank and acceptable condition;
- the final injected gausslet sector is orthonormal;
- remaining true residuals are orthogonal to the injected sector `F`, not only
  to the original `G`;
- duplicate injected directions across owners are merged or dropped by an
  explicit injected-subspace rank rule.

If the injected subspace is too large, nearly singular, or cannot be merged
without unstable rotations, the construction should stop and report the
blocker rather than falling back silently.

## One-Body And Interaction Convention

The proposed convention is:

- exact one-body operators use the true injected/raw representation;
- true residual Gaussians use the existing exact augmented one-body
  transformation;
- injected-sector two-body IDA inherits the original gausslet IDA semantics;
- only true residual directions get residual-GTO/MWG interaction channels.

The inherited IDA treatment for injected directions is an approximation. The
one-body basis changes exactly and one-body operators use the injected/raw
representation exactly; the two-body IDA keeps the original gausslet-sector
IDA treatment for the replaced subspace. This is stable in the limit
`lambda -> 0`: a direction already represented by gausslets should not create
a normalized residual function or a residual MWG density.

Do not give injected functions their own MWG descriptors in the first design.
That would reintroduce the near-zero residual-density problem through another
path.

## Distinct Thresholds

Future work should keep at least these policies separate:

```text
candidate_overlap_cutoff       raw GTO candidate linear dependence
residual_injection_cutoff      optional near-gausslet injection threshold
injected_overlap_cutoff        global duplicate injected-mode merge threshold
residual_occupation_cutoff     true RG retention after injection
identity_atol                  final residual identity validation tolerance
```

These thresholds answer different questions. They should not be collapsed into
one "stability" knob.

## HP-RG-INJECT-AUDIT-01 - Measurement-Only First Step

Status: approved by user direction, measurement-only. This is not production
source authority.

It should use ignored probes only and report, for Cr atom, Cr2 monomer
counterpoise if available, and Cr2:

- raw candidate counts;
- stable owner-local candidate counts after `S_AA` rank cleanup;
- provisional injected count by owner for trial `lambda_inj`;
- globally retained injected count after duplicate merge;
- rank and condition of `B = G' S Y_inj`;
- true RG count by owner;
- residual occupation spectra before and after injection;
- `K_RR` and `H1_RR` low eigenvalues for true RGs;
- projected one-body errors of the injected sector versus the original GTO
  span, separated into `K`, each unit `U_A`, and `H1`;
- owner weights and residual-occupation composition of low modes;
- GTO-span one-body accuracy or mismatch against the non-injected path where
  cheap;
- whether the low two-owner residual ghost sector disappears, shrinks, or
  persists.

Approved surfaces:

- ignored `tmp/work/*.jl` probes only;
- durable text/TSV output under `/Users/srw/dmrgtmp/...` or CR2 run
  directories.

Forbidden:

- production source changes;
- committed tests or fixtures;
- artifact schema/provenance/reader/manifest changes;
- driver changes;
- public API/export changes;
- RG default changes;
- automatic pruning or residual-selection implementation;
- MWG/IDA convention changes;
- dense Vee, full HF, or solver workflow;
- Cr2 full Hamiltonian, Cr2 artifact, or Cr2-specific workflow.

Validation for the audit:

- `git diff --check`;
- package load;
- ignored injection audit probe for the target Cr/Cr2 fixture;
- no full HF and no new Hamiltonian artifact.

Failure rule: if the audit cannot reconstruct the needed owner-local
candidate spans, injected-sector projection `B`, or residual-sector one-body
blocks cheaply from existing construction seams, stop and report the exact
missing reusable seam. Do not add production source instrumentation as part of
this lane.

## Audit Result And Implementation Decision

The first `HP-RG-INJECT-AUDIT-01` probe did not remove the current Cr2 low
two-center residual sector under the tested reconstruction. The best tested
value, `lambda_inj = 1.0e-4`, reduced the severity only modestly:

```text
lambda_inj        injected   true RG count   min K_RR   min H1_RR
0                 0          132             0.428594   -7.349209
1.0e-4            38         100             0.445326   -7.061948
```

The audit also showed that the trial injected sector was numerically healthy:
for the dimer at `lambda_inj = 1.0e-4`, `B` condition was about `1.002` and
the implicit `F' S F` error was about `3.9e-12`. The result should be
interpreted as follows:

- injection did not by itself fix the current Cr2 low-H1 residual sector;
- the first audit `lambda = 0` baseline did not exactly match production RG;
- the surviving low mode after `lambda_inj = 1.0e-4` had no low-occupation
  weight below that threshold, so the remaining issue is not only the
  just-above-cutoff tail;
- nevertheless, RG alone still has the bad singular-complement limit, while
  injection gives near-gausslet GTO directions the correct third fate.

Implementation is therefore approved as a principled construction improvement,
not as a claim that injection alone solves the Cr2 residual-sector safety
problem. A later spectral stop-and-report gate may still be needed.

## HP-RG-INJECT-FN-01 - Default-Off In-Memory Injection Implementation

Status: historical default-off `G`-injection source authority. It is not the
current compact-first implementation target and must not be used to turn on
the existing injection path as-is for the protected-original design below.
Any implementation blurb for the current Cr2 compact-first direction must use
`HP-RG-PROTECT-INJECT-DESIGN-01` as governing design and must name a fresh
source surface or source-amendment authority.

Approved source surface:

- `src/cartesian_residual_gaussians/residual_basis.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`;
- `src/cartesian_residual_gaussians/mwg_interaction.jl`;
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` only
  for narrow internal keyword plumbing, same-construction validation, and
  compatibility wiring.

Approved behavior:

- add an internal `residual_injection_cutoff` option whose default preserves
  current behavior, with `residual_injection_cutoff <= 0` meaning injection
  disabled;
- when injection is disabled, preserve current production residual selection,
  transforms, exact augmented one-body operators, MWG/IDA interaction, H2
  endpoint values, and artifact behavior within roundoff;
- when injection is enabled, form owner-local stable candidate principal modes
  after `S_AA` rank cleanup and classify `y_i = A_tilde v_i`, not raw GTO
  columns;
- globally merge provisional injected modes, drop duplicate injected modes by
  an explicit injected-overlap rank rule, and validate rank/condition of
  `B = G' S Y_inj`;
- construct the injected gausslet sector as the replacement sector
  `F = [Y, G Q_perp]` or an equivalent numerically stable representation;
- residualize true RG candidates against `F`, not against the original `G`;
- keep true RG residual selection owner-local, followed by one final
  inter-owner merge;
- keep injected directions out of residual-GTO/MWG channels;
- transform exact one-body/moment/unit-nuclear operators into the in-memory
  `[F, R]` basis when injection is enabled;
- preserve inherited gausslet-sector IDA for the injected base sector as the
  explicitly documented two-body approximation;
- continue using residual MWG/IDA only for true residual directions;
- carry only compact numerical authority needed for the injected base sector
  and true residual transforms; avoid raw inventories, discarded spectra,
  broad reports, or status payloads.

Approved object/policy facts, if needed:

- `residual_injection_cutoff`;
- candidate-overlap rank threshold values;
- injected-overlap rank threshold value;
- injected dimension and optionally compact injected owner counts;
- numerical authority for the injected base-sector transform or equivalent
  low-rank representation;
- existing `T_G`/`T_A` authority for true residual directions.

Forbidden:

- default-on injection;
- driver input or canonical driver workflow changes;
- public API/export changes;
- artifact schema/provenance/reader/manifest changes;
- writing injection-enabled supplemented artifacts without a later provenance
  amendment;
- full HF, dense Vee, solver workflow, Cr2 full Hamiltonian, Cr2 artifact, or
  Cr2-specific workflow;
- automatic pruning by kinetic or `H1_RR` spectrum;
- kinetic/`H1_RR` spectral-guard implementation;
- MWG descriptors or residual-MWG channels for injected directions;
- global residual selection;
- width/zeta filtering default changes;
- route, shellification, terminal-lowering, raw-block, Qiu-White, or
  Hamiltonian artifact writer rewrites;
- committed tests or fixtures unless separately approved.

Validation for the source pass:

- `git diff --check`;
- package load;
- current H2 residual-GTO/MWG endpoint unchanged with injection disabled;
- ignored default-off replay showing the source path matches current
  production RG counts and residual spectra for the target Cr/Cr2 fixture;
- ignored enabled-injection replay showing finite/symmetric in-memory
  one-body operators, `F' S F`, `F' S R`, `R' S R`, injected counts, true RG
  counts, `B` rank/condition, and `K_RR`/`H1_RR` spectra;
- no full HF and no new Hamiltonian artifact.

Failure rule: if implementing injection requires artifact schema changes,
driver/public API changes, source outside the approved files, raw-block
rewrites, terminal-basis changes, or a broad payload/report framework, make no
source commit and report the blocker. If exact one-body transformation into
`[F, R]` cannot be done without storing an unacceptable dense `nG x nG`
workspace persistently, stop and request a compact-transform design amendment.

## Do Not Confuse

- Candidate GTO overlap rank is not residual occupation.
- Residual occupation is not numerical rank and not residual integral weight.
- Injected directions are not residual Gaussians.
- Injected directions are not appended. In the historical direct-injection
  path they replace a subspace of `G`; in the protected-original compact-main
  design they replace a subspace of `M = [G, R_compact]`.
- Global injection merging is not permission for global residual selection.
- Exact one-body injection is not a new residual-MWG density convention.
- A kinetic or `H1_RR` spectral guard remains a later safety gate, not the
  first automatic pruning rule.

## Source Authority Status

`HP-RG-INJECT-FN-01` approves only the default-off in-memory source lane
above as historical `G`-injection authority. It does not approve the protected
compact-main construction, changing the production default, artifact
provenance, driver workflow, public API, full HF, Cr2 artifact/workflow, or
spectral pruning policy.

## HP-RG-PROTECT-INJECT-DESIGN-01 - Protected-Original Injection Over Compact Main Space

Status: approved design authority only. This is a docs-only amendment before
any implementation blurb. It approves no source edits, no tests, no artifact
schema/provenance changes, no driver input, no public API, and no Cr2
production claim.

### Purpose

The previous injection framing treated injection as replacement inside the
original gausslet sector `G`. That is not the right construction for the
compact-first Cr2 path. The compact-first selector first builds a protected
local correction space from narrow residual Gaussians. Injection should then
ask whether the original Gaussian supplement can be represented by that
improved main space.

Plainly:

```text
build compact RGs first;
define M = [G, R_compact];
inject original Gaussian functions by replacing directions inside M;
do not make broad non-injectable originals into MWG residual Gaussians.
```

### Definitions

- `G`: the original orthonormal terminal gausslet/final-PQS basis.
- `R_compact`: compact or narrow residual-Gaussian functions selected first by
  the existing ordered compact-first selector. This design does not change the
  current ordered selector behavior.
- `M = [G, R_compact]`: the compact main space used as the parent space for
  injection. `M` is the space whose directions may be replaced.
- `A_all`: the original supplement Gaussian candidates, before residualizing
  against `G` or `M`.
- `A_protected`: the original Gaussian functions corresponding to the accepted
  compact/narrow RGs. These are protected originals.
- `A_broad`: remaining original supplement Gaussian candidates after the
  protected originals are identified.
- `Z_protected`: an orthonormal basis for `A_protected` in the original GTO
  overlap metric. It is built without subtracting `M`.
- `Z_broad`: accepted remaining original directions after orthogonalization
  against `Z_protected`, Gaussian Gram cleanup, and representability testing.
- `Z = [Z_protected, Z_broad]`: the full injected original-Gaussian block.
- `B = M' S Z`: the projection of injected originals into the compact main
  space.
- `Q_perp`: an orthonormal complement to `B` inside the coordinate space of
  `M`, satisfying `B' Q_perp = 0`.
- `F = [Z, M Q_perp]`: the injected fixed sector. This is replacement, not
  append.

### Algorithm

1. Build compact/narrow RGs first using the existing ordered compact-first
   selector. Do not change that selector in this design pass.
2. Define `M = [G, R_compact]`.
3. Use the original supplement Gaussians as injection candidates, including
   originals corresponding to accepted compact RGs.
4. Put protected narrow originals first. Orthonormalize them among themselves
   in original GTO overlap. Do not subtract `M` from this protected block.
5. Orthogonalize all remaining original Gaussians against the protected block.
6. Gram-rank-clean the remaining block in its own Gaussian overlap metric,
   using a candidate-overlap rule such as
   `max(candidate_overlap_atol, candidate_overlap_rtol * maxeig)`. This
   removes linearly dependent Gaussian directions only.
7. Test injection representability in the compact main space by forming
   `B = M' S Z`. `B` must be full rank and acceptably conditioned.
8. If representability passes, inject by replacement:

   ```text
   F = [Z, M Q_perp]
   ```

   not by appending `Z` to `M`.
9. Build any remaining MWG residual channels only from compact/local
   directions that remain approved as true residuals. Broad non-injectable
   candidates must not become MWG RGs.

### Gaussian Gram Cleanup Versus Injection Representability

These are different tests:

- Gaussian Gram cleanup asks whether an original Gaussian candidate direction
  is a real independent direction in the supplement overlap metric. Tiny
  eigenvalues here mean raw Gaussian linear dependence or numerical junk.
- Injection representability asks whether a real Gaussian direction is stably
  represented by the compact main space `M`. This is tested by `B = M' S Z`.

A direction can pass Gaussian Gram cleanup and still fail injection
representability. That failure must not be interpreted as an opportunity to
make a broad residual-Gaussian/MWG channel.

### Failure Interpretation

If a good-norm original Gaussian direction is not stably represented by
`M = [G, R_compact]`, stop and report:

```text
insufficient compact main-basis support for desired Gaussian addition
```

For Cr2 `lmax = 2`, this is a useful diagnostic. At small `ns` such as `4` or
`5`, `d`-like original directions may be real Gaussian directions but poorly
represented by the current main basis. The correct action is to report the
failed owner/channel and improve the main gausslet basis, for example by
increasing `ns`, not to force a broad RG/MWG residual.

### Protected-Span Preservation

The invariant is protected-span preservation, not exact column identity.
`Z_protected` should remain the protected narrow original span. A final
well-conditioned Lowdin or inverse-square-root cleanup may make a tiny
rotation, but diagnostics must show that the protected subspace overlap before
and after cleanup is near identity and that broad directions did not
substantially rotate the protected span away.

Required diagnostic:

```text
sigma(proj(final_Z_protected_span, initial_Z_protected_span))
```

or an equivalent principal-angle/subspace-overlap report, plus the final
orthonormality and condition of the cleanup.

### Cr2 `lmax = 2` Diagnostics

A Cr2 `lmax = 2` protected-original injection audit or future implementation
handoff must report:

- `ns`, `lmax`, owners, and candidate labels/channels;
- protected original counts by owner and angular channel;
- Gaussian Gram eigenvalue ranges and discarded Gram-null directions;
- broad remaining counts by owner/channel after protected-block
  orthogonalization;
- failed representability directions by owner/channel, especially `d`-like
  channels;
- rank and condition of `B = M' S Z`, globally and by useful owner/channel
  summaries;
- protected-span preservation before and after final cleanup;
- final `F' S F`, `F' S R`, and `R' S R` errors if an in-memory construction
  is attempted;
- explicit statement that broad non-injectable candidates were not converted
  into MWG residual channels.

### Forbidden In This Design

- source edits;
- public API or driver changes;
- artifact/provenance/schema changes;
- changing current ordered compact-first selector behavior;
- turning on the existing `G`-injection implementation as-is;
- broad non-injectable candidates becoming MWG RGs;
- Cr2 production claims.

## HP-RG-PROTECT-INJECT-FN-01 - Staged Protected-Original Geometry Prototype

Status: approved narrow source authority for an internal, default-off,
in-memory geometry prototype only. This is not public driver/API authority, not
artifact/provenance authority, not a Hamiltonian writer path, not Cr2 HF, and
not a production default.

### Measurement Basis

The staged-filter report in
`docs/src/developer/reports/cr2_staged_subspace_filter_870498b54/` changes the
source direction. Scalar per-mode fake/representability cuts and combined
score prefixes failed because they allowed collectively weak injected
subspaces. The viable geometry applied subspace filters in this order:

```text
1. representability subspace filter in W using B = M' S W;
2. optional localization/shape filter;
3. fake-RDM eigenspace filter in the surviving subspace.
```

The best measurement variant was:

```text
s_cut       = 0.95
shape       = none
occ_cut     = 0.003
broad dim   = 87
Z dim       = 117
B_min       = 0.9934658245
B < 0.99    = 0
fake trace  = 93.3726973285
```

This is geometry authority only. The result does not prove that a production
Hamiltonian, artifact, or HF workflow should use this construction.

### Approved Source Surface

Approved file:

```text
src/cartesian_residual_gaussians/residual_basis.jl
```

Approved work:

- add private helpers for protected-original staged injection geometry;
- reuse the existing ordered compact-first selector to construct
  `R_compact`;
- identify protected original candidate indices corresponding to accepted
  compact RGs;
- build `M = [G, R_compact]` in the existing mixed `(G,A)` coordinate
  convention;
- build the broad original subspace `W` by orthogonalizing remaining original
  supplement Gaussians against protected originals and Gaussian Gram cleaning;
- apply a representability subspace filter using singular values of
  `B = M' S W`, with private parameters such as `s_cut`;
- optionally localize and classify shape for diagnostics only;
- apply a fake-RDM eigenspace filter with private parameter `occ_cut`;
- form the geometry diagnostics for

  ```text
  Z = [Z_protected, Z_broad]
  F = [Z, M Q_perp]
  ```

- report `B` singular values, fake-RDM trace retained/dropped, protected-span
  preservation, `Z' S M Q_perp`, sampled or block-estimated `F' S F`, and
  dropped-direction summaries.

The implementation may add compact internal helper return records inside
`residual_basis.jl` if needed to avoid parsing labels for compact source
indices. These records must stay private and must not become artifact,
manifest, public API, or driver payload shapes.

### Explicitly Not Approved

- public driver/API/input keywords;
- exports or broad module API;
- artifact schema, artifact writing, provenance keys, or reader support;
- exact one-body or IDA/MWG Hamiltonian transformation for the protected
  geometry;
- Cr2 HF, solver work, or production Cr2 claim;
- screened-reference/rho0 work;
- changing the default residual basis;
- enabling this construction through existing `residual_injection_cutoff`;
- converting rejected broad directions into MWG residual channels;
- source files outside `src/cartesian_residual_gaussians/residual_basis.jl`
  without a later amendment.

### Failure Rule

If the geometry prototype cannot reproduce the staged-filter measurement
within roundoff, or if implementation requires operator/Hamiltonian plumbing,
artifact state, public wiring, or source files outside the approved surface,
stop and report the exact missing object. Do not broaden this ID in source.

Line budget: target at most `220` added source lines in
`src/cartesian_residual_gaussians/residual_basis.jl`. If the implementation
needs a larger helper layer, persistent result object surface, or cross-file
plumbing to reproduce the measured geometry, stop and request a follow-up
amendment.

## HP-RG-PROTECT-INJECT-TEST-01 - Staged Geometry Validation

Approved validation for `HP-RG-PROTECT-INJECT-FN-01`:

- `git diff --check`;
- package load;
- H2 residual endpoint/facade smoke showing default behavior unchanged;
- ignored Cr2 staged-geometry probe comparing source-backed geometry against
  the measurement report values for at least the best variant:

  ```text
  s_cut = 0.95, shape = none, occ_cut = 0.003
  ```

- no Cr2 HF;
- no artifact write/readback for protected-original injection;
- no committed tests unless a later source-review pass requests one.

## HP-RG-PROTECT-ONEBODY-AUDIT-01 - Protected Fixed-Sector One-Body Audit

Status: approved measurement-only audit authority. This is not source
implementation authority and does not approve production Hamiltonian
construction.

### Purpose

The source-backed staged protected-original geometry now reproduces the
measured Cr2 target:

```text
Z = [Z_protected, Z_broad]
F = [Z, M Q_perp]
M = [G, R_compact]
```

The next question is whether exact one-body operators can be transformed
consistently into the injected fixed sector `F`, using only existing geometry
and existing exact one-body data available in the Cartesian/Cr2 construction
path. The audit should clarify the operator dataflow before any source helper
or ownership decision is made.

### Answered Design Questions

The first pass is measurement-only:

- use an ignored `tmp/work/*.jl` probe;
- do not approve a source helper yet;
- do not approve changes to `augmented_operators.jl`,
  `residual_basis.jl`, terminal residual code, artifacts, or public wiring.

If a later source lane is justified, likely ownership for exact one-body
transformation is `src/cartesian_residual_gaussians/augmented_operators.jl`,
not `residual_basis.jl`. `residual_basis.jl` remains the geometry owner. This
audit may discover a missing geometry/export seam, but it must report that
seam rather than adding source instrumentation.

The accepted comparison target is not an energy or a production Hamiltonian.
The target is in-memory consistency of exact one-body blocks in the protected
fixed sector:

- `F' S F - I`;
- finite/symmetric `F' K F`;
- finite/symmetric `F' U_A F` for each available nuclear unit block;
- finite/symmetric `F' H1 F`;
- protected original block `H1` before and after replacement/cleanup;
- compact RG block `H1` before and after replacement;
- trace and low-spectrum diagnostics for `K`, each `U_A`, and `H1`;
- low-mode weights in protected original, accepted broad, and
  `M Q_perp` complement pieces;
- protected-span preservation before and after any final cleanup;
- `B = M' S Z` singular values and `Q_perp` orthogonality diagnostics.

Coordinate and second-moment transforms may be included if they are already
available through the same exact one-body dataflow, but they are not required
for the first audit.

Cr2 is the primary target because it is where the staged protected-original
geometry was measured. A bounded H2 or H2+ sanity run is useful if cheap and
already available, but it is not a committed endpoint or a prerequisite for
the Cr2 audit.

### Approved Surfaces

Allowed:

- ignored `tmp/work/*.jl` probes;
- durable measurement output under `/Users/srw/dmrgtmp/...` or a bounded
  report directory if later accepted by manager review;
- consume the current source-backed staged geometry helper;
- consume existing exact one-body matrices or raw blocks already obtainable in
  the Cr2/Cartesian path;
- build in-memory transformed one-body matrices for `F`;
- report symmetry, orthogonality, trace, low-spectrum, protected-span, and
  block-composition diagnostics.

Forbidden:

- production source changes;
- public driver/API/input wiring or exports;
- artifact schema, provenance, writer, reader, or manifest changes;
- exact IDA/MWG interaction transform;
- screened-reference/rho0 work;
- Cr2 HF, solver work, or production Cr2 Hamiltonian claim;
- changing residual defaults;
- changing the staged geometry selector;
- treating rejected broad directions as MWG residuals;
- committed tests or fixtures.

Failure rule: if the audit cannot construct the fixed-sector one-body blocks
from the existing source-backed geometry and available exact one-body data,
stop and report the exact missing reusable seam. Do not add source
instrumentation, artifact fields, public wiring, or a temporary operator
helper under this audit ID.
