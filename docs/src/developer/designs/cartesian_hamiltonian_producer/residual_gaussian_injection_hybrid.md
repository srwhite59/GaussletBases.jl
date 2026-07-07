# Residual Gaussian Injection Hybrid Memo

Status: design memo plus approved measurement-only audit authority under
`HP-RG-INJECT-AUDIT-01`, historical default-off `G`-injection authority under
`HP-RG-INJECT-FN-01`, and the current protected-original compact-main design
authority under `HP-RG-PROTECT-INJECT-DESIGN-01`. The protected-original design
is the current direction for compact-first RG/injection work. Measurement-only
protected fixed-sector one-body and Vee audits are approved under
`HP-RG-PROTECT-ONEBODY-AUDIT-01` and
`HP-RG-PROTECT-VEE-AUDIT-01`; protected-localized EGOI measurement is
approved under `HP-RG-PROTECT-EGOI-AUDIT-01`, and the first retained-GTO
local-product EGOI helper is approved under `HP-RG-PROTECT-EGOI-FN-01`;
same-parent protected-localized ladder transfer is approved as a
measurement-only audit under `HP-RG-PROTECT-LADDER-XFER-AUDIT-01`;
rho0/Galerkin IDA correction measurement is approved under
`HP-RG-RHO0-GAL-AUDIT-01`, with the successor reference-density-matrix target
recorded in
`rho0_reference_density_matrix.md` under `HP-RHO0-REFDENS-AUDIT-01`. This
document does not approve source edits for the rho0/reference-density design,
production defaults, artifact schema changes, driver inputs, public API, full
HF, solver workflow, or Cr2 production workflow.

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

## HP-RG-PROTECT-ONEBODY-FN-01 - Protected Fixed-Sector Exact One-Body Transform

Status: approved narrow internal source authority. This is not public
driver/API authority, not artifact/provenance authority, not IDA/MWG
interaction authority, not Cr2 HF, and not a production default.

### Measurement Basis

The audit recorded in
`docs/src/developer/reports/cr2_protected_onebody_audit_eaf05a38c/` showed
that the source-backed staged geometry can receive exact one-body operators
from existing in-memory data:

```text
F' K F
F' U_A F by nuclear center
F' H1 F
```

The source geometry matched the staged-filter target with `Z = 117`,
`F = 6945`, `B_min = 0.993465824505872`, and `B < 0.99 = 0`. The in-memory
checks reported `F' S F - I` block estimate `1.164e-9`,
`Z' S M Qperp = 9.873e-16`, finite/symmetric one-body blocks, and converged
low `H1_FF` values. The lowest `H1_FF` modes were dominated by protected
originals, with no obvious low-`H1` broad injected mode before any IDA/MWG
design.

### Approved Source Surface

Primary file:

```text
src/cartesian_residual_gaussians/augmented_operators.jl
```

Optional only if transform-ready geometry fields or accessors are missing:

```text
src/cartesian_residual_gaussians/residual_basis.jl
```

`residual_basis.jl` remains the geometry owner. `augmented_operators.jl` owns
the exact one-body transformation once the geometry is available.

### Approved Behavior

- add private helper(s) to transform exact dense one-body matrices/blocks into
  the protected fixed sector `F = [Z, M Qperp]`;
- consume the source-backed protected geometry from `residual_basis.jl`;
- support exact kinetic `K`, per-center uncharged nuclear `U_A`, and assembled
  `H1`;
- produce in-memory dense transformed matrices and diagnostics only;
- keep the path internal/default-off and unreachable from public driver/API or
  artifact writing;
- preserve the existing default Residual Gaussian behavior and existing
  `[G, R]` exact operator transforms when protected-original geometry is not
  explicitly used;
- report or expose internal diagnostics needed for symmetry,
  orthogonality, low-spectrum replay, protected-span behavior, and block
  composition.

The first source lane is a dense in-memory transform lane. It does not approve
a new matrix-vector action framework. If Cr2 replay requires a general action
interface, stop and request a follow-up amendment.

### Explicitly Not Approved

- public driver/API/input keywords or exports;
- artifact schema, provenance, writer, reader, manifest, or sidecar changes;
- exact IDA/MWG interaction transform;
- screened-reference/rho0 work;
- Cr2 HF, solver work, or production Cr2 Hamiltonian claim;
- residual default changes;
- staged geometry selector changes;
- using rejected broad directions as MWG residual channels;
- source files outside the approved surface;
- committed tests or fixtures by default.

### Failure Rule

If exact protected fixed-sector one-body transformation cannot be implemented
as private internal helpers in `augmented_operators.jl`, with only narrow
geometry access in `residual_basis.jl` if needed, stop and report the missing
object. Do not add artifact state, public wiring, route/terminal/raw-block
changes, a matrix-action framework, or IDA/MWG interaction plumbing under this
ID.

Line budget: target at most `180` added source lines across the approved files.

## HP-RG-PROTECT-ONEBODY-TEST-01 - Protected One-Body Transform Validation

Approved validation for `HP-RG-PROTECT-ONEBODY-FN-01`:

- `git diff --check`;
- package load;
- H2 default residual/facade smoke unchanged;
- ignored Cr2 one-body replay reproducing the audit geometry:
  protected `30`, broad `87`, `Z = 117`, `F = 6945`,
  `B_min = 0.993465824505872`, and `B < 0.99 = 0`;
- replay `F' S F - I`, `Z' S M Qperp`, one-body symmetry, trace, and low
  `H1_FF` diagnostics against the audit report within reviewed numerical
  tolerances;
- finite/symmetric dense in-memory `K`, per-center `U_A`, and `H1` transformed
  blocks;
- no protected-original injection artifact write/readback;
- no IDA/MWG interaction transform;
- no Cr2 HF;
- no committed tests unless a later source-review pass requests one.

## HP-RG-PROTECT-VEE-AUDIT-01 - Protected Fixed-Sector Vee Interaction Audit

Status: approved measurement-only audit authority, now interpreted through the
recorded audit results. The direct `C' V C` protected interaction transform
was invalidated and must not be reused. The viable convention from the later
protected-localized probe is angular-gausslet-style: build a protected
localized injected basis `L`, use exact one-body operators in `L`, and inherit
the pre-injection site-order `Vee_M` interaction. This is not source
implementation authority, not source-backed IDA/MWG authority, not artifact
authority, and not production Hamiltonian authority.

### Purpose

After `HP-RG-PROTECT-ONEBODY-FN-01`, exact dense one-body transforms into the
protected fixed sector are source-backed for

```text
F = [Z, M Qperp]
M = [G, R_compact]
```

The next question is whether an in-memory electron-electron interaction
candidate for the same protected-original fixed sector is numerically sane
before any source implementation. The audit should test the interaction
interpretation, block costs, and broad-injection incentives without changing
production Hamiltonian construction.

### Interaction Interpretation

For this audit only:

- the `G`/base part of `M` keeps the current IDA interaction;
- `R_compact` keeps the current compact residual/MWG interaction;
- injected original directions replace a subspace of `M`; they are not
  appended on top of it;
- rejected broad directions do not become residual MWG channels;
- the in-memory Vee candidate is built by transforming the available
  `M`-space interaction through `F = [Z, M Qperp]`.

This is a measurement model for the protected-original interaction problem. It
does not approve a production IDA/MWG convention or a Hamiltonian workflow.

### Approved Surfaces

Allowed:

- ignored `tmp/work/*.jl` probes only;
- outputs under `/Users/srw/dmrgtmp/...`;
- consume current source-backed protected geometry and one-body helpers;
- consume existing in-memory interaction data already available through the
  Cr2/Cartesian path;
- build in-memory transformed Vee candidates for `F`;
- compute diagnostics only;
- optionally run one bounded in-memory Cr2 HF replay in the same measurement
  lane, but only after the Vee diagnostics pass the gate below.

Forbidden:

- tracked source edits;
- public driver/API/input wiring or exports;
- artifact schema, provenance, writer, reader, manifest, or sidecar changes;
- production Hamiltonian workflow;
- source-backed IDA/MWG interaction implementation;
- screened-reference/rho0 work;
- treating Vee scaling as the primary fix;
- treating rejected broad directions as MWG residual channels;
- Cr2 production claims;
- committed tests or fixtures.

### Required Diagnostics

The audit must report:

- geometry dimensions: `G`, `R_compact`, `M`, `Z`, `Qperp`, and `F`;
- singular spectrum of `B = M' S Z`;
- finite and symmetry checks for the Vee candidate;
- block diagnostics by protected-`Z`, broad-`Z`, and `Qperp` sectors;
- diagonal ranges or self-costs for protected-`Z` and broad-`Z` directions;
- low eigenvalues or representative quadratic-form checks;
- broad-`Z` interaction costs compared to compact-RG and default residual
  costs;
- interaction self-costs of low `H1_FF` modes;
- a residual/broad-`Z` occupation incentive proxy;
- comparison to the default bad residual sector, compact-only sector, and the
  protected one-body audit.

### Gate

If the Vee candidate is finite, symmetric, and the broad-`Z` directions are not
anomalously cheap, the same measurement lane may run one bounded in-memory Cr2
HF replay. That replay remains an audit result, not a production workflow or
artifact claim.

If broad-`Z` remains too cheap, Vee has bad low modes, or the transform needs
source/Hamiltonian/artifact plumbing, stop and report the protected-original
interaction design as the blocker. Do not add source instrumentation,
artifact fields, public wiring, Hamiltonian helpers, or source-backed IDA/MWG
plumbing under this audit ID.

### Output Policy

The first pass may use only ignored probes plus `/Users/srw/dmrgtmp/...`
tables. A compact committed report under `docs/src/developer/reports/` is not
required for the initial measurement, but should be added before using the
result to justify source authority.

### Recorded Outcome

The `C' V C` interaction transform for the replacement fixed sector failed
the algebraic null/projected energy checks. A two-index IDA/MWG
density-density matrix cannot be rotated as `V_F = C' V_M C` for arbitrary
`F = [Z, M Qperp]` densities and expected to preserve many-electron energies.
Do not resume that lane or interpret the resulting broad-`Z` occupation as a
basis-fate signal.

The current viable protected interaction baseline is the protected-localized
inherited-site-order convention:

```text
M = [G, R_compact]
L = protected localized injected basis
one-body(L) = exact transformed one-body operators
Vee(L) = inherited pre-injection site-order Vee_M
```

This baseline is judged by bounded physics diagnostics, not by arbitrary
interaction rotation invariance.

## HP-RG-PROTECT-ART-FN-01 - Protected-Localized Hamiltonian Artifact Variant

Status: approved narrow source/artifact authority.

### Goal

Persist the protected-localized injection Hamiltonian as an explicit,
opt-in `.jld2` artifact variant so solver and MP2-NO consumers can resume
from the current Cr2 protected-localized Hamiltonian without rebuilding the
geometry, one-body transform, and inherited-site interaction in memory.

The artifact records the positive angular-gausslet-style convention:

```text
M = [G, R_compact]
Z = protected original injected directions
L = localized protected replacement basis
H1_L = exact one-body Hamiltonian transformed into L
Vee_L = inherited localized-site IDA/MWG interaction matrix in L site order
```

Protected injection is a small localized one-particle improvement. The
IDA/MWG density interaction remains attached to the localized gausslet-like
site basis. This lane does not revive `C' V C`, does not define alternative
interaction rotations, and does not create a general all-injection artifact
contract.

### Approved Source Surface

Primary approved files:

- `src/cartesian_residual_gaussians/augmented_operators.jl`
- `src/cartesian_ida_hamiltonian.jl`

Optional only if directly required to expose already-computed
protected-localized geometry, diagnostics, or producer plumbing:

- `src/cartesian_residual_gaussians/residual_basis.jl`
- `src/cartesian_base_hamiltonian.jl`

No driver or public API surface is approved in this lane. If a later consumer
needs a command-line entry point, design-manager must approve that separately.

### Required Artifact Contract

The artifact must be versioned and convention-identified before any consumer
uses it. Minimal required convention facts:

- artifact variant, for example
  `:protected_localized_injection_hamiltonian`;
- convention ID, for example `:protected_localized_injection_v1`;
- variant/schema version.

Required numerical datasets:

- `H1_L`;
- `Vee_L`;
- `nup`, `ndn`;
- final dimension.

Required provenance and controls:

- source recipe/artifact provenance;
- source commit and current commit;
- public basis controls and geometry inputs;
- basis/injection convention ID;
- compact RG count and selection metadata;
- protected original count;
- broad injected count;
- localized basis ordering.

Required maps and diagnostics:

- sector maps for `G`/base, compact-`R`, protected-`Z`, broad-`Z`, and
  `Qperp`/localized complement;
- representability diagnostics including `B_min`, singular-value thresholds,
  and singular counts;
- orthogonality and localization diagnostics;
- inherited-site interaction diagnostics.

### Readback Checks

Readback must validate:

- recognized artifact variant, convention ID, and version;
- dimension consistency for `H1_L`, `Vee_L`, sector maps, and final dimension;
- finite/symmetric `H1_L`;
- finite/symmetric `Vee_L`;
- `nup`/`ndn` present and consistent with the recorded system;
- sector counts sum to the final dimension;
- source recipe/provenance present;
- `B_min`, singular-count, and orthogonality diagnostics present.

If the existing reader cannot safely distinguish the protected-localized
variant from ordinary Cartesian IDA Hamiltonian artifacts, the implementation
must add the minimal convention/version field and reject unrecognized or
missing conventions before Cr2 consumers use the file. It must not silently
load this variant as a standard PQS/WL/RG artifact.

### Forbidden

- changing default producer behavior;
- replacing existing PQS, WL, or ordinary RG artifact semantics;
- rho0/reference-density correction work;
- broad public workflow, driver flags, or exported API;
- new solver methods;
- paper or production energy claims;
- changes to RG/injection selection policy;
- treating rejected broad directions as MWG residual channels;
- `C' V C` or other alternative interaction rotations;
- artifact schema changes for existing default artifacts;
- Cr2-specific branches.

### Decision Rule

Approve implementation only if it fits as a clearly versioned opt-in artifact
variant of the existing `.jld2` Hamiltonian family. If that is not possible
without broad reader/schema redesign, stop and report the exact reader or
schema blocker rather than creating a compatibility wrapper.

## HP-RG-PROTECT-ART-TEST-01 - Protected Artifact Validation

Status: approved validation gates for `HP-RG-PROTECT-ART-FN-01`.

Approved validation:

- `git diff --check`;
- package load;
- one H or Be small smoke if convenient;
- one bounded Cr2 protected-localized artifact write/readback or
  readback/resume smoke;
- compare loaded `.jld2` `H1_L`, `Vee_L`, dimensions, sector maps, and
  occupations against the current in-memory protected-localized replay;
- confirm unrecognized or missing convention/version fields are rejected;
- no converged Cr2 energy claim required;
- no public workflow, solver method, rho0, or production default validation.

## HP-RG-PROTECT-ARTLOC-FN-01 - Protected Artifact Row-Locality Metadata

Status: approved narrow source/artifact amendment.

### Goal

Add row-locality metadata to the protected-localized Hamiltonian artifact so
Cr2 solver and MP2-NO consumers can make z-ordered working copies without
reconstructing protected-localized geometry in memory or relying on stale
manifest labels.

The canonical matrices remain native-order:

```text
H1_L[native, native]
Vee_L[native, native]
```

Any z-order information is metadata for consumers. It must not silently imply
that `H1_L`, `Vee_L`, sector ranges, or existing sector counts have been
permuted.

### Center Definition

The row centers are diagonal position expectations in the actual
protected-localized working basis. Given inherited main-space position
operators `X_M`, `Y_M`, `Z_M` and the native protected-localized transform
`ML`, define:

```text
center_x = diag(ML' * X_M * ML)
center_y = diag(ML' * Y_M * ML)
center_z = diag(ML' * Z_M * ML)
```

Do not use construction labels, source-region labels, or manifest center
metadata as numerical authority for these centers.

### Approved Source Surface

Primary approved files:

- `src/cartesian_residual_gaussians/augmented_operators.jl`
- `src/cartesian_ida_hamiltonian.jl`

Optional only if directly required to expose already-computed transform or
position data:

- `src/cartesian_residual_gaussians/residual_basis.jl`
- `src/cartesian_base_hamiltonian.jl`

### Required Metadata

Required native-order vectors:

- `center_x`;
- `center_y`;
- `center_z`;
- a per-row sector label or native-sector index, so consumers do not have to
  reinterpret contiguous native sector ranges after z sorting.

Required z-order vectors:

- `z_order_to_native`, where entry `k` is the native row index at z-sorted
  position `k`;
- `native_to_z_order`, the inverse permutation.

The z sort must be deterministic. Sort by `center_z`, then by native row index
as the tie-breaker unless a later source review approves a more specific
localized-order tie rule.

Approved optional diagnostics, when existing second-moment data are already
available without new raw-block or operator construction:

- `spread_x`;
- `spread_y`;
- `spread_z`.

Spreads should be finite, nonnegative diagnostics such as
`sqrt(max(<axis^2> - center_axis^2, 0))`. If second moments are not already
available in the approved path, do not synthesize spreads from labels or add a
new raw-block lane under this ID; record the missing spread diagnostic as out
of scope.

### Readback Checks

Readback must validate:

- locality vector lengths equal the final dimension;
- centers are finite;
- spreads, when present, are finite and nonnegative;
- `z_order_to_native` and `native_to_z_order` are inverse permutations of
  `1:dim`;
- z-order centers are monotone under the recorded permutation, with native
  index tie-breaks;
- per-row sector labels or indices are present and agree with the native
  sector counts;
- native-order `H1_L` and `Vee_L` dimensions and symmetry checks remain the
  authoritative matrix checks.

### Forbidden

- mutating the existing artifact matrix order;
- writing only z-sorted matrices under the existing convention ID;
- reusing native contiguous sector ranges as if they remained contiguous after
  z sorting;
- changing RG/injection selection, localization, or Vee semantics;
- new driver/API/solver workflow;
- new raw-block or position-second-moment construction;
- rho0/reference-density work;
- Cr2 production energy claims.

### Decision Rule

If native position operators or the `ML` transform are not available at the
protected artifact seam, stop and report the exact missing source fact. Do not
fall back to manifest labels or route metadata as numerical center authority.

## HP-RG-PROTECT-ARTLOC-TEST-01 - Row-Locality Validation

Status: approved validation gates for `HP-RG-PROTECT-ARTLOC-FN-01`.

Approved validation:

- `git diff --check`;
- package load;
- one protected-localized artifact write/readback smoke;
- readback confirms center vector lengths, finite centers, permutation inverse
  checks, monotone z order, and sector-label/count consistency;
- if spreads are present, readback confirms finite nonnegative spreads;
- compare stored centers and z permutation against centers recomputed from the
  in-memory `ML' * X_M/Y_M/Z_M * ML` diagonal constructions for the same
  fixture;
- confirm loaded `H1_L`/`Vee_L` still match native-order in-memory replay;
- no Cr2 converged energy claim, public workflow, solver method, rho0, or
  production default validation.

## HP-RG-PROTECT-EGOI-AUDIT-01 - Protected-Localized EGOI Measurement Audit

Status: approved measurement-only audit authority. This is not source
authority, artifact authority, public workflow authority, solver authority, or
a Cr2 production claim.

### Goal

Test whether the existing matrix-level EGOI correction is a better fit for
protected-localized interaction errors than the stalled rho0/reference-density
one-body correction path.

The protected-localized artifact now provides the working objects needed by
the audit:

- `H1_L`;
- `Vee_L`;
- native ordering;
- row-locality metadata;
- sector maps.

The audit may use the existing matrix-level EGOI routines:

- `egoi_target_product_matrix`;
- `egoi_target_coulomb_matrix`;
- `egoi_density_density_correction`;
- `egoi_stationary_hamiltonian_correction`.

### Allowed

- ignored `tmp/work/*.jl` measurement probes only;
- outputs under `/Users/srw/dmrgtmp/...`;
- H, Be, and Be2 first;
- existing EGOI matrix routines;
- reconstruct `Qtarget` from current protected/injection geometry in the
  probe;
- exact Gaussian target Coulomb for the selected target orbitals;
- optional bounded Cr2 diagnostic only after H/Be/Be2 diagnostics look sane.

### Required Diagnostics

The audit must report:

- target definition and target-selection rule;
- `Qtarget` dimension, Gram matrix, and orthogonality diagnostics;
- target representability and projection loss;
- exact target Coulomb construction details;
- EGOI residual before and after correction;
- `DeltaV` max norm, Frobenius norm, and relative size;
- product-matrix singular values and rank;
- corrected `Vee` finite and symmetry checks;
- low Fock spectra before and after correction;
- H one-electron and self-interaction sanity check;
- Be/Be2 corrected behavior compared with the rho0/P0 audit;
- if Cr2 is reached, the exact small-system gate that justified it.

### Forbidden

- tracked source edits;
- artifact/schema/provenance/writer/reader changes;
- EGOI-corrected artifact variants;
- public driver/API/export or solver workflow;
- Cr2 production claims;
- RG/injection selection-policy changes;
- protected-localized artifact convention changes;
- rho0/P0 revival as part of this audit;
- treating rejected broad directions as MWG residual channels.

### Decision Rule

If H/Be/Be2 EGOI reduces target residuals with small or moderate `DeltaV` and
benign Fock behavior, the result may justify a later source/artifact lane for
protected-localized EGOI target metadata and corrected interaction variants.

If EGOI requires large or cancellation-dominated `DeltaV`, creates bad low
modes, or fails the H/Be/Be2 sanity checks, keep it diagnostic-only. Do not
run Cr2 until small cases pass.

## HP-RG-PROTECT-EGOI-FN-01 - Retained-GTO Local-Product EGOI Helper

Status: approved narrow internal source authority.

### Purpose

Turn the successful retained-original-GTO EGOI measurement convention into a
source-backed, in-memory helper without changing public workflow, artifact
semantics, or production defaults.

The target is not broad protected-`Z`, not atom-HF orbitals, not final basis
rows, and not residualized RG functions. The approved first target is:

```text
retained original supplement GTOs mapped from compact retained source indices
owner-balanced retained s1+s2 only
molecular retained-original-GTO target
local M2 mask
symmetric/local-product EGOI formulation
```

### Physical Convention

- Local products on each atom are first-class EGOI products.
- The inter-atom local-product Coulomb block is included in the exact
  target/acceptance metric.
- AB overlap products are not first-class targets by default.
- Long-range or otherwise disallowed `DeltaV` entries remain exactly zero.
- `s3`, `p`, `d`, and broader target classes remain measurement-only until
  separately approved.

### Evidence

- Be2 retained-GTO target ladder:
  `/Users/srw/dmrgtmp/protected_localized_retained_gto_target_ladder_22f051741`.
  `s1+s2` compact targets project at roundoff, reduce residuals by about
  `97.55%`, require small `DeltaV`, and have benign low-Fock shifts.
- Cr2 molecular M2/M3/M4 isolation:
  `/Users/srw/dmrgtmp/cr2_molecular_mask_radius_egoi_22f051741`.
  `s1+s2`/`M2` gives residual reduction about `99.603%`,
  `DeltaV/V` Frobenius about `9.076e-5`, relative `p95` about `2.421e-4`,
  and low-Fock shift about `+2.090e-5` Ha. `M3`/`M4` are not materially
  better.
- Cr2 cross-term scaling audit:
  `/Users/srw/dmrgtmp/cr2_egoi_cross_term_scaling_22f051741`.
  The remaining residual after `s1+s2`/`M2` is the AA-BB local-product
  Coulomb block, `99.09%` diag-diag. AB overlap product blocks are negligible
  and are not promoted as default targets.

### Approved Source Surface

Preferred neutral/internal owner:

- `src/hamiltonian_corrections.jl`

Optional narrow consumer/helper wiring only if required to obtain
protected-localized retained source mapping or transform-ready `Qtarget`:

- `src/cartesian_residual_gaussians/augmented_operators.jl`
- `src/cartesian_residual_gaussians/residual_basis.jl`

No new public export/API is approved. Any new helpers should remain internal
unless a later public workflow lane explicitly approves them.

### Helper Responsibilities

The helper may:

1. build retained original-GTO target metadata from protected geometry/source
   indices;
2. select only the first approved owner-balanced retained `s1+s2` target
   class;
3. build `Qtarget` for those original GTOs in the protected-localized native
   basis;
4. build symmetric local-product target products while excluding AB overlap
   products as first-class products;
5. include the AA-BB local-product Coulomb block in the exact target and
   acceptance metric;
6. build and apply the `M2` local mask:
   `r_ij <= 1.75 * max(ell_i, ell_j)`, where
   `ell_i = max(nearest-neighbor ell_i, core_spacing)`;
7. solve local constrained EGOI using existing matrix-level routines or a
   narrow symmetric/local-product wrapper;
8. return in-memory `DeltaV` and compact diagnostics.

### Required Diagnostics

Return or report compact diagnostics for:

- target labels, source indices, owners, and channels;
- projection loss;
- product counts and symmetric rank/singular values;
- local-product block residuals before and after correction:
  `AA-AA`, `BB-BB`, and `AA-BB`;
- AA-BB diag-diag, diag-offdiag, and offdiag-offdiag split;
- `DeltaV/V` Frobenius norm, max norm, `p95`, and median by variable class;
- saturated variables by class;
- `max_disallowed_delta_v`;
- corrected `Vee` finite/symmetric checks;
- low-Fock shift diagnostic;
- cache/probe parity against the accepted measurement outputs.

### Forbidden

- public API/export/driver workflow;
- artifact/schema/provenance/writer/reader changes;
- corrected protected-localized artifact variants;
- solver/HF/MP2-NO workflow integration;
- Cr2 production energy claims;
- RG/injection selection-policy changes;
- broad protected-`Z` targets;
- atom-HF/P0/rho0 revival;
- AB overlap products as default targets;
- `s3`, `p`, or `d` target promotion;
- committed large Cr2 tests or fixtures.

### Decision Rule

Approve source only while this remains an internal in-memory helper with
compact diagnostics and no artifact/workflow changes. If source work requires
a corrected artifact or solver workflow to be meaningful, stop and keep the
lane measurement-only.

## HP-RG-PROTECT-EGOI-TEST-01 - Retained-GTO EGOI Validation

Status: approved validation gates for `HP-RG-PROTECT-EGOI-FN-01`.

Approved validation:

- package load;
- `git diff --check`;
- H, Be, and Be2 retained-GTO smoke for `s1` and `s1+s2`;
- ignored Cr2 replay matching the accepted measurement:
  - `s1+s2`/`M2` residual reduction about `99.603%`;
  - `DeltaV/V` Frobenius about `9e-5`;
  - relative `p95` `DeltaV` about `2.4e-4`;
  - benign low-Fock shift;
  - `max_disallowed_delta_v = 0`;
- no production Cr2 HF;
- no committed large Cr2 tests or fixtures.

## HP-RG-PROTECT-LADDER-XFER-AUDIT-01 - Same-Parent Ladder Transfer Audit

Status: approved measurement-only audit authority. This is not source
authority, not artifact/schema authority, not public driver/API authority, and
not a Cr2 production claim.

### Purpose

Test whether the current Cr2 protected-localized UHF discrepancy is primarily
basis/contraction/Hamiltonian convergence rather than UHF basin failure. The
starting evidence is a Cr2 UHF state with global spin diagnostics close to the
Yann/Sandeep reference, including `<S^2>` about `4.866` and large AFM local
moments, but energy about `36 mHa` below the cc-pwCV5Z UHF reference.

The audit builds a same-parent protected-localized ladder, for example:

- fixed parent lattice, same supplement, and same Cr2 geometry;
- `ns = 7` protected-localized inherited-site Hamiltonian;
- `ns = 9` protected-localized inherited-site Hamiltonian;
- optional `ns = 11` if affordable.

Then compute exact final-basis cross overlaps:

```text
S_9,7 = <L_ns9 | L_ns7>
S_11,9 = <L_ns11 | L_ns9>   optional
```

and transfer occupied orbitals or density information by cross overlap only:

```text
C_B = S_BA C_A
S_BA = <B | A>
```

The transferred state is evaluated with the target-basis `H1_L` and `Vee_L`.
A few bounded UHF sweeps may run only after trace and orthonormality checks
show that the transfer is meaningful.

### Critical Transfer Convention

Final working bases are intended orthonormal. Transfer must use only the
cross overlap between final bases. Do not use generalized self-overlap
transfer, do not transform source Hamiltonians into the target basis, do not
transform source `Vee`, and do not use `C' V C` or any interaction rotation.
After transfer, all fixed-density and UHF-continuation diagnostics must use
the target protected-localized Hamiltonian.

### Allowed

- ignored `tmp/work` probes;
- outputs under `/Users/srw/dmrgtmp`;
- existing protected-localized inherited-site Hamiltonian construction and
  writer;
- in-memory or ignored sidecar cross-overlap matrices;
- transfer of saved occupied orbitals from one ladder basis to the next;
- fixed-density target-energy evaluation;
- small bounded UHF continuation only if transferred trace and occupied
  overlap checks pass.

### Forbidden

- tracked source edits;
- new public API/export;
- production workflow or driver wiring;
- durable artifact schema changes;
- changes to the protected-localized `Vee` convention;
- transforming source `Vee` into the target basis;
- `C' V C` or any interaction rotation;
- rho0/P0 revival;
- EGOI expansion or corrected artifact behavior;
- Cr2 production claims.

### Required Diagnostics

- exact geometry and shared parent-lattice controls;
- `ns` values and final dimensions;
- protected/localized counts and `B_min` for each basis;
- `H1_L` and `Vee_L` finite/symmetry checks;
- cross-overlap dimensions and singular spectrum;
- transferred electron trace loss and occupied-overlap loss;
- `E_target[P_transferred]` before any sweep;
- returned or recomputed energy if bounded sweeps are run;
- `<S^2>`, local spin diagnostics, and sector occupations before and after
  transfer;
- wall times and output paths.

### Decision Rule

If `ns = 7 -> ns = 9` transfer has small trace/occupied loss and the evaluated
energy moves toward the Yann/Sandeep reference, the discrepancy is likely
basis/contraction/Hamiltonian convergence. If transfer loss is large, the
ladder construction is not comparable and final-basis capture/cross-overlap
diagnostics are needed before physics interpretation. If energy stays too low
after clean transfer and a few bounded sweeps, suspect protected-localized
`Vee`/IDA/EGOI/injection convention accuracy rather than UHF basin failure.

If a reusable source helper or durable artifact sidecar appears necessary, the
audit must report the smallest owner and exact stored fields for a later
source or artifact lane. It must not promote source or artifact changes from
this measurement authority.

## HP-RG-RHO0-GAL-AUDIT-01 - Rho0/Galerkin IDA Correction Audit

Status: approved measurement-only audit authority. This is not source
implementation authority, not source-backed IDA/MWG authority, not artifact
authority, not production Hamiltonian authority, and not a Cr2 production
claim. Later row-gauge audits showed this formulation was algebraically
under-specified as a correction target. Keep this lane as historical
measurement evidence; use `HP-RHO0-REFDENS-AUDIT-01` in
`rho0_reference_density_matrix.md` for the current reference-density-matrix
target.

### Purpose

Test whether a reference-density / Galerkin correction can improve the
existing IDA interaction after the protected-localized injection convention
has removed the broad residual-collapse mechanism. This is an IDA-improvement
lane on top of the sane protected-localized inherited-site baseline. It is not
a basis-fate rule, not a replacement for compact RG/injection selection, and
not permission to revive `C' V C`.

### Baseline

The audit starts from the protected-localized convention:

- build protected-localized injected basis `L`;
- use exact one-body operators in `L`;
- inherit pre-injection site-order `Vee_M` as the IDA/MWG interaction;
- judge by small-system and Cr2 physics diagnostics.

The rho0/Galerkin correction is tested as an additive or replacement
candidate for improving IDA accounting relative to that inherited-site
baseline. The audit must state the exact convention used before reporting
energies or occupations.

### Approved Surfaces

Allowed:

- ignored `tmp/work/*.jl` measurement probes only;
- outputs under `/Users/srw/dmrgtmp/...`;
- in-memory experiments over existing protected-localized geometry and
  one-body/Vee data;
- analytic IDA/Coulomb sanity checks;
- small H, He, and H2 checks;
- Cr2 fixed-density diagnostics;
- one bounded Cr2 HF replay only if static rho0/Galerkin diagnostics are
  sane.

Forbidden:

- tracked source edits;
- public driver/API/input wiring or exports;
- artifact schema, provenance, writer, reader, manifest, or sidecar changes;
- production Hamiltonian workflow;
- `C' V C` interaction transform revival;
- treating rejected broad directions as MWG residual channels;
- treating Vee scaling as the fix;
- screened-reference production claims;
- Cr2 production energy claims;
- publication-scale validation sweeps;
- committed tests or fixtures.

### Required Diagnostics

The audit must report:

- exact rho0 definition and normalization;
- whether rho0 is a spherical one-Gaussian, multi-Gaussian, or fitted atomic
  density;
- Galerkin/reference vector or matrix convention used;
- finite and symmetry checks;
- analytic 1s self-Coulomb checks where applicable;
- H, He, and H2 IDA sanity shifts;
- Cr2 low `H1` / interaction incentive diagnostics;
- fixed-density energy shift on the saved bad density if available;
- bounded HF residual, compact-`R`, and injected-site occupation;
- comparison to the protected-localized inherited-site `Vee_M` baseline.

### Decision Rule

If rho0/Galerkin improves the small IDA sanity checks and Cr2 diagnostics
without reviving broad/residual occupation, the result may justify a later
source-design amendment. That later amendment must name the source owner,
operator convention, artifact/public exclusions, validation, and line budget.

If rho0/Galerkin introduces negative broad residual/interaction modes,
inconsistent energy accounting, or large occupation incentives, stop and
record rho0/Galerkin as the current interaction-design blocker.

This is repo-level algorithm engineering validation. Larger molecule and
convergence sweeps belong to a consumer-oriented workflow after the repo path
is stable.

### Successor Target

The successor lane replaces scalar/row-gauge `rho0` reasoning with a fixed
reference density matrix `P0`. It requires:

```text
Delta_F0_sigma = F_exact0_sigma[P0] - F_app0_sigma[P0]

C0 =
    E_exact0[P0]
  - E_app0[P0]
  - sum_sigma Tr(P0_sigma * Delta_F0_sigma)
```

so that the corrected model satisfies both:

```text
E_new[P0] = E_exact0[P0]
dE_new/dP_sigma at P0 = F_exact0_sigma[P0]
```

Do not continue row-action `(J*w)/w`, `diag(J)`, or scalar `u0/q0` matching as
the acceptance target. Those are diagnostics of different objects, not the
definition of the correction.
