# Residual Gaussian Injection Hybrid Memo

Document role: transitional multi-lane design/history ledger pending subsystem
split. It is not the owner of current status or a single source-authority
contract. Each section's exact registry entry and canonical subsystem document
govern; rejected alternatives and measurements remain here as rationale.

The historical default-off direct-`G` implementation remains a
preservation-only compatibility surface under `HP-RG-INJECT-FN-01`. The
current direction is occupied-first/protected-main injection. Its geometry is
source-backed under `HP-RG-OCC-FIRST-INJECT-FN-01` and
governed by [Occupied-first injection geometry](occupied_first_injection.md).
Protected-original staging, exact one-body transformation, and inherited-site
interaction semantics are governed by
[Protected-localized basis convention](protected_localized_basis.md).
Artifact/locality, retained-GTO EGOI, and ladder facilities have their own
source IDs later in this memo.
`HP-RG-PROTECT-ADDREF-*` governs the prospective combined additive-reference
consumer. This top-level summary does not broaden any of those lanes.

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
Current protected work uses the implemented protected geometry/one-body IDs
and the canonical protected-localized basis contract below; this historical
direct-`G` authority does not broaden them.

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

## Protected Compact-Main Design Contract

`HP-RG-PROTECT-INJECT-DESIGN-01` supplied the compact-first rationale that is
now implemented by the protected geometry and one-body source lanes. Its
canonical numerical contract has moved to
[Protected-localized basis convention](protected_localized_basis.md).

Historical design and measurement evidence remains in manager running-log
Passes 235-253 and in
`docs/src/developer/reports/cr2_staged_subspace_filter_870498b54/`. The durable
rules are replacement over `M = [G, R_compact]`, separate Gaussian Gram and
representability gates, protected-span preservation, and rejection of
unsupported broad directions without creating MWG residual channels.

## Occupied-First Injection History And Contract

`HP-RG-OCC-FIRST-INJECT-AUDIT-01` is completed historical measurement
authority. Passes 323-324 established the occupied-first direction: Be/Ne
occupied subspaces were recovered at roundoff after mandatory inclusion, while
weak optional p-like directions demonstrated why capture-cutoff selection must
remain separate from occupied protection. Numerical evidence remains in the
manager running log.

The implemented geometry and selection contract is now owned by
[Occupied-first injection geometry](occupied_first_injection.md) under
`HP-RG-OCC-FIRST-INJECT-FN-01` and
`HP-RG-OCC-FIRST-INJECT-TEST-01`.

The source-backed helper is not wired into the protected-localized builder and
is not a direct replacement for staged protected-original geometry over
`M = [G, R_compact]`. That composition belongs to
`HP-RG-PROTECT-ADDREF-*` and
[Protected additive atomic reference correction](protected_additive_reference_correction.md).

## Protected Geometry, One-Body, And Interaction History

The implemented staged geometry under `HP-RG-PROTECT-INJECT-FN-01` /
`HP-RG-PROTECT-INJECT-TEST-01` and exact one-body transformation under
`HP-RG-PROTECT-ONEBODY-FN-01` / `HP-RG-PROTECT-ONEBODY-TEST-01` are governed
by [Protected-localized basis convention](protected_localized_basis.md).

The completed one-body audit and source replay remain documented in manager
running-log Passes 254, 255, and 259 and in
`docs/src/developer/reports/cr2_protected_onebody_audit_eaf05a38c/`.

The completed `HP-RG-PROTECT-VEE-AUDIT-01` record is intentionally retained as
negative and positive evidence rather than active measurement authority:

- Pass 269 rejected direct `C' V C` interaction rotation because it fails
  null/projected many-electron energy invariance.
- Pass 270 established localized `L`, exact `H1_L`, and inherited
  pre-injection site-order `Vee_M` as the viable protected-localized
  convention.

Artifact persistence, row-locality metadata, EGOI, ladder facilities, and
rho0/reference-density work remain in their separate sections below and are
not part of this extracted basis contract.

## HP-RG-PROTECT-ART-FN-01 - Protected-Localized Hamiltonian Artifact Variant

Status: implemented.

The artifact identity, native sector/order law, writer/readback behavior, and
exclusions are canonical in
[Protected-localized artifact contract](protected_localized_artifact.md).
The persisted basis numerics remain governed by
[Protected-localized basis convention](protected_localized_basis.md).

## HP-RG-PROTECT-ART-TEST-01 - Protected Artifact Validation

Status: implemented validation contract. Exact validation and historical
evidence are indexed by the
[canonical artifact contract](protected_localized_artifact.md) and registry.

## HP-RG-PROTECT-ARTLOC-FN-01 - Protected Artifact Row-Locality Metadata

Status: implemented.

The native-center calculation, deterministic inverse permutations,
native-sector labels, optional spreads, compatibility behavior, and strict
matrix-order boundary are canonical in
[Protected-localized artifact contract](protected_localized_artifact.md).

## HP-RG-PROTECT-ARTLOC-TEST-01 - Row-Locality Validation

Status: implemented validation contract. Exact validation and historical
evidence are indexed by the
[canonical artifact contract](protected_localized_artifact.md) and registry.

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

## HP-RG-PROTECT-LADDER-BUNDLE-FN-01 - Protected Ladder Bundle Facility

Status: approved opt-in source/artifact authority. This is not driver-default
authority, not solver workflow authority, not a protected-localized
interaction-convention change, and not a Cr2 production claim.

### Purpose

Make same-parent protected-localized Hamiltonian ladders a reusable repo-owned
facility instead of a repeated CR2 one-off. Normal usage should be one
module-qualified call, or a small pair of calls, that builds a directory bundle
containing ladder Hamiltonians, final-basis cross overlaps, optional
transferred-orbital restart data, and human-readable summaries.

The source motivation is the first Cr2 `ns = 7 -> ns = 9` transfer audit:
the transfer was clean enough to make the fixed-density target-energy shift a
useful diagnostic, and repeated protected-localized ladder rebuilds are too
expensive and error-prone to leave as ad hoc consumer scripts.

### Approved Bundle Shape

The preferred durable layout is a directory bundle, not one large JLD2 file:

```text
protected_ladder_bundle/
  manifest.jld2
  members/ns7/protected_localized_hamiltonian.jld2
  members/ns9/protected_localized_hamiltonian.jld2
  transfers/S_ns9_ns7.jld2
  restarts/ns9_from_ns7_occupied_orbitals.jld2   optional
  summaries/ladder_members.tsv
  summaries/transfers.tsv
```

The implementation may choose exact member names, but the manifest must carry
`artifact_kind`, `format_version`, bundle convention ID, member paths, transfer
paths, shared-parent proof, source/current commit facts, and the diagnostics
needed for readback validation.

Bundle members are ordinary protected-localized inherited-site Hamiltonian
artifacts:

- `H1_L` is exact one-body in each protected-localized basis;
- `Vee_L` is the inherited-site IDA/MWG interaction in that basis;
- matrices remain in native protected-localized order;
- existing row-locality metadata and sector maps remain native-order facts.

Transfer sidecars store exact final-basis cross overlaps:

```text
S_BA = <L_B | L_A>
C_B = S_BA C_A
```

Optional restart sidecars may store transferred occupied orbitals in the
target native order with electron counts, source-state provenance, source and
target member IDs, and transfer diagnostics.

### Approved Source Surface

Preferred source owner:

- `src/cartesian_protected_ladder_bundle.jl`;
- `src/GaussletBases.jl` only to include that file.

Approved reuse/wiring surfaces:

- `src/cartesian_ida_hamiltonian.jl` only for protected-localized artifact
  read/write helpers or validation hooks that are genuinely missing;
- `src/cartesian_representation_transfer.jl` for final-basis cross-overlap
  and orbital-transfer helpers;
- `src/cartesian_residual_gaussians/augmented_operators.jl` only if needed to
  expose the in-memory protected-localized representation needed to compute
  cross overlaps without duplicating construction logic.

No package export is approved. A module-qualified internal entry point such as
`GaussletBases.build_protected_localized_ladder_bundle(...)` is allowed if it
stays opt-in and does not alter canonical driver defaults.

### Required Behavior

- build or collect one protected-localized inherited-site Hamiltonian artifact
  per requested ladder member;
- prove shared parent-lattice/geometry/supplement compatibility before writing
  transfer sidecars;
- compute exact final-basis cross overlaps for requested adjacent member
  pairs;
- transfer occupied orbitals or densities using cross overlap only;
- write optional restart sidecars only after trace and orthonormality checks
  pass;
- write manifest/provenance with geometry, parent controls, `ns`,
  `core_spacing`, `basisname`, `lmax`, dimensions, sector counts,
  protected/localized counts, `B_min`, git/source facts, and artifact paths;
- write bounded TSV summaries for human inspection.

### Forbidden

- changing the protected-localized `Vee` convention;
- source-Hamiltonian transforms;
- transforming source `Vee` into target bases;
- `C' V C` or any interaction rotation;
- rho0/P0 revival;
- EGOI expansion or corrected artifact behavior;
- solver/HF/MP2-NO workflow integration;
- default driver behavior changes;
- package exports or broad public API;
- Cr2 production claims.

### Decision Rule

If cross overlaps cannot be computed from existing protected-localized
construction facts without reconstructing large in-memory objects, add only
the smallest source-owned representation seam needed and report it. If the
facility requires changing protected-localized Hamiltonian artifact semantics
or adding new durable fields to those member artifacts, stop unless that exact
field is already approved here. Bundle manifest and sidecar fields are
approved under this lane; existing protected-localized Hamiltonian artifact
semantics are not changed.

## HP-RG-PROTECT-LADDER-BUNDLE-TEST-01 - Protected Ladder Bundle Validation

Status: approved validation gates for `HP-RG-PROTECT-LADDER-BUNDLE-FN-01`.

Approved validation:

- `git diff --check`;
- package load;
- small H, Be, or Be2 bundle/readback smoke if feasible;
- ignored Cr2 `ns = 7 -> ns = 9` bundle validation may verify:
  - shared-parent proof;
  - each member `H1_L`/`Vee_L` finite and symmetric;
  - cross-overlap dimensions and singular spectrum;
  - transferred trace and orthonormality loss;
  - fixed-density target energy reproduces the ladder audit result within
    tolerance;
  - readback roundtrip of manifest, member paths, transfer sidecars, optional
    restart sidecars, and summaries;
- no committed large Cr2 tests or fixtures;
- no production Cr2 HF requirement.

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
