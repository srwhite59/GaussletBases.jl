# Residual Gaussian Injection Hybrid Memo

Status: design memo, not production source authority. This document records
the proposed optional injection-plus-residual construction for near-gausslet
GTO supplement directions. It does not approve source changes, tests, artifact
schema changes, driver inputs, public API, Cr2 workflow, or a production
default.

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

When injection is enabled, `lambda_inj` should be at least the active
`residual_occupation_cutoff`. Otherwise there is an ambiguous band of modes
that are neither injected nor retained as true residuals.

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

Required guards for any future source lane:

- `dim(Y_inj) < nG`;
- the projection of `Y_inj` into the original gausslet span has full rank;
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

This is an approximation for the injected sector, but it is stable in the
limit `lambda -> 0`: a direction already represented by gausslets should not
create a normalized residual function or a residual MWG density.

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

## Measurement-Only First Step

The first approved lane should be measurement-only, likely under a candidate
ID such as:

```text
HP-RG-INJECT-AUDIT-01
```

It should use ignored probes only and report, for Cr atom, Cr2 monomer
counterpoise if available, and Cr2:

- raw candidate counts;
- stable owner-local candidate counts after `S_AA` rank cleanup;
- provisional injected count by owner for trial `lambda_inj`;
- globally retained injected count after duplicate merge;
- true RG count by owner;
- residual occupation spectra before and after injection;
- `K_RR` and `H1_RR` low eigenvalues for true RGs;
- owner weights and residual-occupation composition of low modes;
- GTO-span one-body accuracy or mismatch against the non-injected path where
  cheap;
- whether the low two-owner residual ghost sector disappears, shrinks, or
  persists.

No production source change, committed test, artifact write, full HF, solver,
driver input, public API, MWG convention change, or Cr2 workflow should be
included in the audit lane.

## Do Not Confuse

- Candidate GTO overlap rank is not residual occupation.
- Residual occupation is not numerical rank and not residual integral weight.
- Injected directions are not residual Gaussians.
- Injected directions are not added to `G`; they replace a subspace of `G`.
- Global injection merging is not permission for global residual selection.
- Exact one-body injection is not a new residual-MWG density convention.
- A kinetic or `H1_RR` spectral guard remains a later safety gate, not the
  first automatic pruning rule.

## Source Authority Status

No source authority is granted by this memo. A future docs-only amendment must
approve exact IDs, files, validation, and defaults before implementation.
Likely source ownership, if promoted, belongs in:

```text
src/cartesian_residual_gaussians/residual_basis.jl
```

with the old terminal residual file allowed only for narrow compatibility
keyword plumbing if explicitly approved.
