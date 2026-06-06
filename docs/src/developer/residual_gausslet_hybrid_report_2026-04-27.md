# Residual Gausslet Hybrid Report

Date: 2026-04-27

## Purpose

This note summarizes the first residual-sector study for the experimental
high-order nested shell construction on distorted Cartesian parents.

The main question was:

- Can one keep a structured shell sector that preserves the intended shell
  hierarchy, and then add only a small residual correction sector that repairs
  the completeness defect caused by distortion?

The present report covers the work done after the residual-sector idea was
introduced. It includes:

- the two candidate structured-shell constructions that were tested,
- the residual overlap spectra,
- the distortion-to-identity homotopy check,
- He+ energies for the hybrid structured-plus-residual basis,
- and a first look at the occupancies of the earliest residual modes.

The emphasis is on paper-facing definitions and conclusions rather than code
internals.

## Definitions

### Parent basis

The parent basis is the distorted Cartesian gausslet basis on the full parent
lattice. All overlaps, projections, and He+ Hamiltonians in this report are
formed in that parent space.

### Local full block

For a given side length `n_s`, the local full block is the full tensor-product
space of the retained one-dimensional local contractions on the centered
`n_s x n_s x n_s` cube.

Its dimension is `n_s^3`.

Examples:

- `n_s = 5` gives a full local block of dimension `125`
- `n_s = 6` gives a full local block of dimension `216`
- `n_s = 8` gives a full local block of dimension `512`

### Structured shell

At each outer shell level `R`, the structured shell is the subset of the full
local block lying on the tensor shell, meaning that at least one tensor index
is on the local block boundary.

Its dimension is:

- `n_s^3 - (n_s - 2)^3`

Examples:

- `n_s = 5`: shell dimension `98`
- `n_s = 6`: shell dimension `152`
- `n_s = 8`: shell dimension `296`

### Inner span

At shell level `R`, the inner span is the basis already accepted from strictly
smaller shells.

It is denoted conceptually by `I_{R-1}`.

### Accepted structured shell

The structured shell is first projected against the current inner span and then
orthonormalized in the parent overlap metric.

This produces the accepted shell sector:

- the part of the shell that is genuinely new relative to the interior.

### Full local target

At shell level `R`, the full local target is the entire local full block on the
current cube.

This is the space whose completeness one wants to reproduce.

### Leftover / residual seed

After adding the accepted structured shell to the inner span, the full local
target is projected against the accumulated accepted space.

The remaining columns are the raw residual seed:

- the part of the full local target not captured by shell plus interior.

### Residual Gram matrix

The unnormalized leftover is denoted by `Y`.

Its parent-overlap Gram matrix is

- `G = Y' S Y`

where `S` is the parent overlap matrix.

The eigenvalues of `G` are squared residual singular values:

- `lambda_i = sigma_i^2`

This is why:

- singular-value cutoff `1e-5` corresponds to eigenvalue cutoff `1e-10`
- singular-value cutoff `1e-6` corresponds to eigenvalue cutoff `1e-12`

The present report mostly discusses the eigenvalues `lambda_i`, because those
were the directly measured quantities.

### Residual gausslets

Residual gausslets are obtained by diagonalizing the residual Gram matrix,
keeping the eigenmodes above a chosen cutoff, and orthonormalizing the retained
leftover directions.

These functions are not part of the structured shell sector. They are the
completeness-repair sector.

The name "residual gausslet" in this memo is one-body/completeness shorthand.
It does not mean these residual directions automatically carry ordinary
gausslet IDA, quadrature, or retained-weight semantics. Any interaction use
would need its own operator contract and validation.

## Two different structured-shell constructions

Two different candidates were tested.

### A. Reference-rule shell

This construction builds the shell transform on the undistorted reference
lattice and then applies the same coefficient map to the distorted parent
gausslets.

This was a natural first experiment because it preserves the exact undistorted
subspace identities in reference coordinates.

### B. Distorted-shell hybrid

This construction builds the structured shell directly from the distorted local
cube, then uses the distorted full local cube as the target whose completeness
is to be repaired.

This is the construction most closely aligned with the residual-gausslet idea:

- structured shell from the distorted problem itself,
- residuals only for the part of the distorted full block not already captured.

## Main discriminator: which shell sector gives a small residual?

The crucial first result is that the reference-rule shell is not the right
residualization sector.

In contrast, the distorted-shell hybrid produces small residual norms and exact
zero residuals in the first shell step.

### Reference-rule shell: negative control

Selected largest-shell residual-overlap eigenvalues:

| Case | Largest residual eigenvalue |
|---|---:|
| `11^3`, `n_s = 5`, outer side `11` | `7.140944e-01` |
| `14^3`, `n_s = 4`, outer side `14` | `8.677589e-01` |
| `14^3`, `n_s = 6`, outer side `14` | `9.292356e-01` |
| `14^3`, `n_s = 8`, outer side `14` | `9.345860e-01` |

These values are far too large to represent a small distortion defect. This
construction was therefore treated as the wrong sector for residualization.

### Distorted-shell hybrid: positive result

Selected results for the distorted-shell hybrid:

| Case | First added shell | Largest residual eigenvalue there | Largest outer shell | Largest residual eigenvalue there |
|---|---:|---:|---:|---:|
| `5^3`, `n_s = 3` | side `5` | `5.408435e-30` | side `5` | `5.408435e-30` |
| `11^3`, `n_s = 5` | side `7` | `3.874253e-29` | side `11` | `1.696710e-04` |
| `14^3`, `n_s = 4` | side `6` | `4.757954e-30` | side `14` | `4.092069e-04` |
| `14^3`, `n_s = 6` | side `8` | `3.610439e-26` | side `14` | `9.633565e-05` |
| `14^3`, `n_s = 8` | side `10` | `4.024828e-27` | side `14` | `1.368865e-05` |

Interpretation:

- the first added shell is exact to numerical roundoff in every tested case;
- later shells do produce nonzero residuals;
- but the residual norm remains small;
- and the residual scale decreases as `n_s` increases in the tested `14^3`
  family.

## Detailed distorted-shell residual spectra

The following numbers refer to the distorted-shell hybrid and use the
unnormalized residual Gram matrix `G = Y' S Y`.

### `5^3`, `n_s = 3`

- side `5`
  - trace `2.965546e-29`
  - largest eigenvalue `5.408435e-30`
  - no eigenvalues above `1e-10`

This is the exact-zero residual limit.

### `11^3`, `n_s = 5`

- side `7`
  - trace `4.474057e-28`
  - largest eigenvalue `3.874253e-29`
  - no eigenvalues above `1e-10`

- side `9`
  - trace `2.061642e-03`
  - largest eigenvalue `1.245477e-04`
  - `25` eigenvalues above `1e-10`

- side `11`
  - trace `2.356384e-03`
  - largest eigenvalue `1.696710e-04`
  - `27` eigenvalues above `1e-10`
  - `26` eigenvalues above `1e-6`

### `14^3`, `n_s = 4`

- side `6`
  - trace `7.241825e-29`
  - largest eigenvalue `4.757954e-30`
  - no eigenvalues above `1e-10`

- side `8`
  - trace `2.471244e-03`
  - largest eigenvalue `4.144437e-04`
  - `8` eigenvalues above `1e-10`

- side `10`
  - trace `3.344465e-03`
  - largest eigenvalue `5.149137e-04`
  - `8` eigenvalues above `1e-10`

- side `12`
  - trace `3.290599e-03`
  - largest eigenvalue `4.886153e-04`
  - `8` eigenvalues above `1e-10`

- side `14`
  - trace `2.833810e-03`
  - largest eigenvalue `4.092069e-04`
  - `8` eigenvalues above `1e-10`

### `14^3`, `n_s = 6`

- side `8`
  - trace `7.880211e-26`
  - largest eigenvalue `3.610439e-26`
  - no eigenvalues above `1e-10`

- side `10`
  - trace `2.088708e-03`
  - largest eigenvalue `1.044928e-04`
  - `48` eigenvalues above `1e-10`

- side `12`
  - trace `2.732001e-03`
  - largest eigenvalue `1.104139e-04`
  - `64` eigenvalues above `1e-10`
  - `56` eigenvalues above `1e-6`

- side `14`
  - trace `2.613858e-03`
  - largest eigenvalue `9.633565e-05`
  - `64` eigenvalues above `1e-10`
  - `57` eigenvalues above `1e-6`

### `14^3`, `n_s = 8`

- side `10`
  - trace `4.620826e-26`
  - largest eigenvalue `4.024828e-27`
  - no eigenvalues above `1e-10`

- side `12`
  - trace `4.548931e-04`
  - largest eigenvalue `1.364444e-05`
  - `96` eigenvalues above `1e-10`
  - `60` eigenvalues above `1e-6`

- side `14`
  - trace `6.481622e-04`
  - largest eigenvalue `1.368865e-05`
  - `192` eigenvalues above `1e-10`
  - `126` eigenvalues above `1e-6`

## Distortion-to-identity homotopy

To test whether the residual really behaves as a distortion defect, a one-parameter
homotopy was introduced:

- `t = 0`: identity mapping
- `t = 1`: full distorted mapping

The monitored quantity was the largest-shell distorted-hybrid residual spectrum.

### `11^3`, `n_s = 5`

| `t` | Largest eigenvalue | Trace | Count above `1e-10` | Count above `1e-6` |
|---:|---:|---:|---:|---:|
| `0.00` | `3.590124e-27` | `8.039041e-27` | `0` | `0` |
| `0.25` | `2.016036e-06` | `2.958576e-05` | `27` | `16` |
| `0.50` | `1.647805e-05` | `2.251688e-04` | `27` | `26` |
| `0.75` | `6.229869e-05` | `8.630527e-04` | `27` | `26` |
| `1.00` | `1.696710e-04` | `2.356384e-03` | `27` | `26` |

This is the expected smooth behavior:

- exact zero in the identity limit,
- then monotone growth as distortion is increased.

## He+ energies for the structured-plus-residual hybrid

The next question was whether the residual sector matters for the one-electron
physics, rather than only for abstract subspace completeness.

For each case, three energies were compared:

1. structured shell basis only
2. structured shell plus residual modes above a chosen eigenvalue cutoff
3. null-reduced full-cube union reference on the same parent

The cutoffs tested were:

- `1e-10`
- `1e-8`
- `1e-6`
- `1e-5`

These are eigenvalue cutoffs on `G = Y' S Y`, not singular-value cutoffs.

### `11^3`, `n_s = 5`

- structured only
  - count `419`
  - energy `-1.996519729676`

- full-cube union
  - count `471`
  - energy `-1.997320846514`

- hybrid
  - cutoff `1e-10`: count `471`, energy `-1.997320846514`
  - cutoff `1e-8`: count `470`, energy `-1.997320831810`
  - cutoff `1e-6`: count `469`, energy `-1.997320828459`
  - cutoff `1e-5`: count `469`, energy `-1.997320828459`

Conclusion:

- the residual sector closes essentially the entire structured-to-union gap;
- the result is almost insensitive across the tested cutoffs.

### `14^3`, `n_s = 4`

- structured only
  - count `344`
  - energy `-1.887653192039`

- full-cube union
  - count `376`
  - energy `-1.984437068615`

- hybrid
  - all tested cutoffs gave the same result
  - count `376`
  - energy `-1.984437068615`

Conclusion:

- in this case the residual sector is extremely clean;
- the hybrid exactly reproduces the full-union one-electron energy.

### `14^3`, `n_s = 6`

- structured only
  - count `824`
  - energy `-1.992106669046`

- full-cube union
  - count `999`
  - energy `-1.999190666762`

- hybrid
  - cutoff `1e-10`: count `999`, energy `-1.999190666762`
  - cutoff `1e-8`: count `990`, energy `-1.999190583318`
  - cutoff `1e-6`: count `969`, energy `-1.999163728807`
  - cutoff `1e-5`: count `956`, energy `-1.999163606220`

Conclusion:

- even the coarse `1e-5` eigenvalue cutoff recovers more than `99.6%` of the
  structured-to-union energy gap;
- `1e-10` and `1e-8` are essentially exact on this scale.

### `14^3`, `n_s = 8`

- structured only
  - count `1400`
  - energy `-1.998456850843`

- full-cube union
  - count `1639`
  - energy `-1.999233138926`

- hybrid
  - cutoff `1e-10`: count `1651`, energy `-1.999233137511`
  - cutoff `1e-8`: count `1622`, energy `-1.999233137511`
  - cutoff `1e-6`: count `1556`, energy `-1.999233136617`
  - cutoff `1e-5`: count `1412`, energy `-1.999011785410`

Conclusion:

- `1e-10`, `1e-8`, and `1e-6` all recover essentially the full-union energy;
- `1e-5` is too aggressive in this harder case and only closes about `71.5%` of
  the structured-to-union energy gap.

### Interpretation of the He+ data

The hybrid construction is therefore successful in the following practical
sense:

- the residual sector is small enough in norm to be well controlled;
- adding it recovers nearly all of the one-electron energy lost by keeping only
  the structured shell sector;
- and the needed residual count can still be substantially smaller than the
  full-cube union count, especially once a moderate cutoff is used.

One numerical caveat remains:

- in the `14^3`, `n_s = 8`, `1e-10` case, the hybrid count exceeded the
  orthonormalized full-union count slightly.

This is interpreted as a threshold-level redundancy issue in the shell-by-shell
hybrid construction, not as evidence for a genuinely larger physical subspace,
because the energy matched the full-union value to roundoff.

## First residual modes: `14^3`, `n_s = 8`

Because the residual count may need to remain small in practice, the next step
 was to inspect the earliest residual modes directly.

The first ten retained residual modes all came from the first nonzero residual
shell, which was the `side = 12` shell.

Their eigenvalues were:

1. `1.364444455551e-05`
2. `1.363487797119e-05`
3. `1.363487797119e-05`
4. `1.363487797119e-05`
5. `1.363106695945e-05`
6. `1.363106695945e-05`
7. `7.583188220424e-06`
8. `7.583188220424e-06`
9. `7.583188220424e-06`
10. `7.579328793777e-06`

### Energy lowering from the first ten residual modes

- structured only energy:
  - `-1.998456850843`

- after adding the first residual mode:
  - `-1.998899691714`
  - energy lowering `4.428408712838e-04`

- after adding modes `2` through `10`:
  - the energy stayed at `-1.998899691714` to numerical precision
  - each further incremental change was only roundoff-scale

This means that, in this case, the first residual mode captures essentially all
of the observable one-electron gain available in the first ten residual modes.

### Final occupancies of the first ten residual modes

Occupancies were examined in the final `k = 10` hybrid ground state and are
reported here in descending order, as is customary in orbital-analysis work.

Sorted occupancies:

1. `1.485859848054e-06`
2. `2.711018484582e-27`
3. `1.202701724123e-27`
4. `2.045674522442e-28`
5. `1.774649602894e-28`
6. `5.854855036402e-30`
7. `3.653238334367e-30`
8. `7.477061783768e-31`
9. `4.976518884810e-31`
10. `1.323004552286e-32`

### Interpretation of the occupancy data

The occupancy analysis must be interpreted carefully.

The direct final occupancies are tiny, even for the mode that produced the
entire first energy drop. This means:

- the residual mode is not acting as a heavily occupied extra orbital-like
  direction;
- instead, admitting it changes the variational subspace and allows the later
  structured-shell directions to reorganize more effectively.

So, in this case, energy importance and final occupancy are not strongly
correlated.

The first residual mode is important because it changes the accessible subspace,
not because it ends up carrying a large final occupation weight.

## Current conclusions

1. The residual-gausslet idea is mathematically viable if the structured shell
   is built from the distorted local block itself.

2. The reference-rule shell is not the right correction sector:
   its leftover eigenvalues are too large to be called a small residual.

3. In the distorted-shell hybrid:
   - the first added shell is exact to roundoff in every tested case,
   - later residuals are small in norm,
   - and the residual spectrum goes smoothly to zero as the distortion is
     reduced toward identity.

4. For He+:
   - the hybrid residual sector recovers almost all of the structured-to-full
     energy gap in every tested nontrivial case,
   - and in the easier cases the recovery is almost cutoff-insensitive.

5. For the hardest tested case (`14^3`, `n_s = 8`):
   - a very coarse eigenvalue cutoff `1e-5` is too aggressive,
   - but `1e-6` remains essentially safe.

6. The first-residual-mode analysis suggests that only a very small number of
   residual modes may matter energetically, even when many residual modes remain
   visible in the overlap spectrum.

## Immediate next questions

The next scientifically useful steps are:

1. Repeat the “first residual modes” analysis for:
   - `11^3`, `n_s = 5`
   - `14^3`, `n_s = 6`

2. Determine whether energy importance correlates better with:
   - shell of origin,
   - residual eigenvalue,
   - or coupling to the parent He+ ground state,
   rather than with final occupancy alone.

3. Decide on a practical residual-selection rule:
   - fixed eigenvalue cutoff,
   - fixed singular-value cutoff,
   - capped residual count per shell,
   - or some hybrid criterion.
