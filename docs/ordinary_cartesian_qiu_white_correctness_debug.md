# Qiu-White Reference Correctness Debug Note

This note records the first correctness-focused pass after the Qiu-White
reference route became fast enough to run the He `1s^2` checks cleanly.

The structural and timing work was real and useful, but it also exposed two
separate correctness problems:

- the nearest-center / GGT branch still crashed
- the reported scalar outputs could be physically misleading or outright wrong

The priority in this pass is therefore correctness before any further
optimization.

## What caused the nearest/GGT crash

The nearest-center branch failed in
`_qwrg_interaction_matrix_nearest(...)` with a `DimensionMismatch`.

The cause was a stale row/vector assignment:

- the code wrote `transpose(gausslet_interaction[:, index])` into a vector
  slice
- that produced a shape mismatch once the nearest/GGT branch actually ran

The fix was simply to assign the column vector directly.

## What caused the near-zero example outputs

The alarming near-zero outputs from the public Qiu-White example were not, by
themselves, evidence that the Qiu-White path had collapsed.

That example was still using a deliberately truncated Coulomb expansion
(`nterms = 3`) as a structural smoke check. On that same truncated setup, the
simpler ordinary path also gives very small one-body and interaction scalars.

So the very small example values were mainly a benchmark-definition problem:

- the example was structural
- the `1.25 Eh` He `1s^2` target is the full-expansion physical benchmark

The example output is now labeled more explicitly so that this truncated check
is not confused with the physical benchmark.

## What caused the full-expansion collapse

Once the full light reference case was run, the real correctness issue showed
up.

Two problems were contributing:

1. The active Qiu-White path was still using unsafe proxy data in places where
   the residual-space and MWG moment construction needed more faithful raw
   matrices.
2. The residual-space construction was retaining every tiny positive residual
   overlap eigenvalue as if it were a real residual Gaussian direction.

The first fixes in this pass were:

- use the midpoint-based `ga` overlap and one-body-side cross blocks again for
  the active Qiu-White path
- use exact gausslet-side `x^2` data for MWG moment extraction instead of the
  shared proxy `x^2` bundle entry

Those fixes removed the worst residual-space corruption, but the full light
case still collapsed if every tiny residual-overlap eigenmode was kept.

On the light He case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- `cc-pVTZ`
- full `doacc = false` Coulomb expansion

the residual-overlap eigenvalues were:

`[2.617e-4, 3.963e-4, 1.333e-3, 2.678e-3, 6.912e-3, 1.449e-2]`

Keeping all six of those directions produced a catastrophically unstable
one-body spectrum.

The current reference constructor therefore treats the very small residual
eigenmodes as numerical null modes and keeps only residual-overlap eigenvalues
larger than `1e-2`.

This is a debugging/reference safeguard, not a final scientific claim about
the best residual-space threshold.

## Current light-case status after the fix

On the same full light reference case:

- surrogate MWG path:
  - `E1 = -1.9997595432512447`
  - `⟨Vee⟩ = 0.6597491338297566`
- Qiu-White nearest / GGT:
  - `E1 = -2.5258145386216424`
  - `⟨Vee⟩ = 1.282969072925221`
- Qiu-White MWG:
  - `E1 = -2.5258145386216424`
  - `⟨Vee⟩ = 2.092184338838476`

So the reference path is no longer numerically absurd:

- the nearest/GGT branch runs
- the catastrophic one-body collapse is gone
- the `1s^2` scalar is back in a physically recognizable range

But the path is not yet fully trustworthy:

- the one-body orbital is still too low compared with the hydrogenic `-2`
  target
- MWG still overshoots the `1.25` interaction benchmark badly

That means the next remaining correctness target is narrower than before:

- not residual-space collapse anymore
- but the quantitative accuracy of the Qiu-White raw one-body / MWG
  interaction construction

## Regression guard added

The slow test suite now includes a stronger full-expansion nearest/GGT scalar
guard on the light He case. It does not claim final scientific accuracy, but
it does prevent the earlier catastrophic collapse from passing silently again.
