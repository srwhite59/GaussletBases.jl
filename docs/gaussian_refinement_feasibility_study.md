# Gaussian Refinement Feasibility Study

This note records a narrow numerical feasibility study for the proposed
one-dimensional Gaussian-refinement idea behind the distorted-gausslet PGDG
refinement hierarchy framing.

This is not yet a hierarchy implementation pass.

The immediate question is simpler:

- can one coarse Gaussian be represented very accurately by a uniform array of
  finer Gaussians
- at modest refinement levels such as `1/3`, `1/5`, `1/9`
- with coefficients that stay local and numerically stable

## Scope of this pass

The study is intentionally minimal:

- one-dimensional
- one centered target Gaussian first
- one exact continuous `L2` fit formulation
- a finite symmetric fitting window of fine translated Gaussians

Only after that centered base case looks clear does it make sense to ask about
shifted cases or to connect the result back to the full distorted-gausslet
line.

## Fitting formulation

The study uses the repo's Gaussian convention

`exp(-0.5 * ((x - center) / width)^2)`

with:

- coarse target width fixed to `1`
- fine Gaussian width set equal to the refinement spacing `r`
- fine Gaussian centers on the uniform line `k * r`

The coefficients are obtained from the exact continuous `L2` normal equations

`S c = b`

where:

- `S_ij = <g_i, g_j>`
- `b_i = <g_i, g_target>`

and all overlaps are evaluated analytically from the Gaussian overlap formula.

For each refinement level, the study reports:

- relative `L2` fit error
- max pointwise error on a dense check grid
- conditioning of `S`
- coefficient decay away from the center
- the smallest symmetric fitting window that appears saturated

## Current status

The temporary study script is:

- `tmp/work/gaussian_refinement_feasibility.jl`

It solves the exact continuous `L2` fit for the centered target Gaussian and
then reports a practical saturated fitting window using the criteria:

- relative `L2` error below `1e-6`
- edge coefficient ratio below `1e-5`

## Centered results

Target:

- center `0`
- width `1`

Fine fitting Gaussians:

- width = spacing = refinement level `r`
- centers `k r` on a symmetric finite window

Results:

| refinement | chosen half-window | fine functions | relative `L2` error | max pointwise error | `cond(S)` | edge coefficient ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `1/3` | `±5` | `31` | `4.61e-8` | `4.76e-8` | `7.78e3` | `1.03e-6` |
| `1/5` | `±5` | `51` | `1.40e-7` | `2.18e-7` | `8.93e3` | `4.15e-6` |
| `1/9` | `±5` | `91` | `4.01e-7` | `7.63e-7` | `9.43e3` | `8.59e-6` |
| `1/27` | `±6` | `325` | `0.0` at this precision | `5.50e-9` | `9.65e3` | `5.73e-8` |

The coefficient envelopes are broad but smooth and strongly local in the
practical sense: by the chosen windows, the edge coefficients are already
around `1e-6` to `1e-8` of the peak coefficient.

Representative centered coefficient shells `|k| -> max |c_k|`:

- `1/3`: `0 -> 4.23e-1`, `1 -> 3.98e-1`, `2 -> 3.30e-1`, `3 -> 2.41e-1`, `4 -> 1.56e-1`, `5 -> 8.87e-2`
- `1/5`: `0 -> 4.07e-1`, `1 -> 3.99e-1`, `2 -> 3.75e-1`, `3 -> 3.38e-1`, `4 -> 2.92e-1`, `5 -> 2.42e-1`
- `1/9`: `0 -> 4.01e-1`, `1 -> 3.99e-1`, `2 -> 3.92e-1`, `3 -> 3.79e-1`, `4 -> 3.63e-1`, `5 -> 3.43e-1`

So the finer arrays are not sparse in coefficient count, but they are local,
smooth, and extremely accurate.

## Window scan interpretation

The window scans were useful for judging how much spatial support is really
needed.

For the centered case:

- `1/3` already reaches `4.7e-6` relative `L2` error at `±4`, and `4.6e-8` by `±5`
- `1/5` reaches `2.2e-5` at `±4`, and `1.4e-7` by `±5`
- `1/9` reaches `5.0e-5` at `±4`, and `4.0e-7` by `±5`

So a practical reading is:

- `±5` coarse units is already enough for `1/3`, `1/5`, and `1/9`
- `1/27` wants one more step to `±6` if the same edge-decay criterion is kept

## Shifted check

One cheap richer check was also run:

- refinement `1/9`
- target center `0.25`
- target width `1`

That shifted case still gives:

- relative `L2` error `0.0` at this precision with `±6`
- max pointwise error `6.82e-9`
- condition number `9.50e3`

So the good behavior is not restricted to the perfectly centered case.

## Judgment

This looks numerically promising enough to pursue further.

Main takeaways:

- `1/3` already looks very good in the clean centered test
- `1/5` and `1/9` are even more than accurate enough for this Gaussian-scale
  building block
- the fit systems are not perfectly conditioned, but `cond(S) ~ 1e4` is still
  mild enough for a controlled hierarchy study
- the coefficient profiles remain local and smooth rather than wild or
  alternating

So the basic refinement-hierarchy picture

- `1/3`
- `1/9`
- `1/27`

looks numerically plausible as a systematic proxy hierarchy. The next step
should still stay narrow: connect this building block to the distorted-gausslet
proxy line carefully, rather than jumping straight to a full hierarchy
implementation.
