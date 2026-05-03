# Bounded full-shell basis versus full-block union distorted completeness study

Date: 2026-05-03

This note records a first bounded study of how close the **full-shell basis**
(FSB) is to the **full-block union** (FBU) as distortion is reduced in the
current high-order one-body lane.

## Question

The intended question was:

- how quickly does FSB approach FBU as the mapping distortion scale `s` is
  reduced?

In the undistorted limit, FSB and FBU are expected to agree exactly. The goal
here was to check whether the same is already true, or nearly true, for modest
nonzero distortions in the current same-backend PGDG one-body path.

## Metric definitions

All final working bases were treated as intended orthonormal final bases. In
that setting:

- self-overlaps were used only as diagnostics
- transfer/capture was defined with the cross overlap only

The primary metrics were:

1. He+ energy gap

   `E_FSB - E_FBU`

2. FBU ground-state capture deficiency by FSB

   Let `C_FSB` and `C_FBU` be the final basis coefficient matrices in the
   parent overlap metric `S_parent`. Let `psi_FBU` be the normalized FBU He+
   ground-state coefficient vector in the FBU basis. Then the FSB capture uses
   the cross overlap

   `S_cross = C_FSB' S_parent C_FBU`

   and the deficiency is

   `1 - ||S_cross psi_FBU||^2`

3. Cross-overlap subspace error

   For equal FSB/FBU dimensions, the diagnostic

   `||S_cross' S_cross - I||_inf`

   was reported to show how close the two final subspaces are as orthonormal
   working spaces.

## Study setup

- backend: `:pgdg_localized_experimental`
- target problem: He+
- parent families:
  - `count = 11`, `doside = 5`, sides `5,7,9,11`
  - `count = 13`, `doside = 5`, sides `5,7,9,11,13`
- sweep:
  - identity control using the actual `IdentityMapping()`
  - distorted cases with `AsinhMapping(c = 0.2, s = (s/s0) * sqrt(0.4), tail_spacing = 10.0)`
  - `s/s0 = 1.0, 0.8, 0.6, 0.4, 0.2`

The production high-order stack helper still gates the supported distorted
mapping family to the exact White-Lindsey He baseline map. To avoid widening
that contract just for this study, the sweep driver assembled FSB directly from
the same lower-level physical-block shell helpers already used by the
production route.

Driver:

- [tmp/work/high_order_fsb_fbu_distortion_sweep.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_distortion_sweep.jl)

Saved artifacts from the run:

- [TSV table](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_distortion_sweep_2026-05-03_161107.tsv)
- [text summary](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_distortion_sweep_2026-05-03_161107.txt)

## Results

### Case summary

| case | max `|E_FSB - E_FBU|` | max capture deficiency | max cross-overlap error |
| --- | ---: | ---: | ---: |
| `count=11, doside=5, sides=5:2:11` | `2.60e-14` | `1.22e-15` | `9.67e-13` |
| `count=13, doside=5, sides=5:2:13` | `9.10e-15` | `1.55e-15` | `5.05e-13` |

### Sweep table

| case | mapping | FSB dim | FBU dim | `E_FSB - E_FBU` | capture deficiency | total wall time |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `count11_doside5` | `identity` | `419` | `419` | `0.00e+00` | `0.00e+00` | `2.232 s` |
| `count11_doside5` | `s/s0=1.0` | `419` | `419` | `7.11e-15` | `0.00e+00` | `1.435 s` |
| `count11_doside5` | `s/s0=0.8` | `419` | `419` | `8.88e-16` | `0.00e+00` | `1.217 s` |
| `count11_doside5` | `s/s0=0.6` | `419` | `419` | `-1.11e-14` | `1.22e-15` | `0.980 s` |
| `count11_doside5` | `s/s0=0.4` | `419` | `419` | `1.87e-14` | `0.00e+00` | `0.817 s` |
| `count11_doside5` | `s/s0=0.2` | `419` | `419` | `-2.60e-14` | `0.00e+00` | `1.009 s` |
| `count13_doside5` | `identity` | `517` | `517` | `-5.11e-15` | `0.00e+00` | `2.182 s` |
| `count13_doside5` | `s/s0=1.0` | `517` | `517` | `9.10e-15` | `1.55e-15` | `2.027 s` |
| `count13_doside5` | `s/s0=0.8` | `517` | `517` | `-3.11e-15` | `1.33e-15` | `1.923 s` |
| `count13_doside5` | `s/s0=0.6` | `517` | `517` | `-5.33e-15` | `8.88e-16` | `1.842 s` |
| `count13_doside5` | `s/s0=0.4` | `517` | `517` | `-8.88e-16` | `1.11e-15` | `1.856 s` |
| `count13_doside5` | `s/s0=0.2` | `517` | `517` | `-6.88e-15` | `0.00e+00` | `1.665 s` |

## Interpretation

Within this bounded same-backend PGDG one-body study, FSB and FBU were already
numerically indistinguishable.

The important points are:

- FSB and FBU had the same final dimension in every tested case.
- The He+ energy gap stayed at roundoff level throughout the identity control
  and the full distorted sweep.
- The FBU ground-state capture deficiency by FSB also stayed at roundoff level.
- The cross-overlap subspace error remained around `1e-13` to `1e-12`, which
  is fully consistent with the same final subspace to numerical precision.

So the present bounded study does **not** show a visible trend in which FSB
approaches FBU only as `s` is reduced. In the tested regime, the current
same-backend PGDG one-body path already gives FSB and FBU that are effectively
the same for these He+ completeness measures.

## Performance note

The recent one-body optimizations made this sweep practical. The total per-point
time in the bounded study was roughly `0.8` to `2.2` seconds, with the larger
costs coming from parent projected one-body setup and the FBU/FSB builds, not
from the final He+ eigensolves themselves.

So this kind of bounded FSB-versus-FBU sweep is now feasible without turning
into a long run.

## Ambiguities and next decision

This leaves one scientific ambiguity:

- if the intended FSB-versus-FBU completeness question is asked on the current
  same-backend PGDG one-body route, the answer appears trivial in the bounded
  cases tested here

That means a more informative next study would need one of:

- a larger or different parent/side regime where FSB and FBU actually separate
- a different target metric, such as Gaussian target leftover
- or a clarified cross-backend/comparison contract, if the real question is not
  same-backend algebraic completeness

This note does **not** address Gaussian target leftover yet.
