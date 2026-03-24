# Bond-aligned Diatomic H2 Completion State

This note records the current `H2` diatomic completion state on the Cartesian
QW line.

The landed diatomic milestones are:

- `627ae48` `Add first bond-aligned diatomic geometry path`
- `8ab3b6e` `Hook up bond-aligned diatomic QW consumer`
- `65870c4` `Activate bond-aligned diatomic molecular SP supplement`

Together, those commits mean that the repo now has a real bond-aligned `H2`
line through:

- diatomic coordinate distortion
- diatomic box/fixed-block geometry
- downstream ordinary-QW nearest/GGT consumer hookup
- true molecular `s/p` supplement and nonzero residual-Gaussian sector

The source-of-truth policy pages for that line are:

- [docs/src/algorithms/cartesian_nested_diatomic_coordinate_distortion.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/src/algorithms/cartesian_nested_diatomic_coordinate_distortion.md)
- [docs/src/algorithms/cartesian_nested_diatomic_box_policy.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/src/algorithms/cartesian_nested_diatomic_box_policy.md)

## What Is Truly Landed

What is now truly landed on the main diatomic `H2` line:

- a bond-aligned diatomic distortion path in the ordinary QW reference line
- a first bond-aligned diatomic nested fixed-block construction following the
  written split policy
- a downstream ordinary-QW consumer path reading that diatomic fixed block
- a true two-center molecular `s/p` supplement route used physically by both
  the ordinary and nested diatomic QW paths

In plain terms, the diatomic line is no longer only geometry-on-paper. The
ordinary and nested `H2` paths are both real and physically active.

## H2 Validation

### 1. Distortion-first ordinary diatomic path

Bond-axis mapping family:

- `CombinedInvsqrtMapping`

Transverse mapping family:

- the same `CombinedInvsqrtMapping` family, but single-center at the shared
  transverse projection

Validated local spacing for `R = 1.4`:

- left nucleus: `0.5`
- right nucleus: `0.5`

Representative center pattern for `R = 1.4`, `core_spacing = 0.5`,
`xmax_parallel = 6`, `xmax_transverse = 4`:

- midpoint neighborhood: `[-0.5026381793, 0.0, 0.5026381793]`
- near left nucleus `z = -0.7`:
  `[-1.6802083312, -1.0145329517, -0.5026381793, 0.0]`
- near right nucleus `z = 0.7`:
  `[0.0, 0.5026381793, 1.0145329517, 1.6802083312]`

Representative ordinary diatomic checks:

| Case | Basis size | Overlap error | `E1` | `⟨Vee⟩` |
| --- | ---: | ---: | ---: | ---: |
| `R = 1.4` | `9 × 9 × 13 = 1053` | `6.66e-15` | `-1.2797693345` | `0.7772550166` |
| `R = 2.0` | `9 × 9 × 15 = 1215` | `1.49e-14` | `-1.1005089867` | `0.6598455552` |

### 2. First bond-aligned nested fixed-block geometry

Representative `H2`, `R = 1.4`:

- parent box: `ix = 1:9`, `iy = 1:9`, `iz = 1:13`
- after one shared shell: `ix = 2:8`, `iy = 2:8`, `iz = 2:12`
- midpoint split index on `z`: `7`
- child boxes:
  - left: `(2:8, 2:8, 2:7)`
  - right: `(2:8, 2:8, 8:12)`
- child physical widths:
  - left: `(5.3881, 5.3881, 4.2424)`
  - right: `(5.3881, 5.3881, 3.7391)`

So the first `H2` case does split, and the child boxes remain roughly cubic by
the implemented check.

Nested fixed-block summary:

| Quantity | Value |
| --- | ---: |
| fixed dimension | `577` |
| overlap error | `5.995e-15` |
| min transformed weight | `0.355856` |
| max transformed weight | `5.136620` |
| parent one-body ground | `-1.2797693345` |
| projected nested Rayleigh value | `-1.2555770224` |
| parent `⟨Vee⟩` | `0.7772550166` |
| projected nested `⟨Vee⟩` | `0.7770686844` |
| parent ground-orbital capture | `0.9985616107` |

### 3. Molecular supplement / residual completion

The first molecular completion pass activates a real two-center explicit-shell
supplement and a nonzero residual-Gaussian sector for `H2`.

Representative `H2`, `R = 1.4`, `H cc-pVTZ`, `lmax = 1`:

| Path | Gausslet count | Residual count | Overlap error | `E1` | `⟨Vee⟩` |
| --- | ---: | ---: | ---: | ---: | ---: |
| ordinary diatomic QW | `1053` | `10` | `1.04e-12` | `-1.2837120054` | `0.7802860060` |
| nested diatomic fixed-block QW | `577` | `3` | `2.65e-14` | `-1.2805903152` | `0.7779695961` |

That is the critical evidence that the `H2` residual sector is now genuinely
active rather than an empty placeholder.

## Current Local Packaging Follow-up

One small follow-up is currently present in the local tree but is not part of
the landed milestone commits above:

- the vendored curated `BasisSets` file now includes:
  - `H cc-pVTZ`
  - `H cc-pVQZ`
- the `H2` diatomic supplement tests are therefore repo-self-contained in the
  current local tree rather than depending on `~/BasisSets`

This is a packaging/self-containment improvement, not a change to the diatomic
physics/design state summarized above.

## Current Local Visualization-Debug Consequence

The current local visualization/debug artifacts make one policy correction
clear enough to record before another geometry pass:

- loosening the representative `xz`-plane tolerance from `1e-12` to `1e-5`
  does not change the nested selected count at `R = 1.4`
  (`82 / 580` in both cases)
- the raw source-region geometry on the same plane recovers the full parent
  `xz` slice (`117 / 1053`)
- so the visually sparse nested picture is mainly a compressed-center artifact,
  not a slice-tolerance artifact

Those same artifacts also show that the current midpoint split creates an
avoidable homonuclear asymmetry:

- after one shared shell, the working box is `7 × 7 × 11`
- the current split gives child boxes `7 × 7 × 6` and `7 × 7 × 5`
- so the midpoint row is effectively assigned to one child rather than kept
  shared

That is enough evidence to justify a policy correction on the box-policy page:

- for odd-length homonuclear working intervals, reserve the midpoint as a
  shared `nx × ny × 1` slab
- split the remainder into equal left/right child boxes

The related constant-resolution shell-matching issue is also real, but it is
not yet fixed as a hard rule here. That should wait for explicit 1D `doside` /
`COMX` diagnostics.

## Current Local Shared-Shell Resolution Correction

The subsequent local `doside` / `COMX` trace pass and two-point follow-up on
`H2` now support one narrow default correction on the homonuclear shared shell:

- on symmetric tangential shared-shell intervals, if the provisional face
  retain count is even, reduce it by one
- in the current `H2` line this changes the implicated shared-shell face
  retains from `(4, 3)` to `(3, 3)` on the symmetric tangential directions
- keep that rule confined to the shared shell of the current homonuclear
  bond-aligned route

What that correction achieved at both `R = 1.4` and `R = 2.0`:

- the traced shared-shell tangential contractions regained the near-zero center
- the representative nested `xz` projection became visibly more centered
- the fixed dimension dropped slightly
- the main fixed-only and nearest/GGT diagnostics stayed essentially unchanged

What it does not yet claim:

- it is not a general adaptive local-side-count formula
- it is not yet extended to heteronuclear diatomics
- it is not yet extended to chains
- it is not yet applied to edges or child/core local retains

## Validation

Landed milestone validation:

- `julia --project=. test/runtests.jl`

Current local-tree validation after the `H` packaging follow-up:

- `julia --project=. test/runtests.jl`
- `julia --project=docs docs/make.jl`

The docs build emitted only the usual Documenter no-deploy-environment warning.

## What Remains Out Of Scope

What is still deliberately out of scope for this `H2` completion state:

- heteronuclear diatomics
- linear-chain implementation
- arbitrary non-linear geometries
- a general molecule-wide basis-placement subsystem
- anything beyond the first honest two-center molecular `s/p` supplement route

So the next geometry-family work should extend outward from this bond-aligned
`H2` line rather than reopen its settled distortion or box policy.

## Roadmap Consequence

The main diatomic `H2` line is now complete enough that the next missing piece
is no longer basic geometry or consumer hookup.

The next likely implementation priorities are:

- commit the small vendored-`H` packaging follow-up so the `H2` tests are
  repo-self-contained
- then extend the same distinguished-axis policy to:
  - heteronuclear diatomics
  - linear chains

In short, the bond-aligned `H2` line has crossed from “policy and scaffolding”
to “real landed working path.”
