# Atomic True-SP Support Milestone

This note records the atomic `SP`-support milestone landed in commit
`fa2b194`.

## Milestone Summary

Exact capability gap, in one sentence:

- before this pass, `lmax = 1` was only loader metadata because the active
  ordinary-QW and nested-QW supplement algebra still assumed one centered
  separable 1D supplement channel cubed across `x/y/z`, which cannot honestly
  represent `p_x/p_y/p_z`

Exact fix, in one sentence:

- the active ordinary-QW and nested-QW routes now expand the shared atomic
  named-basis supplement into an explicit atomic-centered 3D Cartesian shell
  supplement carrying `s`, `p_x`, `p_y`, and `p_z`

Clear active-path statement:

- the active atomic ordinary-QW and nested-QW routes now support a true
  atomic-centered 3D Cartesian `SP` supplement

Clear caveat:

- the separate one-dimensional hybrid builder remains intentionally `s`-only
  and still rejects `l > 0`

## Active Representation

The shared loader object remains:

- `LegacyAtomicGaussianSupplement`
- `legacy_atomic_gaussian_supplement(atom, basis_name; lmax = 1)`

The active QW-side supplement object is now the explicit internal 3D shell
representation:

- `_AtomicCartesianShellSupplement3D`

That object expands the named-basis shell metadata into true atomic-centered
Cartesian orbitals:

- `s`
- `p_x`
- `p_y`
- `p_z`

Both ordinary-QW and nested-QW now consume that same 3D supplement physically
when non-`s` shells are present.

## Focused He Checks

Shared supplement for these checks:

- `legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0 or 1)`

### Ordinary QW

| Case | `E1` | `⟨Vee⟩` |
| --- | ---: | ---: |
| `lmax = 0` | `-1.9653729039662613` | `1.2571992738145714` |
| `lmax = 1` | `-1.9550949132829811` | `1.2520606339789417` |

### Nested Shell-plus-core QW

| Case | Fixed dim | Residual count | `E1` | `⟨Vee⟩` |
| --- | ---: | ---: | ---: | ---: |
| `lmax = 0` | `1405` | `2` | `-1.9945837622850549` | `1.247625954403538` |
| `lmax = 1` | `1410` | `7` | `-1.9933438009534143` | `1.2460620688949136` |

These changes are the important evidence that the `p` channel is now physically
active rather than decorative metadata.

## Validation

Validation run for this milestone:

- `julia --project=. test/runtests.jl`
- `julia --project=docs docs/make.jl`

The docs build emitted only the usual Documenter no-deploy-environment warning.

## Current Trust Status

What is now trusted:

- the atomic ordinary-QW and nested-QW basis-design line is now effectively
  finished through true atomic `SP` support
- both active QW routes now use the same real atomic-centered 3D Cartesian
  supplement for `lmax = 1`

What remains intentionally narrower:

- the separate one-dimensional hybrid builder remains honestly `s`-only

## Roadmap Consequence

This is the point where the main roadmap can move from atomic completion to
diatomic design.

Why:

- the atomic nesting line is already settled
- the active atomic hybrid supplement line now has true `SP` support
- the remaining `s`-only limitation is confined to the separate 1D hybrid
  builder rather than the main active ordinary/nested QW line

So this pass closes the last major atomic basis-design gap that still had to be
resolved before opening the next geometry phase.
