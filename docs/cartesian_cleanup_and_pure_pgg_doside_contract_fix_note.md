# Cartesian Cleanup And Pure-PGG `doside` Contract Fix

## Scope

This note records one narrow cleanup and one narrow contract repair:

- the mistaken old one-dimensional COMX-cleaned hybrid ordinary route was
  removed from the supported public `GaussletBases` surface
- the legacy pure-PGG constructor used through
  `PureGaussianGausslet.getNGgaussletonly(...)` was corrected so top-level
  `doside` again controls nesting/shell resolution

The point is not a broad Cartesian redesign. It is to remove a misleading
public route and to make the legacy pure-PGG diagnostics honest again.

## What was removed or quarantined

The old branch centered on `hybrid_mapped_ordinary_basis(...)` was an early
mistake. It is a COMX-cleaned one-dimensional hybrid surrogate, not the
paper-faithful 3D gausslet-plus-Gaussian formulation.

What changed on the repo side:

- `hybrid_mapped_ordinary_basis(...)` and `HybridMappedOrdinaryBasis1D` were
  removed from the exported public surface in `src/GaussletBases.jl`
- the public/reference docs and current example guides no longer present that
  route as a supported workflow
- the legacy examples that still use it were marked as
  legacy/internal experimental regressions rather than current examples
- the old implementation remains in `src/ordinary_hybrid.jl` only as
  quarantined source for surrogate comparison and regression checks

The supported ordinary public story is now:

- mapped ordinary Cartesian backbone work
- plus the separate paper-faithful Qiu-White 3D residual-Gaussian path

not the old one-dimensional COMX-cleaned hybrid route.

## What the actual `doside` bug was

The live `getNGgaussletonly(...)` implementation had drifted so the requested
top-level `doside` no longer controlled local side counts.

The concrete problem was in the current local count selection:

- `finishmain(...)` chose `dosideyz` from the mapped `u`-span, effectively
  `round(uvec[end] - uvec[1] + 1)`
- the main face branches similarly recomputed `dosidex` and `dosidey` from the
  local mapped `u`-span

So density/mapping adaptation replaced the requested nesting knob instead of
refining it. In practice, changing top-level `doside` could leave the basis
unchanged or only weakly coupled to the requested shell resolution.

## Corrected nesting contract

The corrected contract is:

- top-level `doside` is the base requested shell/nesting resolution
- local side counts are derived from that requested `doside`
- rectangular adaptation may scale that requested count by physical aspect
  ratio on a face
- density-based mapping still shapes the local basis functions, but it no
  longer replaces the requested count

Operationally, the fix now does:

- local finish counts use the requested `doside`, capped by the available
  local span
- face counts use the requested `doside`, optionally scaled by in-face
  physical span ratio before capping
- the old density/mapping fit still determines the distortion and basis shape,
  but not the primary shell-count contract

That restores the intended meaning:

- `doside` sets nesting resolution in a scale-independent way
- rectangular adaptation is secondary
- density-based adaptation is subordinate

## Validation

Repo-side cleanup validation:

- public docs/examples/reference pages were updated so the old 1D hybrid route
  is no longer presented as current workflow
- the remaining legacy hybrid examples were explicitly marked as
  legacy/internal

Legacy pure-PGG `doside` validation:

- minimal diagnostic script:
  `tmp/work/pure_pgg_doside_contract_check.jl`
- focused result after the fix on a tiny atomic probe:
  - `doside = 3` gave basis dimension `63`
  - `doside = 5` gave basis dimension `125`

That is the key regression signal: changing top-level `doside` now produces a
meaningfully different nested basis again, instead of being silently replaced
by density/mapping-derived local counts.
