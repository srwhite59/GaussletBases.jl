# GaussletBases.jl — Design And Documentation Authority

This file is now a **stable repo-wide design/contracts note**, not an omnibus
"single evolving design note."

The older root-level `DESIGN.md` had become an early 1D/radial API sketch that
no longer matched the current repo surface. The current authority split is now
explicit.

## Documentation authority

At the repository root:

- `README.md`
  - onboarding and first trustworthy repo overview
- `DESIGN.md`
  - repo-wide design/contracts note and documentation-authority map
- `STATUS.md`
  - current capability / trust matrix
- `ROADMAP.md`
  - strategic next-pressure note, not a schedule

In the rendered docs:

- `docs/src/developer/architecture.md`
  - shortest current architecture map
- `docs/src/explanations/current_atomic_branch.md`
  - current user-facing atomic branch status
- `docs/src/explanations/current_ordinary_branch.md`
  - current user-facing ordinary / Cartesian branch status
- `docs/src/algorithms/*.md`
  - path-specific construction recipes and code pointers

Older flat `docs/*.md` files remain useful as supporting-note history, but they
are not the first place to look for the current public repo story.

## Stable repo-wide design contracts

The current package-level design contracts are:

- GaussletBases is a **basis, quadrature, and operator** package first.
- The mature onboarding path is still the **radial / atomic** line.
- The repo now also has a real **ordinary / Cartesian** surface, including:
  - exact Cartesian basis-to-basis `cross_overlap`
  - exact `basis_projector`
  - exact `transfer_orbitals`
  - one-center nested Cartesian fixed-block routes
  - bond-aligned diatomic nested source / fixed-block / diagnostics / geometry
    payload routes
- The package is **not** yet a broad general molecular solver workflow.
- The old one-dimensional COMX-cleaned hybrid route is **legacy/internal** and
  should not be presented as a supported public workflow.
- Experimental chain, square-lattice, and angular lines are real, but they are
  still narrower and less settled than the radial and current ordinary lines.

## How to read the repo today

For a new reader:

1. trust `README.md` for the first repo overview
2. trust `STATUS.md` for what is mature, experimental, or quarantined
3. trust `ROADMAP.md` for what the project is actually pushing on next
4. trust the rendered `current_*_branch.md` pages for branch-specific current
   status
5. trust algorithm pages for path-specific recipes

That ordering is deliberate. It is meant to remove the old ambiguity about
which top-level note was supposed to be current.

## Historical note

The old long-form `DESIGN.md` API sketch is no longer the current authority.
Its useful ideas have either become real code and docstrings or have been
superseded by:

- the rendered architecture page
- the rendered current-branch pages
- the rendered reference and algorithm pages

Historical flat notes remain in `docs/` as supporting material, but they should
not override the authority split above.
