# Paper-Parity Radial Tail-Extension Milestone

This note records the narrow repo milestone introduced by commit `7c82ff4`.

Suggested framing:

- restore paper-parity radial tail extension to `Rmax`

This is distinct from the earlier prototype milestone `53a16a4`. That earlier
change established the named cached manuscript boundary prototype. This
follow-on milestone restores the missing contract layer that turns that sealed
boundary object into the full old paper-style runtime radial basis.

## Regression before the patch

After the named cached boundary prototype was introduced, rerun-style radial
bases that should have extended to atom-specific `Rmax` could effectively stop
at the sealed boundary prototype alone.

That meant:

- the repo had the correct cached manuscript boundary object
- but the full paper-style runtime basis was not automatically restored
- old extent anchors such as heavy-atom `Rmax = 30.0` were no longer guaranteed
- the sealed `32`-function boundary prototype was therefore not a drop-in
  replacement for the full old paper-style atom basis

## Restored contract

The repo now supports both of the intended paper-parity layers:

- the sealed boundary prototype itself
- the full paper-style radial basis obtained by appending the ordinary
  positive-side `:G10` tail out to atom-specific `Rmax`

This restores the old operational distinction:

- the cached manuscript object is the boundary block
- the full atom-style runtime basis is boundary prototype plus ordinary outer
  tail

## Public interfaces

The restored public entry points are:

- `build_paper_parity_radial_basis(prototype; rmax=..., mapping=..., rmax_count_policy=:legacy_strict_trim)`
- `build_basis(prototype; mapping=..., rmax=..., rmax_count_policy=:legacy_strict_trim)`

The no-`rmax` form:

- `build_basis(prototype; mapping=...)`

still returns only the sealed boundary prototype itself.

## Semantics

The restored full-basis construction is intentionally simple and explicit:

- the cached boundary prototype remains sealed
- the extension step is append-only
- no global relocalization or re-orthogonalization is rerun
- the outer tail is ordinary shifted positive-side `:G10` gausslets
- shorter `Rmax` cases are produced by taking the needed leading prefix under
  the legacy strict-count policy

This means the prototype continues to define the manuscript boundary object,
while the extension layer restores the old paper-style runtime extent contract.

## Restored parity anchors

The restoration was checked against the old paper-style extent anchors.

Direct repo-side anchor checks recovered:

- Ne with `Rmax = 30.0`:
  - `nr = 46`
  - last center `28.34082360930011`
- He with `Rmax = 10.0`:
  - `nr = 30`
  - last center `8.57496365549267`

These are the operational parity checks that matter for rerun drivers using the
paper-style atom basis.

## Why this is distinct from the prototype milestone

- `53a16a4` was the scientific-object milestone for the named cached boundary
  prototype
- `7c82ff4` is the follow-on restoration that makes that prototype usable as
  the front end for the full old paper-style radial extent contract

The prototype milestone settled the boundary object. This milestone restores
how that boundary object is turned into the full atom-style runtime basis for
paper reruns.

## Practical consequence

Rerun drivers should stop using the boundary prototype alone when they want the
full paper-style atomic radial basis.

They should instead use the prototype-plus-tail-to-`Rmax` path:

- `build_paper_parity_radial_basis(...)`
- or `build_basis(prototype; ..., rmax=...)`

## Validation status

Focused coverage for this restored contract lives in `test/runtests.jl`. It
checks:

- Ne-style extension to `Rmax = 30.0`
- He-style shorter extension to `Rmax = 10.0`
- agreement of the restored prototype-plus-tail basis with the corresponding
  direct `RadialBasisSpec(:G10; ...)` construction
- parity of the overload `build_basis(prototype; ..., rmax=...)` with
  `build_paper_parity_radial_basis(...)`

According to the already-reported validation for this commit, the full radial
test group passed:

- `env GAUSSLETBASES_TEST_GROUPS=radial JULIA_DEPOT_PATH=/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/julia_depot julia --project=. test/runtests.jl`

There was also a direct anchor-check command for the Ne/He extent values quoted
above.

This note itself did not rerun those commands in the current turn. So the trust
statement here relies on the already-reported passing full radial test-group
run plus the targeted anchor check, rather than on a fresh rerun tied to
writing this note.
