# Roadmap

GaussletBases is currently strongest as a radial gausslet package: it builds
1D, half-line, and radial bases, constructs explicit quadrature grids, and
forms radial one-body and current two-index IDA-style multipole operators.

This roadmap is intentionally short. It is a guide to likely next work, not a
promise that every item will appear on a fixed schedule.

## Recently completed

- `implemented`: quadrature accuracy profiles
  - `radial_quadrature` now supports `accuracy = :medium`, `:high`, and
    `:veryhigh`
  - `:high` is the normal default
  - `refine` and `quadrature_rmax` remain expert overrides

- `implemented`: broader adaptive quadrature checks
  - the default quadrature path now uses more than overlap alone
  - cheap basis-aware stability checks are part of the stopping logic

## Near-term work

- `planned`: calibrate and tune the new quadrature profiles on representative
  atomic cases
  - validate on `Z = 1`, `2`, and `10`
  - include `0`, `1`, and `2` `XGaussian` cases
  - keep the README workflows directly tested

- `possible`: refine the quadrature profile thresholds and weightings
  - keep the default path conservative without making it unnecessarily heavy
  - improve the balance between overlap checks, center stability, and simple
    operator-weighted checks

- `possible`: revisit offline higher-precision radial basis construction
  - likely as an external or precomputed path
  - not as a normal runtime requirement

## Likely later extensions

- `possible`: precomputed standard radial tables
  - especially recommended `1 XGaussian` and `2 XGaussian` cases
  - keep the general runtime builder as fallback

- `possible`: exact non-diagonal radial electron-electron API
  - separate from the current two-index IDA-style `multipole_matrix`
  - naming and public shape still need design work

- `possible`: more conventional gausslet functionality
  - broader gausslet capabilities beyond the current radial-centered slice
  - the exact scope is intentionally left open for now because this will likely
    be revisited soon

## Longer-term possibilities

- PGDG-related work
- hybrid Gaussian extensions
- interop and export helpers
- downstream workflow integrations

## Open questions

- what convergence contract the public quadrature API should promise
- how much offline or precomputed basis data should ship directly in the package
- which standard `1 XGaussian` and `2 XGaussian` radial tables are stable enough
  to ship directly in the source
- how closely future exact electron-electron APIs should resemble the current
  radial IDA operator surface
