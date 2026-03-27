# Radial Paper-Parity Prototype Milestone

This note records the repository milestone introduced by commit `53a16a4`
(`Add analytic paper-parity radial prototype`).

This is a separate milestone from the earlier paper-parity x-Gaussian-default
cleanup. The earlier change established that the repo default supplement pair
should match the paper contract. This milestone goes further: the repo now owns
the paper-parity radial construction itself as a first-class named prototype and
vendored cache artifact, with a public loader path and strict construction
semantics.

## What changed

Before this milestone, paper parity on the radial side was expressed mainly as:

- repo defaults chosen to match the manuscript supplement pair
- helper-driven reconstruction of the manuscript radial object
- numerical comparison against legacy and manuscript-style rebuilds

After this milestone, the repo owns a first-class paper-parity radial prototype:

- prototype name: `:paper_parity_g10_k6_x2`
- implementation surface: `src/radial_boundary_prototypes.jl`
- vendored cache artifact: `data/radial/paper_parity_g10_k6_x2.jld2`
- public loader: `radial_boundary_prototype(:paper_parity_g10_k6_x2)`
- runtime realization path: `build_basis(prototype; mapping = ...)`

This means the manuscript radial contract is no longer only a generic runtime
default or an external legacy reference. It is now embodied directly inside
`GaussletBases` as a named cached scientific object.

## Fixed prototype contract

The fixed manuscript prototype contract is:

- prototype name: `:paper_parity_g10_k6_x2`
- family: `:G10`
- reference spacing: `1.0`
- odd seed half-width: `L = 24`
- even-tail parameter: `K = 6`
- paper-parity high-precision x-Gaussian pair:
  - `0.09358986806`
  - `0.02357750369`

The manuscript build/evaluation provenance recorded with the prototype is:

- `h = 0.001`
- `sigma = 3`
- `s0 = 6.5`
- `rmax_int = 80`

## Cache and endpoints

The vendored cache artifact stores two related endpoints.

The canonical scientific endpoint is:

- the final manuscript prototype coefficients in the basis of the seed
  gausslets plus the two paper-parity x-Gaussians

The runtime endpoint is:

- the same final prototype expanded into the underlying half-line Gaussian
  primitives plus the two x-Gaussians

This is better than treating the manuscript contract as only a generic runtime
default because it preserves both:

- a compact authoritative representation of the scientific object
- a directly consumable runtime representation for downstream mapped/operator
  work

## Public interface

The public interface for the named prototype is:

- `radial_boundary_prototype(:paper_parity_g10_k6_x2)`
- `build_basis(prototype; mapping = ...)`

This keeps the manuscript prototype distinct from both:

- the raw family-level front door such as `:G10`
- and later atom-specific mapped basis or operator choices

## Strict semantics and trust boundary

This prototype is built with strict manuscript semantics:

- analytic `S` and `X` construction for the undistorted half-line prototype
- high-precision build path
- explicit expected final dimension
- no silent mode loss

If the expected final mode count is not retained, the build must fail rather
than silently dropping a direction.

This trust boundary is specific to the named manuscript prototype. It does not
automatically imply that every generic runtime radial build in the package uses
the same strict contract.

## Validation

Named-prototype coverage lives in `test/runtests.jl` and includes:

- prototype name and cache presence
- fixed contract fields
- checksum stability
- strict dimension / identity / mode-drop diagnostics
- analytic-versus-numerical paper-parity basis agreement
- mapped operator parity on the manuscript contract

The trust story includes explicit checks on:

- expected final dimension
- retained final dimension
- `mode_drop_count`
- overlap identity error
- basis/center/coefficient checksums

Validation status for this milestone:

- radial test group passed:
  - `env GAUSSLETBASES_TEST_GROUPS=radial JULIA_DEPOT_PATH=/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/julia_depot julia --project=. test/runtests.jl`
- docs build:
  - not rerun specifically for this milestone note

## Practical consequence

This milestone settles the paper-parity radial prototype line inside the repo:

- the manuscript radial object now has a first-class named implementation
- the cache artifact is vendored and loadable directly
- paper-parity reconstruction no longer depends on ad hoc numerical rebuilds

This milestone does not settle every later radial engineering choice. In
particular, it does not by itself fix or define:

- downstream mapped/operator quadrature policy
- later multipole-builder stabilization policy
- paper-driver orchestration or figure-generation workflow

For the paper-parity radial prototype line itself, this is the right stopping
point unless the manuscript contract itself changes.
