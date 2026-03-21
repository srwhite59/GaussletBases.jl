# Qiu-White Residual-Gaussian Route

## Pseudocode

1. Build the one-dimensional distorted gausslet bases on the three Cartesian
   axes: `Gx`, `Gy`, `Gz`.
   Code: `src/bases.jl` via `build_basis(...)`

2. Form the full three-dimensional gausslet product basis
   `G = {g_i^x(x) g_j^y(y) g_k^z(z)}`.
   Code: `src/ordinary_cartesian_ida.jl`

3. Form the added three-dimensional Gaussian orbitals `A = {a_I}`, typically
   from a standard basis set.
   Code: `src/legacy_basis_adapter.jl`

4. Define the residual Gaussians by orthogonalizing the added 3D Gaussian
   orbitals to the full 3D gausslet space:
   `a'_I = (1 - P_G) a_I`,
   then orthonormalize the surviving `{a'_I}` among themselves to obtain the
   residual-Gaussian set `R = {r_A}`.
   Code: `src/ordinary_qiu_white_rg.jl`

5. Define the final orthonormal hybrid basis as `B = G ∪ R`.
   Code: `src/ordinary_qiu_white_rg.jl`

6. Construct the one-particle matrices exactly by Galerkin matrix elements in
   the raw gausslet-plus-GTO space:
   `S_raw`, `T_raw`, `Vnuc_raw`, ...
   Then transform them into the final basis `B` using the linear
   transformation from the raw basis to the orthonormal basis `{G, R}`.
   Code: `src/ordinary_qiu_white_rg.jl`

7. Construct the gausslet-gausslet part of the two-electron interaction using
   the integral-diagonal approximation (IDA).
   Code: `src/ordinary_qiu_white_rg.jl`

8. For all interaction terms involving one or more residual Gaussians, use a
   residual-Gaussian approximation.
   Code: `src/ordinary_qiu_white_rg.jl`

   8a. In the nearest-center / GGT version, assign each residual Gaussian to
   the nearest/core gausslet and use the corresponding diagonal interaction
   data.

   8b. In the matched-width-Gaussian (MWG) version, compute the exact first and
   second moments of each residual Gaussian:
   `<x>`, `<x^2>`, `<y>`, `<y^2>`, `<z>`, `<z^2>`.

   8c. From those moments, define an effective separable 3D Gaussian for each
   residual Gaussian, with matched centers and widths along `x`, `y`, `z`.

   8d. Use those effective Gaussian orbitals to evaluate the diagonal
   RG-gausslet and RG-RG two-electron terms, keeping the interaction in the
   same two-index integral-diagonal-approximation (IDA) form used for the
   gausslet channel. This is the form used in the work of Qiu and White,
   although the Qiu-White paper was unfortunately ambiguous on this point.

9. Return the final hybrid Hamiltonian in the basis `B`.
   Code: `src/ordinary_qiu_white_rg.jl`

## References

- Y. Qiu and S. R. White, "Hybrid gausslet/Gaussian basis sets"
- Residual-Gaussian definition and raw-space one-body transformation:
  Sec. III
- Residual-Gaussian GGT and MWG interaction approximations: Sec. IV

## What This Builds

This algorithm builds the paper-faithful Qiu-White hybrid basis and Hamiltonian
route:

- a fixed 3D Cartesian gausslet product basis
- plus orthonormalized 3D residual Gaussians
- with exact one-particle matrices in that final basis
- and a diagonal / two-index interaction representation in which the
  residual-Gaussian terms are approximated by GGT or MWG

The key object is the final basis `B = G ∪ R`, where `R` is defined in the
orthogonal complement of the full 3D gausslet space.

## Where the diagonal approximation enters

The diagonal approximation enters only in the two-electron terms involving
residual Gaussians.

- The gausslet-gausslet two-electron part uses the ordinary gausslet IDA.
- Terms involving residual Gaussians are approximated into the same two-index
  diagonal / IDA form by GGT or MWG.
- The one-particle matrices are not approximated by GGT or MWG.

The very approximate residual-Gaussian treatment is justified by the very low
occupancies of the RGs.

## Current Repo Status

The repo now has a separate paper-faithful reference path for this algorithm.

Current split:

- `src/ordinary_qiu_white_rg.jl` implements the reference path described on
  this page:
  - full 3D gausslet product basis first
  - 3D residual Gaussians orthogonalized to that full space
  - exact contracted 1D raw-space blocks assembled into the 3D one-body and
    moment matrices, then transformed into the final basis
  - RG interaction terms kept in the same two-index IDA form
- the public hybrid path in `src/ordinary_hybrid.jl` is still the later
  COMX/localized hybrid route
- the current `:residual_gaussian_nearest` and `:residual_gaussian_mwg`
  treatments on `HybridMappedOrdinaryBasis1D` remain useful surrogate or
  comparison paths, but they are not the paper-faithful Qiu-White route

This page is still the source-of-truth algorithm for the Qiu-White reference
implementation.

## Implementation Notes

For this algorithm family, code comments should track the pseudocode steps.

Recommended comment style:

```julia
# Alg QW-RG step 4: Define residual Gaussians by orthogonalizing 3D GTOs
# to the full 3D gausslet space.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
```

Guidelines:

- use the same step number as in this page
- keep the comment wording close to the pseudocode wording
- include the docs path exactly
- place the comment at the code block that implements that step
- if the implementation is only partial, say so nearby in the code or docstring
