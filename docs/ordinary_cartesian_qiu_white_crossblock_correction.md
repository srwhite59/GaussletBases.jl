# Qiu-White Cross-Block Correction

This note records the correction made to the separate Qiu-White residual-
Gaussian reference path.

The first reference implementation was already moving in the right conceptual
direction:

- it built the full 3D gausslet product basis first
- it defined true 3D residual Gaussians against that full space
- it kept residual-Gaussian interaction terms in the same two-index
  integral-diagonal approximation (IDA) form

But one important part of the implementation had still gone too numerical.

The gausslet-Gaussian cross blocks were being built by midpoint-grid
cross-quadrature and dense kernel application, even though the surrounding
ordinary code already had the machinery to build the corresponding raw-space
one-dimensional contracted blocks directly from the primitive layer.

The correction in this pass is therefore:

- keep the separable one-dimensional / Kronecker-style algorithm structure
- but move the raw gausslet-Gaussian cross blocks onto the contracted
  one-dimensional raw-space route
- and assemble the resulting 3D matrices from those 1D blocks directly, using
  explicit factorized loops rather than `kron(...)` hotspots

The intended exact raw-space pieces are now taken from the same contracted 1D
machinery wherever the primitive operators support it:

- overlap
- kinetic
- position and second moments
- Gaussian-factor blocks for the one-body nuclear expansion
- pair-factor blocks for the two-index IDA interaction assembly

This correction matters for two reasons:

- it matches step 6 of the Qiu-White algorithm page more faithfully
- it should also move runtime closer to the legacy style, where the main work
  is in 1D contracted blocks plus structured 3D assembly, not basis-level
  midpoint cross-quadrature

The algorithm source of truth remains:

- [`docs/src/algorithms/qiu_white_residual_gaussian_route.md`](src/algorithms/qiu_white_residual_gaussian_route.md)

The current COMX/localized hybrid path is still left in place as a separate
later route. This pass only corrects the Qiu-White reference implementation.
