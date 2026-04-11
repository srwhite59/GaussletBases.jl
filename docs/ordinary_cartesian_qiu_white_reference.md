# Ordinary Cartesian Qiu-White Reference Path

This note records the purpose of the current pass.

The ordinary Cartesian code had accumulated a mixed hybrid story:

- a current COMX/localized hybrid path aimed at the later White-Lindsey-style
  ordinary branch
- plus residual-Gaussian interaction experiments that were still partly built
  in a one-dimensional framework

That was the wrong mixture for the original Qiu-White residual-Gaussian route.

The present pass therefore adds a separate paper-faithful reference path while
leaving the older COMX/localized surrogate path in the tree only for legacy
regression and comparison work.

There is now one important follow-on correction to that first implementation:

- the reference path keeps the one-dimensional separated structure
- but the gausslet-Gaussian cross blocks are no longer built by midpoint-grid
  cross-quadrature
- they now come from the contracted raw-space 1D matrix route used by the
  ordinary code where the primitive operators support it

That correction is recorded in:

- [`ordinary_cartesian_qiu_white_crossblock_correction.md`](ordinary_cartesian_qiu_white_crossblock_correction.md)

The reference algorithm is:

- [`docs/src/algorithms/qiu_white_residual_gaussian_route.md`](src/algorithms/qiu_white_residual_gaussian_route.md)

The important implementation split is now explicit:

- the older `hybrid_mapped_ordinary_basis(...)` route remains only as a
  quarantined legacy/internal surrogate path
- the new Qiu-White reference route builds the full 3D gausslet product basis
  first, then orthogonalizes the added 3D Gaussian orbitals to that full 3D
  space

The scientific goal of this pass is narrow:

- get the Qiu-White route correct in source, docs, examples, and tests
- keep RG terms in the same two-index integral-diagonal approximation (IDA)
  form used for the gausslet channel
- recheck the He `1s^2` scalar in the friendly paper-like backbone regime

This pass should behave materially better than the current surrogate MWG path
if the earlier large RG effects were mostly caused by the incorrect mixture of:

- one-dimensional residual construction
- COMX-cleaned final hybrid bases
- transfer of RG interaction data into a later localized basis

What this pass does not do:

- it does not delete the legacy COMX/localized surrogate path from source yet
- it does not open an ordinary He solver
- it does not broaden the Gaussian adapter into a larger basis subsystem
