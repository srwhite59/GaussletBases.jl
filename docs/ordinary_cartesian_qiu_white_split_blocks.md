# Qiu-White Split-Block 1D Construction

This note records the next correction to the separate Qiu-White residual-
Gaussian reference path.

The timing pass changed the diagnosis.

After removing the earlier midpoint-grid cross-block mistake, the main
bottleneck was still not the later dense 3D assembly. It was the supposed
"exact 1D raw-block build," because the mixed raw one-dimensional layer was
still selecting `_NumericalPrimitiveMatrixBackend` for the combined
distorted-gausslet-plus-Gaussian primitive set.

That mixed-layer abstraction was therefore the next thing to remove.

The replacement in this pass follows the legacy split-block pattern more
closely:

- gausslet-gausslet (`gg`) 1D blocks from the existing gausslet side
- Gaussian-Gaussian (`aa`) 1D blocks analytically
- gausslet-Gaussian (`ga`) 1D cross blocks from a dedicated cross-block route
- then assemble the raw-space block matrices from those pieces

The required block families are:

- overlap `gg / ga / aa`
- kinetic `gg / ga / aa`
- position `gg / ga / aa`
- `x^2` `gg / ga / aa`
- Gaussian-factor `gg / ga / aa`
- pair-factor `gg / ga / aa`

This is closer to the legacy Qiu-White-style structure because it avoids the
one bad abstraction point:

- no combined distorted-gausslet-plus-Gaussian primitive layer
- no mixed-layer primitive-backend selection dominating phase 1

It is still a reference path, so dense 3D assembly remains in place for now.
Only after this split-block correction is it meaningful to judge how much of
the remaining cost is intrinsic dense reference work and how much is still an
implementation artifact.
