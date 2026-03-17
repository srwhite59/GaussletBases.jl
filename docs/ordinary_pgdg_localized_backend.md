# Ordinary Mapped PGDG: Localized Backend

This note records the next ordinary-branch step after the static Cartesian
IDA ingredient layer.

## 1. What the next issue is

The encouraging result from the first ordinary Cartesian IDA pass is that the
separable `Vee` assembly is already in reasonable shape in the mild mapped
regime.

The larger remaining issue is the one-body / overlap side:

- the raw refined PGDG proxy is still not orthonormal
- the raw product-basis overlap is therefore not yet close enough to identity
- the branch is not yet clean enough to serve as a solver-facing basis path

So the next task is not more work on `Vee` assembly.

It is making the experimental PGDG ordinary branch into a more faithful
orthonormal localized basis path.

## 2. What is not the right main answer

A post hoc three-dimensional Lowdin cleanup of the final Cartesian product
basis is not the right main answer.

That would blur the localized diagonal-interaction structure we are trying to
preserve, and it would postpone the real basis-construction issue until after
the product basis is already formed.

## 3. What the right next move is

The cleanup should happen in one dimension first:

1. build the refined analytic primitive proxy
2. clean up overlap / orthogonalize
3. localize with COMX-style position diagonalization
4. fix ordering and signs deterministically
5. only then build the three-dimensional Cartesian product basis

That gives a localized orthonormal one-dimensional backend, and the 3D product
basis inherits that structure rather than being repaired after the fact.

## 4. Backend roles after this pass

The ordinary mapped branch should now be read as:

- `:numerical_reference`
  the trusted validation route
- `:pgdg_experimental`
  the refined pre-COMX analytic proxy path
- `:pgdg_localized_experimental`
  the candidate solver-ready experimental one-body path

The localized backend is still experimental. It is only the next step toward a
clean ordinary solver-facing basis route.

## 5. What the comparison should focus on

For the localized backend, the comparison priorities are:

- one-dimensional overlap error
- three-dimensional product overlap error
- one-body `H1` agreement
- ordinary Cartesian IDA `Vee` agreement

Strong-distortion cases remain stress tests. The main target regime is still
the mild-to-moderate mapped regime where the earlier distortion study already
showed the PGDG route is faithful enough in the White-Lindsey sense.
