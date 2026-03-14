# Terminology guide

GaussletBases uses a few words that may not be standard in every physics or chemistry group. This page gives short plain-language definitions.

The point is not to force a special vocabulary on you. The point is to make the package easier to read once or twice, and then let the words fade into the background.

## Gausslet

A gausslet is a localized basis function built from a sum of Gaussians.

In this package, gausslets are meant to combine several useful features at once:

* they are smooth
* they can be made orthonormal
* they are local in space
* they can be paired with explicit quadrature grids

For most users, it is enough to think:

**a gausslet is a compact basis function built from Gaussian pieces**

## Basis recipe (`BasisSpec`)

A `BasisSpec` is just a recipe for how to build a basis.

For example, `RadialBasisSpec(...)` is not yet the basis itself. It is the set of instructions used to build one.

This distinction is useful because:

* the recipe is compact and easy to save or reason about
* the actual basis contains the built functions and shared coefficient data

## Mapping

A mapping controls how resolution is distributed in physical space.

For radial work, the mapping lets you place more resolution near the nucleus without forcing the entire outer region to use that same fine spacing.

In practice, you usually meet this through something like:

```julia
map = AsinhMapping(c = s / (2Z), s = s)
```

If you are new to the package, it is fine to think:

**the mapping tells the basis where to spend its resolution**

## Radial basis

A radial basis is the set of basis functions used for the reduced radial function

```math
u(r) = r R(r).
```

In this package, the radial basis is represented by `RadialBasis`, and an individual radial basis function can be accessed as `rb[i]`.

## Quadrature grid

A quadrature grid is the numerical integration grid used to evaluate integrals.

In GaussletBases, the quadrature grid is separate from the basis.

That is important enough to repeat:

**the basis is not the quadrature grid**

You build the basis with `build_basis(...)`, and you build the radial quadrature grid separately with `radial_quadrature(...)`.

## Stencil

`stencil(f)` returns the exact Gaussian expansion used to build the function `f`.

This word is not standard chemistry language, so if it feels unfamiliar, that is normal.

You do not need it for ordinary use. It becomes useful when:

* you want to inspect how a function is built
* you want to work with the common Gaussian layer shared by a whole basis
* you are doing method development

## Primitive

A primitive in this package is one of the underlying Gaussian-type building blocks from which a function or basis is assembled.

For a whole basis, you can inspect the common primitive layer with:

```julia
P = primitives(rb)
C = stencil_matrix(rb)
```

Here:

* `P` is the list of shared building blocks
* `C` is the matrix of coefficients that combines them into each basis function

If you are not doing low-level method work, you can ignore this at first.

## `reference_center`, `center`, and `moment_center`

These names are similar, so it helps to separate them.

### `reference_center(f)`

The center before the coordinate mapping is applied.

### `center(f)`

The physical-space center after mapping.

### `moment_center(f, grid)`

The first-moment center computed numerically on a supplied quadrature grid.

For radial functions near the origin, these do not have to agree exactly.

## `D`

`D` is the package’s aggregate center-mismatch measure used in the radial diagnostics.

Very roughly, it measures how much the nominal basis centers differ from the first-moment centers computed on the quadrature grid.

You do not need to memorize the exact formula to use it well. In practice, it is a diagnostic that helps you notice when the near-origin behavior may need attention.

## `multipole_matrix`

In the current v0 package, `multipole_matrix` means the supported two-index IDA-style radial multipole matrix.

It does **not** mean an exact four-index electron-electron tensor.

That larger exact object is a separate future extension.

## `atomic_operators`

`atomic_operators(rb, grid; Z, lmax)` is a convenience constructor that builds a bundle of radial operator matrices.

It gives you:

* `ops.overlap`
* `ops.kinetic`
* `ops.nuclear`
* `centrifugal(ops, l)`
* `multipole(ops, L)`

This is helpful when you want a compact “give me the radial pieces” interface without rebuilding each matrix yourself.

## Final advice

If you are new to the package, you do not need to master all of these terms at once.

A good order is:

1. learn what a radial basis is
2. learn why the quadrature grid is separate
3. learn the recommended atomic starting setup
4. learn hydrogen
5. only then worry about stencils and primitives if you need them

That is enough to get started productively.
