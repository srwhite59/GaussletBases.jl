# Ordinary Cartesian hydrogen through Coulomb expansion

This note fixes the starting point for the ordinary-gausslet hydrogen line.

The right first mainline is **not**:

- build a full 3D grid first
- sample `-Z/r` directly on that grid
- treat the grid as the central object

The better first path is:

1. expand `1/r` as a short sum of Gaussians
2. exploit the separability of each Gaussian term
3. build the needed one-dimensional Gaussianized operator factors
4. assemble the 3D hydrogen operator from those separable pieces

## 1. Why this is the cleaner first path

For ordinary Cartesian gausslets, the attractive feature is not a 3D quadrature
grid by itself. It is the combination of:

- a short Gaussian expansion of `1/r`
- tensor-product separability
- one-dimensional operator building blocks

If

```text
1 / r ≈ sum_k c_k exp(-zeta_k r^2)
```

then for a nucleus at the origin,

```text
exp(-zeta_k (x^2 + y^2 + z^2))
= exp(-zeta_k x^2) exp(-zeta_k y^2) exp(-zeta_k z^2).
```

So the 3D nuclear-attraction operator becomes a sum of products of 1D matrices.

That is the natural first ordinary-gausslet hydrogen path.

## 2. Where this comes from in the legacy code

The main legacy sources are:

- `ShortGaucoulomb.jl`
- `Gaucoulomb.jl`
- `Gausslet1dints.jl`

The most important conceptual source is `ShortGaucoulomb.jl`.

That file gives the clean deterministic construction

- fixed uniform sampling in mapped coordinate `u`
- fixed sinh map
- direct coefficients and exponents
- no runtime fitting

`Gaucoulomb.jl` contains the longer historical context and additional fits, but
for the present repository stage it is not the best source to port first.

`Gausslet1dints.jl` then shows how the old code used mapped uniform-in-`u`
integration to build 1D matrices such as overlap, kinetic, and simple
multiplicative operators.

## 3. What the first port should be

The first ordinary Cartesian hydrogen layer should therefore be:

- a deterministic Coulomb Gaussian expansion object
- a one-dimensional Gaussianized operator builder

The key one-dimensional matrix is

```text
<phi_mu | exp(-zeta (x - X)^2) | phi_nu>.
```

Once that exists, the 3D nuclear operator follows from tensor products.

## 4. How this fits the current repository structure

The present repository already has:

- an explicit primitive layer
- contraction from primitives to basis functions
- basis representations

So the first ordinary Cartesian hydrogen path should use the same layering:

1. build the Gaussianized factor matrix on the primitive layer
2. contract it upward to the working ordinary basis
3. assemble the 3D Hamiltonian from the contracted 1D pieces

This keeps the ordinary line aligned with the rest of the package architecture.

## 5. What the mapped-`u` grid is still good for

A mapped-`u` grid is still useful.

It is just **not** the best first mainline for ordinary Cartesian hydrogen.

It remains valuable as:

- a debug path
- a comparison path
- a fallback for more general external potentials that do not have the same
  separable Gaussian form

This is where the logic from files such as:

- `Gausslet1dints.jl`
- `Grid1dNs.jl`
- `Grid1ds.jl`
- `Other/DistortedGaussiansv1.jl`

still matters.

## 6. Bottom line

For the first ordinary-gausslet hydrogen path:

- use Coulomb expansion first
- use separability first
- build 1D Gaussianized pieces first
- treat a mapped-`u` grid as a later debug and general-potential path

That is the cleaner and more historically faithful starting point.
