# Ordinary Cartesian IDA Ingredients

This note records the next ordinary-branch milestone after the mapped
Cartesian hydrogen path.

## 1. What this pass is building

The target here is the static He-style interacting ingredient for the ordinary
Cartesian branch.

That means:

- one mapped full-line basis on each Cartesian axis
- one-body ordinary Cartesian ingredients
- a separable ordinary Cartesian IDA electron-electron matrix

It does **not** mean a He solver yet.

This pass stops at the static object:

- overlap
- one-body Hamiltonian
- two-index IDA interaction matrix
- explicit product-orbital indexing

## 2. Why this is the right next ordinary milestone

The ordinary branch now has a clean one-body hydrogen path built on:

- global coordinate maps on each axis
- Coulomb expansion first
- the experimental PGDG-style analytic backend in the mild/moderate regime
- the numerical mapped route retained as the reference path

The natural follow-on is therefore the static two-electron ingredient for
ordinary Cartesian He-style work, not a solver framework.

## 3. What the implementation route is

This pass keeps the same separable Coulomb-expansion architecture.

For

```text
1 / r12 ≈ Σ_k c_k exp(-ζ_k |r1 - r2|^2),
```

each term separates into Cartesian factors:

```text
exp(-ζ_k (x1 - x2)^2)
exp(-ζ_k (y1 - y2)^2)
exp(-ζ_k (z1 - z2)^2).
```

So the ordinary Cartesian IDA matrix is assembled from one-dimensional pair
factors, not from a brute-force 6D grid.

## 4. Relation to the legacy ordinary branch

This follows the same structural idea as the older ordinary code:

- build one-dimensional pair-kernel factors for each Coulomb Gaussian term
- transform them to the working one-dimensional basis
- divide by the one-dimensional basis weights
- assemble the final three-dimensional IDA matrix from separable products

That keeps the new layer close to the historical ordinary-branch algebra while
still using the current primitive/contraction architecture.

## 5. Backend choice

For the mapped ordinary branch, the backend choice remains explicit:

- `:pgdg_experimental`
- `:numerical_reference`

The experimental PGDG-style analytic backend is the preferred implementation
route in the mild-to-moderate distortion regime.

The numerical mapped route remains the validation and reference path where that
comparison is still practical.

## 6. What this pass is not doing

This is still narrow.

It is not:

- a He solve
- a HF/FCI/DMRG layer
- a brute-force 3D or 6D grid framework
- a broad public PGDG surface

It is only the static ordinary Cartesian He-style ingredient layer that a
later solver could build on.
