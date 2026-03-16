# Internal `ShortGaucoulomb` port

This note records the implementation choice for the Coulomb Gaussian expansion
used by the ordinary Cartesian hydrogen line.

## 1. Decision

The legacy `ShortGaucoulomb.jl` file should be treated as the source-of-truth
implementation for this part of the package.

It is:

- short
- deterministic
- already well documented
- scientifically clear about the mapped-`u` construction

So the right move is not to keep rewriting the same logic in a more "modern"
style. The right move is to bring that construction into the repository almost verbatim
as an internal implementation.

## 2. What is ported nearly unchanged

The internal implementation should preserve, as closely as practical:

- the sinh-mapping construction
- the `doacc` split between the low-cost and high-accuracy presets
- the direct coefficient/exponent generation
- the explanatory comments about what the mapped-`u` grid is doing

This keeps the numerical source clear and makes it obvious that the current
repository is using the same tested construction as the legacy code.

## 3. What is adapted lightly

Only small adaptations are needed for the present repository:

- remove the standalone module wrapper
- fit the implementation into the current package namespace
- wrap the output in `CoulombGaussianExpansion`
- keep the public entry point as `coulomb_gaussian_expansion(...)`

So the legacy construction becomes the internal engine, while the current repo
keeps the simpler package-facing object and naming.

## 4. What does not change publicly

The public story should remain:

- `CoulombGaussianExpansion`
- `coulomb_gaussian_expansion(...)`
- `gaussian_factor_matrix(...)`

In other words, this pass is an internal cleanup and source-of-truth
clarification. It is not a public API redesign.

## 5. Why this is worthwhile

This gives the repository:

- a trusted and well documented Coulomb expansion source
- less risk of silent drift from the legacy implementation
- a cleaner foundation for the separable Cartesian hydrogen path

That is the right balance here: preserve the tested scientific core, while
keeping the current repository structure around it.
