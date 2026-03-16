# A minimal UHF layer on top of the atomic IDA approximation

This note defines the smallest self-consistent mean-field layer that is worth
adding on top of the current atomic IDA ingredients.

The goal is narrow:

- keep the present atomic IDA approximation explicit
- choose a physically clear spin convention
- add the smallest useful UHF fixed-point kernel

This is **not** a broad HF framework. It is only the first self-consistent
consumer of the current one-body, direct, and exchange pieces.

## 1. What the state variables are

The basic UHF state is:

- `density_alpha`
- `density_beta`

Both are spatial one-particle density matrices in the current channel-major
orbital basis used by `AtomicIDAOperators`.

Because the radial basis and the `Y_{lm}` channels are orthonormal in the
intended numerical sense, this is an ordinary orthonormal-orbital mean-field
problem. No generalized overlap machinery is needed in the solver path.

## 2. The Fock matrices

The spin-aware Fock matrices are

```math
F^\alpha = h + J[\rho^\alpha + \rho^\beta] - K[\rho^\alpha]
```

```math
F^\beta = h + J[\rho^\alpha + \rho^\beta] - K[\rho^\beta].
```

So:

- the direct term depends on the total density
- the exchange term depends only on the same-spin density

This is the cleanest first physical convention. Closed-shell RHF-like behavior
is then the special case

```math
\rho^\alpha = \rho^\beta.
```

## 3. Orbital occupations

The first UHF layer keeps occupations completely explicit:

- `nalpha`
- `nbeta`

At each step, the lowest `nalpha` eigenvectors of `F^\alpha` and the lowest
`nbeta` eigenvectors of `F^\beta` are occupied.

That is enough for the present He-like tests and keeps the implementation easy
to inspect.

## 4. The UHF total energy

The total energy is evaluated as

```math
E = \mathrm{tr}[h(\rho^\alpha + \rho^\beta)] +
    \frac{1}{2}\mathrm{tr}[J[\rho^\alpha + \rho^\beta](\rho^\alpha + \rho^\beta)] -
    \frac{1}{2}\mathrm{tr}[K[\rho^\alpha]\rho^\alpha] -
    \frac{1}{2}\mathrm{tr}[K[\rho^\beta]\rho^\beta].
```

This is still the current atomic IDA model:

- `direct_matrix(...)` uses the present atomic IDA/local-diagonal approximation
- `exchange_matrix(...)` uses full radial-pair density blocks
- but the radial interaction data are still the current two-index radial
  multipole tables

So this is not a fully general four-index Coulomb UHF layer.

## 5. The iteration

The first SCF kernel should be a simple fixed-point iteration:

1. start from a one-body guess
2. build `F^\alpha` and `F^\beta`
3. diagonalize them
4. rebuild `\rho^\alpha` and `\rho^\beta`
5. mix with simple damping if needed

That is enough to validate the current mean-field structure.

There is no need yet for:

- DIIS
- a broad SCF framework
- open-ended occupation control

## 6. Why this is the right next step

The repository now has:

- a corrected radial primitive/contraction substrate
- an explicit atomic `(l,m)` layer
- a sectorized angular backend
- direct and exchange mean-field pieces

The next scientifically meaningful question is whether those ingredients
support a small self-consistent He-like UHF calculation.

That is exactly what this minimal layer is meant to answer.
