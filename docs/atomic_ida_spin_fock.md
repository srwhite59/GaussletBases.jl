# Spin conventions for the atomic IDA Fock layer

This note fixes the density and spin convention before any SCF loop is added.

That issue matters because the algebraic helper

```julia
fock_matrix(ops, density)
```

is useful, but by itself it does **not** define a physically complete HF update
step unless the meaning of `density` is pinned down.

## 1. The convention chosen first

The first physically explicit convention should be UHF-style.

That means the inputs are:

- `density_alpha`
- `density_beta`

both expressed in the same channel-major spatial-orbital basis already used by
the atomic IDA layer.

## 2. The formulas

For the present atomic IDA approximation, the first spin-aware Fock matrices are

```math
F^\alpha = h + J[\rho^\alpha + \rho^\beta] - K[\rho^\alpha]
```

```math
F^\beta = h + J[\rho^\alpha + \rho^\beta] - K[\rho^\beta].
```

So:

- the direct/Hartree term uses the total density
- the exchange term uses only the same-spin density

This is the cleanest first HF convention because RHF then becomes the special
case

```math
\rho^\alpha = \rho^\beta.
```

## 3. How the existing pieces fit into this

The current pieces now have clear meanings:

- `direct_matrix(ops, density)`
  - the direct term for the present atomic IDA/local-diagonal approximation
- `exchange_matrix(ops, density)`
  - the exchange term for the present atomic IDA model
  - it uses full radial-pair density blocks
  - but still with the current two-index radial multipole approximation

Neither one should be described as a fully general four-index Coulomb
contraction.

## 4. What role the old `fock_matrix(ops, density)` still has

The existing one-density helper can still remain, but only with a narrow
meaning:

- a simple algebraic helper
- a spinless-model helper
- or a same-spin model helper

It should **not** silently become the basis of a physical HF update step.

That is why the spin-aware layer should be explicit:

- `fock_matrix_alpha(ops, density_alpha, density_beta)`
- `fock_matrix_beta(ops, density_alpha, density_beta)`

## 5. What this pass is not

This pass is still not:

- a full RHF implementation
- a full UHF workflow
- a complete SCF driver

It only fixes the spin convention so that the next SCF-like step can be built
on a physically explicit foundation.
