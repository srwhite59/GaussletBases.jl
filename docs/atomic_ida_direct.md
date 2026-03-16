# Sector-aware direct Coulomb on top of the atomic IDA layer

This note describes the first small physical consumer of the sectorized atomic
angular backend.

The goal is not yet Hartree-Fock, exchange, or a larger atomic workflow.

The goal is simpler:

**given a spatial one-particle density matrix, build the corresponding direct
Hartree one-body term cleanly on top of `AtomicIDAOperators`**

## 1. What density object this pass takes

The input is a spatial one-particle density matrix in the current channel-major
orbital basis exposed by:

```julia
orbitals(ops)
```

So the density lives on the same orbital ordering already used by the one-body
and IDA layers.

For this first direct/Hartree pass, the important point is that the current
local diagonal approximation only uses the **radial-diagonal** density blocks.

In other words, if the orbital labels are written as `(p, α)` with:

- `p` = radial basis index
- `α` = angular channel index

then this direct term uses

```math
\rho_{(p,\alpha),(p,\beta)}
```

but ignores

```math
\rho_{(p,\alpha),(q,\beta)} \quad \text{for } p \ne q.
```

That matches the local-DA Hartree structure used in the older code.

## 2. What effective output it builds

The output is a one-body direct matrix in the same spatial-orbital basis.

It is therefore something that can later be added to the one-electron atomic
Hamiltonian as an effective Coulomb term.

Within the present approximation, the output is block diagonal in the radial
index:

```math
J_{(p,\alpha),(q,\alpha')} = 0 \qquad \text{for } p \ne q.
```

Each radial block is built by combining:

- the radial multipole tables
- the angular Gaunt data
- the supplied radial-diagonal density blocks

## 3. Why the sectorized angular structure matters here

The old dense route would conceptually write

```math
J_{(p,\alpha),(p,\alpha')}
= \sum_{q,L,\beta,\beta'}
M^{(L)}_{p q}\,
Q_L(\alpha,\alpha',\beta,\beta')\,
\rho_{(q,\beta),(q,\beta')}.
```

That is correct, but it expands immediately into the dense four-index angular
kernel.

The better internal route is to keep the contraction in a more factorized
angular form:

1. group the angular data by fixed multipole projection `M`
2. contract the density against the corresponding `G_{L,-M}` slices
3. rebuild the output block from the matching `G_{L,M}` slices

This uses the conserved angular structure directly, instead of routing the main
computation through dense `angular_kernel(...)`.

## 4. Relation to the dense `angular_kernel` path

The dense `angular_kernel(ops, L)` view is still useful for:

- tiny reference comparisons
- debugging
- examples
- inspection of the small-`lmax` structure

But it should no longer be the main implementation path for consumers.

For the direct/Hartree contraction, the real working path should be:

- sparse/block Gaunt data
- sector-aware angular contraction
- dense `angular_kernel` only as a derived check

## 5. Why this is the right first application

This is the right next step before exchange or larger workflows for four
reasons.

First, it is the first real physical payoff from the new angular preparation.

Second, it uses the same orbital ordering and radial multipole structure already
established by `AtomicIDAOperators`.

Third, it stays much smaller than a whole Hartree-Fock or UHF layer.

Fourth, it gives the cleanest bridge to later work:

- exchange
- self-consistent atomic mean-field steps
- larger interacting applications

## 6. What this pass is not

This pass is not:

- a wholesale `VeeYlm` port
- a Hartree-Fock workflow
- an exchange implementation
- a larger solver layer

It is only the first clean direct/Hartree contraction built on the corrected
atomic IDA architecture.
