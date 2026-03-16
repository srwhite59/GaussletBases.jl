# A tiny two-electron consumer of the atomic IDA layer

This note describes the first genuinely interacting calculation built on top of the present atomic IDA ingredients.

The scope is intentionally narrow.

The goal is **not** to build a general many-electron framework. The goal is to answer one cleaner question:

**does the present `AtomicIDAOperators` object support a real two-electron calculation in a way that is explicit, inspectable, and scientifically sensible?**

## 1. What this layer is

The first consumer is a tiny `1 up, 1 down` application layer built on top of:

1. the radial one-body substrate
2. the explicit `(l,m)` channel layer
3. the radial multipole tables
4. the angular kernels built from Gaunt factors

This gives the smallest real interacting problem that still uses the present IDA factorization.

## 2. Why start with `Nup = 1`, `Ndn = 1`

This is the narrowest useful interacting atomic test.

It has three advantages:

- the spatial basis can still be inspected orbital by orbital
- the full two-electron space is small enough for dense exact diagonalization
- the construction exercises the present interacting ingredients directly, without first building a broader fermion or solver framework

So this is a validation layer, not yet a workflow layer.

## 3. What basis is used

The one-electron orbitals are the channel-major spatial orbitals already exposed by `AtomicIDAOperators`.

The two-electron basis is then the product basis

```math
| p \uparrow, q \downarrow \rangle,
```

with:

- the spin-up orbital index running slowest
- the spin-down orbital index running fastest

That makes the indexing explicit and keeps the formulas simple.

## 4. One-body part

For the present atomic IDA line, the intended basis is already orthonormal in the scientific sense:

- the radial gausslet basis is orthonormal to numerical precision
- the `Y_{lm}` channel functions are orthonormal
- therefore the orbital basis is orthonormal to numerical precision
- and the `1 up, 1 down` product basis is orthonormal to numerical precision

So the overlap matrices remain useful diagnostics, but they are not the organizing principle of the solver.

If `h` is the one-electron Hamiltonian matrix in the spatial-orbital basis, then the `1 up, 1 down` product basis carries

```math
H^{(1)} = h \otimes I + I \otimes h.
```

The corresponding overlap diagnostics are still worth checking:

- orbital overlap, which should be close to `I`
- two-electron product overlap, which should be close to `I`

But the actual solve should be treated as an ordinary Hermitian eigenproblem.

## 5. Two-body part

The interaction uses the present IDA factorization in its intended two-index form:

- radial two-index multipole matrices
- angular kernels `Q_L`
- a density-density interaction in the orbital basis

For two spatial orbitals `i = (p,\alpha)` and `j = (q,\beta)`, the interaction
strength is

```math
V_{ij} =
\sum_L M^{(L)}_{p q}\; Q_L(\alpha,\alpha,\beta,\beta).
```

This is a two-index orbital interaction, not an exact four-index electron-electron tensor.

For the `1 up, 1 down` product basis, the two-body part is therefore diagonal:

```math
H^{(2)} | i \uparrow, j \downarrow \rangle
= V_{ij} | i \uparrow, j \downarrow \rangle.
```

That is the intended meaning of the present IDA layer, and it is the one that
matches the older density-density Hartree-Fock-style use of the same Coulomb
object.

## 6. What this object should contain

A good first object should make visible:

- the orbital list
- the two-electron product-state list
- the orbital and two-electron overlap diagnostics
- the one-body contribution
- the two-body contribution
- the full dense Hamiltonian

That is enough to validate the assembly and run a tiny dense solve.

For somewhat larger toy models, a standard Hermitian Lanczos solve is also natural.

## 7. What this layer is not

This is not yet:

- a general FCI layer
- an HF layer
- a DMRG layer
- a performance-oriented many-electron implementation

It is only the first clean interacting consumer of `AtomicIDAOperators`.

## 8. Why this is the right next step

Once this tiny two-electron layer works, the repository will have shown a complete line from:

- radial basis construction
- radial operators
- explicit angular channels
- static interacting IDA ingredients
- an actual interacting calculation

Only after that does it make sense to consider broader many-electron infrastructure.
