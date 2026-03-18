> **Status:** supporting note. For the current atomic branch, read
> `docs/current_atomic_branch.md` first. For this mean-field note chain, use
> `docs/atomic_mean_field_supporting_notes.md`.

# Sector-aware exchange on top of the atomic IDA layer

This note describes the matching exchange/Fock-style one-body term for the
present atomic IDA approximation.

The goal is still narrow.

This is **not** a Hartree-Fock workflow, not a full SCF loop, and not a
wholesale port of the older atomic mean-field code.

The goal is only:

**given a spatial one-particle density matrix, build the corresponding exchange
term consistently with the current atomic IDA structure**

## 1. What “exchange” means here

In the present atomic IDA setting, exchange means the effective one-body term
obtained by contracting the same radial multipole tables and angular kernels
against the off-diagonal orbital density structure of one spin sector.

So this is the partner of the new direct/Hartree builder:

- `direct_matrix(...)` gives the Hartree-style term
- `exchange_matrix(...)` gives the Fock-style exchange term

Both live in the same spatial-orbital basis.

## 2. What density object it takes

The input is again a spatial one-particle density matrix in the current
channel-major orbital basis:

```julia
orbitals(ops)
```

The important difference from the direct term is:

- the direct term uses only **radial-diagonal** density blocks
- the exchange term uses the full **radial-pair** density blocks

That is, exchange depends on

```math
\rho_{(p,\alpha),(q,\beta)}
```

with both `p = q` and `p \ne q`.

## 3. How it differs from the direct builder

The direct term has the form

```math
J_{(p,\alpha),(p,\alpha')}
= \sum_{q,L,\beta,\beta'}
M^{(L)}_{p q}\,
Q_L(\alpha,\alpha',\beta,\beta')\,
\rho_{(q,\beta),(q,\beta')}.
```

So it is block diagonal in the radial index on the output side.

The exchange term instead has the form

```math
K_{(p,\alpha),(q,\beta)}
= \sum_{L,\alpha',\beta'}
M^{(L)}_{p q}\,
Q_L(\alpha,\alpha',\beta,\beta')\,
\rho_{(p,\alpha'),(q,\beta')}.
```

So exchange keeps the full radial-pair structure `(p,q)` and therefore
produces a full two-site one-body matrix.

## 4. How the sectorized angular representation is used

The key angular fact is the same as before:

```math
m_\alpha + m_\beta = m_{\alpha'} + m_{\beta'}
```

for any nonzero kernel element.

So the current sectorized angular structure can be used directly:

1. group channel pairs by conserved `m` sum
2. gather the density block for one radial pair `(p,q)` into those sectors
3. apply the per-sector angular kernels
4. scale by the radial multipole value `M^{(L)}_{p q}`
5. scatter the result back into the orbital matrix

That means the real implementation path can stay on:

- sparse/block Gaunt data
- sectorized pair kernels
- dense `angular_kernel(...)` only for tiny reference comparisons

## 5. What approximation is being made

This exchange term is still tied to the present atomic IDA/local-diagonal
approximation.

That means:

- the radial information still comes from the current two-index multipole
  tables in `AtomicIDAOperators`
- this is **not** a fully general four-index Coulomb contraction
- this is **not** yet the full atomic exchange story beyond the present IDA
  approximation

So the result is the exchange term for the current atomic IDA model, not a
claim of complete generality.

## 6. Why this is the right next step

This is the right step after the direct term because it completes the first
small mean-field-style pair:

- direct
- exchange

without requiring a full SCF driver.

Once both terms exist cleanly, the next natural small extension is a tiny
Fock-builder helper on top of:

```math
F = h + J - K,
```

still without committing to a full mean-field workflow.
