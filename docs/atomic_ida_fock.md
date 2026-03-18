> **Status:** supporting note. For the current atomic branch, read
> `docs/current_atomic_branch.md` first. For this mean-field note chain, use
> `docs/atomic_mean_field_supporting_notes.md`.

# A tiny Fock-style helper on top of the atomic IDA layer

This note describes the next small step after the direct and exchange builders.

The goal is narrow:

**assemble the effective one-body matrix**

```math
F = h + J - K
```

for the current atomic IDA approximation.

This is not a full SCF framework. It is only a small helper that combines
pieces the repository already has.

## 1. What the helper combines

The current atomic IDA layer already provides three pieces:

- the one-body atomic Hamiltonian `h`
- the direct/Hartree term `J`
- the exchange term `K`

So the natural next small helper is simply the Fock-style matrix built from
those pieces.

## 2. What it takes as input

The input remains the same spatial one-particle density matrix used by
`direct_matrix(...)` and `exchange_matrix(...)`.

That density lives in the current channel-major spatial-orbital basis.

## 3. What it does not do

This helper does **not**:

- choose occupations
- diagonalize and update orbitals
- manage mixing
- run a self-consistent loop

It only assembles the effective one-body matrix for a supplied density.

## 4. Why this is useful

This helper is worthwhile because it gives the first clean mean-field-style
combination layer without committing the repository to a larger SCF workflow.

It also makes the structure explicit:

- the one-body term comes from the present radial-plus-angular atomic layer
- the direct term comes from the sector-aware Hartree contraction
- the exchange term comes from the sector-aware exchange contraction

## 5. Approximation level

The approximation level remains exactly the same as in the underlying pieces:

- `direct_matrix(...)` uses the current atomic IDA/local-diagonal approximation
- `exchange_matrix(...)` uses the full radial-pair density blocks, but still
  with the current two-index radial multipole approximation
- `fock_matrix(...)` is therefore still a Fock-style matrix for the present
  atomic IDA model, not a fully general four-index Coulomb mean-field object

## 6. What the smallest next step would be

Once this helper exists, the smallest possible SCF-like step would be:

1. choose a trial density
2. build `F = h + J - K`
3. diagonalize `F`
4. rebuild the density from the occupied orbitals

That is now a natural next possibility.

But unless it falls out almost for free, a full SCF loop is still a separate
step and should remain separate from this small helper pass.
