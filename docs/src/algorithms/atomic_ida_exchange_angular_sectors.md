# Atomic IDA Exchange Angular-Sector Rule

This page records the current angular-sector contract for
`exchange_matrix(::AtomicIDAOperators, ...)`. It is a narrow algorithm page for
the present atomic IDA approximation. It is not a full four-index Coulomb
algorithm and not a complete Hartree-Fock workflow.

## Pseudocode

1. Work in the channel-major spatial orbital basis returned by:

   ```julia
   orbitals(ops)
   ```

   An orbital is indexed by a radial basis index and an angular channel.

2. Accept a one-spin spatial density matrix `D`.
   Unlike the direct/Hartree builder, exchange must use radial-pair density
   blocks, not only radial-diagonal blocks.

3. For every radial pair `(p, q)`, use the radial multipole value `M[p, q, L]`.
   The exchange output keeps the full radial-pair structure:

   ```math
   K_{(p,\alpha),(q,\beta)}
   =
   \sum_{L,\alpha',\beta'}
   M^{(L)}_{p q}\,
   Q_L(\alpha,\alpha',\beta',\beta)\,
   D_{(p,\alpha'),(q,\beta')}.
   ```

4. Use the exchanged angular channel ordering. The nonzero angular condition
   is:

   ```math
   m_\alpha - m_\beta = m_{\alpha'} - m_{\beta'}.
   ```

   Equivalently, exchange groups angular pairs by `m_left - m_right`. It must
   not reuse the direct-term sector rule in a way that keeps only the `m = 0`
   exchange partner.

5. For each angular sector:
   - gather the density entry `D[(p, alpha_prime), (q, beta_prime)]`
   - read the exchanged kernel element `Q_L(alpha, alpha_prime, beta_prime, beta)`
   - scale by the radial multipole value
   - scatter into `K[(p, alpha), (q, beta)]`

6. Use dense angular kernels only for small oracle comparisons. The intended
   production path is the sectorized angular representation.

## Code Pointers

- Public exchange entry point:
  - `src/atomic_ida_fock.jl:exchange_matrix`
- Sectorized angular kernels:
  - `src/atomic_angular_coulomb.jl`
- Atomic IDA operator payload:
  - `src/atomic_ida.jl:AtomicIDAOperators`
- Current examples:
  - `examples/20_atomic_ida_exchange.jl`

## Validation Contract

The minimal rotational-degeneracy check is:

- direct terms are rotationally degenerate across all `m`
- exchange for a p shell finds all 3 angular partners
- exchange for a d shell finds all 5 angular partners

If exchange finds only one partner for p or d, the implementation is using the
wrong angular sector rule.

## Scope

This algorithm remains inside the present IDA/local-diagonal atomic
interaction model:

- the radial data are two-index multipole tables from `AtomicIDAOperators`
- the output is an effective one-body exchange matrix
- this is not a general Gaussian-style `(pq|rs)` backend
- this is not an SCF policy by itself

The Fock helpers build on this convention as:

```math
F = h + J - K.
```

## References

- Supporting historical note: `docs/atomic_ida_exchange.md`
