# R3 Same-Construction Base Reuse

Status: implemented internal orchestration and completed validation contract.

This page is the canonical contract for:

- `HP-R3BASE-FN-01`;
- `HP-R3BASE-TEST-01`;
- `HP-R3BASE-DRV-WIRE-01`;
- `HP-R3BASE-DRV-TEST-01`.

The facility avoids rebuilding exact terminal `G-G` kinetic and by-center unit
nuclear blocks when a supplemented Residual Gaussian construction already owns
the matching base Hamiltonian. It is reuse of trusted matrices, not another
operator algorithm or cache layer.

## Implemented Surfaces

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_base_hamiltonian.jl
bin/cartesian_ham_builder.jl
```

The low-level keyword inputs are:

```text
base_kinetic
base_unit_nuclear
```

They are consumed by
`pqs_terminal_residual_gto_augmented_products(...)` and
`pqs_terminal_residual_gto_augmented_unit_nuclear(...)`. The staged base
producer and `cartesian_residual_gto_mwg_hamiltonian(...)` pass
`base_ham.kinetic` and `base_ham.nuclear_attraction_unit_by_center` from their
own construction. The protected-ladder composer is another internal
same-construction consumer.

Commit `b9ad881df` implemented source reuse. Commit `720912ca4` wired the
canonical driver. Manager-log Passes 109 and 110 preserve the accepted parity,
allocation, and driver evidence.

## Trust Contract

The reused base Hamiltonian, terminal realization, parent bundles, supplement,
and residual augmentation must come from one construction. In the implemented
facade and driver paths, this is guaranteed by the local object graph: the
caller constructs the base once and immediately passes its matrices with the
corresponding terminal objects.

The reuse seam does not persist a same-construction proof. It has no matrix
fingerprint, center metadata, provenance payload, or runtime comparison with a
fresh fallback. Callers must not pass unrelated matrices merely because their
dimensions agree.

Before reuse, source performs these local checks:

- `base_kinetic` has final-basis square dimensions;
- `base_unit_nuclear` has one matrix per atom location;
- every unit-nuclear matrix has final-basis square dimensions;
- atom-location and nuclear-charge counts agree;
- the residual contract matches the terminal basis and atom locations;
- the carried Coulomb expansion matches the PGDG exponent sequence.

The center ordering of `base_unit_nuclear` is the same atom-location ordering
used to construct the base Hamiltonian and augmented operator call.

## Reuse And Fallback

With trusted inputs:

- `base_kinetic` supplies the exact `G-G` block for the augmented kinetic
  transform;
- `base_unit_nuclear[A]` supplies the exact uncharged `G-G` block for center
  `A` before the `G-A` and `A-A` blocks are transformed.

When either keyword is omitted, the corresponding terminal matrix is
recomputed by the implemented exact fallback. The fallback contracts remain in
[terminal G-G products](r3_terminal_gg_product_matrices.md) and
[unit-nuclear Gaussian sums](r3_unit_nuclear_ugg_gaussian_sum.md).

The canonical driver passes the trusted matrices only in its supported
supplemented path. It does not add inputs, hooks, timing labels, stage fields,
or artifact data. Base-only driver behavior is unchanged.

## Validation

The maintained focused source gate is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

Its direct low-level construction exercises fallback; its supported facade
construction exercises reuse. It checks base-block equality, exact augmented
operator finiteness and symmetry, the endpoint, and artifact readback.

Driver validation is accepted evidence from manager-log Pass 110. There is no
dedicated committed driver test for these two keyword handoffs.

## Boundaries

This contract does not authorize:

- public inputs, exports, or API changes;
- persistent cache, provider, provenance, or stage objects;
- report, status, timing, metadata, or artifact fields;
- source-Hamiltonian or interaction transformation;
- terminal product, Gaussian-sum, raw-block, residual, IDA, or MWG changes;
- solver behavior or Cr2 workflow claims.

The optimization remains a same-construction internal call-site decision.
