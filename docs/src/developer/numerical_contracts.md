# Numerical Contracts

This page records small internal engineering policies that are easy for
generated code to miss. It is developer-facing rather than part of the public
manual.

## Orthonormal Blocks

When a block is constructed to be orthonormal, the intended contract is:

1. build it to be orthonormal in the relevant metric;
2. check that the resulting overlap is `I +` small unavoidable Float64 noise;
3. then treat the block as orthonormal.

Examples include finalized PGDG/COMX-cleaned blocks, nested fixed blocks, and
other internal blocks explicitly constructed as orthonormal bases.

Do not propagate a near-identity overlap matrix as meaningful working data. In
particular:

- do not store `S = I + eps` by default when `eps` is only Float64 residue;
- do not build downstream logic that repeatedly consults such an `S`;
- do not interpret `1e-12` to `1e-14` nonorthogonality as a structural feature.

Use self-overlaps only for construction, validation, assertions, or one final
cleanup when the deviation is too large. For transfer between final orthonormal
working bases, use the cross overlap only. Do not turn that path into a
generalized-overlap formulation.

For decomposed White--Lindsey plus a GTO supplement, the raw combined Galerkin
generalized solve is diagnostic. The working path projects the supplement into
the residual space, orthonormalizes retained residual directions, and solves an
ordinary Hermitian problem in the resulting final basis.

## Nested Fixed-Block Kinetic

`_NestedFixedBlock3D.kinetic` follows the nested packet contract:

- it is the kinetic matrix carried by the assembled nested packet;
- it is the kinetic payload downstream nested operators use;
- it is not automatically interchangeable with a later contraction of a
  separately assembled ordinary parent kinetic.

For current nested diatomic routes, these can differ measurably even on the same
final basis dimension.

## One-Body Reassembly

`assembled_one_body_hamiltonian(...)` reassembles from the payload actually
stored:

- `kinetic_one_body`;
- `nuclear_one_body_by_center`;
- requested `nuclear_charges`.

Compare reassembled one-body matrices only with operators built under the same
stored kinetic convention. Do not use this helper to claim equality between
routes that disagree about their kinetic payload.

## Current Cartesian Route Boundary

The active PQS/White--Lindsey producer realizes terminal blocks on disjoint
owned support and assembles exact one-body and IDA operators through current
terminal owners. Retained-unit and transform-contract records may feed that
construction, but the former unit-pair index, pair-operator-plan, source-safe
term, bridge/readiness, and pair-block-materialization ladder was
[retired](designs/cartesian_hamiltonian_producer/cartesian_pair_planning_materialization_retirement.md).

Old decomposed-WL pair-count ladders, timing tables, local block placement
results, and He/H2 acceptance numbers were development evidence for that
abandoned architecture. They are available in Git history, not current
numerical contracts.

The continuing rules are:

- shellification owns disjoint support;
- terminal realization owns shell-local rank, Lowdin, and sign conventions;
- cross-block overlap is structurally zero, but cross-block operators need not
  be zero;
- exact one-body and IDA matrices use the realized terminal basis directly;
- source/support summaries are diagnostics, never numerical authority;
- retained density normalization uses final function integrals at the explicit
  IDA boundary;
- uncharged nuclear matrices remain separated by center until physical charges
  are applied;
- no retired pair status, index table, placement plan, or materialization flag
  should be restored without a new source-backed consumer and docs-only
  authority.

The mapped-COMX axis-transform helper remains live at its existing path. Its
presence inside the reduced `CartesianPairBlockMaterialization` module is a
narrow ownership detail, not a surviving pair-materialization contract.
