# Dictionary For The Current Cartesian Route

This developer note defines current CPB, shellification, lowering, retained,
and terminal vocabulary. It is navigation, not source authority. Numerical
contracts live under the
[Cartesian Hamiltonian producer design index](designs/cartesian_hamiltonian_producer/README.md).

The central distinction is:

```text
shellification chooses owned geometry
lowering chooses a representation recipe
terminal realization builds basis functions
operator assembly builds H1 and IDA in that realized basis
```

The old unit-pair, pair-operator-plan, and pair-block-materialization ladder was
[retired](designs/cartesian_hamiltonian_producer/cartesian_pair_planning_materialization_retirement.md).
Its types and status vocabulary are historical, not the intended architecture.

## Current Route Map

```text
mapped parent lattice
-> common shellification and disjoint owned support
-> route-specific terminal lowering
-> retained-unit and transform contracts
-> PQS or White--Lindsey terminal realization
-> exact terminal one-body and IDA assembly
-> optional residual/supplement composition
```

Mapped-COMX is one source-axis option within PQS realization. Its axis-fact
helper remains at
`src/cartesian/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl`; the
module path does not make the retired pair ladder live.

## Geometry Terms

| Term | Current meaning | Not the same as |
| --- | --- | --- |
| **Parent lattice** | Full mapped Cartesian grid and its one-dimensional axis data. | A retained or final basis. |
| **Parent row/site** | One product-grid function or coefficient row. | A final contracted function. |
| **Coordinate Product Box (CPB)** | `I_x x I_y x I_z`, with singleton intervals allowed. | A general shell difference or arbitrary set of rows. |
| **Owned support** | Disjoint parent rows assigned to one terminal region. | A source CPB, which may overlap another source CPB. |
| **Shellification** | Route-independent assignment of parent rows to cores, atom-local shells, shared shells, slabs, and mismatch pieces. | Contraction, COMX, or operator assembly. |
| **Terminal region** | A geometry leaf requiring no further ownership subdivision. | A terminal basis block; “terminal” here is geometric. |
| **Direct core** | Nucleus-centered identity sector retained without shell contraction. | A PQS source box or a globally orthogonalized core. |
| **Slab** | Thin face-like terminal support used near molecular interfaces. | A complete shell or arbitrary bond-axis panel. |

A shell such as `B_outer \ B_inner` is owned support, not a CPB. It may be
described by CPB relations while remaining a distinct object.

## Lowering And Realization Terms

| Term | Current meaning | Guardrail |
| --- | --- | --- |
| **Lowering** | Selects source CPBs and a route-family recipe for one terminal region. | It need not materialize coefficients. |
| **Lowering contract** | Stable record tying owned support to the chosen source and recipe. | Algorithmic data must not hide in metadata. |
| **Retained unit** | Column-owning planned piece derived from one selected lowering contract. | It is not a unit-pair or global matrix. |
| **Transform contract** | Describes how source modes become retained columns. | It is not itself the coefficient matrix. |
| **Terminal realization** | Builds final support-local coefficients, checks rank/metric, canonicalizes signs, and appends a terminal block. | No previous-block projection or global Lowdin repair. |
| **Terminal basis block** | Realized functions supported on one authoritative owned row set, with a native column range. | Cross-block operators may still be nonzero. |
| **Materialization** | Construction of actual coefficients, operators, or files. | A status record saying “ready” is not materialization. |

Current retained-unit and transform-contract objects remain live route inputs.
The retirement removes only the unused pair inventory and materialization
successors.

## PQS Terms

| Term | Meaning |
| --- | --- |
| **PQS** | Projected q-Shell construction from a filled product source and boundary product-mode selection. |
| **Source shape** | Route-authoritative source dimensions, commonly `(q,q,L)` for a complete molecular shell. |
| **Boundary product mode** | A product mode whose local index is first or last on at least one axis. |
| **Owned-support restriction** | Restriction of candidate source columns to the terminal region's exact support rows. |
| **Shell-local Lowdin** | Symmetric orthogonalization of the retained Gram matrix on those owned rows only. |
| **Mapped-COMX source span** | Optional protected-polynomial plus mapped-Chebyshev axis span entering the existing COMX cleanup. |

For a rectangular `(q,q,L)` product source, the boundary-mode count is:

```text
q^2 L - (q - 2)^2 (L - 2)
```

That count is a consequence of boundary-mode selection, not a license to build
PQS by subtracting a previously contracted interior cube.

## White-Lindsey Terms

| Term | Meaning |
| --- | --- |
| **Boundary stratum** | Face, edge, or corner CPB derived from common complete-shell support. |
| **White-Lindsey lowering** | Route-specific contraction over those strata after common shellification. |
| **Matched aspect source shape** | Shared physical `(q,q,L)` outer shape whose family-specific contractions have equal aggregate shell dimension. |
| **Thin slab contraction** | Compact face-like realization for shared thin slab support. |

White-Lindsey and PQS share parent geometry and eligible shell support. Their
contractions differ only after those common objects exist.

## Operator Terms

| Term | Meaning | Guardrail |
| --- | --- | --- |
| **Raw axis factor** | One-dimensional overlap, kinetic, moment, nuclear, or Coulomb factor. | Preserve the coefficient/exponent convention. |
| **Terminal pair block** | Operator block between two realized terminal blocks. | It is assembled by current terminal operator owners, not the retired pair-plan ladder. |
| **Exact one-body** | Galerkin kinetic plus separated uncharged nuclear matrices, combined with physical charges at Hamiltonian assembly. | `U_A = -1/r_A`; do not bake in `Z_A`. |
| **IDA interaction** | Two-index density-density interaction in the same localized terminal basis. | Do not rotate it as `C' V C`. |
| **Residual augmentation** | Optional Gaussian complement plus exact augmented one-body and MWG interaction blocks. | It is downstream of the terminal basis. |

Disjoint terminal support makes cross-block overlap structurally zero. It does
not make kinetic, nuclear, one-body, or IDA blocks diagonal.

## Status And Provenance Terms

| Term | Use |
| --- | --- |
| **Inventory** | A vector-backed collection of actual route records. Avoid variable-size tuple types. |
| **Plan/contract** | An algorithmic object consumed by a later stage. It should have typed fields and a compact fingerprint. |
| **Summary** | Disposable diagnostic view. It must not become an algorithmic data bus. |
| **Provenance** | Stable construction identity, policy, fingerprints, and source labels. |
| **Oracle** | Independent or bounded reference path, not automatically production authority. |
| **Historical** | Evidence retained for archaeology; it grants no work. |

Avoid resurrecting retired terms such as `UnitPairPlan`,
`PairOperatorPlanInventory`, `PairBlockMaterializationPlan`, source-shell
bridge/readiness statuses, or false materialization field clouds. If a future
physics target genuinely needs a pair abstraction, design it from the current
consumer and authorize it separately.

## Short Reading Rule

For a numerical task, read the assigned machine-authority record and its
canonical contract. Use this dictionary only to decode route vocabulary. Old
plans and the retirement ledger are task-gated history, not implementation
instructions.
