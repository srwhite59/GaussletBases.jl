# PQS Source-Box Operator Framework

This note is the compact current framework for the projected q-shell (PQS)
source-box lane. It is a developer contract for private/internal
implementation work, not a public API, artifact schema, or production-science
claim.

Current Cartesian Hamiltonian producer authority still lives under
`docs/src/developer/designs/cartesian_hamiltonian_producer/`. If this note and
that authority disagree, stop and resolve the documentation conflict before
changing code.

## Current Contract

PQS is source-box first:

```text
filled raw product source CPB
-> one-dimensional COMX/product source transforms
-> boundary product-mode retained rule
-> source-space operator factors / pair plans
-> optional shell realization or support-row adapter
-> final retained blocks and diagnostics
```

The raw product-box stage and shell-realization stage are distinct. Boundary
COMX/product-mode selection is the retained-source authority. Shell projection,
support-row contraction, and Lowdin cleanup are later realization, debug, or
oracle machinery. They must not be used as raw operator authority.

The coordinate-product vocabulary is defined in
[`cartesian_coordinate_product_box_contract.md`](cartesian_coordinate_product_box_contract.md).
A filled source box is a CPB. A shell support such as `B_outer \ B_inner` is
owned support, not a CPB. Filled source boxes, owned shell support rows, and
final retained functions are different objects even when their counts happen
to be related.

## Raw Product Source Stage

A PQS raw product source starts from a filled local source CPB with total
source-mode dimensions. Source dimensions are not "retain plus two" interior
counts; for a `q x q x L` source they are exactly `(q, q, L)` after the chosen
source-span policy is applied.

The retained rule selects boundary product modes from the full product source:

```text
orthogonal 1D COMX/product transforms
-> full q x q x L product modes
-> columns whose local mode index is first or last on at least one axis
```

For rectangular examples, the useful count picture is:

```text
dim PQS(q, L) = q^2 L - (q - 2)^2 (L - 2)
```

For `q = 5`, this is `16L + 18`; cubic `PQS(5, 5)` has retained count `98`.
This is only a compact count example. The construction rule is boundary
product-mode selection from the full source product block, not subtraction of
a previously contracted interior cube.

Raw source-box operator references are built in product-source mode space:

```text
O_boundary = P_boundary' * O_product_box * P_boundary
```

No shell-row projection, Lowdin cleanup, final retained weights, density
normalization, Hamiltonian export, or artifact writing belongs in this raw
source-box reference stage.

## Source Pair Operators

Retained source blocks should be ordinary transformed raw-source blocks:

```text
O_retained = T_left' * O_raw_product(left, right) * T_right
```

This framing should cover direct slabs, rectangular source boxes, PQS shells,
PQS/product mixed blocks, and later supplement/source adapters without adding a
new special Hamiltonian case for each pair family.

Safe one-body source terms are overlap, position, `x2`, and kinetic. For
electron-electron or IDA-style terms, keep the Gaussian expansion/factor index
as the short inner reduction and carry raw pair numerators through
source-transform and realization steps before any reviewed final-basis density
normalization.

PGDG analytic pair-factor terms remain the required production path for raw
pair numerators unless an explicit debug/reference/numerical path is selected.
Do not silently replace the analytic PGDG path with a convenience numerical
quadrature path.

## Shell Realization

Shell realization takes selected source-box boundary modes and realizes them
on owned shell support rows:

```text
selected full source-box boundary modes
-> restrict rows to support.support_indices / support.support_states
-> shell-local Gram on owned support
-> symmetric shell-local Lowdin
-> sign canonicalization / block assembly
```

Lowdin is shell-local. Previous-block projection, recursive projection,
projection against locked prior spans, support growth onto previous regions,
and global final-basis Lowdin repair are forbidden in the normal PQS terminal
route.

Support-row contraction is useful for final basis construction, selected
oracle comparisons, and diagnostics. It is not the source-box operator
algorithm. Cross-block overlap is structurally zero for disjoint owned support;
cross-block kinetic, nuclear, and IDA terms may still be nonzero and are
assembled over terminal block pairs.

## Weights

Keep three weight roles separate:

- raw/source quadrature weights attached to source points or source factors;
- retained-column diagnostic weights produced during provisional transforms;
- final retained-function integrals used by final-basis density/IDA logic.

Final PQS IDA weights are final retained-function integrals. They are not
provisional retained-column diagnostic weights, not squared normalization
constants, and not a license to divide inside raw source-box construction.
Angular GTO supplements and other non-quadrature final functions must not be
treated as positive-weight quadrature carriers.

## Current Module Boundaries

`CartesianRawProductSources` owns source CPBs, source-mode dimensions,
deterministic source-mode ordering, metadata-only axis transform facts, and
source-mode boundary selector facts tied to raw source ordering. It does not
own retained-rule policy, shell projection, Lowdin cleanup, final retained
units, IDA weights, pair blocks, Hamiltonians, exports, or artifacts.

`CartesianPairOperatorPlans` consumes retained-unit pairs plus retained-unit
transform contracts. It must not infer realization paths directly from
retained-unit kinds.

`CartesianPairBlockMaterialization` may materialize safe source-space one-body
blocks and adapter summaries. Metadata-only bridge/readiness summaries are not
route authority and must not grow into Hamiltonian, export, artifact, Coulomb,
IDA/MWG, or public-driver surfaces.

`pqs_terminal_basis_realization` is the active terminal realization surface for
PQS shells. It realizes boundary source modes on owned support with shell-local
Lowdin and construction checks. Older complete-core/shell or source-shell
helpers are donor/oracle/reference material unless the current design docs say
otherwise.

## Common Shell Geometry

PQS and White-Lindsey share first-step shellification geometry: direct core
regions, shell regions, owned support rows, ordering, and coverage are
route-family-free. The family split begins after common shell records exist:

```text
PQS: common shell support + full source CPB -> retained source-box modes
WL:  common shell support -> faces/edges/corners/strata -> 1D contractions
```

Thin slabs are face-like compact slab objects shared by PQS and
White-Lindsey. They are not real shell regions, direct identity sectors, or a
PQS-specific shell projection path.

Mapped-COMX source spans are source-span options at the existing doside/COMX
seam. They may change source-axis transform facts and provenance, but must not
create a second COMX route, Hamiltonian branch, artifact schema, or raw-source
builder.

## Provenance And Artifacts

This framework does not define artifact schemas. Artifact schema authority
lives in the Cartesian Hamiltonian producer design docs.

The durable provenance rule inherited by artifact/manifest work is that basis
identity is a status-bearing construction label and source/fixed-column
provenance fact, not `center_xyz`. Centers are representative metadata only.
Do not infer route identity, basis identity, or retained-rule authority from
center coordinates.

## Guardrails

- Do not use shell-row projection or Lowdin as raw source-box operator
  authority.
- Do not define PQS by subtracting a contracted inner cube or by stitching
  panel/endcap leftovers in the regular shell body.
- Do not preserve fake-PQS/source-backed WL/QW retained transforms as
  independent PQS authority.
- Do not compare gausslet-only PQS values to supplemented WL/QW references as
  same-basis evidence.
- Do not divide by retained-column diagnostic weights inside raw source-box
  construction.
- Do not use source/support metadata summaries as numerical authority.
- Do not promote private H1, H1-J, RHF, supplement preflight, CR2, export, or
  public API readiness without a separate manager/design decision.

## Retired Notes

The older raw-product retained-transform policy was folded into this compact
framework and deleted as duplicated transition narrative.
`projected_q_shell_policy.md` is retained only as a short compatibility pointer
for historical inbound links; the current PQS framework is this file plus the
current Cartesian Hamiltonian producer design authority.
