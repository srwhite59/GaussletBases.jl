# Dictionary for the CPB / shellification / retained-unit route

This developer note is a dictionary for Cartesian route reports, staged helper
objects, and current CPB/shellification terminology. Its purpose is to make
report fields readable. The most important rule is:

```text
shellification chooses owned geometry;
lowering chooses how that geometry will be represented;
construction/materialization builds actual transforms, basis functions, and operator blocks.
```

The [Cartesian coordinate product box contract](cartesian_coordinate_product_box_contract.md)
defines the layer split as:

```text
shellification -> ShellificationRegion / owned shell support
lowering       -> CPBs plus a lowering recipe
construction   -> intermediate and final retained spaces
pair planning  -> pairs of final retained units
```



---

## Core route stages

| Term                               | Meaning                                                                                                                                                                                                                                        | Not the same as                                                                                                                        |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| **Parent lattice**                 | The full Cartesian grid of parent gausslet sites. All later regions are subsets or derived coordinate boxes inside this lattice.                                                                                                               | A retained basis. The parent lattice is the large starting space.                                                                      |
| **Parent site / parent row**       | One grid point or row in the parent lattice.                                                                                                                                                                                                   | A final contracted basis function.                                                                                                     |
| **Shellification**                 | The geometry step that partitions parent sites into **owned support regions**: atom-local shells, shared molecular shells, cores, midpoint slabs, mismatch slabs, etc.                                                                         | Lowering, contraction, operator construction. Shellification should not choose COMX transforms or face/edge/corner contractions.       |
| **Atom-growth**                    | A shellification policy: grow atom-local regions outward from nuclei, then handle contact/midpoint/shared molecular regions.                                                                                                                   | The whole route. It should not name the low-order or PQS method.                                                                       |
| **Terminal shellification region** | A leaf region output by shellification. “Terminal” means no further shellification subdivision is needed before lowering.                                                                                                                      | “Final basis unit.” Terminal means final **for geometry**, not final for the basis.                                                    |
| **Lowering**                       | The step that takes a terminal owned region and chooses the source CPBs and recipe that will later produce retained functions. For LW, lowering chooses boundary-stratum CPBs; for PQS, it chooses a filled source CPB and boundary-mode rule. | Contraction. Lowering may plan a contraction, but it does not necessarily build the coefficient matrix.                                |
| **Construction**                   | The step that actually builds retained spaces, coefficient maps, transforms, shell realizations, etc.                                                                                                                                          | Shellification. Construction uses the shellification/lowering plan.                                                                    |
| **Materialization**                | Actually building concrete numerical objects: coefficient matrices, basis blocks, operator blocks, Hamiltonian matrices, or files.                                                                                                             | Metadata planning. A plan can exist without being materialized.                                                                        |
| **Final retained unit**            | A column-owning retained object used by pair planning and Hamiltonian assembly. It should link back to its owned support, source CPB, lowering recipe, intermediate retained space, and shell realization.                                     | Terminal shellification region. A terminal region becomes a final retained unit only after lowering/construction metadata is attached. |

### Current Cartesian metadata chain

The current staged metadata chain is:

```text
CartesianShellification
-> CartesianTerminalLowering
-> CartesianRetainedUnits
-> CartesianRetainedUnitTransformContracts
-> CartesianUnitPairs
-> CartesianPairOperatorPlans
-> CartesianPairBlockMaterialization preflight
-> direct/direct one-body local pair-block pilot
-> PQS/PQS raw source-space safe-term pilot
-> PQS source-to-shell-realization bridge summaries
-> metadata-only final PQS pair-block readiness summaries
-> LW boundary-stratum adapter preflight
-> local old-kernel-backed LW unit coefficient maps
-> local LW one-body pair-block pilots and summaries
-> private mixed one-body batch consumer and compact summary
-> future full route final pair-block materialization
-> future assembly
```

`CartesianRawProductSources` is a side module used by the PQS source-space
pilot. It records source CPBs, source-mode dimensions, deterministic
source-mode ordering, and metadata-only axis transform facts. It is not a
retained-rule, shell-realization, or pair-block module.

Plainly:

```text
geometry ownership
-> source boxes and lowering recipes
-> final retained-unit records
-> per-unit transform/realization contracts
-> retained-unit pairs
-> pair-operator construction plans
-> pair-block materialization readiness
-> direct/direct one-body pair blocks for the current pilot
-> PQS/PQS raw source-space one-body blocks for the current source pilot
-> metadata-only bridge summaries for later PQS shell realization
-> metadata-only readiness summaries for later final PQS pair blocks
-> local LW boundary-stratum one-body pair blocks for the current adapter pilot
-> private mixed one-body consumer over existing local one-body families
-> broader final pair blocks and assembly later
```

`CartesianPairBlockMaterialization` is now the preflight layer after
`CartesianPairOperatorPlans`. It has also started numerical local pair-block
pilots: direct/direct final local one-body terms, PQS/PQS raw source-space safe
one-body terms, and local White--Lindsey boundary-stratum one-body adapter
blocks for overlap, position, x2, and kinetic. The PQS pilot is not a final
shell-realized PQS pair-block path. Its bridge summaries record source
term/status, source-mode facts, transform/source contract keys, realization
paths, status/blocker counts, and nonmaterialized final flags for single
results or batches. They do not build shell projection, Lowdin, final PQS
retained blocks, Hamiltonian data, exports, artifacts, IDA/MWG data, or
Coulomb blocks. The LW pilot is local boundary-stratum adapter materialization
only; it is not full White--Lindsey route assembly, Coulomb/IDA, a Hamiltonian
bundle, export, or artifact writer. The local PQS and LW selector/summary
helpers are private implementation metadata, not route API.

The private mixed one-body consumer in `CartesianPairBlockMaterialization`
wraps a `PairBlockMaterializationPlan`, one safe one-body term, and
caller-supplied factor/provider facts. It dispatches only to existing local
family selectors: direct/direct local blocks, PQS/PQS raw source-space blocks,
and White--Lindsey boundary-stratum local adapter blocks. It returns a local
batch result plus compact summary. It does not wire into route drivers, build
global operators or Hamiltonians, construct Coulomb, change IDA/MWG semantics,
export artifacts, build PQS shell projection/Lowdin data, or assemble a full
White--Lindsey route.

The private mixed one-body block-set consumer keeps results term-separated.
For routing/reporting, use `_one_body_pair_block_set_view(consumption)` as the
compact status/count view; it is not a matrix or batch-result container. Term
accessors retrieve materialized result records or skipped records for one
requested term. Pair-key accessors retrieve materialized result records or
skipped records for one retained-unit pair key across materialized terms.
`_one_body_pair_block_lookup(consumption, term, pair_key)` is the explicit
exact lookup for one term and one pair key. It may expose one matrix-bearing
`result` field for a materialized hit, or one `skipped_record` field for a
skipped hit. These helpers are retrieval and reporting conveniences only. They
do not sum terms, assemble local/global operators, build Hamiltonians, add
Coulomb, change IDA/MWG semantics, wire route drivers, export artifacts, build
PQS shell projection/Lowdin data, or assemble a full White--Lindsey route.

The local final-readiness helper
`pqs_source_pair_final_block_readiness_summary(bridge_summary)` consumes single
or batch PQS source shell-realization bridge summaries and reports whether a
future final retained PQS pair block could be attempted. Current summaries are
blocked by `:shell_realization_not_materialized`, and blocked bridge summaries
propagate their blockers.

For White--Lindsey boundary-stratum retained-unit pairs,
`CartesianPairBlockMaterialization` still recognizes the pair-operator
`:white_lindsey_boundary_stratum_adapter_path` as
`:white_lindsey_boundary_stratum_adapter_preflight`. That preflight remains the
metadata classification boundary. Behind it, local old-kernel-backed unit
coefficient maps now exist for facet/face, edge, and corner strata, local
pair-level coefficient gathering exists, and local one-body adapter blocks now
exist for overlap, position_x/y/z, x2_x/y/z, and kinetic. These are adapter
pilot blocks only: they are not full route/operator assembly and they build no
Hamiltonians, exports, artifacts, IDA/MWG data, or Coulomb.

The companion internal helper
`white_lindsey_boundary_stratum_adapter_summary(record)` consumes those
preflight records and records old-kernel reuse guidance: facet/face strata
point to `_nested_face_product`, edge strata to `_nested_edge_product`, corner
strata to `_nested_corner_piece`, and facet/edge side helpers to
`_nested_doside_1d`. For batch or plan inputs, the summary reports
`reuse_metadata_available_count` and `reuse_metadata_blocked_count`; these
counts describe reuse metadata availability, not route assembly readiness.

The unit-level descriptor
`white_lindsey_boundary_stratum_unit_adapter_descriptor(unit)` records compact
source-CPB and kernel-input facts for one LW boundary-stratum retained unit:
stratum kind, source CPB role/codimension, active/fixed axis metadata, and the
planned old kernel/helper symbols. The pair-level descriptor
`white_lindsey_boundary_stratum_pair_adapter_descriptor(record[, unit_pair])`
classifies adapter-preflight pairs as facet/facet, facet/edge, facet/corner,
edge/edge, edge/corner, or corner/corner. The old-seed oracle helper
`white_lindsey_materialized_seed_oracle_summary(...)` returns compact
validation counts/ranges/roles and one-body matrix availability metadata only.
It is not route authority and not adapter authority.

The materialized LW adapter helpers are local and explicit:
`white_lindsey_boundary_stratum_unit_coefficients(...)`,
`white_lindsey_boundary_stratum_pair_unit_coefficients(...)`,
`white_lindsey_boundary_stratum_one_body_block(...)`,
`white_lindsey_boundary_stratum_one_body_blocks(...)`, and
`white_lindsey_boundary_stratum_one_body_adapter_summary(...)`. They reuse old
kernels as adapter inputs, not as route authority. The old-seed oracle helper
remains validation-only and not adapter authority.

Current focused oracle validation covers selected representative local
White--Lindsey boundary-stratum pairs, not an exhaustive enumeration. The
adapter fixture maps the x-low `yz` facet to old seed face index 5, retained
range `162:170`; the x-high/y-low `z` edge to old seed edge index 11, retained
range `210:212`; and the adjacent x-high/y-low/z-high corner to old seed
corner index 6, retained range `221:221`. Representative local pairs now
validated against old fixed-block slices are facet/edge, facet/facet,
edge/edge, edge/corner, and corner/corner. The validated one-body terms are
overlap, position_x/y/z, x2_x/y/z, and kinetic. Corner unit coefficients expose
parent support indices when `parent_dims` and fixed coordinates are available;
without `parent_dims`, the support-local metadata path remains explicit and no
parent support index is guessed. This is selected local adapter validation,
not full White--Lindsey route validation, not route-driver wiring, and not
Coulomb, IDA/MWG, Hamiltonian/export, artifact, or production dense-parent
fallback readiness.

The focused oracle-comparison test exposes an aggregate validation coverage
summary for this checkpoint: 5 selected pair families x 8 one-body terms = 40
value comparisons, with zero blocked or metadata-shape-only comparisons in the
current focused coverage. Treat that test as a focused validation gate for LW
oracle coverage, not as a tiny per-pass smoke test.

Routine validation should use the smallest test that validates the edit. For
the private mixed one-body consumer, the tiny synthetic smoke test is the
default per-pass check. The generic record/plan mixed one-body tests are
contract tests for semantic changes or closeout. The White--Lindsey focused
mixed tests use real local adapter fixtures and are boundary tests. The
White--Lindsey oracle comparison remains a gate-only validation, not a casual
smoke test.

### Why “lowering”?

“Lowering” is compiler-style language: it means translating a higher-level geometric object into a lower-level construction representation. Here:

```text
terminal shellification region
-> source CPBs plus lowering recipe
```

For example:

```text
complete shell
-> LW: facet/edge/corner CPBs
```

or:

```text
complete shell
-> PQS: filled source CPB plus boundary COMX-product retained rule
```

It is **not** the same as contraction. “Contraction” should be reserved for an actual linear combination of basis functions, usually represented by a coefficient matrix. Lowering can say which contraction should eventually be built, but lowering itself may remain metadata-only.

A clearer report label might be:

```text
lowering_recipe_status = :planned_not_materialized
```

rather than just:

```text
lowering_status = :planned_not_lowered
```

---

## Geometry terms

| Term                              | Meaning                                                                                                                              | Notes                                                                                                                            |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------- |
| **Coordinate Product Box / CPB**  | An axis-aligned product of coordinate intervals, with singleton intervals allowed: `I_x × I_y × I_z`.                                | A CPB can be a filled box, face/facet, edge, corner, or slab.                                                                    |
| **Filled CPB / Volume CPB**       | A codimension-0 CPB, such as `I_x × I_y × I_z`, with all intervals non-singleton.                                                    | PQS uses filled source CPBs.                                                                                                     |
| **Facet CPB / Face CPB**          | A codimension-1 boundary CPB of a shell, such as `{i_x} × I_y × I_z`.                                                                | “Facet” is more mathematical; “face” matches existing code vocabulary.                                                           |
| **Edge CPB**                      | A codimension-2 boundary CPB, such as `I_x × {i_y} × {i_z}`.                                                                         | Used by low-order LW shell lowering.                                                                                             |
| **Corner CPB**                    | A codimension-3 CPB, such as `{i_x} × {i_y} × {i_z}`.                                                                                | Used by low-order LW shell lowering.                                                                                             |
| **Slab CPB**                      | A codimension-1 CPB that need not be a shell boundary face.                                                                          | Midpoint slabs, contact slabs, and mismatch slabs are slab CPBs.                                                                 |
| **Source CPB**                    | A CPB used as the source domain for constructing retained functions or operator blocks.                                              | In PQS, the source CPB is a filled box even though the owned support is a shell.                                                 |
| **Support CPB**                   | A CPB used as an owned or diagnostic support piece.                                                                                  | Support CPBs describe rows/sites, not necessarily source-mode spaces.                                                            |
| **Shell**                         | Usually `B_outer \ B_inner`: the region between an outer box and an inner exclusion box.                                             | A shell is **not** itself a CPB. The CPB note explicitly says `B_outer \ B_inner` is owned support or shell support, not a CPB.  |
| **Complete shell**                | A shell where the inner box is formed by removing exactly one low and one high boundary layer from every axis of the outer box.      | This is the shell type that decomposes cleanly into 6 facets, 12 edges, and 8 corners.                                           |
| **Owned support**                 | The parent sites assigned to a shellification region.                                                                                | Owned support may be a shell, a CPB, or a union of CPBs.                                                                         |
| **Outer box**                     | The outer CPB that bounds a shell or region.                                                                                         | For a complete shell, `owned_support = outer_box \ inner_box`.                                                                   |
| **Inner exclusion box**           | The inner CPB removed from an outer box to form shell support.                                                                       | Not present for direct slabs/cores.                                                                                              |
| **Core / direct core**            | A direct atom-local central region.                                                                                                  | Usually a filled CPB and often handled by direct or identity-like lowering.                                                      |
| **Midpoint slab / contact slab**  | A direct CPB region between two atom-local regions when a central gap remains.                                                       | Usually lowered as a direct slab, not as a shell.                                                                                |
| **Outer mismatch slab**           | A CPB piece used to fill leftover parent-boundary support when the regular shell growth does not exactly end at the parent boundary. | Should be explicit, not hidden as a fake complete shell.                                                                         |
| **Central distorted product box** | A central gap region large enough to require a product-box-style treatment rather than one-site direct slabs.                        | Not identity, not a shell, and not LW face/edge/corner lowering.                                                                 |

---

## Shellification terms

| Term                                 | Meaning                                                                                                 | Recommended usage                                                                                 |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- |
| **Terminal shellification scaffold** | The metadata object holding terminal shellification regions before unit lowering.                       | Good term. “Scaffold” means planning object, not materialized basis.                              |
| **Terminal region**                  | A shellification leaf: atom-local shell, shared shell, core, slab, mismatch piece, central product box. | Prefer “terminal region” over “terminal unit.”                                                    |
| **Terminal unit inventory**          | A metadata inventory made from terminal regions, describing what retained units will later be built.    | This is okay, but define it carefully: it is not yet a materialized retained basis.               |
| **Aggregate atom box**               | A temporary planning container such as “left atom box” or “right atom box.”                             | Avoid letting this cross the shellification boundary. It should not become a final retained unit. |
| **Nested atom box**                  | Better plain-language name for an atom-local planning container containing its own shells/core.         | Should be opened into terminal regions before lowering.                                           |
| **Atom-local shell**                 | A terminal shell around one atom.                                                                       | LW lowers it into facets/edges/corners; PQS lowers it through a filled source CPB.                |
| **Shared molecular shell**           | A terminal shell surrounding both atoms.                                                                | Same lowering alternatives as atom-local shell, but different role.                               |
| **Central gap**                      | The space between grown atom-local boxes before shared shells begin.                                    | May become midpoint slabs or a central distorted product box.                                     |

Important rule:

```text
Shellification may use aggregate atom boxes internally,
but terminal shellification output should not contain left_atom_box/right_atom_box
as final units.
```

---

## Lowering and construction terms

| Term                                      | Meaning                                                                                                                                                       | Notes                                                                                                     |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| **Lowering recipe / lowering family**     | The chosen method for turning a terminal region into source CPBs and retained-space construction metadata.                                                    | Examples: `:white_lindsey_boundary_strata`, `:pqs_filled_source_cpb`, `:direct_identity_cpb`.             |
| **LW / White–Lindsey low-order lowering** | The low-order route where a complete shell is decomposed into boundary CPBs: facets, edges, corners.                                                          | The complete-shell enumeration is explicitly described in the CPB note.                                   |
| **Boundary-stratum CPB lowering**         | The LW shell lowering: shell → facet CPBs + edge CPBs + corner CPBs.                                                                                          | This is LW-specific. It should not become the universal lowering model.                                   |
| **Facet/edge/corner decomposition**       | The 6/12/8 CPB split of a complete shell.                                                                                                                     | For LW only, or for diagnostics. PQS should not use this as its primary construction.                     |
| **Direct CPB lowering**                   | A direct/core/slab region uses its own CPB as the source/support object.                                                                                      | This is the closest thing to identity-like lowering.                                                      |
| **PQS lowering**                          | A shell region chooses a filled source CPB, forms product-box modes, selects boundary COMX-product modes, then optionally realizes those modes on shell rows. | Not identity. The source geometry is simple, but the retained rule and shell realization are nontrivial.  |
| **Source-mode dimensions**                | The dimensions of the compact source-mode space on a source CPB, such as `(q, q, L)`.                                                                         | Not necessarily the same as physical parent-site dimensions.                                              |
| **Intermediate retained space**           | The retained space after applying a retained rule in source space, before final shell realization.                                                            | PQS boundary product modes live here.                                                                     |
| **Retained rule**                         | A rule saying which source-space modes or columns are kept.                                                                                                   | Examples: identity selector, product/doside transform, PQS boundary-mode selection.                       |
| **Shell realization**                     | The map from intermediate retained space to final shell-supported representation.                                                                             | For PQS, this is shell projection plus Lowdin cleanup.                                                    |
| **Lowdin cleanup**                        | Orthogonalization/cleanup after projecting selected PQS source modes onto shell rows.                                                                         | Belongs to shell realization, not raw product-box operator construction.                                  |
| **Contraction**                           | Actual linear combination of basis functions, usually via coefficient matrices.                                                                               | More specific than lowering. Lowering can plan a contraction; construction/materialization builds it.     |
| **Coefficient map / coefficient matrix**  | Numerical object mapping source/support functions into retained functions.                                                                                    | This is materialized construction data, not shellification.                                               |

---

## PQS-specific terms

| Term                                     | Meaning                                                                                                                    | Notes                                                                                 |
| ---------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| **PQS**                                  | Projected q-Shell. A high-order shell construction based on filled product boxes and boundary COMX-product mode selection. | Replaces stitched endcap/panel/annulus mental models.                                 |
| **q**                                    | The source-mode order/side size used in COMX/source construction.                                                          | For cubic PQS, source dims are `(q, q, q)`; for rectangular shells often `(q, q, L)`. |
| **Filled source CPB**                    | The full coordinate-product box used by PQS before shell projection.                                                       | Not the shell itself.                                                                 |
| **Boundary COMX-product mode selection** | Retained rule selecting product modes whose local mode index is on a boundary in at least one axis.                        | This gives the PQS count such as `q^2 L - (q-2)^2 (L-2)`.                             |
| **Raw product-box operator**             | Operator built in the source product-box mode space using 1D factors.                                                      | Should be built before shell projection/Lowdin.                                       |
| **Shell-row projection**                 | Restricting or projecting selected source modes onto the owned shell support rows.                                         | Part of shell realization.                                                            |
| **PQS realization**                      | Projection of selected product-box modes to shell rows plus Lowdin cleanup.                                                | Not the raw operator rule.                                                            |

Plain answer to your question:

```text
For PQS, the source CPB is simple, but the lowering is not identity.
```

The PQS source geometry is “one filled box,” but the actual retained space comes from boundary COMX-product selection and later shell realization.

---

## Pair/operator terms

| Term                                | Meaning                                                                                               | Notes                                                            |
| ----------------------------------- | ----------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| **Unit inventory**                  | A list of retained-unit records. At metadata stage, it describes units planned from terminal regions. | It does not necessarily contain numerical basis functions yet.   |
| **Final retained unit**             | A unit ready to be used for pair planning.                                                            | It should not be an aggregate atom box.                          |
| **Retained-unit transform contract** | Metadata saying how one retained unit will later get from source rows or modes to final retained columns. It records transform and realization paths without building matrices. | It is not a coefficient matrix, COMX transform, Lowdin matrix, or operator block. |
| **Unit pair**                       | A pair of final retained units whose operator block must eventually be built.                         | For symmetric matrices, only upper-triangular pairs are counted. |
| **Upper-triangular pair inventory** | All pairs `(i, j)` with `i ≤ j`. If there are `N` units, count is `N(N+1)/2`.                         | Pure bookkeeping.                                                |
| **Pair operator plan**              | Metadata describing how an operator block between two units should be built.                          | It reads transform and realization paths from retained-unit transform contracts, not directly from retained-unit kinds. It may be ready, blocked, adapter-only, or not materialized. |
| **Source operator block**           | Operator block built between source CPBs/intermediate retained spaces.                                | For PQS, this is the natural first numerical block.              |
| **Final pair block**                | Operator block between final retained units after any realization/transform maps.                     | This is what assembly eventually places into the global matrix.  |
| **Pair-block materialization**      | The step that preflights and then builds concrete pair blocks from pair-operator plans.                | Current numerical pilots include direct/direct local one-body blocks, PQS/PQS raw source-space one-body blocks, and local LW boundary-stratum one-body adapter blocks; not Coulomb/IDA, Hamiltonian assembly, or artifact export. |
| **Mixed one-body pair-block consumer** | Private `CartesianPairBlockMaterialization` helper over a `PairBlockMaterializationPlan`, safe one-body term or term set, and caller-supplied factors/provider facts. | Dispatches to existing direct/direct, PQS/PQS source-space, and LW boundary-stratum local one-body selectors and returns term-separated local results plus compact views/accessors; not route-driver wiring, local/global assembly, Coulomb, IDA/MWG, PQS shell projection/Lowdin, Hamiltonian/export, artifact writing, or full LW route assembly. |
| **Mixed one-body block-set view/accessors** | Private status/count view and explicit retrieval helpers for a mixed one-body block-set consumption result. | `_one_body_pair_block_set_view` is compact and matrix-free; term and pair-key accessors can explicitly retrieve result or skip records; `_one_body_pair_block_lookup` is the exact term/pair accessor that may expose one matrix-bearing `result`. None of these helpers perform assembly or route adoption. |
| **Pair operator block**             | Numerical block for one pair of final retained units and one or more operator terms.                  | Not yet built when report says metadata-only.                    |
| **PQS source-space block**          | A raw source-mode block for a PQS/PQS safe one-body term, built from caller-supplied 1D source factors. | Not a final shell-realized PQS pair block.                       |
| **PQS source shell-realization bridge summary** | Metadata-only summary describing how a PQS source-space block or batch can later be consumed by shell realization. | It records keys, source-mode facts, statuses, blockers, paths, and flags; it builds no shell projection, Lowdin, final pair block, Hamiltonian, export, artifact, IDA/MWG data, or Coulomb. |
| **PQS final pair-block readiness summary** | Metadata-only summary over a single or batch PQS source shell-realization bridge summary. | Reports whether a future final retained PQS pair block could be attempted; currently blocks on `:shell_realization_not_materialized` and builds no shell projection, Lowdin, final block, Hamiltonian, export, artifact, IDA/MWG data, or Coulomb. |
| **LW boundary-stratum adapter preflight** | Pair-block materialization classification for `:white_lindsey_boundary_stratum_adapter_path`. | Uses `:white_lindsey_boundary_stratum_adapter_preflight` as the adapter boundary. Local unit coefficient maps and one-body adapter blocks now exist behind this boundary, but full route/operator assembly, Coulomb, IDA/MWG, Hamiltonian export, and artifacts remain unavailable. |
| **LW adapter reuse summary** | Internal metadata helper `white_lindsey_boundary_stratum_adapter_summary(record)` over LW adapter preflight records. | Records reuse targets: facet/face -> `_nested_face_product`, edge -> `_nested_edge_product`, corner -> `_nested_corner_piece`, facet/edge side helper -> `_nested_doside_1d`; batch summaries expose reuse metadata counts only, not full route assembly readiness. |
| **LW unit adapter descriptor** | Internal metadata helper `white_lindsey_boundary_stratum_unit_adapter_descriptor(unit)` over one LW boundary-stratum retained unit. | Records compact source-CPB/kernel-input facts before coefficient materialization. |
| **LW pair adapter descriptor** | Internal metadata helper `white_lindsey_boundary_stratum_pair_adapter_descriptor(record[, unit_pair])` over one LW adapter-preflight pair. | Classifies facet/facet, facet/edge, facet/corner, edge/edge, edge/corner, and corner/corner pairs; it does not build pair blocks or call old kernels. |
| **LW boundary-stratum unit coefficients** | Local adapter helper `white_lindsey_boundary_stratum_unit_coefficients(...)` for facet/face, edge, and corner strata. | Builds old-kernel-backed unit coefficient maps as local adapter inputs only; not route-global state, Hamiltonian data, export, artifact, IDA/MWG data, or Coulomb. |
| **LW one-body adapter block** | Local adapter helper `white_lindsey_boundary_stratum_one_body_block(...)` and batch helper `white_lindsey_boundary_stratum_one_body_blocks(...)` for overlap, position_x/y/z, x2_x/y/z, and kinetic. | Produces local/source/final pair-block pilot results from prepared unit coefficients; not full route/operator assembly, Coulomb, IDA/MWG, Hamiltonian assembly, export, artifact, or production dense-parent fallback. |
| **LW one-body adapter summary** | Compact readiness helper `white_lindsey_boundary_stratum_one_body_adapter_summary(...)`. | Reports supported terms/strata, local-only coefficient-map scope, materialized/skipped batch counts, and false global-path flags. |
| **LW seed oracle summary** | Internal validation helper `white_lindsey_materialized_seed_oracle_summary(...)`. | Summarizes old seed counts, ranges, roles, packet/operator availability, and one-body matrix dimensions for validation only; it is not route or adapter authority. |
| **PQS source safe-term descriptor** | Private local metadata in `CartesianPairBlockMaterialization` for supported PQS source safe one-body terms. | Selector/code-organization helper only, not public API or route behavior. |
| **Assembly**                        | Placing pair blocks into full retained operator/Hamiltonian matrices.                                 | Comes after pair-block construction.                             |
| **Hamiltonian matrix / Ham bundle** | Final or export-ready operator/Hamiltonian data.                                                      | Much later than shellification and lowering.                     |

Important rule:

```text
PairOperatorPlans reads transform and realization paths from
RetainedUnitTransformContracts, not directly from retained-unit kinds.
```

Retained-unit kinds can still help classify pair families and source-operator
paths. They should not be the authority for how a retained unit is realized in
the final pair block.

---

## Metadata/status terms

| Term             | Meaning                                                                                          | Suggested clearer phrase                             |
| ---------------- | ------------------------------------------------------------------------------------------------ | ---------------------------------------------------- |
| **Scaffold**     | Planning object that carries geometry/metadata but not numerical matrices.                       | “metadata scaffold”                                  |
| **Sidecar**      | Additional typed object attached to an existing staged named tuple without replacing old fields. | “typed sidecar”                                      |
| **Inventory**    | A list/table of planned objects: regions, units, pairs, etc.                                     | “metadata inventory” if not materialized.            |
| **Readiness**    | Whether enough metadata exists to safely start the next step.                                    | “precondition status”                                |
| **Preflight**    | Checks before building something expensive or numerical.                                         | “pre-materialization check”                          |
| **Blocked**      | A required input/design condition is missing.                                                    | Not an error by itself.                              |
| **Materialized** | Numerical object has actually been built.                                                        | Strong word: do not use for metadata-only records.   |
| **Adapter**      | Temporary bridge from old/current implementation to new contract.                                | Should not be mistaken for the final algorithm.      |
| **Oracle**       | Reference/checking path used to verify another method.                                           | Not production algorithm unless explicitly promoted. |
| **Fallback**     | Less-preferred path used when primary path is unavailable.                                       | Should be labeled clearly.                           |

---

## Terms to avoid or rename

| Current / confusing term                        | Problem                                                         | Better wording                                                                    |
| ----------------------------------------------- | --------------------------------------------------------------- | --------------------------------------------------------------------------------- |
| **Aggregate subtree**                           | Too abstract; unclear to humans.                                | `nested atom-local box` or `compound atom-local planning container`               |
| **Terminal unit**                               | Sounds like final basis unit.                                   | `terminal shellification region` or `terminal region inventory`                   |
| **Lowered**                                     | Can sound like contracted/materialized.                         | `lowering recipe chosen`, `lowering metadata planned`, or `lowering materialized` |
| **Final unit** when referring to geometry       | “Final” should be reserved for retained basis/assembly objects. | Use `terminal region` for geometry.                                               |
| **Source box** without saying CPB or mode space | Ambiguous between physical grid support and source-mode space.  | `source CPB` for geometry; `source-mode space` for modes.                         |
| **Materialization ready**                       | Can sound like materialized.                                    | `ready to begin materialization`, `not materialized yet`                          |

---

## Recommended report wording

Instead of:

```text
lowering_status = :planned_not_lowered
```

use:

```text
lowering_recipe_status = :planned_not_materialized
```

Instead of:

```text
terminal_shellification_unit_inventory
```

use:

```text
terminal_region_unit_inventory
```

or:

```text
unit_inventory_from_terminal_regions
```

Instead of:

```text
aggregate_subtree_operator_plan_required
```

use:

```text
atom_local_box_not_terminal_region
```

or, after the new shellification fix:

```text
blocked_pending_terminal_atom_shell_units
```

Instead of:

```text
PQS source CPB identity
```

say:

```text
PQS filled source CPB selected; boundary-mode retained rule not materialized
```

---

## Short conceptual map

```text
Parent lattice
    full grid

Shellification
    partitions parent sites into terminal owned regions

Terminal shellification region
    geometry leaf: shell, core, slab, mismatch piece

CPB
    coordinate-product source/support object

Lowering
    chooses CPBs and recipe for a terminal region

Construction
    builds retained spaces, transforms, shell realizations

Final retained unit
    column-owning unit used for pair planning

Pair planning
    pairs final retained units

Pair operator block
    numerical operator block for one pair

Assembly
    global retained operator/Hamiltonian matrices
```

The two most important definitions are:

```text
terminal = final for shellification geometry only
final = final retained basis/unit level
```

and:

```text
lowering = choosing/organizing the construction recipe
contraction = actual linear-combination coefficient construction
```
