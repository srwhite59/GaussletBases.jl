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
-> future source/final pair-block materialization
-> future assembly
```

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
-> broader pair blocks and assembly later
```

`CartesianPairBlockMaterialization` is now the preflight layer after
`CartesianPairOperatorPlans`. It has also started numerical local pair-block
pilots, but only for direct/direct one-body terms: overlap, position, x2, and
kinetic. It is not yet a PQS block path, White--Lindsey block path,
Coulomb/IDA path, full operator assembly, Hamiltonian bundle, export, or
artifact writer.

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
| **Pair-block materialization**      | The step that preflights and then builds concrete pair blocks from pair-operator plans.                | Current numerical pilots are direct/direct one-body local pair blocks only; not PQS/LW blocks, Coulomb/IDA, Hamiltonian assembly, or artifact export. |
| **Pair operator block**             | Numerical block for one pair of final retained units and one or more operator terms.                  | Not yet built when report says metadata-only.                    |
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
