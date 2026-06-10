# Cartesian Route Retirement Ledger

This developer ledger tracks old Cartesian, White--Lindsey, and PQS route
surfaces while the implementation moves toward the module spine:

```text
CartesianCPB
-> CartesianShellification
-> CartesianTerminalLowering
-> CartesianRetainedUnits
-> CartesianRetainedUnitTransformContracts
-> CartesianUnitPairs
-> CartesianPairOperatorPlans
-> CartesianPairBlockMaterialization
```

Use this file to answer: "Which old authority surface is this pass moving
toward adapter, oracle, or retired?" The ledger is not a public API contract.
It is a cleanup map for preventing old route glue, report aliases, and oracle
paths from becoming new route authority.

| Old surface / file | Current role | Replacement module/stage | Equivalence test | Status | Deletion condition | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| `src/pqs_source_box_route_driver_helpers.jl` | Legacy route-driver glue and report compatibility aliases. | Module-owned route state/summaries and future route adapter over the module spine. | Focused driver/report compatibility checks comparing compact structured summaries to the still-required report aliases. | compatibility | Final reports and downstream consumers no longer depend on flat scalar aliases; route code consumes module-owned objects directly. | No new route concepts should be introduced here. |
| `src/cartesian_pair_block_materialization/white_lindsey_*` and old White--Lindsey nested kernels such as `_nested_doside_1d`, `_nested_face_product`, `_nested_edge_product`, `_nested_corner_piece` | Old-kernel reuse behind the White--Lindsey boundary-stratum adapter. | `CartesianPairBlockMaterialization` White--Lindsey boundary-stratum adapter. | Focused local adapter equivalence tests against selected old-kernel/fixed-block slices for retained-unit coefficients and low-order one-body blocks. | adapter/oracle | New spine can produce and validate all low-order one-body blocks and global overlap/one-body placement without treating old kernels as route authority. | Old kernels may feed adapters or oracle comparisons; they should not decide route structure. |
| Old White--Lindsey materialized seed / fixed-block safe one-body matrices | Oracle/reference for selected safe one-body global matrix slices. | `CartesianPairBlockMaterialization` route-global safe one-body path: route-local one-body block collection -> term-specific placement plan -> term-specific dense global retained matrix. | `test/nested/cartesian_pair_block_global_overlap_oracle_runtests.jl` covers selected old-oracle overlap slices through `White-Lindsey local adapter overlap blocks -> local one-body block collection -> overlap placement plan -> global overlap matrix`: face/face `162:170, 162:170`, face/edge `162:170, 210:212`, and transpose-filled edge/face `210:212, 162:170`. Observed max absolute errors were `1.7636953211319766e-15`, `2.6343504968563524e-46`, and `2.6343504968563524e-46`, respectively. `test/nested/cartesian_pair_block_global_kinetic_oracle_runtests.jl` covers selected White--Lindsey global kinetic equivalence through `White-Lindsey local adapter kinetic blocks -> local one-body block collection -> kinetic placement plan -> global kinetic matrix` on the same slices. Observed kinetic max absolute errors were `3.4638958368304884e-14`, `3.1907967999947086e-30`, and `3.1907967999947086e-30`, respectively. `test/nested/cartesian_pair_block_global_position_oracle_runtests.jl` covers selected global `position_x`, `position_y`, and `position_z` equivalence through `White-Lindsey local adapter position blocks -> local one-body block collection -> position placement plan -> global position matrix` on the same slices. Observed position_x max absolute errors were `1.7435358081325719e-15`, `9.857551802004026e-46`, and `9.857551802004026e-46`; position_y errors were `8.459873058838228e-16`, `8.758115402030107e-46`, and `8.758115402030107e-46`; position_z errors were `8.474501317837309e-16`, `4.926439913641935e-47`, and `4.926439913641935e-47`. `test/nested/cartesian_pair_block_global_x2_oracle_runtests.jl` covers selected global `x2_x`, `x2_y`, and `x2_z` equivalence through `White-Lindsey local adapter x2 blocks -> local one-body block collection -> x2 placement plan -> global x2 matrix` on the same slices. Observed x2_x max absolute errors were `1.844473375912524e-15`, `6.0243961376664915e-34`, and `6.0243961376664915e-34`; x2_y errors were `4.541938536395539e-16`, `1.7964170116877891e-32`, and `1.7964170116877891e-32`; x2_z errors were `4.577863961878137e-16`, `2.7369110631344083e-47`, and `2.7369110631344083e-47`. Synthetic global retained matrix pilots now exist for `:overlap`, `:kinetic`, `:position_x`, `:position_y`, `:position_z`, `:x2_x`, `:x2_y`, and `:x2_z`. | adapter/oracle-only; partially replaced for selected face/face and face/edge safe one-body slices | Broader retained-unit inventory equivalence, broader local/global block equivalence, and route-driver adoption must land before deleting old authority. | This is selected-slice overlap, kinetic, position, and x2 equivalence plus synthetic safe one-body global assembly coverage, not full-route replacement. Keep the old path as validation authority. No full White--Lindsey route assembly, driver integration, Hamiltonian assembly, term summing, Coulomb, IDA/MWG, PQS projection/Lowdin, exports, or artifacts are claimed here. |
| Raw PQS support-row contraction / shell-row fallback paths | Debug/oracle only. | Raw product source-space blocks plus shell-realization bridge. | PQS source-space safe-term blocks plus explicit projection/Lowdin realization match the support-row or shell-row oracle cases. | oracle | PQS source-space blocks plus explicit projection/Lowdin realization match oracle cases. | Must not become the production PQS algorithm. Shell/support rows are compatibility and debug surfaces. |
| Flat `terminal_shellification_*` / `route_core_typed_pair_operator_*` report aliases | Final report compatibility. | Compact summaries from `terminal_route_state`, `pair_operator_summary`, and module summaries. | Report compatibility checks comparing alias values to compact module summaries, without copying large staged metadata objects. | compatibility | Reports, tests, and users consume structured summaries rather than flat aliases. | Do not add new flat driver fields. Derive required aliases at the report boundary only. |
| Local one-body block collection | New local result organizer, not old code. | Future global placement/assembly layer. | Focused tests that local entries, skipped records, summaries, and materialization flags match the mixed one-body block-set consumer. | active | Not a deletion target; becomes typed if it crosses a stable module boundary. | It organizes local pair-block results and skips only. It does not place global matrices, sum terms, build Hamiltonians, export artifacts, or change IDA/MWG semantics. |

## Rules For Future Passes

- New work should move one ledger row toward adapter, oracle, or retired.
- Do not add new flat driver fields.
- Do not add metadata-only readiness layers unless they replace a duplicated
  surface, feed overlap/global placement, or retire an old route surface.
- Old code may be used as adapter or oracle, not as new route authority.
- Deletion requires an explicit equivalence test and a stated deletion
  condition.

## Development Steering Corrections

The Cartesian route cleanup has exposed two repo-wide steering corrections.

First, new route concepts should be organized as module-owned objects with
clear dependency direction, not as flat helper and metadata sprawl. When a
concept crosses a stage boundary, it should usually have a compact record,
summary, and owning module.

Second, long-term tests should preferentially be scientific or workflow
acceptance gates plus compact module-contract tests. Development-era helper
tests, exhaustive metadata checks, and route-shadow scaffolding should be
deleted or quarantined once they no longer protect the intended algorithm.

This policy also counters an observed default agent bias: agents tend to add
code, tests, adapters, and compatibility layers more readily than they delete
them. In this repo, extra source, tests, docs, metadata, compatibility glue,
and adapter code all have carrying cost. New artifacts should earn that cost
by protecting a live contract, improving clarity, reducing duplication,
improving performance, or enabling a current workflow.

Scientific acceptance checks should force the intended algorithm to exist. The
H/H2+ acceptance checkpoint is the current example: the full-window CPB
workaround produced plausible energies, but it was not a decomposed
White-Lindsey route and was retired rather than kept as an active passing
test. The active readiness audit now materializes decomposed route-global
overlap, kinetic, one-center electron-nuclear by-center, and a charge-applied
one-electron Hamiltonian. It blocks at the solve boundary because the current
decomposed overlap metric is singular:
`:decomposed_wl_overlap_metric_not_positive_definite`.

## Current Replacement Pressure

- Selected global overlap, kinetic, position_x, position_y, position_z, x2_x,
  x2_y, and x2_z equivalence are now covered for White--Lindsey face/face and
  face/edge slices.
- Synthetic global retained matrix pilots now exist for the safe one-body
  terms `:overlap`, `:kinetic`, `:position_x`, `:position_y`, `:position_z`,
  `:x2_x`, `:x2_y`, and `:x2_z`.
- A private route-shaped global one-body adapter now composes
  `PairBlockMaterializationPlan` inputs through the local block collection,
  term-specific placement plans, and term-specific dense global retained matrix
  pilots for the individual safe one-body terms `:overlap`, `:kinetic`,
  `:position_x`, `:position_y`, `:position_z`, `:x2_x`, `:x2_y`, and `:x2_z`.
  This is module-level adapter coverage, not route-driver wiring, term summing,
  one-body Hamiltonian construction, or Hamiltonian assembly.
- A narrow decomposed WL one-electron Hamiltonian assembly helper now consumes
  route-global kinetic and separated route-global electron-nuclear by-center
  matrices. By-center matrices remain uncharged and center-separated; recorded
  nuclear charges and center summation are applied only at Hamiltonian assembly.
  The current H atom audit reaches this assembled Hamiltonian but still does not
  produce an accepted scientific energy because the decomposed overlap solve
  metric is singular.
- A private route-shaped safe one-body matrix-set adapter now wraps individual
  `route_global_one_body_matrix(...; term = ...)` calls and returns
  term-separated result objects plus a compact summary. It does not copy dense
  matrices out of term results, allows partial success with blocked per-term
  results for unsupported terms, and is still not term summing, a one-body
  Hamiltonian object, or Hamiltonian assembly. PQS source-space records remain
  blocked inside term results until shell projection/Lowdin realization exists.
- A private overlap-only route-state adapter now accepts structured state that
  already carries a `PairBlockMaterializationPlan` and delegates to the existing
  route-global overlap matrix path. This is replacement pressure only: it is
  not production route-driver wiring and does not retire the old route by
  itself.
- A private overlap-only driver-facing hook,
  `driver_global_overlap_result`, delegates structured state carrying
  `pair_block_materialization_plan` to the route-state overlap adapter. It
  blocks with `:missing_pair_block_materialization_plan` when no structured
  pair-block materialization plan is present. It is not production driver
  wiring and does not retire the old route by itself.
- A private opt-in overlap driver option,
  `private_global_overlap_requested`, carries
  `private_global_overlap_result` plus `private_global_overlap_summary` from the
  materialization stage. It returns the structured overlap result when enough
  structured state, global dimension, and overlap inputs are present, or a
  compact blocked summary otherwise. This remains nonproduction and off by
  default; it is not old-route retirement. The override example
  `examples/private_global_overlap_option.jl` documents option shape only, and
  focused tests validate that option plumbing without running the full builder
  route.
- Private overlap input facts now label the 1D factor space and convention
  explicitly, following the legacy PGDG lesson that overlap factors are
  axis-space objects with a specific contraction convention.
- Real dry-run overlap fingerprints currently split by parent-axis construction
  mode. The manual-count dry report keeps `probe_parent_axis_construction =
  false`, derives `global_dimension = 221` from retained-dimension
  compatibility and `parent_axis_counts = (5, 5, 5)` from route-axis counts, but
  remains blocked on `:missing_parent_axis_bundle_overlap_factors`. The
  probe-enabled dry report carries a structured
  `route_materializer_payload.parent_axis_bundle_object`; private overlap input
  facts can read that object as
  `:route_materializer_payload_parent_axis_bundle_object`, yielding available
  overlap facts with factor provenance
  `:parent_axis_bundle_pgdg_intermediate` /
  `:axis_bundle_one_body_overlap`. This is a fingerprint of carried facts only,
  not overlap driver adoption or route retirement.
- The local CPB overlap provider chain is implemented through
  `parent_overlap_axis_factor_packet -> cpb_interval_pair ->
  cpb_overlap_axis_blocks -> cpb_overlap_dense_block ->
  cpb_local_overlap_block_record -> cpb_local_overlap_block_collection`. The
  real probe-enabled report can build one local overlap collection from the
  `(:product, :product)` source pair with dense local block shape `(25, 25)`;
  the placement-candidate probe uses `report.retained_dimension` as its
  retained/global dimension source. The private local collection adapter
  recognizes the collection as structured local overlap source data, but global
  overlap remains blocked on `:missing_placement_or_retained_transform`.
  Placement candidates now split partial placeholder facts into available and
  missing requirements, and even all placeholder facts remain blocked on
  `:placement_not_implemented`. This milestone claims no retained transforms,
  no placement plan, no global matrix materialization, no route-global overlap
  stage adoption, and no kinetic, position, x2, Coulomb, Hamiltonian, IDA/MWG,
  PQS projection/Lowdin, export, or artifact work.
- A real-report overlap placement source audit now shows that the
  probe-enabled report carries a structured retained dimension source through
  `report.retained_units`; `report.retained_dimension` is not the only
  dimension source in that fixture. The same audit still finds no structured
  retained transforms, no retained column ranges for the `(:product, :product)`
  local overlap source pair, no reviewed overlap placement plan, and no overlap
  accumulation rule. With only `:local_cpb_overlap_collection` and
  `:global_dimension` available, the private placement candidate and placement
  plan skeleton remain blocked on
  `:missing_placement_or_retained_transform`; missing requirements remain
  `:missing_retained_transform`, `:missing_left_column_range`,
  `:missing_right_column_range`, `:missing_placement_plan`, and
  `:missing_accumulation_rule`. The next implementation boundary is therefore
  structured retained-transform and column-range carry design, not numerical
  placement. This claims no placement, retained transforms, global matrix
  materialization, route-global overlap adoption, kinetic, position, x2,
  Coulomb, Hamiltonian, IDA/MWG, PQS projection/Lowdin, export, or artifact
  work.
- The overlap placement metadata boundary is now explicit. The provider layer
  owns `CPBRetainedTransformCarry`, `CPBSourcePairPlacementRange`, and
  `CPBOverlapPlacementFacts`; the private placement skeleton can consume
  `CPBOverlapPlacementFacts` directly. The real probe report negative
  fingerprint builds placement facts from the real local collection without
  placeholder transforms or ranges and still reports only
  `available_requirements = (:local_cpb_overlap_collection,)`. Missing
  requirements are `:missing_retained_transform`,
  `:missing_left_column_range`, `:missing_right_column_range`,
  `:missing_global_dimension`, `:missing_placement_plan`, and
  `:missing_accumulation_rule`. Left and right CPB summaries are preserved
  through placement facts into the private skeleton. Missing placement ranges
  also report `:missing_global_dimension`, and non-matrix retained-transform
  references block with `:unsupported_retained_transform_reference`. This still
  claims no transform application, placement, global overlap accumulation,
  route adoption, kinetic, position, x2, Coulomb, Hamiltonian, IDA/MWG, PQS
  projection/Lowdin, export, or artifact work.
- The reviewed overlap placement plan object now exists as metadata-only
  contract data. It owns placement plan kind, accumulation rule, symmetry or
  transpose policy, duplicate record policy, accepted block keys and record
  inventory, required global dimension source, status/blocker, and route/global
  nonclaim flags. `CPBOverlapPlacementFacts` and the private placement
  skeleton can use it instead of placeholder `placement_plan` and
  `accumulation_rule` values. A real-report fingerprint can supply this
  reviewed plan, making `:placement_plan` and `:accumulation_rule` available
  while real retained transforms, source-pair column ranges, and
  placement-derived global dimension remain missing. It still does not apply
  transforms or assemble a matrix. Non-reviewed `placement_plan` values are
  now labeled as `:placeholder_placement_plan_compatibility`, preserving
  transition compatibility without treating them as reviewed placement
  authority.
- `CPBOverlapPlacementFacts` now validates the reviewed plan's accepted record
  inventory against local collection block keys. A matching real-report
  fingerprint accepts the `(:product, :product)` local overlap source pair;
  unaccepted records block with `:unaccepted_overlap_placement_record`, and
  duplicate block keys under `:reject_duplicate_block_keys` block with
  `:duplicate_overlap_placement_record`. It also validates local ordering
  contracts for available records and required global dimension source for
  available placement ranges; missing ranges still report missing range and
  dimension requirements rather than a source mismatch. This is compact
  metadata validation only; it does not apply transforms, place blocks,
  assemble global overlap, or adopt route-global overlap.
- The overlap placement metadata gates are now sufficient for a strategic
  implementation fork rather than more status layers. The next step should be
  either a real-source carry audit/design for retained transforms and
  source-pair column ranges, or, if no credible real source is nearby, a tiny
  synthetic provider-level numerical overlap placement pilot using fully
  reviewed metadata. The preferred order is real-source carry first, then a
  deliberately synthetic pilot if the audit does not find usable route facts
  quickly. A broader retained-unit metadata audit remains an option if the
  existing route metadata is ambiguous. None of these options is route-driver
  adoption by itself.
- The first real-source carry audit for the probe-enabled `(:product,
  :product)` local overlap source found source support count, local ordering,
  reviewed-plan metadata, accumulation rule, and retained/global dimension
  source, but no structured retained-transform source and no source-pair
  retained column ranges. Current blockers are
  `:missing_retained_transform_source` and
  `:missing_source_pair_column_ranges`. Without a newly identified real carry
  source, the next implementation should be the tiny synthetic provider-level
  numerical overlap placement pilot, not another metadata layer. This still
  claims no transform application, placement, global matrix assembly, or
  route-global overlap adoption.
- A tiny synthetic provider-level numerical overlap placement pilot now exists
  for one fully reviewed local dense overlap block. It applies
  `T_left' * O_cpb * T_right` and accumulates into an explicitly supplied dense
  provider-level global overlap matrix range. The pilot is gated by available
  transform carries, placement range, reviewed placement facts, reviewed
  placement plan status, explicit accumulation rule, explicit symmetry policy,
  and duplicate-record policy. It is synthetic/provider-level only: no real
  route/report placement, private driver global overlap adoption, route-global
  overlap availability, Hamiltonian assembly, or new physics path is claimed.
  Provider-level materialization is recorded with provider-level summary fields;
  unqualified and route-global matrix materialization flags remain false.
- The synthetic provider-level overlap placement pilot also has a small
  collection wrapper. It matches local overlap records by `block_key`, requires
  dense local source blocks plus reviewed placement facts, and adds each
  retained block into one provider-level target matrix under the explicit
  accumulation rule. It still does not infer transpose/symmetry fill, consume
  real route placement facts, wire the private driver, or claim route-global
  overlap availability.
- The provider-level overlap placement pilot is now a milestone, not a prompt
  for more synthetic feature work. The implemented synthetic boundary includes
  one-block placement, collection placement, additive accumulation under
  `:add_explicit_blocks_into_ranges`, reviewed metadata gates, collection/facts
  alignment, compact missing transform/range block-key coverage, provider-level
  dimension consistency, and compact summaries with provider-level matrix flags
  distinct from route/global flags. Remaining route blockers are real
  retained-transform sources, real source-pair retained column ranges,
  route-owned placement-plan/range carry, route-global overlap adoption, and
  driver wiring. Prefer architectural review or a broader retained-unit /
  terminal-route audit before any further implementation; do not keep extending
  synthetic placement features without a reviewed plan for real transforms and
  ranges.
- Architecture pivot: the overlap pilot is better evidence for a CPB-local
  operator-block layer than for route-global placement machinery. The parent /
  axis-factor layer owns universal 1D factors. The CPB operator layer fills
  rectangular local blocks for source CPB pairs and may materialize dense
  provider-level local matrices while remaining realization-neutral and
  route-neutral. The realization layer then differs by route: White-Lindsey can
  consume CPB blocks nearly directly, while PQS applies local retained
  transforms such as `O_retained = T_left' * O_cpb * T_right`. The route/global
  layer should assign retained/global ranges and accumulate global matrices
  later. Overlap is the current pilot term; kinetic, position, and x2 should
  become analogous one-body CPB operator blocks later, and Coulomb/two-body
  work should remain CPB-local but may require pair-pair or factorized records.
  Do not add more placement metadata layers before a reviewed WL/PQS
  realization design identifies real transform and range ownership. The next
  likely implementation should be an overlap-only cleanup/renaming pass or a
  small generic `CPBOperatorBlock` design, not more placement code.
- A generic CPB-local axis-product operator block primitive now exists. It
  materializes one dense local product-space block from prepared 1D axis
  operators `(x = Ox, y = Oy, z = Oz)` using x-slowest, z-fastest ordering and
  compact route-neutral metadata. Overlap is a thin wrapper over that primitive:
  parent overlap packet -> CPB interval pair -> overlap axis blocks ->
  axis-product block with `term = :overlap`. This is still CPB-local operator
  construction, not WL/PQS realization and not route/global placement. Future
  kinetic, position, and x2 CPB one-body blocks should be sums of axis-product
  terms such as `Kx Sy Sz + Sx Ky Sz + Sx Sy Kz`, `Xx Sy Sz`, and
  `Sx X2y Sz`; inactive directions should use explicit overlap factors, not
  ambiguous identity labels.
- The CPB operator-block layer now distinguishes kernel families. Axis-product
  and sum-of-axis-products kernels cover simple separable one-body terms such
  as overlap, kinetic pieces, position, and x2. Electron-nuclear and
  electron-electron Coulomb-family blocks should use specialized CPB-local
  Gaussian-sum kernels with the Gaussian expansion index as an inner loop, not
  an outer loop over repeated calls to `cpb_axis_product_operator_block`. The
  shared abstraction should eventually be compact result records and metadata,
  not necessarily a single universal kernel implementation.
- A small `cpb_sum_of_axis_products_operator_block` primitive now materializes
  dense CPB-local sums of prepared axis-product terms with compact per-term
  summaries. It is for simple separable one-body terms only; no production
  kinetic, position, x2, Coulomb-family kernel, WL/PQS realization, or
  route/global placement was added with it.
- Parent axis bundle one-body factor support has been audited and added
  narrowly where the bundle already owns structured 1D matrices. The audited
  real-bundle paths are `axis.pgdg_intermediate.kinetic`,
  `axis.pgdg_intermediate.position`, and `axis.pgdg_intermediate.x2`; structured
  top-level fallbacks `axis.kinetic`, `axis.position`, and `axis.x2` are also
  accepted. The parent axis factor packet now carries optional `kinetic_1d`,
  `position_1d`, and `x2_1d` with matrix/count validation and compact category
  metadata. CPB-local wrappers materialize kinetic as
  `Kx Sy Sz + Sx Ky Sz + Sx Sy Kz`, position as `Xi` on one active axis with
  overlap factors on inactive axes, and x2 similarly. This does not add
  Coulomb-family kernels, WL/PQS realization, retained transforms,
  route/global placement, driver wiring, Hamiltonian assembly, IDA/MWG, PQS
  projection/Lowdin, exports, or artifacts.
- Coulomb-family CPB kernels have been audited against
  `PureGaussianGausslet.jl` and should not be forced through the simple
  axis-product one-body primitive. Existing separable kernels
  `cpb_axis_product_operator_block` and
  `cpb_sum_of_axis_products_operator_block` cover overlap, kinetic, position,
  x2, and similar one-body terms. Electron-nuclear Coulomb needs a separate
  CPB-local Gaussian-sum one-body kernel over left/right CPBs, parent
  Coulomb/nuclear factors, nucleus data, and Gaussian expansion data.
  Electron-electron Coulomb needs a separate CPB-local two-body or pair-pair
  Gaussian-sum kernel. In both cases the Gaussian expansion index must be an
  inner optimized loop or contraction inside the kernel, following the
  `getAtomPot`, `gethamsNoPot`, and `addinGaussians` pattern of precomputed or
  transformed 1D Gaussian factors followed by local alpha summation. Do not
  allocate one dense CPB product matrix per Gaussian term unless a later
  performance review justifies it. No Coulomb production code, parent Coulomb
  packet structs, WL/PQS realization, route/global placement, driver wiring, or
  Hamiltonian assembly was added with this audit.
- The current parent axis-bundle Coulomb source audit found structured
  ingredients but no reviewed parent Coulomb factor packet. Axis bundles carry
  centered Gaussian one-body terms at
  `axis.pgdg_intermediate.gaussian_factor_terms` and
  `axis.pgdg_intermediate.gaussian_factors`, pair Gaussian ingredients at
  `axis.pgdg_intermediate.pair_factor_terms`,
  `axis.pgdg_intermediate.pair_factors`, and
  `axis.pgdg_intermediate.pair_factor_terms_raw`, and expansion exponents at
  `axis.pgdg_intermediate.exponents` / `axis.exponents`. They do not carry
  expansion coefficients or a full `CoulombGaussianExpansion` object. Real
  reports carry center data separately through `report.nuclear_charges`,
  `report.atom_locations`, and `report.center_table`, while
  `report.route_materializer_payload.parent_axis_bundle_object` carries the
  axis bundle. One-body electron-nuclear factors are not currently carried as
  per-nucleus Gaussian axis terms on the parent bundle; existing nuclear
  helpers build those per-center term tables on demand. Electron-electron
  1D pair factors are present as ingredients, but not as a CPB pair-pair
  Coulomb source record. No numerical Coulomb kernel, packet struct,
  WL/PQS realization, route placement, driver wiring, Hamiltonian assembly,
  IDA/MWG change, PQS projection/Lowdin, export, or artifact was added.
- The current Cartesian/White-Lindsey Coulomb port map identifies the working
  code that should guide the CPB-local adaptation. Electron-nuclear construction
  currently lives mainly in `src/ordinary_qw_raw_blocks.jl` through
  `_qwrg_diatomic_factor_term_cache`,
  `_qwrg_contracted_nuclear_axis_term_tables`,
  `_qwrg_fill_direct_contracted_nuclear_matrix!`,
  `_qwrg_fill_staged_nuclear_submatrix!`,
  `_qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center`,
  `_qwrg_bond_aligned_staged_by_center_nuclear_one_body_by_center`, and
  `_qwrg_diatomic_nuclear_one_body_by_center`, with final/by-center packaging
  in `src/ordinary_qw_operator_assembly.jl`. Electron-electron construction
  currently uses `_qwrg_gausslet_interaction_matrix`,
  `_qwrg_diatomic_interaction_matrix`, `_qwrg_fixed_block_interaction_matrix`,
  `_qwrg_interaction_matrix_nearest`, `_qwrg_mwg_interaction_components`, and
  `_qwrg_bond_aligned_molecular_interaction_matrix` to produce WL/QW
  two-index density-density interaction matrices. The smallest CPB boundary is
  not full Hamiltonian export: first adapt the existing parent/axis sources into
  a compact Coulomb packet/source summary, then port the staged nuclear
  submatrix fill to a CPB-local by-center one-body kernel and the pair-factor
  interaction fill to a CPB pair-pair density-density record. Existing global
  Hamiltonian/export tests remain oracle surfaces only.
- A metadata-only parent Coulomb source adapter now exists as
  `CartesianParentGaussletBases.parent_coulomb_axis_source_summary`. It consumes
  existing structured inputs: parent axis bundle, `CoulombGaussianExpansion`,
  and optional center metadata from route/parent/QW objects. It reports compact
  availability and source paths for expansion coefficients/exponents,
  Gaussian one-body ingredients, electron-electron pair-axis ingredients, and
  center metadata. The current expected fingerprint remains partial:
  electron-electron pair-axis ingredients can be available, but per-center
  electron-nuclear axis term tables and CPB pair-pair Coulomb source records
  are not built or carried by this adapter. No numerical Coulomb block, CPB
  Coulomb kernel, WL/PQS realization, route/global placement, Hamiltonian
  assembly, IDA/MWG/PQS semantic change, export, or artifact was added.
- The mixed gausslet/GTO supplement port map is now documented as a CPB-local
  operator-layer adaptation of existing polynomial-Gaussian code, not a formula
  rewrite. The reusable source surfaces are
  `_cartesian_basis_supplement_axis_primitive_cross`,
  `_cartesian_basis_supplement_cross`,
  `_qwrg_atomic_axis_cross_data`, `_qwrg_atomic_axis_aa_data`,
  `_qwrg_cartesian_shell_cross_moment_blocks_3d`,
  `_qwrg_cartesian_shell_self_moment_blocks_3d`,
  `_qwrg_atomic_axis_factor_cross_data`,
  `_qwrg_atomic_axis_factor_aa_data`,
  `_qwrg_diatomic_cartesian_shell_overlap_blocks_3d`, and
  `_qwrg_diatomic_cartesian_shell_blocks_3d`, with
  `GaussianAnalyticIntegrals` as the polynomial-Gaussian formula source. The
  first mixed overlap test now exists as a tiny provider-level CPB-local pilot:
  `CPBMixedGTOLocalOverlapBlock` / `cpb_mixed_gto_overlap_block(parent, cpb,
  orbital)` reuses `_cartesian_basis_supplement_axis_primitive_cross` and
  compares local CPB rows against `_cartesian_basis_supplement_cross`. It is
  one supplement orbital only. The simple CPB-local mixed one-body pilots now
  cover all `position_*` and `x2_*` coordinate moments through
  `CPBMixedGTOCoordinateMomentLocalBlock` plus `cpb_mixed_gto_position_operator_block`
  and `cpb_mixed_gto_x2_operator_block`, and kinetic through
  `CPBMixedGTOKineticLocalBlock` / `cpb_mixed_gto_kinetic_operator_block`. They
  reuse existing QW polynomial-Gaussian axis-integral wrappers and compare
  local rows against `_qwrg_cartesian_shell_cross_moment_blocks_3d` fields
  `position_*_ga`, `x2_*_ga`, and `kinetic_ga`. No mixed GTO nuclear/Coulomb
  term, route/global placement, Hamiltonian assembly, IDA/MWG/PQS semantic
  change, export, or artifact was added.
- The whole-supplement source boundary for CPB-local GTO work now has a focused
  fingerprint over existing objects and the modern conversion surface. Vendored
  `data/legacy/BasisSets` H/He cc-pVTZ fixtures load through
  `LegacyAtomicGaussianSupplement`,
  `LegacyBondAlignedDiatomicGaussianSupplement`, and
  `LegacyBondAlignedHeteronuclearGaussianSupplement`, then convert through
  `basis_representation(...)` into
  `CartesianGaussianShellSupplementRepresentation3D`. The fingerprint records
  source kind, atom/basis or nuclei metadata, `lmax`, contracted/uncontracted
  convention, orbital labels, angular powers, centers, exponents,
  coefficients, and `:axiswise_normalized_cartesian_gaussian` primitive
  normalization.
- Whole-supplement CPB-local mixed GTO one-body blocks now consume
  `CartesianGaussianShellSupplementRepresentation3D` directly through the
  existing overlap, position, x2, and kinetic wrapper names. The returned
  `CPBMixedGTOSupplementLocalBlock` stores a dense provider-level local matrix
  shaped as CPB support rows by GTO orbital columns, plus compact
  whole-supplement and per-orbital summaries. Tests compare overlap,
  representative position/x2 axes, and kinetic against existing Cartesian/QW
  oracle rows for contracted and uncontracted source counts. This is not yet a
  GTO/GTO one-body Galerkin record, mixed nuclear/Coulomb term, route/global
  placement, Hamiltonian assembly, IDA/MWG/PQS semantic change, export, or
  artifact.
- GTO/GTO provider-level one-body Galerkin records now exist through
  `CPBGTOSupplementOneBodyBlock` and wrappers for overlap, kinetic, position,
  and x2. Overlap reuses `_cartesian_supplement_cross_overlap`; kinetic,
  position, and x2 reuse `_qwrg_cartesian_shell_self_moment_blocks_3d`. Tests
  cover contracted He cc-pVTZ `lmax = 0` overlap, kinetic, all position axes,
  all x2 axes, and an uncontracted overlap count. This is a supplement-sector
  one-body provider record only; no mixed nuclear/Coulomb term, route/global
  placement, Hamiltonian assembly, IDA/MWG/PQS semantic change, export, or
  artifact was added.
- The next overlap implementation boundary is no longer additional placement
  fingerprinting. First decide the CPB operator-block and WL/PQS realization
  design: what local block objects exist, how White-Lindsey consumes them, how
  PQS obtains real retained transforms, and where source-pair ranges are owned.
  Only after that should a route/global placement plan define column/range
  placement, global dimension, accumulation rule, ordering, and blockers.
- Selected White--Lindsey old-oracle equivalence currently covers overlap,
  kinetic, position_x, position_y, position_z, x2_x, x2_y, and x2_z.
- `src/cartesian_pair_block_materialization/one_body_global_matrix_helpers.jl`
  owns only behavior-neutral symmetric placement validation/insertion. The
  term-specific global matrix files still own result statuses, blockers,
  object kinds, and materialization flags.
- Route-global adapter validation should stay focused:
  `cartesian_pair_block_route_global_one_body_adapter_runtests.jl` is the
  broader individual-term adapter contract at about 40--45 seconds, and
  `cartesian_pair_block_route_global_matrix_set_smoke_runtests.jl` is the
  focused matrix-set smoke at about 40 seconds after slimming. Do not add these
  casually to broader default runners; old-oracle tests remain gate/oracle
  tests, not routine per-pass checks.
- The next replacement slice should decide whether to broaden selected
  safe-term oracle coverage or move toward carefully scoped driver wiring.
- The current global one-body pilots do not assemble Hamiltonians, build
  one-body Hamiltonian objects, sum terms, build Coulomb, touch IDA/MWG, build
  density-density, nuclear attraction, or Gaussian local terms, perform PQS
  projection/Lowdin, wire route drivers, export artifacts, or assemble a full
  White--Lindsey route.

## CartesianContractedParentMetrics CPB-Layer Retirement Audit

Audit target: `src/CartesianContractedParentMetrics.jl` as of the CPB-local
operator layer pivot. This file is about 19.8k lines and currently mixes five
different roles: core retained metric packet construction, source-box operator
kernels, oracle comparisons, route diagnostics/report writers, and private
shadow/fallback transition paths. No code was deleted in this audit.

Classification labels:

- `KEEP_CORE_METRICS`: retained overlap, weights, centers, first moments, and
  retained position packet construction.
- `MOVE_TO_CPB_PROVIDER`: local source-box/operator block kernels now better
  expressed as CPB-local operator blocks.
- `KEEP_AS_ORACLE_REFERENCE`: comparison helpers useful for testing CPB blocks
  or route adapters.
- `MOVE_TO_ROUTE_DIAGNOSTICS`: readiness, route-smoke, schema, report, and
  inventory diagnostics.
- `RETIRE_OR_QUARANTINE`: shadow, fallback, and private transition paths
  superseded by the CPB operator layer once tests are moved.

Classified regions by line/function family:

| Line range | Function family | Classification | Notes |
| --- | --- | --- | --- |
| `65-313` | `CartesianContractedParentMetricPacket3D`, `_AxisMetricData1D`, coefficient-entry helpers, `_contract_pair_matrix`, `_contract_linear_vector`, `_contracted_first_moments` | `KEEP_CORE_METRICS` | This is the retained metric packet core for overlap, weights, centers, and position moments. |
| `315-655` | `_metric_dispatch_*`, `_pqs_product_mixed_block_policy`, `_contracted_parent_metric_dispatch_shadow_plan` | `MOVE_TO_ROUTE_DIAGNOSTICS`; `RETIRE_OR_QUARANTINE` for the shadow wrapper | These are readiness/dispatch diagnostics, not CPB kernels. |
| `657-809` | staged-axis interval/count/projection helpers and `_product_doside_retained_unit_plan` | `KEEP_AS_ORACLE_REFERENCE` | Useful retained transform/range facts for old product/doside units; do not make them new provider authority. |
| `811-1221` | product/doside source-box pair plans, safe one-body factor blocks, reference blocks, and `_product_doside_source_box_shadow_blocks` | `MOVE_TO_CPB_PROVIDER`; `RETIRE_OR_QUARANTINE` for shadow layout | Safe one-body local block math maps to CPB axis-product and sum-of-axis-products blocks. |
| `1222-1705` | local Gaussian-sum and density-density product/doside source-box helpers | `MOVE_TO_CPB_PROVIDER` | CPB-local Coulomb-family kernels or precursors. |
| `1716-2026` | IDA pair-factor provenance and raw-weight to density-normalized conversion | `KEEP_AS_ORACLE_REFERENCE` / future CPB source adapter | Keep as convention oracle until CPB source summaries carry the same facts. |
| `2028-2177` | centered local Gaussian term-table construction | `MOVE_TO_CPB_PROVIDER` | Belongs with CPB-local electron-nuclear/electron-electron Gaussian-sum kernels. |
| `2178-2294` | product/doside retained low-order, separable-sum, kinetic retained blocks | `KEEP_AS_ORACLE_REFERENCE` | Useful retained-block oracle while CPB blocks and realization transforms are validated. |
| `2295-2580` | PQS/product support-local low-order and kinetic reference blocks | `KEEP_AS_ORACLE_REFERENCE` | Compares factored PQS/product blocks against support-local oracle entries. |
| `2599-3352` | raw product-box mode matrices, structural plans, operator factors, shell realization plans | `MOVE_TO_ROUTE_DIAGNOSTICS` | Route-state/readiness facts; future home is route-spine modules. |
| `3364-4557` | PQS/product source-box safe one-body, Gaussian-sum, density-density, and by-center nuclear helpers | `MOVE_TO_CPB_PROVIDER`; selected outputs `KEEP_AS_ORACLE_REFERENCE` | Local operator kernels belong to CPB provider/result records. |
| `4559-5722` | PQS/PQS source-box safe one-body, Gaussian-sum, density-density, and by-center nuclear helpers | `MOVE_TO_CPB_PROVIDER`; selected outputs `KEEP_AS_ORACLE_REFERENCE` | Same as PQS/product, but for PQS/PQS raw-box pairs. |
| `5722-6170` | PQS/product source-box shadow blocks and all-pairs inventory | `RETIRE_OR_QUARANTINE` | Private shadow layout; replace with CPB provider collections plus explicit route inventory. |
| `6172-6793` | PQS/PQS/product route descriptor, raw-box geometry facts, raw-box route producer | `MOVE_TO_ROUTE_DIAGNOSTICS` / `RETIRE_OR_QUARANTINE` | Fixture route producer and descriptor diagnostics. |
| `6795-8042` | route descriptor diagnostics, route fact audits, contact-cap and mismatch fixture construction | `MOVE_TO_ROUTE_DIAGNOSTICS` | Readiness and fixture diagnostics. |
| `8042-8897` | contact-cap, outer-mismatch, and atom-box safe-term operator comparisons | `KEEP_AS_ORACLE_REFERENCE` | Old retained/support path comparisons, not provider-layer authority. |
| `8913-12153` | current-route inventories, shell-realization transform facts, source metadata exports, safe-term matrices, authority comparisons | `MOVE_TO_ROUTE_DIAGNOSTICS`; selected comparison payloads `KEEP_AS_ORACLE_REFERENCE` | Route/report inventory and authority comparison plumbing. |
| `12155-12617` | route descriptor diagnostic and three-unit safe-term shadow consumer | `RETIRE_OR_QUARANTINE` | Private route-shaped shadow consumer. |
| `12619-16016` | density-density and nuclear-attraction route producers/consumers | `MOVE_TO_CPB_PROVIDER` for pair kernels; `MOVE_TO_ROUTE_DIAGNOSTICS` for route-shaped assembly | Local pair blocks belong to CPB provider Coulomb-family kernels; retained route assembly is diagnostic. |
| `16018-18187` | component route smoke, summaries, report adapters, sidecar/schema writers | `MOVE_TO_ROUTE_DIAGNOSTICS` | Smoke reports and schema/report writers. |
| `18187-18877` | raw-product reference blocks, fallback staged separable/kinetic/metric/x2 blocks, staged linear vectors | `KEEP_AS_ORACLE_REFERENCE`; fallback entry points `RETIRE_OR_QUARANTINE` | Useful dense/support-local oracle contractions. |
| `18892-19255` | packet-build safe-field shadows, `_cartesian_packet_build_source_safe_field_shadow` | `RETIRE_OR_QUARANTINE` | Explicit shadow/source-field transition path. |
| `19256-19749` | packet diagnostics, `_metric_packet_from_matrices`, `_support_local_metric_packet`, `_resolved_payload_product_staged_metric_packet`, public packet constructors, dense reference oracle | `KEEP_CORE_METRICS`; dense reference constructor `KEEP_AS_ORACLE_REFERENCE` | This is the CCPM core API surface; dense parent reference remains tiny-test oracle only. |

Tests that still call private CCPM helpers:

- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
  is the main private-helper dependency. It calls raw-box plan producers,
  source-box reference/shadow helpers, density-density and nuclear source-box
  helpers, retained low-order/kinetic helpers, route-shaped consumers,
  descriptor diagnostics, and staged entry helpers.
- `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
  calls current-route inventory, retained pair inventory, current-route
  safe-term matrices, contact-cap / outer-mismatch / atom-box comparison
  helpers, route fact audit, shell-realization transform fact, and authority
  comparison helpers.
- `test/nested/bond_aligned_diatomic_endcap_panel_shared_shell_source_policy_runtests.jl`
  calls packet-build shadow helpers, staged-unit entries, retained kinetic
  shadow matrices, resolved-payload product-staged metric packet, and dispatch
  shadow plan helpers.

Recommended first cleanup target:

Shrink the product/product safe one-body source-box test burden instead of
moving its full transitional architecture into new scaffolding. The CPB
provider layer already owns active coverage for overlap, kinetic, position,
x2, axis-product, and sum-of-axis-products blocks. The old
`_product_doside_source_box_pair_plan(...)`,
`_product_doside_source_box_block_from_factors(...)`,
`_product_doside_source_box_reference_block(...)`, and
`_product_doside_source_box_shadow_blocks(...)` family should remain only a
private oracle / route-shadow bridge until the route diagnostics are retired.
Do not add new routine tests that preserve every old diagnostic flag.

First retirement cut:

- The exhaustive product/product source-box pair-plan assertions and
  shadow-component oracle comparisons were removed from
  `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`.
  The slow integration file now keeps only a small legacy smoke that
  `_product_doside_source_box_shadow_blocks(...)` returns finite blocks with
  expected terms, ranges, and shapes while route diagnostics still call it.
- Detailed safe one-body correctness is intentionally left to focused CPB
  provider tests for axis-product, sum-of-axis-products, overlap, kinetic,
  position, and x2 blocks.
- `src/CartesianContractedParentMetrics.jl` marks the product/product
  source-box family as a private oracle and route-shadow bridge. New safe
  one-body product-box checks should use CPB provider axis-product or
  sum-of-axis-products blocks directly.
- Deletion is still blocked by route-shadow and diagnostics paths in
  `CartesianContractedParentMetrics.jl` that call
  `_product_doside_source_box_reference_block(...)` internally, especially
  `_product_doside_source_box_shadow_blocks(...)` and downstream current-route
  inventory/oracle comparisons. The next deletion-oriented cut should replace
  those route-shadow component calls with CPB provider blocks or quarantine the
  whole shadow route as oracle-only.

Second deletion pass:

- The two-unit PQS/product shadow helper family was deleted:
  the raw-plan helper, descriptor overload, and private all-pairs inventory.
  At that checkpoint, the three-unit
  `_pqs_pqs_product_source_box_shadow_blocks(...)` and
  `_pqs_pqs_product_source_box_all_pairs_inventory(...)` family still existed.
- The matching slow-test helper in
  `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl` was
  removed with its call site.
- Removed assertions that preserved old all-pairs inventory vocabulary,
  reference-helper names, dense-raw-storage flags, retained-weight flags, and
  explicit sub-block equality against development-era reference helpers. Those
  details are transitional architecture, not active scientific contract.
- No current source or test caller requires the deleted two-unit helper.

Third deletion-oriented pass:

- The three-unit `_pqs_pqs_product_source_box_shadow_blocks(...)` and
  `_pqs_pqs_product_source_box_all_pairs_inventory(...)` family was not
  deleted in that pass because `_pqs_pqs_product_raw_box_route_producer(...)`
  still built the inventory for descriptor validation.
- The slow integration tests for this family were reduced to smoke coverage:
  expected ranges, retained dimension, finite block payloads, and route-shaped
  descriptor checks. Detailed all-pairs inventory vocabulary,
  helper-name assertions, reference sub-block equality, and repeated metadata
  nonclaims were removed.
- `src/CartesianContractedParentMetrics.jl` marks this family as
  route-diagnostic / oracle-only. It should not be extended as a production
  placement or provider-layer contract.

Fourth caller-driven audit:

- `_pqs_product_source_box_reference_blocks_from_pair_plan(...)` and
  `_pqs_product_source_box_reference_blocks(...)` were `KEEP_ORACLE_ONLY` at
  that point. Their exact source callers were inside the then-live
  `_pqs_pqs_product_source_box_shadow_blocks(...)`.
- `_pqs_pqs_source_box_reference_blocks_from_pair_plan(...)` and
  `_pqs_pqs_source_box_reference_blocks(...)` were also `KEEP_ORACLE_ONLY`.
  Exact source callers were the PQS/PQS reference wrapper and the then-live
  `_pqs_pqs_product_source_box_shadow_blocks(...)`.
- `_pqs_pqs_product_raw_box_route_producer(...)` is `SHRINK_TEST_ONLY`, not a
  deletion target in this pass. Exact source callers are
  `_pqs_pqs_product_raw_box_route_from_geometry_facts(...)` and the
  density-density route producer path. Tests should not preserve every raw-box
  helper flag as contract.
- `_pqs_pqs_product_route_shaped_safe_term_consumer(...)` was `DELETE_NOW`.
  Exact caller audit found no source caller outside its own definition; only
  slow-test route-shadow assertions exercised it. The source function and those
  slow-test call sites were deleted rather than migrated.
- The slow PQS local-layer integration test was reduced again: the
  PQS/product mixed source-box check no longer keeps the unused shared raw-plan
  scaffold, no longer exercises the single-term reference wrapper, and no
  longer asserts unsupported-term helper vocabulary. It now keeps only compact
  pair-plan shape/status checks plus finite numerical equivalence for the
  active multi-term oracle. CPB provider tests remain responsible for detailed
  local one-body operator correctness.

Fifth batched caller-driven sweep:

| Helper family | Caller classification | Action |
| --- | --- | --- |
| `_product_doside_source_box_reference_block(...)` | `KEEP_ORACLE_ONLY` / `MOVE_LATER_TO_CPB` | Kept. Source callers remain in product/doside shadow/reference comparisons, route authority comparisons, and product-product entries in `_pqs_pqs_product_source_box_all_pairs_inventory(...)`. |
| `_pqs_product_source_box_reference_block(...)` | `DELETE_NOW` | Deleted. Exact caller audit showed only slow-test/docs references; the surviving multi-term PQS/product oracle covers the remaining narrow private check. |
| `_pqs_product_source_box_reference_blocks(...)` | `KEEP_ORACLE_ONLY` | Kept. Source caller remains the PQS/product wrapper path and product entries in the all-pairs inventory metadata. |
| `_pqs_product_source_box_reference_blocks_from_pair_plan(...)` | `KEEP_ORACLE_ONLY` | Kept. Source caller remains `_pqs_product_source_box_reference_blocks(...)`. |
| `_pqs_pqs_source_box_reference_block(...)` | `DELETE_NOW` | Deleted. Exact caller audit showed only slow-test/docs references; the multi-term PQS/PQS oracle remains. |
| `_pqs_pqs_source_box_reference_blocks(...)` | `KEEP_ORACLE_ONLY` | Kept. Source caller remains its wrapper/from-pair-plan path and PQS/PQS entries in the all-pairs inventory metadata. |
| `_pqs_pqs_source_box_reference_blocks_from_pair_plan(...)` | `KEEP_ORACLE_ONLY` | Kept. Source caller remains `_pqs_pqs_source_box_reference_blocks(...)`. |
| `_pqs_pqs_product_source_box_shadow_blocks(...)` | `DELETE_NOW` | Deleted. After the route-shaped safe-term consumer deletion, exact caller audit showed no source caller and only slow-test/docs pressure. |
| `_pqs_pqs_product_source_box_all_pairs_inventory(...)` | `SHRINK_TEST_ONLY` / `KEEP_ORACLE_ONLY` | Kept. `_pqs_pqs_product_raw_box_route_producer(...)` still uses it for descriptor inventory validation. |
| `_pqs_pqs_product_raw_box_route_producer(...)` | `SHRINK_TEST_ONLY` | Kept. Source callers remain `_pqs_pqs_product_raw_box_route_from_geometry_facts(...)` and the density-density route producer. |
| `_pqs_pqs_product_raw_box_density_density_route_producer(...)` | `MOVE_LATER_TO_CPB` | Kept. Source callers remain the IDA-provenance adapter and component route smoke path. |
| `_pqs_pqs_product_nuclear_attraction_pair_block(...)` | `MOVE_LATER_TO_CPB` | Kept. Source caller remains `_pqs_pqs_product_route_shaped_nuclear_attraction_by_center(...)`. |
| `_pqs_pqs_product_route_shaped_nuclear_attraction_by_center(...)` | `MOVE_LATER_TO_CPB` | Kept. Source caller remains the component route smoke path. |

Batch actions:

- Deleted the single-term PQS/product and PQS/PQS source-box reference wrappers.
  Tests now use the surviving multi-term helpers where a narrow oracle remains.
- Deleted the three-unit safe-term shadow block helper and removed the slow-test
  helper/call sites that only preserved complete route-shaped safe-term
  matrices. The raw-box route producer still checks descriptor and inventory
  facts through `_pqs_pqs_product_source_box_all_pairs_inventory(...)`.
- Updated developer docs so the deleted shadow/helper families are described
  as retirement history, not present-tense provider or route contracts.

Sixth Coulomb-family route-shaped sweep:

| Helper family | Caller classification | Action |
| --- | --- | --- |
| `_pqs_pqs_product_raw_box_density_density_route_producer(...)` | `SHRINK_TEST_ONLY` / `MOVE_LATER_TO_CPB` | Kept. Source callers remain the IDA-provenance adapter and component route smoke path. Slow tests now keep only finite route-shaped matrix smoke and raw-vs-normalized agreement. |
| `_pqs_pqs_product_density_density_route_ranges(...)` | `KEEP_ORACLE_ONLY` | Kept. Source caller remains `_pqs_pqs_product_route_shaped_density_density_consumer(...)`. Tests no longer assert route-range implementation vocabulary directly. |
| `_pqs_pqs_product_density_density_all_pairs_inventory(...)` | `SHRINK_TEST_ONLY` / `MOVE_LATER_TO_CPB` | Kept. Source caller remains `_pqs_pqs_product_route_shaped_density_density_consumer(...)`. Slow tests no longer preserve all-pairs inventory object names or helper-family maps. |
| `_pqs_pqs_product_density_density_pair_block(...)` | `MOVE_LATER_TO_CPB` | Kept. Source caller remains the route-shaped density-density consumer. CPB-local electron-electron provider tests now own detailed local pair-factor numerical coverage. |
| `_pqs_pqs_product_route_shaped_density_density_consumer(...)` | `SHRINK_TEST_ONLY` / `MOVE_LATER_TO_CPB` | Kept. Source caller remains `_pqs_pqs_product_raw_box_density_density_route_producer(...)`. Slow tests now assert only retained shape, finite output, pair counts, symmetry, and raw-weighted equivalence. |
| `_pqs_pqs_product_nuclear_attraction_pair_block(...)` | `MOVE_LATER_TO_CPB` | Kept. Source caller remains `_pqs_pqs_product_route_shaped_nuclear_attraction_by_center(...)`. CPB-local electron-nuclear provider tests own detailed by-center local block coverage. |
| `_pqs_pqs_product_nuclear_attraction_all_pairs_inventory(...)` | `SHRINK_TEST_ONLY` / `MOVE_LATER_TO_CPB` | Kept. Source caller remains `_pqs_pqs_product_route_shaped_nuclear_attraction_by_center(...)`. Slow tests no longer assert all-pairs inventory helper-family vocabulary. |
| `_pqs_pqs_product_route_shaped_nuclear_attraction_by_center(...)` | `SHRINK_TEST_ONLY` / `MOVE_LATER_TO_CPB` | Kept. Source caller remains the component route smoke path. Slow tests now keep by-center separation, finite matrices, center labels/charges, and summed-center smoke only. |

Batch actions:

- No Coulomb-family source helper was deleted in this sweep because each
  candidate still has a source caller inside the remaining CCPM IDA provenance
  or component route smoke path.
- The slow PQS local-layer integration test no longer reconstructs route-shaped
  density-density or nuclear matrices from every old source-box pair helper.
  It also no longer asserts all-pairs inventory object names, helper maps, or
  repeated diagnostic nonclaims for these Coulomb-family paths.
- Detailed local Coulomb numerical correctness is intentionally owned by the
  CPB provider electron-electron and electron-nuclear local block tests. The
  CCPM route-shaped Coulomb helpers remain retirement-bound diagnostics/oracles
  until a real CPB route adaptation replaces their remaining source callers.
