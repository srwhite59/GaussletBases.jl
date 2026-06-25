# Current Cartesian Hamiltonian Producer Authority

This page is the compact live status page for future agents. It is not the
algorithm manual. Read it first to understand what is implemented and where the
current authority lives.

Normal startup reading:

- `README.md` for orientation;
- `current.md` for live status;
- `registry.md` for approved IDs, ownership, files, and function surfaces;
- `invariants.md` for architecture-wide guardrails;
- `residual_gaussian_domain_module.md` for the canonical Residual Gaussian
  algorithm contract;
- `residual_gaussian_orthogonality_robustness.md` for the narrow final
  residual identity-check robustness lane;
- `cartesian_gaussian_raw_blocks_nuclear.md` for the neutral uncharged nuclear
  raw-block owner;
- `cartesian_gaussian_raw_blocks_non_nuclear.md` for the neutral
  overlap/kinetic/moment raw-block owner;
- `r3_terminal_gg_product_matrices.md` for the narrow R3/RG terminal `G-G`
  product-matrix optimization lane;
- `r3_remaining_exact_operator_allocation_audit.md` for the measurement-only
  decision on remaining exact-operator allocation after terminal `G-G`
  workspace reuse;
- `r3_unit_nuclear_ugg_gaussian_sum.md` for the narrow terminal final-basis
  unit-nuclear `U_GG` Gaussian-sum allocation lane;
- `r1_one_center_base_atoms.md` for explicit origin-centered all-electron
  one-center base atoms beyond H;
- `cartesian_driver_usability_workflow.md` for the compact artifact-producing
  canonical driver lane;
- `cartesian_driver_atom_workflow.md` for explicit origin-centered base atom
  driver inputs;
- `r3_homonuclear_diatomic_supplemented_workflow.md` for the explicit
  homonuclear z-axis diatomic supplemented facade/driver relaxation;
- `white_lindsey_terminal_basis_realization.md` for the narrow terminal-basis
  seam required by the `nesting = :wl` construction family;
- `nesting_supplement_composition_plan.md` for the target 2 x 2 x 2
  composition matrix over geometry, nesting, and supplement state;
- `cartesian_hamiltonian_artifact_manifest.md` for compact Hamiltonian
  artifact sidecar groups and recipe provenance;
- `route_inventory_type_surface_cleanup.md` for the first route-inventory
  type-surface cleanup lane;
- `raw_product_source_mode_inventory_cleanup.md` for the raw product
  source-mode inventory cleanup lane;
- `contract_plan_vector_cleanup.md` for terminal-lowering and retained-unit
  transform contract-plan vector cleanup;
- `route_stage_type_surface_cleanup.md` for Be2 q5 compile-attributed
  route/stage compatibility-inventory cleanup;
- `route_stage_carrier_cleanup.md` for the post-cleanup route/stage carrier
  and plan-signature cleanup lane;
- `complete_core_shell_rhf_retirement.md` for the narrow stale
  complete-core-shell RHF payload-stack deletion lane;
- `route_driver_materialization_retirement.md` for retiring the old
  route-driver materialization/report/save wrapper workflow and stale tool/test
  pressure;
- `docs/src/developer/algorithm_implementation_index.md` for existing kernels
  and donor paths.

Historical design and review material remains under `history/`, `reviews/`,
and the R3 amendment pages. Those files are evidence and rationale, not normal
startup authority when they conflict with the compact current files.

## Live Status

The internal base PQS Hamiltonian lane is implemented for origin-centered H and
Cartesian z-axis H2. `HP-R1-ATOM-*` approves a narrow source relaxation for
explicit origin-centered all-electron one-center atoms beyond H through the
same base facade. This is internal numerical producer authority, not broad
public API polish.

Implemented base path:

- Slice A terminal basis realization over owned terminal supports;
- Slice B final-basis one-body assembly;
- Slice C localized IDA matrix assembly and existing
  `CartesianIDAHamiltonian` construction;
- Slice D base PQS materialization handoff;
- R1 public base facade for the approved H/H2 scope and fixed
  `producer_provenance/` artifact group.
- R1 one-center base atom relaxation for explicit origin-centered all-electron
  atoms through the existing base facade; current committed endpoint remains H.
- `HP-COMP-NS-*` approves public size-parameter normalization: `ns` is the
  durable user-facing cube/source/nesting size, while route-local `q` is
  derived as `q = ns` for `nesting = :pqs` or `q = ns - 2` for
  `nesting = :wl`.
- `HP-COMP-WLNS-*` approves the narrow WL z-axis diatomic `ns` contract
  cleanup: normalized `nesting = :wl`, `Natom = 2`, `ns < 4` should reject
  early in `src/cartesian_base_hamiltonian.jl`, and WL diatomic retained
  support may saturate across `ns` ranges when physical parent extent
  dominates.

Implemented Residual Gaussian path:

- Residual Gaussian basis construction, exact augmented one-body/moment
  operators, and residual MWG/IDA interaction live in
  `src/cartesian_residual_gaussians/`;
- compatibility entry points and artifact/facade hooks remain in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- the H2 owner-local residual-GTO/MWG endpoint has augmented dimension `489`
  and lowest-orbital IDA self-Coulomb `0.4574265214362075` within `1.0e-10`;
- Be2 remains an internal/performance-supported usability proxy, not a
  committed public gate;
- `HP-R3U-ZDI-*` relaxes the supported supplemented molecule scope to explicit
  homonuclear z-axis diatomics with no element-specific defaults or
  Cr2-specific branch;
- Cr2 may be used only as an explicit generic homonuclear z-axis ignored/user
  stress or usability run after H2/Be2 validation, not as a committed gate or
  special workflow.

Approved Residual Gaussian robustness lane:

- `HP-RG-ORTHO-FN-01` approves only final residual normalization/validation
  robustness in `src/cartesian_residual_gaussians/residual_basis.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for narrow internal keyword plumbing if needed;
- the target is strict N2 q5 p10 at `core_spacing = 0.042857`, where
  `G' S R` passes, owner metrics are full rank, retained counts are `9,9`, the
  merge metric is healthy, and only `R' S R - I` slightly exceeds the old
  absolute tolerance;
- the lane may use a symmetric final-overlap validation and a combined
  absolute/relative final identity check, but it must not change residual
  selection semantics, occupation cutoff, merge failure rule, MWG/IDA, raw
  blocks, artifacts, driver workflow, or public API.
- `HP-RG-IDTOL-FN-01` approves only changing the default final residual
  `R' S R` identity validation tolerance to `1.0e-8` in the same RG owner,
  with optional narrow compatibility keyword plumbing in
  `pqs_terminal_residual_gto.jl`. The default
  `residual_occupation_cutoff` remains `1.0e-8`; width/zeta filtering remains
  explicit and user-controlled; owner-local metric checks, final merge metric
  checks, and `G' S R` orthogonality checks remain active.
- The Be atom cc-pV5Z `lmax = 1` evidence for `HP-RG-IDTOL-*` is a tiny final
  identity overshoot: `21` retained residual directions, minimum occupation
  `6.151e-6`, final merge condition `1.0`, `max |G' S R| = 1.776e-14`, and
  `max |R' S R - I| = 2.183e-10` against an old allowed error of about
  `2.000e-10`.

Approved stale complete-core-shell RHF retirement:

- `HP-RETIRE-CCS-RHF-FN-01` approves only removing the include for
  `src/pqs_multilayer_complete_core_shell_rhf.jl` from `src/GaussletBases.jl`
  and deleting that RHF payload-stack file;
- the path is stale route-era workflow machinery with no live source/bin/test/
  tool caller found outside the file itself, while current CR2-facing work
  consumes canonical driver `CartesianIDAHamiltonian` artifacts;
- do not add replacements, adapters, compatibility wrappers, new status or
  payload objects, driver changes, artifact changes, route/shellification/
  terminal-lowering/raw-block/RG/MWG/IDA changes, or Cr2 workflow;
- do not change the older complete-core-shell H1/final-basis files or
  source-box materialization under this ID.

Approved route-driver materialization/report/save retirement:

- `HP-RETIRE-DRV-MAT-FN-01` approves only removing the old
  `cartesian_materialization`, `cartesian_print_summary`,
  `cartesian_print_details`, `cartesian_save`, and matching underscored
  route-driver materialization/report/save helpers;
- `HP-RETIRE-DRV-MAT-TOOL-01` approves only deleting or quarantining old tools
  that exist to drive that retired wrapper workflow;
- `HP-RETIRE-DRV-MAT-DOC-01` approves only active docs/index cleanup that stops
  presenting the wrapper workflow as canonical or active authority;
- `HP-RETIRE-DRV-MAT-TEST-01` approves only focused live-reference scans,
  canonical base/supplemented artifact smokes, unchanged H2 RG endpoint, and
  removal/update of stale docs-policy wrapper assertions;
- the canonical staged driver and current staged producer/artifact path must
  remain unchanged.

Approved neutral Cartesian Gaussian raw-block owner:

- `src/cartesian_gaussian_raw_blocks/` is approved for exact uncharged
  by-center Cartesian Gaussian nuclear `G-A` and `A-A` raw blocks;
- `HP-CGRB-FN-02` approves only a neutral-kernel reorganization around unique
  one-dimensional supplement axis families and term-first table reuse;
- `HP-CGRB-NN-*` approves only non-nuclear overlap, kinetic, coordinate
  moment, and second-moment `G-A`/`A-A` raw blocks;
- `src/cartesian_gaussian_axis_integrals.jl` remains optional support only for
  `HP-CGAI-FN-01` helper work needed by that family-reuse kernel;
- the owner may be consumed by Residual Gaussian and Qiu-White code after
  behavior-preserving parity;
- it does not own terminal projection, residual Gaussian transforms,
  Qiu-White route objects, final-basis `G-G` product-matrix optimization,
  caches, reports, artifacts, or public API.

Approved R3/RG terminal `G-G` product-matrix lane:

- `HP-R3GG-FN-01` approves only `K_GG`, coordinate moment `G-G`, and
  second-moment `G-G` product matrices used by
  `pqs_terminal_residual_gto_augmented_operators(...)`;
- owner module is `CartesianFinalBasisRealization`;
- approved files are
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` and,
  only if needed for small internal terminal-product workspace/helper reuse,
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`;
- it does not approve `G-A`/`A-A` raw-block changes, unit-nuclear Gaussian-sum
  work, IDA/MWG changes, route setup, public API, artifacts, metadata/status
  fields, persistent caches, or Cr2 workflow.

Approved remaining exact-operator allocation decision:

- after `954c86cd`, terminal `G-G` product workspace is considered crossed for
  the current Cr2 q4 proxy;
- `HP-R3REM-AUDIT-01` attributed the largest in-wrapper remaining allocation to
  unit-nuclear `U_GG` factor lookup plus Gaussian-sum construction;
- `HP-R3UN-FN-01` now approves only that narrow terminal final-basis
  unit-nuclear `U_GG` Gaussian-sum allocation lane under
  `CartesianFinalBasisRealization`;
- `HP-R3BASE-FN-01` approves same-construction reuse of already-built base
  final-basis `K_GG` and unit `U_GG[A]` blocks in supplemented exact augmented
  operators, with dimension and center-count validation and no metadata or
  provenance payload;
- `HP-R3BASE-DRV-WIRE-01` approves only the canonical driver supplemented-mode
  call-site wiring that passes `base_ham.kinetic` as `base_kinetic` and
  `base_ham.nuclear_attraction_unit_by_center` as `base_unit_nuclear`;
- the driver wiring must not change public inputs, hooks, timing labels,
  visible stage sequence, artifact schema, or driver contract;
- route/stage setup, raw-block setup, neutral raw-block kernels,
  residual/MWG/IDA changes, public workflow, and Cr2 facade/artifact work
  remain unapproved.

Approved compact Hamiltonian artifact manifest lane:

- `HP-HAM-MANIFEST-FN-01` approves only sidecar JLD2 groups for existing
  `CartesianIDAHamiltonian{Float64}` artifact files;
- approved source files are `src/cartesian_base_hamiltonian.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, and
  `src/cartesian_ida_hamiltonian.jl` only for a small unexported sidecar writer
  helper if needed;
- `hamiltonian_manifest/` follows the earlier PQS fixed-column/source-mode
  provenance model: basis identity is a status-bearing construction label, not
  a representative center;
- `hamiltonian_manifest/final_basis_labels/` records one row per matrix-order
  final basis column with sector, unit/source labels, shell/ray/radial status,
  representative center metadata, owner nucleus index where meaningful,
  locality/freezing labels, and supplement angular labels where available;
- optional `final_basis_source_relations/`, `source_shells/`, and
  `source_modes/` subgroups may be written only for native construction facts;
- `recipe_provenance/` records the validated public construction recipe and
  base/residual/augmented dimensions, including the public `nesting`
  construction family and a truthful route label derived from `(input.kind,
  input.nesting)`;
- `HP-NEST-ART-FN-01` approves only the narrow source cleanup that records
  `nesting` in base and recipe provenance, writes `:one_center_wl_base` rather
  than a PQS-oriented label for WL one-center artifacts, and rejects
  supplemented `nesting = :wl` before expensive base-stage construction;
- existing matrix keys and `read_cartesian_ida_hamiltonian` behavior must not
  change;
- no public reader API, driver public input change, route report/status
  payload, dense transform, raw inventory, inferred shell/ray labels, solver
  field, Cr2-specific field, or committed Cr2 fixture is approved.
- `HP-HAM-MANIFEST-SRC-FN-01` separately approves a compact construction-native
  source-mode provenance seam from terminal lowering / retained-unit /
  raw-product source plans to the base working basis manifest context, only for
  optional `source_shells/`, `source_modes/`, native source relations, and
  native final-label improvements.
- the source-mode seam must not serialize coefficients, dense transforms,
  `T_G`, `T_A`, raw inventories, route reports, allocation probes, diagnostic
  payloads, or repo-chosen ray/radial labels.

Approved canonical driver usability lane:

- `HP-DRV-FILE-01` approves only `bin/cartesian_ham_builder.jl`;
- `HP-DRV-FN-01` approves a compact functional driver workflow with visible
  defaults, one optional trusted project input file, command-line `key=value`
  overrides, visible public `system`/`basis`/optional `supplement` contract
  construction, coarse user-facing timing/summary printing, artifact write,
  optional readback check, and only the compact hooks `check_file`,
  `print_contract`, `print_timing`, and `expected_dimension`;
- `HP-DRV-STAGE-FN-01` approves a narrow non-exported, non-underscored staged
  producer surface so the driver can execute visible physics-level stages and
  operator-class timings: base working basis / terminal realization, product
  / moment operators, unit-nuclear attraction, electron-electron / IDA or
  residual-MWG interaction, base Hamiltonian assembly, Gaussian supplement,
  residual augmentation, and supplemented Hamiltonian assembly; this must be
  separate named stage functions, not one all-in-one replacement wrapper;
- approved source files for `HP-DRV-STAGE-FN-01` are
  `src/cartesian_base_hamiltonian.jl`,
  `src/pqs_source_box_low_order_materialization.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, only
  for behavior-preserving operator-class stage factoring;
- `HP-DRV-STAGE-WIRE-01` allows the canonical driver to call that staged
  surface as visible top-level stage calls and time/print those workflow stages
  without calling underscored package helpers;
- `HP-DRV-NEST-FN-01` and `HP-DRV-NEST-WIRE-01` approve one visible
  construction-family input, `nesting = :pqs` or `nesting = :wl`, in the
  canonical driver and narrow base-facade plumbing; `:pqs` maps to the existing
  PQS source-box route and remains the default, while `:wl` maps only to the
  existing White-Lindsey low-order route;
- `nesting` is not a diagnostic route switch and must not expose route
  skeletons, retained rules, raw-block switches, stop-after controls,
  diagnostics, route reports, or route-stage labels;
- supplemented `nesting = :wl` must be rejected clearly unless it is already
  valid through the existing supported supplemented facade/staged path; adding
  broad supplemented White-Lindsey behavior requires a separate amendment;
- artifact provenance must record the public `nesting` input truthfully under
  `HP-NEST-ART-FN-01`; a WL artifact must not be labeled as PQS merely because
  helper code remains PQS-named;
- the target producer contract is the three-choice composition recorded in
  `nesting_supplement_composition_plan.md`: geometry (`atom` or z-axis
  diatomic), nesting (`:pqs` or `:wl`), and supplement state (`off` or `on`);
  this is planning authority except for explicitly promoted cells, and
  unsupported cells must keep clear rejection until separately approved;
- the driver may call only approved base, staged, and supported supplemented
  producer surfaces and the approved artifact writer/readback;
- `basisname = nothing` selects base mode; `basisname !== nothing` selects a
  supported supplemented mode; `HP-COMP-SUPPATOM-FN-01` separately approves
  relaxing the old supplemented `Natom == 2` guard so `Natom == 1` and
  `Natom == 2` can use the same supplemented staged path;
- `padding` is physical box padding around the atom or two nuclei and maps
  internally to the existing public facade fields;
- route diagnostics, stop-after internals, ladder probes, raw-block switches,
  underscored package helper calls, status/report/payload fields, allocation
  probes, artifact schema dumps, solver work, public API/export changes,
  artifact schema changes, old route-stage choreography, and Cr2-specific
  workflow remain unapproved in the canonical driver.

Approved White-Lindsey terminal-basis seam:

- `HP-WLTERM-FN-01` approves only terminal-basis realization for the existing
  `:white_lindsey_low_order` route family, returning the same
  `CartesianTerminalBasisRealization` consumed by the staged Hamiltonian path;
- `HP-WLTERM-WIRE-01` approves only the route-helper seam that lets a WL route
  with sufficient native terminal lowering/retained records produce that
  terminal basis instead of being discarded by the PQS-only guard;
- approved files are `src/pqs_source_box_route_driver_helpers.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  optional
  `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
  and the include in
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`;
- this is not approval to adapt the old WL H1/H1+J materialization path,
  change route skeleton semantics, change shellification/retained-selection
  behavior, add diagnostics, change artifacts, or add supplemented WL behavior.

Approved first composition lane:

- `HP-COMP-WLDIAT-FN-01` approves the native WL z-axis diatomic base path for
  `Natom = 2`, `nesting = :wl`, `basisname = nothing`;
- the implementation must produce native WL diatomic terminal records and then
  use the same `CartesianTerminalBasisRealization`, staged base Hamiltonian
  construction, writer, and reader path as the PQS producer;
- `HP-WLDIAT-COMPACT-FN-01` separately approves the compact retained-basis
  correction for WL z-axis diatomics: the current elongated shared-shell
  boundary-stratum identity realization is classified as a mechanical endpoint,
  not the intended production compact WL retained basis;
- each WL diatomic unit must carry or realize compact retained columns from
  products of one-dimensional contractions, not full-support identity retention;
- identity realization remains valid only for true direct/core identity units,
  not for WL boundary-stratum retained units;
- deleted WL coefficient helpers are historical donor/reference material only:
  source work may re-express the essential CPB-local product-of-1D coefficient
  primitive behind the current terminal-basis boundary, but must not revive the
  route-global WL stack, reports, adapters, or H1/H1+J materialization;
- `HP-WLDIAT-PARITY-FN-01` separately approves the boundary-stratum retained
  count parity cleanup: odd-side enforcement belongs to nucleus-centered
  core/contact blocks, not to boundary shell strata; public WL `ns = 4`
  therefore targets `4^3 - 2^3 = 56` boundary columns rather than the inherited
  symmetric-odd donor result `26`;
- approved source files are the narrow diatomic complete-core-shell,
  terminal-shellification, terminal-lowering, route-helper, terminal-basis,
  final-basis include, and base-facade surfaces named in `registry.md`;
- `:z_axis_diatomic_wl_base` is approved only as a truthful route-provenance
  value under existing artifact keys, not as an artifact schema change;
- no driver special cases, fake compactness by dropping rows, old WL H1/H1+J
  materialization, RG/MWG/supplement work, route diagnostics, public
  API/export changes, committed tests, or Cr2 workflow is approved.
- `HP-COMP-BASEDIAT-FN-01` separately approves only
  `src/cartesian_base_hamiltonian.jl` to relax base two-center validation from
  H2-only to explicit homonuclear z-axis all-electron diatomics;
- the base diatomic relaxation requires equal symbols, equal finite positive
  integer-valued charges, two finite distinct z-axis centers, and neutral
  `nup + ndn == sum(charges)`;
- it preserves both `nesting = :pqs` and `nesting = :wl`, but it does not
  approve route/shellification/terminal-lowering changes, supplement/RG/MWG
  work, driver changes, element tables, heteronuclear support, translated or
  non-z-axis geometry, committed tests, or Cr2 workflow.
- `HP-COMP-SUPPWL-FN-01` approves the supplemented White-Lindsey z-axis
  diatomic composition lane for `Natom = 2`, `basisname !== nothing`, and
  `nesting = :wl`;
- approved source surface is `src/cartesian_base_hamiltonian.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  optional only for a direct RG/MWG compatibility genericity blocker;
- the implementation may remove the two early supplemented-WL blockers only if
  the existing Residual Gaussian/MWG path works with the WL
  `CartesianTerminalBasisRealization`;
- it must preserve the supplement contract, residual selection, exact
  augmented operators, MWG/IDA convention, base K/U reuse, artifact keys,
  manifest/provenance, driver inputs, and stage labels;
- it does not approve driver changes, supplemented atoms, route
  skeleton/shellification/terminal lowering changes, raw-block changes,
  residual-selection changes, MWG/IDA convention changes, artifact schema or
  reader changes, public API/export changes, old WL H1/H1+J materialization,
  committed tests, or Cr2 workflow.
- `HP-COMP-SUPPATOM-FN-01` approves the supplemented one-center atom
  composition lane for `Natom = 1`, `basisname !== nothing`, and
  `nesting = :pqs` or `nesting = :wl`;
- approved source surface is `src/cartesian_base_hamiltonian.jl` and
  `bin/cartesian_ham_builder.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  optional only for a direct one-owner RG/MWG genericity blocker;
- the implementation must use the existing base atom validation, terminal
  basis construction, residual Gaussian augmentation, exact augmented
  operators, residual MWG/IDA interaction, base K/U reuse, assembly, writer,
  readback, manifest, and provenance;
- it may select `legacy_atomic_gaussian_supplement(...)` for one-center inputs
  and keep the existing diatomic supplement loader for two-center inputs;
- it may relax only the canonical driver's supplemented `Natom == 2` guard and
  must otherwise preserve driver public inputs, ordering, hooks, spacing,
  stage labels, and artifact contract;
- it does not approve a separate atom-only Hamiltonian builder, new driver
  inputs, route switches, diagnostics, stop-after controls, route
  skeleton/shellification/terminal lowering changes, raw-block changes,
  residual-selection changes, MWG/IDA convention changes, artifact schema or
  reader changes, public API/export changes, solver/ECP work,
  heteronuclear/general geometry, translated atoms, committed tests, or Cr2
  workflow.

Composition lane status:

- the explicit `atom | z-axis diatomic`, `:pqs | :wl`, and
  `supplement = off | on` composition lanes now all have approved
  implementation authority under the current origin-centered atom and
  homonuclear z-axis diatomic geometry constraints;
- `HP-COMP-ATOMBOX-FN-01` separately approves one-center atom parent sizing in
  `src/cartesian_base_hamiltonian.jl`: public `basis.radius` is the physical
  box extent authority, parent axis counts must depend on radius plus
  `core_spacing` / existing spacing policy, and the fixed `2*q + 1` atom
  parent side count artifact must be removed;
- `HP-COMP-NS-FN-01` separately approves public size naming: `ns` is the
  requested source/nesting size, and route-local `q` is derived from `ns` and
  `nesting`;
- this atom box lane preserves origin-centered atom validation, explicit
  charge/electron-count validation, `nesting = :pqs` and `nesting = :wl`,
  supplemented atoms, artifact keys, manifest/provenance, and canonical driver
  inputs;
- it does not approve driver changes, route-family switches, raw-block
  changes, residual-selection changes, MWG/IDA convention changes, artifact
  schema or reader changes, public API/export changes, solver/ECP work,
  diagnostics/status/report payloads, committed tests, Cr2-specific workflow,
  translated atoms, non-origin atom support, element lookup/default tables,
  broad parent-construction rewrites, or diatomic sizing changes;
- deferred geometry, solver, ECP, public export, and Cr2-specific work still
  need later docs-only amendments before implementation may begin.

Approved R1 one-center base atom relaxation:

- `HP-R1-ATOM-FN-01` relaxes `cartesian_base_hamiltonian(...)` from
  one-center H-only validation to explicit origin-centered all-electron
  one-center atoms;
- `HP-R1-ATOM-WIRE-01` maps the public nuclear charge to the existing
  White-Lindsey atomic mapping `Z` and maps resolved public `core_spacing` to
  the private `parent_mapping_d`; public `d` is deprecated and, if
  temporarily accepted, must equal resolved `core_spacing`;
- `HP-R1-CORE-FN-01` freezes `core_spacing` as the single public
  near-nucleus physical scale while keeping `reference_spacing`, tail spacing,
  and box padding separate; White-Lindsey `Z` behavior is a mapping-shape rule
  or preset rule, not a second public knob; visible driver/project defaults
  such as `core_spacing = 0.3` are explicit resolved inputs and may be
  overridden for quick tests; routine correctness-test scalars must be tied to
  their exact override inputs, not described as physics-default results;
- atoms and diatomics must share the same producer workflow after the narrow
  geometry/shellification differences; no atom-only Hamiltonian builder,
  materialization path, route-stage object, report/status payload, or artifact
  shape is approved;
- `HP-R1-ATOM-TEST-01` keeps the committed H endpoint as the regression gate
  and allows ignored/user-run Be or Cr atom artifact checks;
- the R1 atom IDs do not approve translated atoms, ECP, solver workflow,
  artifact schema change, element lookup/default table, public API redesign,
  or source files outside `src/cartesian_base_hamiltonian.jl`;
- supplemented atoms are separately governed by `HP-COMP-SUPPATOM-*`.

Approved route-recipe cleanup:

- `HP-ROUTE-RECIPE-FN-01` allows `cartesian_recipe(...)` in
  `src/pqs_source_box_route_driver_helpers.jl` to construct only the
  `route_family`-selected subrecipe and to set the inactive subrecipe to
  `nothing` where needed for compatibility;
- `:pqs_source_box` route inputs must no longer require inactive
  `white_lindsey_*` fields, and `_cartesian_base_route(kind)` in
  `src/cartesian_base_hamiltonian.jl` may delete those unused fields;
- explicit `:white_lindsey_low_order` route support remains live and must keep
  using its WL route fields where selected;
- no canonical-driver changes, numerical kernel changes, shellification or
  terminal-lowering changes, materialization/artifact schema changes,
  report/status expansion, WL materialization deletion, new tests, or Cr2 run
  are approved by this cleanup.

Approved route-inventory type-surface cleanup:

- `HP-ROUTE-INV-FN-01` approves only
  `src/pqs_source_box_route_driver_helpers.jl`;
- the target is removal of runtime-keyed retained-unit inventory
  `NamedTuple{unit_keys}` shapes and runtime-keyed `pair_family_counts =
  NamedTuple{families}(...)`;
- approved replacements are vector-backed records/tables, stable dictionaries,
  or helper accessors with stable concrete types;
- public input `NamedTuple`s, fixed `NTuple{3,Int}` coordinates/dimensions,
  artifact sidecar tables, `RawProductBoxPlan.source_mode_indices`, terminal
  lowering contract tuples, and retained-unit transform-contract tuples remain
  out of scope.

Approved raw product source-mode inventory cleanup:

- `HP-RAW-SRCMODE-FN-01` approves vector-backed source-mode inventory cleanup
  in `src/cartesian_raw_product_sources/records.jl`,
  `src/cartesian_raw_product_sources/source_mode_indices.jl`, and
  `src/cartesian_raw_product_sources/summaries.jl`;
- narrow consumer wiring is approved only in
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  and `src/cartesian_base_hamiltonian.jl`;
- the target is removal of `RawProductBoxPlan.source_mode_indices` and
  `source_mode_column_indices` variable-length tuple storage, while preserving
  deterministic source-mode ordering, fixed `NTuple{3,Int}` mode coordinates,
  retained-rule behavior, and manifest source-mode/relation facts;
- accessor semantics mean the same ordered facts and column associations, not
  preserving a variable-length `Tuple` concrete return type;
- terminal-lowering contract tuple cleanup and broader source-box/pair-block
  rewrites remain out of scope.

Approved contract-plan vector cleanup:

- `HP-CONTRACT-VEC-FN-01` approves vector-backed plan-level contract
  inventories in `src/cartesian_terminal_lowering/contracts.jl`,
  `src/cartesian_terminal_lowering/selection.jl`,
  `src/cartesian_terminal_lowering/summaries.jl`,
  `src/cartesian_retained_unit_transform_contracts/records.jl`,
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`, and
  `src/cartesian_retained_unit_transform_contracts/summaries.jl`;
- narrow consumer wiring is approved only if needed in
  `src/pqs_source_box_route_driver_helpers.jl`,
  `src/cartesian_base_hamiltonian.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`;
- the target is removal of variable-length tuple storage for
  `TerminalLoweringPlan.available_contracts`, `TerminalLoweringPlan.contracts`,
  and `RetainedUnitTransformContractPlan.contracts`;
- accessors `available_contracts(plan)`, `selected_contracts(plan)`,
  `contracts(plan)`, and `transform_contracts(plan)` must preserve ordered
  behavior and semantics without preserving variable-length `Tuple` concrete
  field types;
- `source_cpbs::Tuple{Vararg{CoordinateProductBox}}` remains out of scope.

Approved route/stage type-surface cleanup:

- `HP-ROUTE-STAGE-TYPE-FN-01` approves only
  `src/pqs_source_box_route_driver_helpers.jl` and
  `src/cartesian_terminal_shellification_geometry.jl`;
- approved targets are Be2 q5 compile-attributed compatibility inventories and
  wide route/stage return surfaces around
  `_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`,
  `cartesian_units`,
  `_pqs_source_box_route_driver_transform_stage_low_order_summary`,
  `cartesian_transforms`, and
  `_cartesian_terminal_shellification_region_unit_inventory`;
- stale compatibility inventories may be deleted, and remaining runtime-sized
  `NamedTuple` / `Tuple` carriers may be replaced with vector-backed compact
  internal objects, stable dictionaries, accessors, or smaller summaries;
- H2 base/supplemented artifact behavior, terminal shellification/lowering
  order, artifact/manifest schema, public driver contract, route semantics, and
  numerical matrices must remain unchanged;
- driver changes, artifact/manifest changes, numerical/raw-block/RG/MWG/IDA
  changes, route diagnostic/status/report expansion, committed tests,
  PackageCompiler/PrecompileTools/sysimage work, and Cr2 workflow remain
  unapproved.

Approved route/stage carrier cleanup:

- `HP-ROUTE-STAGE-CARRIER-FN-01` approves only
  `src/pqs_source_box_route_driver_helpers.jl` and
  `src/pqs_source_box_diatomic_complete_core_shell.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  optional only where directly required to slim terminal realization plan
  carriers in the approved path;
- approved targets are broad stage signatures and plan carriers around
  `cartesian_shells`, `cartesian_units`, `cartesian_transforms`, terminal
  topology support-region planning, terminal retained-rule planning, and
  directly required terminal realization plan carriers;
- giant shellification, route-skeleton, support-plan, retained-rule-plan, and
  terminal-plan `NamedTuple` / tuple shapes may stop crossing stage function
  signatures; necessary facts should move to compact typed/vector-backed
  carriers, smaller summaries, accessors, or local recomputation from canonical
  objects;
- route skeleton construction semantics and
  `src/pqs_source_box_route_driver_skeletons.jl` remain out of scope;
- H2 base/supplemented artifact behavior, H2 R3 endpoint if touched,
  terminal support/shellification/lowering order, public driver contract,
  artifact/manifest schema, route semantics, and numerical matrices must remain
  unchanged;
- driver changes, artifact/manifest changes, numerical/raw-block/RG/MWG/IDA
  changes, route diagnostic/status/report expansion, committed tests,
  PackageCompiler/PrecompileTools/sysimage work, and Cr2 workflow remain
  unapproved.

Approved canonical driver atom workflow:

- `HP-DRV-ATOM-FN-01` approves explicit one-center atom input normalization in
  `bin/cartesian_ham_builder.jl` for `mode = :base`;
- `HP-DRV-ATOM-WIRE-01` lets the driver call the existing
  `cartesian_base_hamiltonian(system; basis, hamfile)` facade for
  origin-centered base atom construction where the base facade already
  supports the requested atom;
- current validation remains origin-centered H; these driver atom IDs do not
  approve broader base atoms, translated atoms, element tables, ECP, solver
  workflow, public API/export changes, or artifact schema changes;
- supplemented atom driver wiring is separately governed by
  `HP-COMP-SUPPATOM-*`.
- `HP-DRV-ATOM-CLEAN-01` approves only removing the stale hidden
  `d = core_spacing` atom-basis field from the canonical driver; public
  inputs, defaults, overrides, hooks, timing labels, stage sequence, artifact
  schema, and driver contract must remain unchanged.

Approved homonuclear z-axis diatomic supplemented workflow:

- `HP-R3U-ZDI-FN-01` relaxes the non-exported supplemented facade from
  hardcoded H2/Be2 checks to explicit homonuclear two-center z-axis diatomic
  validation;
- `HP-R3U-ZDI-WIRE-01` lets canonical driver supplemented mode call the
  supported facade instead of package-internal helper paths;
- required inputs remain explicit: symbols, charges, spin sectors, geometry,
  base basis controls, supplement labels, and optional trusted `basisfile`;
- no heteronuclear, non-z-axis, ECP, solver, public export, artifact schema,
  route diagnostic, or Cr2-specific branch is approved.

## Residual Gaussian Authority

`residual_gaussian_domain_module.md` is the canonical current RG algorithm
contract. Do not duplicate the full residual-selection, exact-transform, or
MWG-interaction algorithm in this file.

Essential live guardrails:

- residual Gaussian directions are selected separately on each atom and then
  merged once;
- residual occupation is a physical retained-direction measure, not numerical
  rank;
- exact augmented one-body/moment transforms are not the MWG approximation;
- MWG descriptors are rotation-dependent and must be computed from the final
  merged residual basis;
- RG does not own artifact writing, artifact provenance, basis loading, facade
  parsing, parent lattice construction, terminal topology, raw Gaussian block
  formula ownership, or public exports.
- final residual identity validation may use the approved
  `HP-RG-ORTHO-FN-01` absolute/relative check only for small floating-point
  overshoots after owner-local selection and a healthy final merge; it is not
  permission to relax residual selection or retain low-occupation directions.
- `HP-RG-IDTOL-FN-01` sets the default final residual identity validation
  tolerance to `1.0e-8`. This is a final validation/cleanup tolerance only; it
  does not change the residual occupation cutoff, width/zeta filtering
  defaults, owner grouping, merge failure rules, or MWG/IDA conventions.

Exact Cartesian Gaussian raw blocks are separate neutral kernel authority:
uncharged by-center nuclear blocks under
`cartesian_gaussian_raw_blocks_nuclear.md`, and non-nuclear
overlap/kinetic/moment blocks under
`cartesian_gaussian_raw_blocks_non_nuclear.md`. RG consumes those blocks as
exact operator inputs; RG does not own their raw construction.

The old R3-A/R3-B/R3-C labels are implementation-history labels. The current
domain names are:

- `build_residual_gaussian_basis`;
- `transform_augmented_operator`;
- `moment_matched_gaussians`;
- `assemble_residual_ida_interaction`.

## Current Boundaries

Base pair/assembly role:

- Cartesian pair terms remain a structural precursor for the base route;
- `cartesian_assembly` remains topology/report context only;
- production Hamiltonian construction must use realized terminal basis and
  final-basis matrix kernels, not source-plan/report reconstruction.

Artifact ownership:

- base Hamiltonian files use the existing `CartesianIDAHamiltonian` artifact
  shape plus the fixed `producer_provenance/` group;
- supplemented R3/RG writes may add only the approved compact
  `supplement_provenance/` group through the existing helper outside the RG
  module;
- `read_cartesian_ida_hamiltonian` ignores provenance and reads matrices
  normally.

Deferred lanes:

- broad driver diagnostics and public-driver polish beyond the approved
  compact artifact-producing workflow;
- supplemented atom Hamiltonians and translated atoms;
- Cr2-specific workflow, committed Cr2 gates, and Cr2 support decisions beyond
  the generic explicit homonuclear z-axis path;
- non-base/supplement public workflow;
- ECP/EGOI/RHF/solver/HamV6 work;
- artifact/public API decisions beyond the approved compact provenance groups.

Before implementation, confirm the approved ID and source surface in
`registry.md`; if the needed surface is absent, do a docs-only amendment first.
