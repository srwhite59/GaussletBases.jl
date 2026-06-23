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
- route/stage setup, raw-block setup, neutral raw-block kernels,
  residual/MWG/IDA changes, public workflow, and Cr2 facade/artifact work
  remain unapproved.

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
- the driver may call only approved base, staged, and supported supplemented
  producer surfaces and the approved artifact writer/readback;
- `basisname = nothing` selects base mode; `basisname !== nothing` selects
  supported supplemented diatomic mode and must reject `Natom == 1`;
- `padding` is physical box padding around the atom or two nuclei and maps
  internally to the existing public facade fields;
- route diagnostics, stop-after internals, ladder probes, raw-block switches,
  underscored package helper calls, status/report/payload fields, allocation
  probes, artifact schema dumps, solver work, public API/export changes,
  artifact schema changes, supplemented atoms, old route-stage choreography,
  and Cr2-specific workflow remain unapproved in the canonical driver.

Approved R1 one-center base atom relaxation:

- `HP-R1-ATOM-FN-01` relaxes `cartesian_base_hamiltonian(...)` from
  one-center H-only validation to explicit origin-centered all-electron
  one-center atoms;
- `HP-R1-ATOM-WIRE-01` maps the public nuclear charge to the existing
  White-Lindsey atomic mapping `Z` while keeping `d`, `core_spacing`, and
  `reference_spacing` separate;
- atoms and diatomics must share the same producer workflow after the narrow
  geometry/shellification differences; no atom-only Hamiltonian builder,
  materialization path, route-stage object, report/status payload, or artifact
  shape is approved;
- `HP-R1-ATOM-TEST-01` keeps the committed H endpoint as the regression gate
  and allows ignored/user-run Be or Cr atom artifact checks;
- no translated atoms, supplemented atoms, ECP, solver workflow, artifact
  schema change, element lookup/default table, public API redesign, or source
  file outside `src/cartesian_base_hamiltonian.jl` is approved.

Approved canonical driver atom workflow:

- `HP-DRV-ATOM-FN-01` approves explicit one-center atom input normalization in
  `bin/cartesian_ham_builder.jl` for `mode = :base`;
- `HP-DRV-ATOM-WIRE-01` lets the driver call the existing
  `cartesian_base_hamiltonian(system; basis, hamfile)` facade for
  origin-centered base atom construction where the base facade already
  supports the requested atom;
- current validation remains origin-centered H; this does not approve broader
  base atoms, translated atoms, supplemented atoms, element tables, ECP,
  solver workflow, public API/export changes, or artifact schema changes.

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
