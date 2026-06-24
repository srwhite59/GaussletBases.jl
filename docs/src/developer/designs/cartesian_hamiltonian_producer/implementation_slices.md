# Implementation Slices

## Slice A — Terminal Basis Realization

Status: implemented.

Implemented boundary:

- `CartesianTerminalBasisBlock`
- `CartesianTerminalBasisRealization`
- `pqs_terminal_basis_realization(...)`
- terminal-basis wiring from `cartesian_transforms`

Corrected construction contract:

```text
full source-box product modes
-> boundary product-mode column selection
-> restrict rows to support.support_indices / support.support_states
-> shell-local Gram on that owned support
-> symmetric shell-local Lowdin
-> final sign canonicalization
-> append block with unchanged owned support
```

Previous-block projection, recursive projection, and effective-support growth
onto previous terminal regions are rejected. Cross-block overlap is structurally
zero because parent gausslet rows are orthonormal and terminal regions own
disjoint parent rows. A nonzero structural overlap means duplicated support
rows, incorrect row restriction, wrong support ownership, or an indexing error;
it is not a physical residual to compute or repair.

Cross-block kinetic, nuclear-attraction, and IDA interactions may still be
nonzero and remain assembled over terminal block pairs.

Validation gates used:

- one-center atomic terminal basis through the shared entry point;
- H2 terminal basis dimension `471`;
- Cr2 terminal basis dimension `4291` during design/implementation validation;
- positive final integrals;
- structurally disjoint terminal supports and shell-local identity overlaps;
- no global Lowdin.

Deletion obligations completed:

- terminal source-realization preflight path removed;
- source-plan summary mirrors removed;
- old blocked-route smoke assertions reduced.

Deferred gates:

- larger performance/stress work is not part of Slice A.

Implementation correction required:

- remove `_subtract_previous` and recursive projection;
- make `_shell_seed` use `support.support_indices` /
  `support.support_states`, not all `outer_box` rows;
- remove `projection_atol` plumbing when it has no remaining construction use;
- replace production cross-overlap audit / `max_cross_overlap` with structural
  checks: every block support equals its authoritative terminal support,
  terminal support sets are pairwise disjoint, and each shell-local overlap is
  identity;
- H2 realized block supports should return to local counts such as `(275, 362,
  578)`, not cumulative supports like `(275, 637, 1215)`;
- rerun Slice A/B/C and R1 validation after source correction.

The recursive-projection Be2 measurements are superseded. Post-`d2bf139c` Be2
measurements are the valid optimization baseline.

## Slice B — Final-Basis One-Body Operators

Status: implemented.

Implemented boundary:

- `assemble_terminal_product_operator!(...)`
- file-local term-first Gaussian-sum nuclear attraction helper

Validation gates used:

- one-center H one-body baseline around `-0.49855234726272035`;
- H2 one-body lowest energy `-0.79460371733658908`;
- light separated-diatomic N2 validation for topology/performance;
- finite/symmetric `K` and `U_A` matrices;
- 64 MiB local workspace cap.

Deletion obligations completed:

- dense direct identity allocation removed from terminal overlap helper;
- no production global support one-body operator was introduced.

Deferred gates:

- Cr2 one-body stress/performance validation.

## Slice C1 — Localized IDA Matrix

Status: implemented.

Implemented boundary:

- `assemble_terminal_ida_interaction!(...)`

Validation gates used:

- H2 one-body energy remains at the Slice B value;
- final IDA weights positive and finite;
- `electron_electron_ida` finite and symmetric;
- reviewed H2 self-Coulomb `0.4569117646737212`.

Deletion obligations completed:

- no CPBM route revival;
- no IDA payload/cache/status framework.

Deferred gates:

- non-base/supplement IDA variants.

## Slice C2 — CartesianIDAHamiltonian Construction

Status: implemented.

Implemented boundary:

- existing `CartesianIDAHamiltonian` construction from `K`, uncharged
  by-center `U_A`, localized IDA `V`, charges, positions, and electron counts.

Validation gates used:

- H2 final dimension `471`;
- constructed type `CartesianIDAHamiltonian{Float64}`;
- `one_body_hamiltonian(ham)` reproduces direct Slice B H1;
- H2 one-body lowest `-0.79460371733658908`;
- self-Coulomb from `ham.electron_electron_ida`
  `0.4569117646737212`.

Deferred gates:

- public workflow polish and larger molecule stress.

## Slice D — Base Materialization Handoff

Status: HP-WIRE-02 approved and implemented.

Implemented boundary:

```julia
cartesian_materialization(report, terminal_basis_realization, materialization_inputs)
```

Validation gates used:

- no-request path returns `nothing`;
- requested H2 path returns `CartesianIDAHamiltonian{Float64}`;
- final dimension `471`;
- H2 one-body lowest `-0.79460371733658908`;
- H2 self-Coulomb `0.4569117646737212`;
- optional artifact write/readback with one-body delta `0.0`.

Deletion obligations completed:

- materialization wrapper/result-kind path removed for the PQS base handoff;
- blocked physical-gausslet source-plan/report mirrors removed;
- physical-gausslet target/supplement payload chain removed;
- H2 and Cr2 stage probes no longer treat blocked source-plan summaries as
  route authority.

Deferred gates:

- public-driver polish;
- public export or driver workflow beyond the approved non-exported R3
  usability facade;
- Residual Gaussian migration cleanup: residual basis construction, exact
  augmented operators, and residual MWG interaction now live under
  `src/cartesian_residual_gaussians/`; keep deleting old R3 wrappers when live
  callers move;
- compact supplemented artifact writing and facade parsing remain outside
  `CartesianResidualGaussians` as terminal/facade workflow glue unless a later
  amendment identifies a real duplication or consumer need;
- Cr2 stress/performance and any full Cr2 supplemented Hamiltonian run;
- measurement-only Cr2-readiness forecasting, consumer workflow beyond the
  approved H2/Be2 R3 usability facade, and basis/supplement-realism decisions
  until separately approved;
- non-base Hamiltonian variants.

## Residual Gaussian Module Migration

Status: residual-basis construction, exact augmented operator transformation,
and matched-width Gaussian residual interaction are migrated to the
`CartesianResidualGaussians` domain module.

Canonical algorithm authority:

- `residual_gaussian_domain_module.md`

Current compatibility boundary:

- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` remains
  the terminal/facade compatibility file for existing R3 entry points, artifact
  writing, and facade wiring;
- moved physics helpers should not be reintroduced there;
- RG does not own `supplement_provenance/`, JLD2 artifact workflow, facade
  input parsing, basis loading, parent lattice construction, or public exports.

Completed cleanup:

- residual-basis helpers moved to `residual_basis.jl`;
- exact `[G,A] -> [G,R]` operator transform moved to
  `augmented_operators.jl`;
- residual MWG descriptor and interaction assembly moved to
  `mwg_interaction.jl`;
- standalone R3 test no longer depends on old test-only R3-B wrapper names.

Remaining cleanup:

- delete any future compatibility wrapper once the exact live caller moves;
- keep artifact and facade hooks outside RG unless a separate design amendment
  approves moving them.

## Compact Hamiltonian Artifact Manifest

Status: approved for implementation under `HP-HAM-MANIFEST-FN-01` and
`HP-HAM-MANIFEST-TEST-01`. The source-mode provenance seam is approved
separately under `HP-HAM-MANIFEST-SRC-FN-01` and
`HP-HAM-MANIFEST-SRC-TEST-01`.

Approved boundary:

- JLD2 sidecar groups for existing `CartesianIDAHamiltonian{Float64}` artifact
  files;
- source file `src/cartesian_base_hamiltonian.jl`;
- source file `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- source file `src/cartesian_ida_hamiltonian.jl` only for a small unexported
  sidecar writer/helper if needed;
- no changes to Hamiltonian matrix keys or `read_cartesian_ida_hamiltonian`
  behavior.

Approved sidecar groups:

- `hamiltonian_manifest/` for source/final-column provenance modeled on the
  prior PQS fixed-column sidecar contract;
- `hamiltonian_manifest/final_basis_labels/` for exact matrix-order final
  basis row labels, sectors, unit/source labels, shell/ray/radial status,
  representative center metadata, owner nucleus indices, locality/freezing
  labels, supplement labels, and angular powers where available;
- optional `hamiltonian_manifest/final_basis_source_relations/`,
  `hamiltonian_manifest/source_shells/`, and
  `hamiltonian_manifest/source_modes/` groups only for construction-native
  relation/source-mode facts;
- `recipe_provenance/` for the validated public recipe: `q`, `core_spacing`,
  padding-derived extents, route, parent axis counts, atom symbols/charges/
  locations, `nup`/`ndn`, supplement label/file/options, and
  base/residual/augmented dimensions.

Required implementation rule:

- basis-function identity is a construction label with status, not a center;
- centers are representative metadata with explicit definition/status;
- derive center conventions only from existing terminal basis blocks, parent
  axes, residual metadata, and augmented moment/MWG descriptors;
- do not infer source-box, shell, ray, radial, or relation labels from centers,
  nearest-grid snapping, support order, support indices, or raw-to-final
  support;
- use explicit `:unavailable` or `:mixed` status when a native construction
  label is not available;
- stop and report if this cannot be done without adding algorithmic metadata;
- do not serialize `T_G`, `T_A`, dense transforms, coefficients, raw
  inventories, dense moments, allocation probes, route reports, or status
  payloads.

Source-mode provenance seam:

- carry one compact construction-native provenance object from terminal
  lowering / retained-unit / raw-product source planning to the base working
  basis manifest context;
- approved source files are terminal lowering contracts/region contracts, raw
  product source records/source-mode ordering, retained-unit records/lowering,
  retained-unit transform-contract records/unit contracts, terminal basis
  realization, and `src/cartesian_base_hamiltonian.jl`;
- the preferred carrier is a `source_mode_provenance` field on the internal
  `cartesian_base_working_basis(...)` result;
- adding one optional source-mode provenance field to
  `CartesianTerminalBasisRealization` is allowed only if it avoids duplication
  or loss of terminal construction ordering;
- populate optional `source_shells/`, `source_modes/`, and native
  `final_basis_source_relations/` only from construction-native facts;
- improve `final_basis_labels/` only where the label is native and directly
  tied to the row's unit/source mode;
- do not add ray/cone/radial labels unless already natively defined.

Padding scope:

- this lane records current recipe provenance only;
- one-center atom padding remains a separate contract clarification and must
  not be bundled into this source pass;
- diatomic padding-derived extents remain active through the existing facade.

Validation gates:

- `git diff --check`;
- package load;
- H atom or H2 base artifact write/readback;
- H2 supplemented artifact write/readback;
- direct JLD2 checks for approved final-basis label keys, dimensions,
  status-bearing unavailable/mixed labels, no inferred labels, and recipe
  values;
- optional Be2 supplemented artifact manifest inspection if practical;
- no Cr2 run.

Forbidden:

- driver public input changes, public reader API/export, matrix-key changes,
  artifact schema dumps in the driver, solver-specific fields,
  CR2-consumer-specific fields, Cr2-specific fields, committed Cr2 fixtures,
  Cr2-specific branches, route reports/status payloads, persistent caches, new
  algorithmic metadata, coefficients, dense transforms, `T_G`, `T_A`, raw
  inventories, diagnostic payloads, or source files outside the approved
  surfaces.

Line budget:

- at most `150` added `src` lines;
- for `HP-HAM-MANIFEST-SRC-FN-01`, at most `180` added `src` lines;
- no new committed test, tool, or input-fixture file.

## Cartesian Gaussian Raw Blocks - Nuclear Slice

Status: approved for implementation; not part of R3/RG public workflow.

Approved boundary:

- neutral owner files under `src/cartesian_gaussian_raw_blocks/`;
- exact uncharged by-center Cartesian Gaussian nuclear raw blocks:
  parent-supplement `G-A` and supplement-supplement `A-A`;
- analytic one-dimensional nuclear factor construction, unique nuclear
  coordinate reuse, upper-triangular `A-A` assembly with mirroring,
  function-local scratch reuse, and term-first contraction.

Implementation sequence:

1. Extract the current Residual Gaussian/Qiu-White nuclear `G-A`/`A-A`
   behavior into the neutral owner without changing numerical conventions.
2. Rewire Residual Gaussian and Qiu-White callers to consume the neutral
   kernel.
3. Delete the duplicate route-local nuclear loops once parity is established.
4. Optimize allocation inside the neutral owner only after extraction parity.

Validation gates:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged as an ignored performance/usability
  measurement if needed;
- Cr2 q4 exact nuclear blocks match the current implementation at roundoff as
  an ignored measurement only;
- a small Qiu-White nuclear parity fixture passes.

Forbidden in this slice:

- overlap, kinetic, coordinate moments, second moments;
- pair factors or matched-width Gaussian interaction;
- terminal projection or residual-basis transformation;
- Qiu-White route objects or parent construction;
- persistent caches, metadata, reports, status fields, artifacts, public API,
  Cr2 facade support, or Cr2 artifact workflow.

Follow-on low-level optimization approval:

- `HP-CGRB-FN-02` approves only
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl` for reorganizing the
  neutral nuclear kernel around unique one-dimensional supplement axis
  families, integer orbital-to-family maps, unique `G-A`/`A-A` table keys, and
  term-first table reuse;
- `HP-CGAI-FN-01` is superseded as a performance endpoint and remains only an
  optional helper surface in `src/cartesian_gaussian_axis_integrals.jl` for a
  specialized nonallocating nuclear-factor scalar integral or tiny table-fill
  helper if needed by `HP-CGRB-FN-02`;
- the bottleneck addressed is repeated rebuilding of equivalent
  one-dimensional Gaussian axis families through flattened 3D orbital loops,
  plus allocation inside scalar nuclear factor integrals;
- independent contraction of x/y/z axis tables into separate scalar
  contractions is forbidden; the kernel must preserve
  `sum_pq c_p c_q Ix[p,q] Iy[p,q] Iz[p,q]`;
- the later source pass must preserve H2/Be2/Qiu-White/Cr2 nuclear parity,
  report unique family/table/scalar-call count reductions, and substantially
  reduce Cr2 q4 nuclear raw-block allocation from the `44552.840 MiB`
  baseline without adding caches, metadata, route objects, reports, artifacts,
  public API, or semantic changes.

## Cartesian Gaussian Raw Blocks - Non-Nuclear Slice

Status: approved for implementation under `HP-CGRB-NN-*`; not part of
`HP-CGRB-FN-02`.

Approved boundary:

- neutral owner file `src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`;
- exact non-nuclear Cartesian Gaussian raw blocks:
  overlap, kinetic, coordinate moments `x`/`y`/`z`, and second moments
  `x^2`/`y^2`/`z^2`;
- parent-supplement `G-A` and supplement-supplement `A-A` only;
- analytic one-dimensional table construction, unique supplement axis-family
  reuse, canonical `A-A` family-pair table keys, upper-triangular `A-A`
  assembly with mirroring, function-local scratch reuse, and coupled
  product-axis contraction.

Implementation sequence:

1. Extract the current Qiu-White non-nuclear `G-A`/`A-A` behavior into the
   neutral owner without changing numerical conventions.
2. Rewire Residual Gaussian exact-operator construction and residual mixed
   overlap setup to consume the neutral output where they currently use the
   Qiu-White donor organization.
3. Rewire Qiu-White consumers to the neutral output.
4. Delete duplicate route-local non-nuclear loops once parity is established.
5. Optimize allocation inside the neutral owner only after extraction parity.

Validation gates:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged as ignored validation if needed;
- Cr2 q4 non-nuclear `G-A`/`A-A` overlap, kinetic, coordinate moment, and
  second-moment blocks match current implementation at roundoff as ignored
  measurement only;
- residual setup mixed overlap `X` matches the current construction at
  roundoff;
- a small Qiu-White non-nuclear parity fixture passes.

Forbidden in this slice:

- nuclear raw-block changes;
- final-basis `G-G` product-matrix optimization;
- terminal projection or residual-basis transformation;
- Qiu-White semantic changes or route objects;
- pair factors, MWG interaction, parent construction, persistent caches,
  metadata, reports, status fields, artifacts, public API, Cr2 facade support,
  or Cr2 artifact workflow.

## R3 Terminal G-G Product-Matrix Optimization

Status: approved for implementation under `HP-R3GG-FN-01`; not part of the
`HP-CGRB-NN-*` raw-block lane.

Approved boundary:

- owner module `CartesianFinalBasisRealization`;
- source files `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  and, only if needed for small internal terminal-product workspace/helper
  reuse, `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`;
- terminal final-basis `G-G` product matrices used by
  `pqs_terminal_residual_gto_augmented_operators(...)`: kinetic, coordinate
  moments, and second moments.

Allowed source shapes:

- accumulate kinetic-axis product contributions into one destination;
- reuse same-construction base Hamiltonian kinetic `G-G` when available and
  validated equal;
- build/transform coordinate and second-moment product matrices axis by axis;
- share function-local terminal-product scratch across consecutive product
  assemblies;
- delete or simplify `_r3a_product_matrix(...)` if replaced.

Validation gates:

- H2 Residual Gaussian endpoint unchanged;
- Be2 usability/performance measurement unchanged except allowed
  timing/allocation improvement;
- Cr2 q4 `K_GG`, coordinate moment, and second-moment `G-G` products match the
  current construction at roundoff as ignored validation;
- Cr2 exact-operator allocation remeasured after parity.

Forbidden in this slice:

- `G-A`/`A-A` raw-block changes;
- nuclear raw-block or unit-nuclear Gaussian-sum changes;
- IDA/MWG, residual selection/orientation/transform changes;
- terminal basis realization, parent construction, Qiu-White semantics, route
  setup, persistent caches, metadata, reports, status fields, artifacts,
  public API, Cr2 facade support, or Cr2 artifact workflow.

Line budget:

- at most `100` added `src` lines;
- no new committed test file in the first source pass;
- stop for a new amendment if a broad product framework, persistent workspace,
  or files outside the approved source surfaces are required.

## R3 Remaining Exact-Operator Allocation Audit

Status: completed measurement-only attribution under `HP-R3REM-AUDIT-01`.

Decision:

- after `954c86cd`, terminal final-basis `G-G` product workspace is crossed for
  the current Cr2 q4 exact augmented-operator proxy;
- the wrapper now measures about `5.8389s / 4605.517 MiB`;
- remaining allocation is reported mainly outside the nine terminal product
  buffers, including unit-nuclear `U_GG` Gaussian-sum work and route/raw-block
  setup;
- those buckets are not covered by `HP-R3GG-FN-01`.

Follow-up audit result:

- Cr2 q4 wrapper total: `5.7739s / 4680.627 MiB`;
- largest in-wrapper owner: unit-nuclear `U_GG` factor lookup plus
  Gaussian-sum construction at `2.0447s / 1856.819 MiB`;
- non-target or crossed buckets: neutral non-nuclear raw blocks, neutral
  nuclear raw blocks, terminal `G-G` kinetic/moment products, and augmented
  nuclear transforms.

Approved audit boundary:

- use ignored `tmp/work` probes only;
- separate allocation and timing for terminal `G-G` products, unit-nuclear
  `U_GG` Gaussian-sum construction, neutral raw blocks, exact augmented nuclear
  transforms, route/stage setup, raw-block setup, and audit/replay overhead;
- report the next proposed owner and source surface only after attribution is
  clear.

Forbidden:

- no `src`, `test`, `tools`, or `bin` edits;
- no unit-nuclear, route/raw-block setup, raw-block, residual Gaussian,
  IDA/MWG, Qiu-White, parent, terminal-basis, public API, artifact, metadata,
  report/status/payload, persistent cache/workspace, Cr2 facade, or Cr2
  artifact source work.

Exit rule:

- approve no implementation handoff unless a later docs-only amendment names
  the exact ID, owner/files/functions, forbidden surfaces, validation gates,
  line budget, deletion/simplification expectation, and failure rule.

The later amendment is now `HP-R3UN-FN-01` / `HP-R3UN-TEST-01` for the narrow
unit-nuclear `U_GG` lane only.

## R3 Unit-Nuclear U_GG Gaussian-Sum Optimization

Status: approved for implementation under `HP-R3UN-FN-01`; not part of the
`HP-R3GG-*` product-matrix lane or `HP-CGRB-*` raw-block lanes.

Approved boundary:

- owner module `CartesianFinalBasisRealization`;
- source files `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`
  and, only for narrow caller wiring if needed,
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- target functions `_accumulate_terminal_gaussian_sum!` and
  `_terminal_gaussian_sum_action`;
- terminal final-basis unit-nuclear `U_GG` Gaussian-sum construction only.

Allowed source shapes:

- function-local scratch/workspace reuse across Gaussian-sum terms and center
  calls;
- in-place accumulation into caller destinations;
- allocation reduction in factor lookup and terminal Gaussian-sum action;
- small internal scratch arguments or file-local helpers with no persistent
  state;
- deletion/simplification of obsolete allocation-heavy Gaussian-sum helper code
  after parity.

Validation gates:

- H2 Residual Gaussian endpoint unchanged;
- Be2 facade/readback unchanged except allowed timing/allocation improvement;
- Cr2 exact-operator audit reports before/after unit-nuclear `U_GG`
  allocation and total wrapper allocation;
- Cr2 `U_GG` block replay parity and exact-operator parity at roundoff;
- exact operators finite and symmetric.

Forbidden:

- neutral raw-block changes, terminal kinetic/moment `G-G` changes, residual
  Gaussian algorithm changes, exact augmented transform semantic changes,
  IDA/MWG changes, Qiu-White semantic changes, route/stage setup, raw-block
  setup, parent construction, terminal basis realization, persistent
  caches/workspaces, broad Gaussian-sum framework, metadata, reports, status
  fields, payload objects, artifacts, public API, committed tests, Cr2 facade
  support, or Cr2 artifact workflow.

Line budget:

- at most `100` added `src` lines;
- no new committed test file;
- stop for a new amendment if the source pass needs a persistent cache, broad
  Gaussian-sum framework, files outside the approved surfaces, or source edits
  outside terminal unit-nuclear `U_GG`.

## R3 Same-Construction Base K/U Reuse

Status: approved for implementation under `HP-R3BASE-FN-01` and
`HP-R3BASE-TEST-01`.

Decision:

- the supported supplemented construction already builds base product and unit
  nuclear blocks before constructing the base Hamiltonian;
- exact augmented operator construction may reuse those same-construction base
  `K_GG` and unit `U_GG[A]` blocks instead of recomputing them;
- prior replay evidence found exact operator delta `0.0` and reduced exact
  augmented-operator replay to `0.8620s / 1237.136 MiB`.

Approved boundary:

- owner module `CartesianFinalBasisRealization` plus narrow caller wiring;
- source files `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  and `src/cartesian_base_hamiltonian.jl`;
- target functions or wrappers around
  `pqs_terminal_residual_gto_augmented_products(...)`,
  `pqs_terminal_residual_gto_augmented_unit_nuclear(...)`, and
  `cartesian_residual_gto_mwg_hamiltonian(...)` / staged helpers.

Allowed source shapes:

- pass `base_ham.kinetic` as trusted same-construction `K_GG` into augmented
  product construction;
- pass `base_ham.nuclear_attraction_unit_by_center` as trusted
  same-construction unit `U_GG[A]` blocks into augmented unit-nuclear
  construction;
- validate matrix dimensions and center count before reuse;
- preserve existing recomputation behavior when trusted base blocks are not
  supplied;
- keep trust local to the same `cartesian_base_working_basis(...)` construction
  path; do not create provenance payloads or metadata proofs.

Validation gates:

- `git diff --check`;
- package load;
- H2 R3 endpoint unchanged;
- Be2 supplemented facade/readback unchanged except allowed
  timing/allocation improvement;
- Cr2 exact-operator attribution audit or focused ignored replay showing base
  `K_GG` / unit `U_GG[A]` reuse parity and allocation effect;
- final exact operators finite and symmetric;
- no Cr2 artifact/workflow.

Forbidden:

- public API/export changes, canonical driver changes, raw-block changes,
  residual selection/orientation/transform changes, MWG/IDA convention changes,
  terminal product or Gaussian-sum kernel rewrites, persistent cache/workspace
  objects, metadata/status/report/artifact schema fields, route/stage setup
  cleanup, committed tests, or Cr2 workflow.

Line budget:

- target under `100` added `src` lines;
- stop without a source commit if same-construction trust cannot be guaranteed
  by local call shape and dimension/center validation, or if implementation
  needs public payloads, metadata, or stage objects.

## R3 Driver Call-Site Base K/U Reuse Wiring

Status: approved for implementation under `HP-R3BASE-DRV-WIRE-01` and
`HP-R3BASE-DRV-TEST-01`.

Approved boundary:

- source file `bin/cartesian_ham_builder.jl`;
- supplemented mode only;
- call sites for `cartesian_residual_gto_augmented_products(...)` and
  `cartesian_residual_gto_augmented_unit_nuclear(...)`.

Allowed source shapes:

- pass `base_ham.kinetic` as `base_kinetic` to augmented products;
- pass `base_ham.nuclear_attraction_unit_by_center` as `base_unit_nuclear` to
  augmented unit-nuclear construction;
- leave public inputs, hooks, timing labels, visible stage sequence, artifact
  schema, and driver contract unchanged.

Validation gates:

- `git diff --check`;
- package load;
- H2 supplemented driver artifact/readback;
- Be2 supplemented driver artifact/readback if practical;
- no Cr2 run.

Forbidden:

- source/kernel changes, diagnostics, new hooks, new timing labels, public
  input changes, visible stage-sequence changes, artifact schema changes,
  committed tests/fixtures, Cr2 workflow, or files outside
  `bin/cartesian_ham_builder.jl`.

Failure rule:

- if this needs any visible driver contract change, stop and report.

## Canonical Cartesian Driver Usability

Status: approved for implementation under `HP-DRV-FILE-01`,
`HP-DRV-FN-01`, `HP-DRV-STAGE-FN-01`, `HP-DRV-STAGE-WIRE-01`,
`HP-DRV-STAGE-TEST-01`, and `HP-DRV-TEST-01`.

Approved boundary:

- source file `bin/cartesian_ham_builder.jl`;
- source file `src/cartesian_base_hamiltonian.jl` for the driver-facing staged
  producer surface;
- source files `src/pqs_source_box_low_order_materialization.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` only for
  behavior-preserving physical operator-class stage factoring;
- compact driver invocation
  `julia --project=. bin/cartesian_ham_builder.jl [input.jl] [key=value ...]`;
- visible editable defaults, one optional trusted project input file,
  command-line overrides, visible public contract construction, compact run
  summary, visible physics-level construction stages, coarse phase timing,
  artifact write, and optional readback check.

Allowed workflow:

- construct public `system`, `basis`, and optional `supplement` objects before
  calling a facade or staged producer surface;
- call approved base and staged producer surfaces for base H/H2;
- call the approved non-exported residual-GTO/MWG usability facade or the
  staged producer surface for supported supplemented H2 and
  internal/performance-supported Be2;
- call the staged producer surface as separate visible top-level stage calls,
  not as one all-in-one staged replacement wrapper;
- time and print visible physics-level stages: base working basis / terminal
  realization, product/moment operators, unit-nuclear attraction operators,
  electron-electron / IDA or residual-MWG interactions, base Hamiltonian
  assembly, Gaussian supplement, residual augmentation, supplemented
  Hamiltonian assembly, and artifact write/check;
- write existing `CartesianIDAHamiltonian` artifacts with approved provenance
  groups;
- print user-facing summaries and timing.

Approved run-level hooks:

- `check_file`;
- `print_contract`;
- `print_timing`;
- `expected_dimension`.

`basisname = nothing` selects base mode. `basisname !== nothing` selects
supported supplemented diatomic mode, is the visible supplement basis label,
and must reject `Natom == 1`. `padding` is physical box padding: it maps to
one-center `radius` for atoms and to z-axis diatomic facade extents around the
two nuclei.

Forbidden:

- private route-stage controls, stop-after internals, ladder probes, stage
  markers, fixture hacks, diagnostic knobs, underscored package helper calls,
  raw-block provider switches, report/status/payload dumps, metadata clouds,
  allocation probes, per-kernel timing frameworks, benchmark harness behavior,
  solver/RHF/ECP/EGOI/HamV6,
  private contract construction, artifact schema dumps, public API/export
  changes, artifact schema changes, committed tests, committed input fixtures,
  supplemented atoms, old route-stage choreography, Cr2-specific driver runs,
  or Cr2-specific workflow support. Generic explicit homonuclear z-axis Cr2
  stress through
  `HP-R3U-ZDI-WIRE-01` is separate ignored/user-run validation authority.

Validation gates:

- package load;
- public contract print/check output for at least one base run when touched;
- visible base-stage timing/summary for H atom or H2 base construction when
  staged wiring changes;
- visible supplemented-stage timing/summary for H2 supplemented construction
  when staged wiring changes;
- H atom base driver artifact write/readback under `HP-DRV-ATOM-TEST-01`;
- H2 base driver artifact write/readback;
- H2 supplemented driver artifact write/readback;
- optional ignored Be2 usability run for supplemented-mode changes.

Line budget:

- at most `150` added `bin` lines;
- at most `200` added `src` lines across the approved staged-driver source
  files;
- no new committed test or tool file;
- stop for a new amendment if a parser framework, source files outside the
  canonical driver and staged producer owner, route-stage diagnostics,
  raw-block changes, kernel rewrites, status/report/payload expansion,
  artifact schema changes, public API/export changes, or Cr2-specific workflow
  support are required.

## Canonical Driver Atom Workflow

Status: approved for implementation under `HP-DRV-ATOM-FN-01`,
`HP-DRV-ATOM-WIRE-01`, and `HP-DRV-ATOM-TEST-01`.

Approved boundary:

- source file `bin/cartesian_ham_builder.jl`;
- explicit one-center atom input in `mode = :base`;
- origin-centered atoms only;
- base Hamiltonian artifact writing through the existing
  `cartesian_base_hamiltonian` facade where that facade already supports the
  requested atom.

Allowed source shapes:

- normalize explicit `system` fields: `atom_symbols`, `nuclear_charges`,
  `atom_locations`, `nup`, and `ndn`;
- require one center at `(0.0, 0.0, 0.0)`;
- validate neutral all-electron count from the explicit nuclear charge;
- pass visible one-center basis fields, including `core_spacing`, to the base
  facade;
- allow visible, easily edited driver/project defaults such as
  `core_spacing = 0.3` and template `padding`, while treating them as explicit
  resolved input values that may be overridden for quick tests;
- routine correctness tests may override driver physics defaults, but any
  asserted scalar must be tied to the exact test input and not described as a
  physics-default result;
- keep any example/default atom input explicit, not inferred from element
  tables.

Decision:

- base atom driver output only is approved;
- current validation remains origin-centered H;
- supplemented atom Hamiltonians remain candidate-only.

Forbidden:

- source edits outside `bin/cartesian_ham_builder.jl`, except under separate
  `HP-R1-ATOM-*` authority;
- supplemented atom Hamiltonians;
- translated atoms;
- broader base atom support beyond the existing base facade;
- element lookup/default tables, ECP, pseudopotentials, solver/RHF workflow,
  public API/export changes, artifact schema changes, route diagnostics,
  metadata/status/report fields, package-internal helper composition from the
  driver, committed atom fixture files, committed tests, or new tool files.

Validation gates:

- H atom base driver artifact write/readback with explicit system and basis;
- optional ignored negative checks for non-origin atom input, nonneutral
  electron count, mismatched temporary `d` if accepted, or unsupported atom
  input;
- no supplemented atom or translated-atom validation.

Line budget:

- at most `80` added `bin` lines;
- no new committed test, tool, or input-fixture file;
- stop for a new amendment if the implementation needs any forbidden surface.

## Canonical Driver Atom Hidden-d Cleanup

Status: approved for implementation under `HP-DRV-ATOM-CLEAN-01`.

Approved boundary:

- source file `bin/cartesian_ham_builder.jl`;
- one-center atom basis construction only.

Allowed source shape:

- remove the hidden `d = vars[:core_spacing]` field from the one-center atom
  `basis` construction;
- keep visible atom basis fields `q`, `core_spacing`, `radius`, and existing
  optional public fields unchanged.

Validation gates:

- `git diff --check`;
- package load;
- H atom base driver artifact/readback;
- H2 base or supplemented driver smoke only if the changed code path shares the
  touched construction.

Forbidden:

- public input changes, default changes, override changes, hook changes, timing
  label changes, visible stage-sequence changes, artifact schema changes,
  diagnostics, source/kernel changes, committed tests/fixtures, Cr2 workflow,
  old `:white_lindsey_low_order` retirement, test/tool route-input cleanup, or
  files outside `bin/cartesian_ham_builder.jl`.

Failure rule:

- if removing the hidden `d` field requires any visible driver contract change
  or producer/source change, make no source commit and report the blocker.

## R1 One-Center Base Atoms

Status: approved for implementation under `HP-R1-ATOM-FN-01`,
`HP-R1-ATOM-WIRE-01`, and `HP-R1-ATOM-TEST-01`.

Approved boundary:

- source file `src/cartesian_base_hamiltonian.jl`;
- existing public call shape
  `cartesian_base_hamiltonian(system; basis, hamfile = nothing)`;
- explicit origin-centered one-center all-electron atoms only;
- current committed regression remains the origin-centered H endpoint.

Allowed source shapes:

- replace H-specific one-center validation with explicit one-center atom
  validation;
- require explicit symbol, nuclear charge, origin location, `nup`, `ndn`, and
  one-center basis controls;
- derive neutral all-electron count from the explicit nuclear charge, not an
  element table;
- map the explicit nuclear charge to private White-Lindsey atomic mapping `Z`;
- map resolved public `core_spacing` to private `parent_mapping_d`;
- keep `reference_spacing`, `tail_spacing`, and box/domain controls separate
  from `core_spacing`;
- reuse existing `HP-R1-ART-01` provenance keys and
  `route = :one_center_pqs_base`.

Shared workflow requirement:

- atoms and diatomics must share the same producer workflow after
  geometry/shellification normalization;
- do not add an atom-only Hamiltonian builder, atom materialization path,
  atom route-stage object, atom report/status/payload object, or atom-specific
  artifact shape;
- if an atom cannot pass through the shared terminal-basis, one-body, IDA,
  Hamiltonian-construction, and writer machinery, stop and report the missing
  shared seam.

Validation gates:

- H public facade endpoint unchanged;
- optional ignored/user-run Be or Cr atom artifact write/readback with explicit
  charge, spin sectors, origin geometry, and basis controls;
- finite/symmetric `K`, unit `U_A`, and IDA `V` for ignored/user-run non-H
  atom checks;
- unsupported translated atoms, noninteger/nonpositive charge, nonneutral
  electron count, mismatched temporary `d`, and element-table/default requests
  throw clear errors where practical.

Forbidden:

- source edits outside `src/cartesian_base_hamiltonian.jl`, private
  materialization-owner edits, atom-only materialization, supplemented atom
  Hamiltonians, translated atoms, ECP/pseudopotentials, solver/RHF workflow,
  public API/export changes, artifact schema changes, element lookup/default
  tables, route diagnostics, metadata/status/report fields, committed non-H
  atom fixtures, committed tests, or driver changes.

Line budget:

- at most `80` added `src` lines;
- no new committed test, tool, driver, or input-fixture file;
- stop for a new amendment if the implementation needs any forbidden surface.

## Route Recipe Family Selection Cleanup

Status: approved for implementation under `HP-ROUTE-RECIPE-FN-01` and
`HP-ROUTE-RECIPE-TEST-01`.

Approved boundary:

- source files `src/pqs_source_box_route_driver_helpers.jl` and
  `src/cartesian_base_hamiltonian.jl`;
- target function `cartesian_recipe(route_inputs)`;
- target base helper `_cartesian_base_route(kind)`;
- existing tests only if they directly construct route inputs that need
  inactive-family fields removed; no new committed tests.

Allowed source shapes:

- validate `route_family` and build only the selected family subrecipe;
- for `:pqs_source_box`, require only source-box recipe fields and set inactive
  `white_lindsey` data to `nothing` or an equivalent non-authoritative empty
  value;
- for `:white_lindsey_low_order`, preserve explicit WL route support and build
  the selected WL subrecipe from existing WL route fields;
- remove unused `white_lindsey_*` fields from `_cartesian_base_route(kind)`,
  because the live base producer route is PQS-only;
- keep existing precomposed recipe compatibility only if a live caller still
  needs it.

Must delete or simplify:

- the PQS base producer should no longer carry inactive
  `white_lindsey_route_shape`, `white_lindsey_mapping_rule`,
  `white_lindsey_nesting_rule`, `white_lindsey_retained_rule`, or
  `white_lindsey_operator_rule` fields.

Validation gates:

- `git diff --check`;
- package load;
- H atom/base artifact readback;
- H2 base artifact readback;
- compact H2 supplemented facade or driver path;
- focused explicit `:white_lindsey_low_order` route recipe smoke if practical,
  or report exact live test/tool callers that block further WL cleanup;
- no Cr2 run.

Forbidden:

- canonical driver changes, numerical kernel changes, terminal lowering policy
  changes, shellification behavior changes, materialization/artifact schema
  changes, route-stage diagnostics, status/report expansion, deletion of WL
  materialization, new route-stage objects, new payload/cache structs, new
  committed tests, or source files outside the approved surfaces.

Failure rule:

- if `cartesian_recipe(...)` cannot be made family-selective without broader
  route-driver, report, materialization, or stage-object changes, make no
  source commit and report the blocker.

Line budget:

- at most `80` added `src` lines, with net simplification expected;
- no new committed test file.

## Route Inventory Type-Surface Cleanup

Status: first lane approved for implementation under `HP-ROUTE-INV-FN-01` and
`HP-ROUTE-INV-TEST-01`.

Approved boundary:

- source file `src/pqs_source_box_route_driver_helpers.jl` only;
- remove runtime-keyed retained-unit inventory `NamedTuple{unit_keys}` shapes;
- remove runtime-keyed `pair_family_counts = NamedTuple{families}(...)`;
- update same-file callers to use vector-backed records/tables, stable
  dictionaries, or storage-hidden helper accessors.

Deferred from `HP-ROUTE-INV-*`:

- `TerminalLoweringPlan.available_contracts` and `contracts`;
- `RetainedUnitTransformContractPlan.contracts`;
- public input `NamedTuple`s;
- small fixed `NTuple{3,Int}` coordinate/dimension values;
- artifact sidecar tables.

Raw product source-mode inventory cleanup is separately approved under
`HP-RAW-SRCMODE-FN-01` / `HP-RAW-SRCMODE-TEST-01`.
Contract-plan tuple cleanup is separately approved under
`HP-CONTRACT-VEC-FN-01` / `HP-CONTRACT-VEC-TEST-01`.

Validation gates:

- `git diff --check`;
- package load;
- H2 base artifact write/readback;
- H2 supplemented artifact write/readback or canonical driver path;
- focused search for absence of `NamedTuple{unit_keys}` and
  `NamedTuple{families}` route inventories;
- no Cr2 run.

Forbidden:

- source files outside the approved file, numerical kernel changes, route
  recipe behavior changes, shellification, terminal lowering, terminal basis,
  Residual Gaussian, raw product source, raw-block changes, canonical driver
  changes, Hamiltonian object changes, matrix-key changes, reader changes,
  artifact schema changes, public API/export changes, report/status/payload
  expansion, compatibility adapters, new committed tests, Cr2 runs, or
  Cr2-specific workflow.

Line budget:

- at most `120` added `src` lines, with net simplification expected;
- no new committed test, tool, benchmark, or input-fixture file.

## Raw Product Source-Mode Inventory Cleanup

Status: approved for implementation under `HP-RAW-SRCMODE-FN-01` and
`HP-RAW-SRCMODE-TEST-01`.

Approved boundary:

- `src/cartesian_raw_product_sources/records.jl`;
- `src/cartesian_raw_product_sources/source_mode_indices.jl`;
- `src/cartesian_raw_product_sources/summaries.jl`;
- narrow consumer wiring only as needed in
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  and `src/cartesian_base_hamiltonian.jl`.

Approved changes:

- replace `RawProductBoxPlan.source_mode_indices::Tuple{Vararg{NTuple{3,Int}}}`
  with vector-backed storage;
- replace `RawProductBoxPlan.source_mode_column_indices::Tuple{Vararg{Int}}`
  with vector-backed storage, or remove it when it is only `1:count` and
  accessors supply the same ordered column numbers;
- preserve deterministic source-mode ordering, fixed `NTuple{3,Int}` mode
  coordinates, existing source-mode facts, retained-rule behavior, and manifest
  source-mode/source-relation output;
- preserve semantic access through stable accessors without preserving
  variable-length `Tuple` concrete return types.

Validation gates:

- `git diff --check`;
- package load;
- H2 base artifact write/readback;
- H2 supplemented artifact write/readback;
- H2 R3 endpoint;
- focused raw-product source order and retained-rule parity;
- manifest source-mode and final-basis source-relation inspection;
- focused search confirming `RawProductBoxPlan` no longer stores source-mode
  inventories as `Tuple{Vararg{...}}`;
- no Cr2 run.

Forbidden:

- source files outside the approved surfaces, terminal-lowering contract tuple
  cleanup, retained-unit transform-contract tuple cleanup beyond narrow caller
  wiring, broad pair-block/source-box rewrites, public input `NamedTuple`
  changes, fixed `NTuple{3,Int}` coordinate/dimension changes, numerical
  kernel changes, route semantic changes, terminal lowering, shellification,
  terminal basis, Residual Gaussian, raw Gaussian block, IDA/MWG, Qiu-White
  semantic changes, canonical driver changes, Hamiltonian object changes,
  matrix-key changes, reader changes, artifact schema changes, public
  API/export changes, report/status/payload expansion, persistent caches,
  compatibility adapters preserving the old tuple-backed shape, new committed
  tests, Cr2 runs, or Cr2-specific workflow.

Line budget:

- at most `150` added `src` lines, with net simplification expected;
- no new committed test, tool, benchmark, or input-fixture file.

Failure rule:

- if vectorizing the raw product plan forces broad pair-block/source-box
  rewrites or source files outside the approved surfaces, stop and report the
  exact callers/blockers instead of adding compatibility layers.

## Contract-Plan Vector Cleanup

Status: approved for implementation under `HP-CONTRACT-VEC-FN-01` and
`HP-CONTRACT-VEC-TEST-01`.

Approved boundary:

- `src/cartesian_terminal_lowering/contracts.jl`;
- `src/cartesian_terminal_lowering/selection.jl`;
- `src/cartesian_terminal_lowering/summaries.jl`;
- `src/cartesian_retained_unit_transform_contracts/records.jl`;
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`;
- `src/cartesian_retained_unit_transform_contracts/summaries.jl`;
- narrow consumer wiring only as needed in
  `src/pqs_source_box_route_driver_helpers.jl`,
  `src/cartesian_base_hamiltonian.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

Approved changes:

- replace
  `TerminalLoweringPlan.available_contracts::Tuple{Vararg{TerminalLoweringContract}}`
  with vector-backed storage;
- replace
  `TerminalLoweringPlan.contracts::Tuple{Vararg{TerminalLoweringContract}}`
  with vector-backed storage;
- replace
  `RetainedUnitTransformContractPlan.contracts::Tuple{Vararg{RetainedUnitTransformContract}}`
  with vector-backed storage;
- preserve `available_contracts(plan)`, `selected_contracts(plan)`,
  `contracts(plan)`, and `transform_contracts(plan)`;
- preserve iteration order, selected-contract semantics, transform-contract
  semantics, summaries, and existing behavior without preserving
  variable-length `Tuple` concrete field types.

Validation gates:

- `git diff --check`;
- package load;
- H2 base artifact write/readback;
- H2 supplemented artifact write/readback;
- H2 R3 endpoint;
- focused terminal-lowering contract order parity;
- focused retained-unit transform-contract order parity;
- focused search confirming the targeted plan inventories no longer store
  contracts as `Tuple{Vararg{...}}`;
- no Cr2 run.

Forbidden:

- `source_cpbs::Tuple{Vararg{CoordinateProductBox}}` changes, source files
  outside the approved surfaces, raw product source-mode changes, retained-unit
  route inventory changes, public input `NamedTuple` changes, fixed
  coordinate/product-box value-object changes, numerical kernel changes, route
  semantic changes, shellification behavior changes, raw Gaussian block,
  Residual Gaussian, IDA/MWG, Qiu-White semantic changes, canonical driver
  changes, Hamiltonian object changes, matrix-key changes, reader changes,
  artifact/manifest schema changes, public API/export changes,
  report/status/payload expansion, persistent caches, compatibility adapters
  preserving the old tuple-backed plan field types, new committed tests, Cr2
  runs, or Cr2-specific workflow.

Line budget:

- at most `150` added `src` lines, with net simplification expected;
- no new committed test, tool, benchmark, or input-fixture file.

Failure rule:

- if vectorizing the plan inventories requires broad route/stage rewrites,
  source files outside the approved surfaces, public API or artifact changes,
  numerical changes, or compatibility layers preserving the old tuple-backed
  plan field types, stop and report the exact callers/blockers.

## Homonuclear Z-Axis Diatomic Supplemented Workflow

Status: approved for implementation under `HP-R3U-ZDI-FN-01`,
`HP-R3U-ZDI-WIRE-01`, and `HP-R3U-ZDI-TEST-01`.

Approved boundary:

- source file `src/cartesian_base_hamiltonian.jl` for the non-exported
  supplemented facade;
- source file `bin/cartesian_ham_builder.jl` for canonical driver
  supplemented-mode wiring;
- explicit homonuclear two-center z-axis diatomics only.

Allowed source shapes:

- replace H/Be-specific guards in `_cartesian_r3_diatomic_inputs(...)` or the
  local equivalent with explicit homonuclear z-axis validation;
- require explicit system, base basis, supplement, and optional `basisfile`
  inputs;
- let the canonical driver call the supported supplemented facade;
- keep Cr2 as an optional ignored/user-run stress through the generic path.

Validation gates:

- H2 supplemented facade/driver artifact write/readback unchanged;
- Be2 supplemented facade/driver artifact write/readback unchanged;
- optional ignored/user-run Cr2 stress/usability run after H2/Be2 pass.

Forbidden:

- element-specific defaults, Cr2-specific branch, committed Cr2 fixture or
  gate, heteronuclear support, non-z-axis/general orientation support, charged
  systems, ECP, solver/RHF, public export/API redesign, artifact schema
  changes, route diagnostics, metadata/status/report fields, package-internal
  helper composition from the driver, new route objects, committed tests, or
  new source/tool/input-fixture files.

Line budget:

- at most `100` added `src`/`bin` lines total;
- net simplification expected where hardcoded H/Be checks are removed;
- stop for a new amendment if the implementation needs any forbidden surface.
