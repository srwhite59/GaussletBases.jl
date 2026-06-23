# Registry

Only entries marked approved/implemented authorize work on the exact surface
they describe. Measurement-only entries do not authorize production source
edits. Candidate or rejected entries do not authorize implementation.

## Approved And Implemented

### HP-OBJ-01 — `CartesianTerminalBasisBlock`

Owner: `CartesianFinalBasisRealization`.

Exact fields:

```julia
struct CartesianTerminalBasisBlock
    unit_key::Symbol
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    coefficients::Union{Nothing,Matrix{Float64}}
    column_range::UnitRange{Int}
end
```

`coefficients === nothing` means direct identity on the listed support rows.
Otherwise `coefficients` is support-local rows by retained columns.
`support_indices` and `support_states` are the authoritative owned terminal
rows for the block. They must not mean a post-projection enlarged support, and
PQS shell realization must not grow a block onto previous terminal regions.
Terminal blocks are block-sparse by representation; direct blocks remain
implicit identity on their owned support.

### HP-OBJ-02 — `CartesianTerminalBasisRealization`

Owner: `CartesianFinalBasisRealization`.

Exact fields:

```julia
struct CartesianTerminalBasisRealization
    blocks::Vector{CartesianTerminalBasisBlock}
    final_dimension::Int
    max_cross_overlap::Float64
end
```

The object does not store parent bundles, stage objects, summaries, metadata,
global coefficients, or global self-overlap matrices. Its blocks are
represented on disjoint owned terminal regions. This block-sparse basis
representation does not imply block-diagonal kinetic, nuclear, or IDA
operators. Overlap between distinct owned-support blocks is structurally zero.

The implemented `max_cross_overlap` field is legacy implementation debt after
the structural-support correction. It must not be treated as a physical
cross-block residual or construction repair signal. Source cleanup should
replace production cross-overlap audit plumbing with structural support checks
under a separate implementation handoff.

### HP-FILE-01 — terminal realization file

Approved file:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

### HP-FN-00 — block-local terminal shell realization

Approved as a file-local/internal helper under HP-FILE-01. It realizes one PQS
terminal shell by:

```text
full source-box product modes
-> boundary product-mode column selection
-> restrict rows to support.support_indices / support.support_states
-> shell-local Gram on that owned support
-> symmetric shell-local Lowdin
-> final sign canonicalization
-> append block with unchanged owned support
```

Previous-block projection, recursive projection, projection basis, and
effective-support growth are forbidden. Full source-box modes are used only to
generate boundary product-mode columns before row restriction to the
shell-owned support.

### HP-FN-01 — terminal basis realizer

Approved internal surface:

```julia
pqs_terminal_basis_realization(
    support_records,
    retained_records,
    transform_contracts,
    bundles;
    identity_atol = 1.0e-8,
    cross_atol = 1.0e-8,
    weight_atol = 1.0e-14,
)
```

Returns `CartesianTerminalBasisRealization` on success. It owns direct-sector
checks, PQS shell realization, positive final-integral sign canonicalization,
and construction of `CartesianTerminalBasisBlock`.

### HP-FN-02 — structural terminal support checks

Corrected design authority: parent gausslet rows are orthonormal to machine
precision and terminal regions own disjoint parent rows. Therefore block-local
terminal basis supports are structurally orthogonal across blocks, and
cross-block overlap is zero by construction.

Production validation should check:

- every block support equals its authoritative terminal support;
- terminal support sets are pairwise disjoint;
- each shell-local overlap is identity within tolerance.

A nonzero structural overlap means duplicated support rows, incorrect row
restriction, wrong support ownership, or an indexing error. It is not a
physical residual to compute or repair, and production code must not mix
coefficients into previous supports to reduce it. This correction does not
approve source cleanup by itself; removing `max_cross_overlap` and replacing
the audit plumbing needs a separate implementation blurb or approved cleanup
surface.

### HP-WIRE-01 — terminal-basis stage integration

Approved owner: `cartesian_transforms`.

It connects supported PQS terminal plans to
`pqs_terminal_basis_realization(...)` from typed terminal support, retained, and
transform records. It must serve one-center atomic, contact-core diatomic, and
separated diatomic terminal plans through the same entry point.

### HP-FN-03 — blockwise one-body assembly

Approved file:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
```

Approved public internal helper:

```julia
assemble_terminal_product_operator!(
    destination,
    basis::CartesianTerminalBasisRealization,
    axis_x,
    axis_y,
    axis_z;
    scale = 1.0,
)
```

Approved file-local Gaussian-sum helper:

```julia
_accumulate_terminal_gaussian_sum!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    factors_x,
    factors_y,
    factors_z;
    scale = -1.0,
)
```

No `K`/`U_A` payload, stage field, report object, persistent cache, or
orchestration API is approved by HP-FN-03.

### HP-FN-04 — localized IDA assembly

Approved file:

```text
src/cartesian_final_basis_realization/pqs_terminal_ida.jl
```

Approved function:

```julia
assemble_terminal_ida_interaction!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    raw_pair_terms_x,
    raw_pair_terms_y,
    raw_pair_terms_z,
    weights_x,
    weights_y,
    weights_z;
    weight_atol = 1.0e-12,
    symmetry_atol = 1.0e-10,
)
```

This is Slice C1 only: it produces final-basis `electron_electron_ida`. It does
not authorize Hamiltonian construction, route wiring, artifacts, or a pair
payload/cache.

### HP-FN-05 — final Hamiltonian construction

Approved as the narrow Slice C2 construction boundary for the existing
`CartesianIDAHamiltonian`.

Conceptual boundary:

```julia
build_cartesian_ida_hamiltonian(
    kinetic,
    nuclear_attraction_unit_by_center,
    electron_electron_ida,
    nup,
    ndn;
    nuclear_charges,
    nuclear_positions,
)::CartesianIDAHamiltonian{Float64}
```

Implementation may call the existing `CartesianIDAHamiltonian(...)`
constructor directly if no helper is needed.

### HP-WIRE-02 — direct materialization Hamiltonian handoff

Approved and implemented Slice D boundary:

```julia
cartesian_materialization(
    report,
    terminal_basis_realization,
    materialization_inputs,
)::Union{Nothing,CartesianIDAHamiltonian{Float64}}
```

The call site passes `transforms.terminal_basis_realization` directly. The
terminal basis must not be embedded in `cartesian_report`, reconstructed from
summaries, or recovered by passing the full `transforms` stage.

Return contract:

- no request returns `nothing`;
- requested base PQS materialization returns `CartesianIDAHamiltonian{Float64}`;
- `save_ham_artifact = true` writes with `write_cartesian_ida_hamiltonian` and
  still returns the same Hamiltonian;
- no materialization wrapper, `result_kind`, `materialized`, status mirror, or
  `ida_hamiltonian` field is approved.

## Approved For R1 Implementation

These entries authorize only the R1 public base producer scope recorded in
`r1_public_base_producer.md`. They do not approve broad driver polish,
additional routes, new artifact shapes beyond the approved `HP-R1-ART-01`
`producer_provenance/` keys in the final Hamiltonian file, solver work,
supplements, corrections, or status/report/payload expansion.

### HP-R1-FILE-01 — public base producer source file

Approved file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved owner: top-level `GaussletBases` public API.

The only new export approved by R1 is `cartesian_base_hamiltonian` from
`src/GaussletBases.jl`.

### HP-R1-FN-01 — public base Hamiltonian producer facade

Approved public call shape:

```julia
cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

Scope:

- origin-centered base H and Cartesian z-axis aligned base H2 first;
- plain `NamedTuple` input groups only;
- no public `method`, `route`, or output group in R1;
- `n_s`, bond length, and private H2 radius are derived internally;
- one-center H requires explicit public `d` with no default;
- one-center H maps `reference_spacing` and `d` separately:
  `spacing_inputs.reference_spacing = basis.reference_spacing`,
  `parent_inputs.parent_mapping_rule = :white_lindsey_atomic_mapping`, and
  `parent_inputs.parent_mapping_d = basis.d`;
- the reviewed H baseline uses explicit public `d = 0.3` and
  `reference_spacing = 1.0`;
- z-axis H2 rejects public `d` because it uses the multicenter mapping
  contract;
- x/y-aligned diatomics, shifted-parallel diatomics, generally oriented
  molecules, translation, and rotation are deferred;
- center-sized public collections must be vectors or other `AbstractVector`
  values, not variable-size tuples;
- unknown public input keys throw `ArgumentError`;
- scalar inputs must be positive and finite where applicable;
- symbols, charges, coordinates, and electron counts must match the approved
  H/H2 scope;
- return the existing `CartesianIDAHamiltonian{Float64}` directly;
- no wrapper, payload, status object, report mirror, or new artifact shape
  except the approved `HP-R1-ART-01` `producer_provenance/` keys in the final
  Hamiltonian file;
- non-`nothing` `hamfile` writes with existing
  `write_cartesian_ida_hamiltonian`; production does not automatically
  read back the artifact.

### HP-R1-WIRE-01 — report-free base producer wiring

Approved wiring for the R1 public facade:

```text
system / specification
-> parent and route geometry
-> terminal basis realization
-> Hamiltonian production
-> optional existing Hamiltonian artifact
```

The recommended base-public path must not require `cartesian_pair_terms` or
`cartesian_assembly`. Existing stages may remain temporarily for legacy
script/report compatibility, but this R1 authority does not approve adding a
new base-route consumer to either stage.

Approved private shared constructor seam:

```julia
_cartesian_base_ida_hamiltonian(
    terminal_basis_realization,
    parent_axis_bundle_object,
    atom_locations::Vector{NTuple{3,Float64}},
    nuclear_charges::Vector{Float64},
    nup::Int,
    ndn::Int,
)::CartesianIDAHamiltonian{Float64}
```

Approved owner file:

```text
src/pqs_source_box_low_order_materialization.jl
```

This PQS-named owner is acceptable for the explicitly PQS-only R1 migration.
It is not permanent method-neutral ownership for R2.

Allowed caller files:

- `src/pqs_source_box_low_order_materialization.jl`
- `src/cartesian_base_hamiltonian.jl`

The existing report-bound materialization path and the new public facade should
call this same private constructor. The implementation must remove or bypass
the current report dependency where `cartesian_assembly` exists chiefly to
supply `route_skeleton` and a low-order shellification summary to
`cartesian_report`, and delete the duplicated report-bound numerical
orchestration it replaces. It must not fabricate a private report object,
duplicate the Hamiltonian builder in the new public file, leave two parallel
base Hamiltonian construction paths, or replace the dependency with a new
report field cloud, status payload, or metadata-carried numerical data.

### HP-R1-ART-01 — public base producer artifact provenance

Approved artifact extension for R1 public facade writes only. When
`hamfile !== nothing`, the final Hamiltonian file may add the following JLD2
keys under `producer_provenance/`:

```text
provenance_version
producer
route
q
core_spacing
reference_spacing
tail_spacing
parent_axis_family
parent_axis_counts
mapping_kind
mapping_d
radius
xmax_parallel
xmax_transverse
atom_symbols
nuclear_charges
atom_locations
nup
ndn
final_dimension
```

Exact values and route-specific `nothing` fields are defined in
`r1_public_base_producer.md`. For one-center H, the provenance must record
`reference_spacing = 1.0`, `mapping_kind = :white_lindsey_atomic_mapping`, and
`mapping_d = 0.3` for the reviewed endpoint. Existing
`read_cartesian_ida_hamiltonian` must continue reading the Hamiltonian matrices
while ignoring these extra keys. This ID does not approve a separate manifest,
provenance file, wrapper object, public provenance reader, status/report
payload, or use of provenance as staged algorithm data after initial
lattice/parent construction.

### HP-R1-TEST-01 — public base producer endpoint test/example

Approved committed validation surface for R1.

Approved standalone integration gate:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

Invocation:

```text
julia --project=. test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

The test/example should exercise the public facade for one-center H and
z-axis H2, verify `CartesianIDAHamiltonian{Float64}` output, validate the
reviewed H baseline using explicit public `d = 0.3` with
`reference_spacing = 1.0` and H2 endpoint facts, validate unknown-key and
malformed input errors including missing H `d` and rejected H2 `d`, validate
x/y-aligned, shifted-parallel, and generally oriented H2 rejection before
expensive construction, and validate existing Hamiltonian artifact
write/readback plus `HP-R1-ART-01` provenance using `mktempdir()`. It is a
standalone integration/endpoint gate, not ordinary tiny unit coverage, and this
ID does not approve adding it to `test/runtests.jl`. It must not assert private
route-stage fields, report mirrors, status/blocker symbols, terminal role
vocabulary, or pair inventories.

## Approved For R1 One-Center Base Atoms

This section approves only the explicit origin-centered one-center
all-electron atom relaxation recorded in `r1_one_center_base_atoms.md`. It
extends the existing base facade scope without adding a new public function,
new export, new artifact schema, new route vocabulary, or supplemented atom
authority.

### HP-R1-ATOM-FN-01 — explicit one-center all-electron base atom facade

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- accept exactly one origin-centered atom in the existing
  `cartesian_base_hamiltonian(system; basis, hamfile)` call shape;
- require explicit vector-valued `atom_symbols`, `nuclear_charges`,
  `atom_locations`, and explicit integer `nup`, `ndn`;
- require finite positive integer-valued nuclear charge supplied by the caller;
- require neutral all-electron count
  `nup + ndn == round(Int, only(nuclear_charges))`;
- treat the atom symbol as provenance/user labeling only, not as a source of
  charge, spin, basis, or ECP defaults;
- keep required one-center basis fields `q`, `core_spacing`, `radius`, and
  explicit public `d`.

This ID does not approve translated atoms, element lookup/default tables,
inferred charge or spin, ECP, solver workflow, supplemented atoms, public API
redesign, or new artifact fields.

### HP-R1-ATOM-WIRE-01 — one-center atom shared workflow wiring

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- map public `only(system.nuclear_charges)` to the existing private
  White-Lindsey atomic mapping `Z`;
- keep `d`, `core_spacing`, and `reference_spacing` independent;
- feed atom geometry/shellification normalization into the same terminal-basis,
  one-body, IDA, `CartesianIDAHamiltonian`, artifact-writing, and provenance
  machinery used by the base producer;
- preserve existing `HP-R1-ART-01` `producer_provenance/` keys with
  `route = :one_center_pqs_base`.

Atoms and diatomics must share the same producer workflow after the narrow
geometry/shellification differences. This ID does not approve an atom-only
Hamiltonian builder, parallel atom materialization path, atom route-stage
object, atom report/status/payload object, or metadata/provenance carrier used
as algorithmic data.

Line budget for `HP-R1-ATOM-FN-01` plus `HP-R1-ATOM-WIRE-01`: at most `80`
added `src` lines. If implementation needs source edits outside
`src/cartesian_base_hamiltonian.jl`, changes to private materialization
owners, atom-only materialization, new artifact keys, translated atoms,
supplemented atoms, ECP behavior, solver workflow, element lookup/default
tables, committed fixtures/tests, or route/report/status/payload expansion,
stop and request a new docs-only amendment.

### HP-R1-ATOM-TEST-01 — one-center base atom validation

Approved validation:

- existing origin-centered H public facade endpoint remains unchanged,
  including the reviewed `d = 0.3`, `reference_spacing = 1.0` baseline;
- optional ignored/user-run Be or Cr one-center base atom artifact
  write/readback using explicit charge, spin sectors, origin geometry, and
  basis controls;
- finite/symmetric `K`, unit `U_A`, and IDA `V` for ignored/user-run non-H
  atom checks;
- clear `ArgumentError` for translated atom input, missing `d`, noninteger or
  nonpositive charge, nonneutral electron count, or element-table/default
  requests where practical.

No new committed test file, committed non-H atom fixture, public non-H
reference scalar, solver run, supplemented atom endpoint, ECP gate,
translated-atom gate, or driver change is approved by this ID.

## Approved For R3/RG Implementation

The R3 labels remain approved compatibility and endpoint-history IDs. Current
Residual Gaussian algorithm authority lives in
`residual_gaussian_domain_module.md`; do not copy that algorithm into this
registry. This section records approved IDs, source owners, function surfaces,
artifact keys, and validation gates.

### R3 Compatibility Boundary

Approved R3 owner/path for compatibility entry points, artifact writing, and the
non-exported usability facade:

```text
Owner module: CartesianFinalBasisRealization
Source file: src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

R3 compatibility surfaces must delegate current residual-basis, exact-operator,
and residual-interaction physics to the Residual Gaussian module where those
module functions exist. They must not preserve retired global raw-candidate
selection, old density gauges, or duplicate RG algorithms.

### HP-R3-OBJ-01 — residual-GTO augmentation object

Approved historical/object authority for the first H2 residual-GTO endpoint.
Current domain object fields and semantics are recorded under `HP-RG-OBJ-01` and
`residual_gaussian_domain_module.md`. Compatibility names may remain only where
live callers still use them.

### HP-R3-FN-01 — residual-basis construction

Approved historical R3-A basis-construction surface. Current production logic is
`build_residual_gaussian_basis(...)` under `HP-RG-FN-01`.

Binding guardrails: residual basis directions are selected separately on each
physical owner atom; residual occupation is not numerical rank; owner-local
sectors are merged once; global raw-candidate Lowdin and global raw-column
pivoted-Cholesky selection are not approved current algorithms.

### HP-R3-FN-02 — exact augmented one-body and moment assembly

Approved historical R3-A exact-operator surface. Current production logic is
`transform_augmented_operator(...)` under `HP-RG-FN-02`.

The exact transformed operators are kinetic, every uncharged by-center nuclear
attraction, `x`, `y`, `z`, `x^2`, `y^2`, and `z^2`. This exact transformation is
not the MWG approximation and must not be replaced by moment-matched interaction
logic.

### HP-R3-FN-03 — residual MWG/IDA and in-memory Hamiltonian

Approved R3-B compatibility entry point:

```text
pqs_terminal_residual_gto_augmented_hamiltonian(...)
```

Output is the existing `CartesianIDAHamiltonian{Float64}`. Current residual MWG
descriptor and interaction math is owned by `moment_matched_gaussians(...)` and
`assemble_residual_ida_interaction(...)` under `HP-RG-FN-03` and `HP-RG-FN-04`.

The accepted H2 owner-local endpoint has augmented dimension `489` and
lowest-orbital IDA self-Coulomb `0.4574265214362075` within `1.0e-10`. Older
R3-B scalars from global-selection or retired density-gauge diagnostics are
historical evidence only and are not active targets.

### HP-R3-ART-01 — compact supplemented artifact provenance

Approved R3-C internal artifact helper remains outside the RG module. It may
write an existing `CartesianIDAHamiltonian` file and add compact
`supplement_provenance/` provenance. RG does not own artifact writing,
artifact schema, JLD2 workflow, or provenance readers.

Approved source owner/path:

```text
Owner module: CartesianFinalBasisRealization
Source file: src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

Approved `supplement_provenance/` keys:

- `provenance_version`;
- `producer = :cartesian_residual_gto_mwg_augmentation`;
- `supplement_policy = :mwg_residual_gto`;
- `basis_by_center`;
- `lmax`;
- `uncontracted`;
- `width_filtering`;
- `candidate_count`;
- `owner_counts`;
- `base_dimension`;
- `residual_dimension`;
- `augmented_dimension`;
- `augmented_basis_order = :base_then_residual`;
- `residual_basis_convention = :owner_local_residual_occupation_final_merge_lowdin`;
- `rank_rule` / owner-local selection rule;
- `occupation_cutoff = 1.0e-8`;
- `tau_neg_abs`, `tau_neg_rel`;
- `tau_merge_abs = 1.0e-12`, `tau_merge_rel = 1.0e-12`;
- `mwg_convention_version`;
- `mwg_convention = :separable_moment_matched_density_normalized`;
- `one_body_source`;
- `interaction_source = :weight_aware_residual_mwg_ida_blocks`;
- compact validation labels and H2 reference value when supplied.

Do not serialize full residual bases, dense moments, `T_G`, `T_A`, MWG centers,
MWG widths, or broad construction inputs in this compact group.

### HP-R3U-FILE-01 — supplemented workflow source and validation files

Approved non-exported usability owner:

```text
src/cartesian_base_hamiltonian.jl
```

Allowed companion surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

No public export, new source file, new committed test file, driver/bin/tool
workflow, report/status/payload object, or artifact shape beyond
`supplement_provenance/` is approved.

### HP-R3U-FN-01 — non-exported supplemented Hamiltonian facade

Approved internal call shape:

```julia
cartesian_residual_gto_mwg_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    supplement::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

Original first systems were z-axis H2 and z-axis Be2. `HP-R3U-ZDI-FN-01`
relaxes that guard to explicit homonuclear two-center z-axis diatomics. Be2
remains an internal/performance-supported proxy. Cr2 is permitted only as an
explicit generic homonuclear z-axis ignored/user-run stress or usability case
after H2/Be2 validation. Heteronuclear systems, non-z-axis or arbitrary
orientations, charged systems, ECP inputs, solver/RHF workflow, public export,
and Cr2-specific branches remain unapproved.

### HP-R3U-WIRE-01 — base-to-RG same-construction workflow

Approved wiring:

```text
validated system/basis/supplement spec
-> R1-style/base normalization and base stages
-> base CartesianIDAHamiltonian plus same-construction terminal basis/bundles
-> legacy named-basis supplement loading
-> basis_representation(supplement)
-> RG/R3 same-construction augmented Hamiltonian path
-> optional R3-C artifact writer
-> CartesianIDAHamiltonian{Float64}
```

The base Hamiltonian, terminal basis realization, and parent axis bundle must
come from the same base construction call. The facade must not expose terminal
basis realizations, bundles, residual objects, augmented-operator objects, MWG
descriptors, pair factors, or provenance payloads.

### HP-R3-TEST-01 / HP-R3U-TEST-01 — standalone H2 endpoint gate

Approved standalone validation file:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

Invocation:

```text
julia --project=. test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

The gate covers the H2 residual-GTO/MWG endpoint family, including residual
basis checks, exact augmented one-body/moment checks, independent weight-aware
`V_GM` comparison, in-memory Hamiltonian checks, and the non-exported usability
facade/artifact readback section. It is not approved for `test/runtests.jl`,
Be2 committed validation, Cr2 validation, or private route/status assertions.

## Approved For Residual Gaussian Domain Migration

Current Residual Gaussian domain algorithm authority is
`residual_gaussian_domain_module.md`. This registry records the approved module
files and function surfaces only.

The Residual Gaussian module does not own compact supplemented artifact writing
or `supplement_provenance/`. Artifact/facade hooks remain outside RG unless a
later amendment names a real duplication or consumer reason.

### HP-RG-FILE-01 — Residual Gaussian module files

Approved internal module and files:

```text
src/cartesian_residual_gaussians/CartesianResidualGaussians.jl
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_residual_gaussians/augmented_operators.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
```

No public export is approved.

### HP-RG-OBJ-01 — residual Gaussian basis object

Approved domain object: a numerical residual Gaussian basis object carrying base
dimension, candidate count, residual dimension, candidate owner indices,
residual source owner indices, owner retained counts, residual occupations,
cutoff/tolerance policy, selection/orientation/sign rules, and final
`T_G::Matrix{Float64}` / `T_A::Matrix{Float64}` transforms.

It must not be a status/result payload and must not carry route metadata,
report fields, status flags, artifact data, MWG descriptors, or public API
state.

### HP-RG-FN-01 — residual Gaussian basis construction

Approved production name:

```julia
build_residual_gaussian_basis(...)
```

Candidate owner indices are required. The current algorithm is the owner-local
residual occupation construction defined in
`residual_gaussian_domain_module.md`.

### HP-RG-FN-02 — exact augmented operator transformation

Approved production name:

```julia
transform_augmented_operator(...)
```

This owns exact `[G,A] -> [G,R]` transformation for kinetic, uncharged
by-center nuclear attraction, and first/second Cartesian moment matrices.

### HP-RG-FN-03 — moment-matched Gaussian descriptors

Approved production name:

```julia
moment_matched_gaussians(...)
```

Descriptors are computed from final merged residual functions and exact moment
matrices. They are residual-containing interaction descriptors only, not exact
residual-GTO Coulomb integrals.

### HP-RG-FN-04 — residual IDA interaction assembly

Approved production name:

```julia
assemble_residual_ida_interaction(...)
```

This owns residual-containing MWG/IDA blocks `V_GM` and `V_MM`, combined with
unchanged base `V_GG`. `V_GM` uses weight-aware final-basis density
normalization for PQS shell blocks.

### HP-RG-WIRE-01 — migration from terminal residual file

`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` may keep
small compatibility wrappers for live callers and artifact/facade hooks. Moved
physics helpers should be deleted from that file once callers use RG-domain
helpers directly.

### HP-RG-TEST-01 — migration validation

Approved validation is the existing standalone H2 residual-GTO/MWG endpoint
with augmented dimension `489`, self-Coulomb `0.4574265214362075`, exact
one-body/moment checks, independent weight-aware `V_GM` check, and optional
ignored Be2 usability/performance measurement when a source pass changes the
interaction path or facade wiring.

No new committed test file, Be2 committed gate, or Cr2 full
Hamiltonian/artifact/facade validation is approved.

## Approved For Cartesian Gaussian Raw-Block Nuclear Owner

This section approves only the neutral uncharged by-center nuclear raw-block
slice recorded in `cartesian_gaussian_raw_blocks_nuclear.md`. It does not
approve a broad Gaussian raw-block framework.

### HP-CGRB-FILE-01 — neutral Cartesian Gaussian raw-block module files

Approved internal module and files:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

`src/GaussletBases.jl` may add only the internal include needed to load the
module, with include-order changes limited to immediate dependency/caller
needs. No public export is approved.

### HP-CGRB-FN-01 — exact uncharged Gaussian nuclear raw blocks

Approved internal kernel family:

```text
cartesian Gaussian nuclear raw blocks by center
```

The exact Julia names may follow local module style, but they must be neutral:
no `r3`, `residual`, `qw`, route-stage, report, or status vocabulary.

Approved numerical outputs:

- parent-supplement `G-A` matrices by nuclear center;
- supplement-supplement `A-A` matrices by nuclear center;
- uncharged unit attraction convention, `U_A = -1/r_A`.

Approved construction details:

- analytic one-dimensional nuclear factor construction;
- reuse across unique nuclear coordinates;
- upper-triangular `A-A` assembly and mirroring;
- function-local scratch/workspace reuse;
- term-first contraction over the Gaussian expansion.

The kernel must not apply physical nuclear charges, perform terminal
projection, transform into residual bases, create overlap/kinetic/moment
blocks, assemble Hamiltonians, or create persistent caches/bundles.

### HP-CGRB-FN-02 — nuclear one-dimensional axis-family reuse

Approved source owner:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

Approved optimization target: reorganize exact uncharged nuclear raw-block
construction around unique supplement one-dimensional axis families rather
than flattened 3D orbital labels.

Approved concepts:

- unique supplement axis-family inventory independent of flattened 3D orbital
  labels;
- integer map `orbital_axis_family[orbital, axis] -> family_id`;
- unique `G-A` table keys
  `(axis, supplement_axis_family, nuclear_axis_coordinate)`;
- unique `A-A` table keys
  `(axis, canonical(left_family, right_family), nuclear_axis_coordinate)`;
- transpose/orientation flags when canonical `A-A` family order is reversed;
- term-first filling of each required one-dimensional table at most once per
  Coulomb Gaussian term;
- reuse of those tables across all 3D orbitals and orbital pairs that reference
  the same axis families;
- coupled primitive-pair contraction
  `sum_pq c_p c_q Ix[p,q] Iy[p,q] Iz[p,q]`;
- function-local workspaces and integer lookup plans only.

Independent contraction of x/y/z axis tables into separate scalar contractions
is forbidden. The kernel must not introduce persistent caches, metadata,
status/report fields, route objects, payload structs, public API/export,
artifact changes, Residual Gaussian algorithm changes, Qiu-White route
semantic changes, overlap/kinetic/moment migration, Cr2 facade support, or Cr2
artifact workflow.

### HP-CGAI-FN-01 — optional Cartesian Gaussian axis helper

Approved source owner:

```text
src/cartesian_gaussian_axis_integrals.jl
```

Status: superseded as a performance endpoint. It remains an optional helper
surface only if needed by `HP-CGRB-FN-02`.

Optional internal helper concept:

```julia
_cartesian_gaussian_axis_integral_table!(
    destination,
    left_exponents,
    left_centers,
    left_powers,
    left_prefactors,
    right_exponents,
    right_centers,
    right_powers,
    right_prefactors,
    term;
    factor_exponent = 0.0,
    factor_center = 0.0,
)
```

The helper fills an already allocated destination matrix and must return the
same values as `_cartesian_gaussian_axis_integral_table(...)` without
allocating the result matrix. The existing scalar
`_cartesian_gaussian_axis_integral(...)` behavior is unchanged. The allocating
helper may delegate to the in-place helper if that preserves behavior cleanly.
The same owner may also add a specialized nonallocating nuclear-factor scalar
integral for the `:factor` term if needed by the `HP-CGRB-FN-02` family-reuse
kernel.

Allowed consumer surface:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

That file may consume the optional helper only inside the exact uncharged
by-center Gaussian nuclear raw-block construction under `HP-CGRB-FN-02`. Do
not treat result-matrix allocation removal as the accepted Cr2 optimization
target. No public API, export, persistent cache, raw-block payload,
metadata/status/report field, artifact change, route object, Residual Gaussian
algorithm change, Qiu-White route semantic change, overlap/kinetic/moment
migration, Cr2 facade, or Cr2 artifact workflow is approved by this ID.

### HP-CGRB-WIRE-01 — Residual Gaussian and Qiu-White rewiring

Approved caller rewiring surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/ordinary_qw_raw_blocks.jl
src/ordinary_qw_operator_assembly.jl
```

The implementation sequence is binding:

1. Extract current nuclear `G-A`/`A-A` behavior preserving conventions.
2. Rewire Residual Gaussian and Qiu-White callers to the neutral kernel.
3. Delete duplicate route-local nuclear loops once parity is established.
4. Optimize allocation inside the neutral owner only after extraction parity.

No Qiu-White route objects, Residual Gaussian selection logic, augmented
operator transforms, terminal projection, parent construction, artifact
workflow, report/status/payload fields, public API, or Cr2 facade/artifact
workflow may be added under this ID.

### HP-CGRB-TEST-01 — nuclear extraction validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged, as ignored validation if needed;
- Cr2 q4 exact nuclear `G-A`/`A-A` blocks match the current implementation at
  roundoff, as ignored measurement only;
- one small Qiu-White nuclear parity fixture.

Approved committed standalone parity file, if needed:

```text
test/nested/cartesian_gaussian_raw_blocks_nuclear_runtests.jl
```

This file must not be added to `test/runtests.jl` without a later amendment.
It must not validate Cr2 workflow, artifacts, public API, route statuses,
report mirrors, payload fields, or Residual Gaussian internals.

## Approved For Cartesian Gaussian Raw-Block Non-Nuclear Owner

This section approves only the non-nuclear raw-block slice recorded in
`cartesian_gaussian_raw_blocks_non_nuclear.md`. It extends the existing neutral
Cartesian Gaussian raw-block owner under new `HP-CGRB-NN-*` IDs. It must not be
implemented under `HP-CGRB-FN-02`.

### HP-CGRB-NN-FILE-01 — non-nuclear raw-block file

Approved owner file:

```text
src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl
```

Allowed module plumbing:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
```

Only the include needed to load `non_nuclear_blocks.jl` is approved there.
Root include changes in `src/GaussletBases.jl` are not expected and are not
approved unless a later amendment identifies a real include-order blocker.

### HP-CGRB-NN-FN-01 — exact non-nuclear Gaussian raw blocks

Approved internal kernel family:

```text
cartesian Gaussian non-nuclear raw blocks
```

The exact Julia names may follow local module style, but they must be neutral:
no `r3`, `residual`, `qw`, route-stage, report, or status vocabulary.

Approved numerical outputs:

- parent-supplement `G-A` overlap, kinetic, coordinate moments, and second
  moments;
- supplement-supplement `A-A` overlap, kinetic, coordinate moments, and second
  moments.

Approved construction details:

- analytic one-dimensional table construction;
- unique supplement axis-family reuse;
- canonical `A-A` family-pair table keys and orientation handling;
- upper-triangular `A-A` assembly and mirroring;
- function-local scratch/workspace reuse;
- coupled product-axis contraction preserving existing Qiu-White values;
- reuse of once-built overlap `G-A` for residual setup mixed overlap
  `X = G' S A` and exact augmented-operator assembly when both are built in
  the same local construction call.

The kernel may return a compact fixed-field internal result containing only the
approved raw matrices. It must not be a status object, route stage, report
payload, metadata carrier, persistent cache, broad provider bundle, or artifact
data.

### HP-CGRB-NN-WIRE-01 — Residual Gaussian and Qiu-White non-nuclear rewiring

Approved caller rewiring surfaces:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/ordinary_qw_raw_blocks.jl
src/ordinary_qw_operator_assembly.jl
```

The implementation sequence is binding:

1. Extract current Qiu-White non-nuclear `G-A`/`A-A` behavior preserving
   conventions.
2. Rewire Residual Gaussian exact-operator construction and residual mixed
   overlap setup to consume the neutral output.
3. Rewire Qiu-White consumers to the neutral output.
4. Delete duplicate route-local non-nuclear loops once parity is established.
5. Optimize allocation inside the neutral owner only after extraction parity.

No nuclear raw-block changes, final-basis `G-G` product-matrix optimization,
terminal projection, Residual Gaussian algorithm changes, augmented-operator
transform changes, Qiu-White semantic changes, Qiu-White route objects, parent
construction, persistent cache, report/status/payload fields, public API,
artifact workflow, or Cr2 facade/artifact workflow may be added under this ID.

### HP-CGRB-NN-TEST-01 — non-nuclear extraction validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged, as ignored validation if needed;
- Cr2 q4 non-nuclear `G-A`/`A-A` overlap, kinetic, coordinate moment, and
  second-moment blocks match the current implementation at roundoff, as
  ignored measurement only;
- residual setup mixed overlap `X` matches the current construction at
  roundoff;
- one small Qiu-White non-nuclear parity fixture.

Approved committed standalone parity file, if needed:

```text
test/nested/cartesian_gaussian_raw_blocks_non_nuclear_runtests.jl
```

This file must not be added to `test/runtests.jl` without a later amendment.
It must not validate Cr2 workflow, artifacts, public API, route statuses,
report mirrors, payload fields, or Residual Gaussian internals.

## Approved For R3 Terminal G-G Product Matrices

This section approves only the terminal final-basis `G-G` product-matrix
optimization recorded in `r3_terminal_gg_product_matrices.md`. It is owned by
`CartesianFinalBasisRealization`, not by `CartesianGaussianRawBlocks`.

### HP-R3GG-FN-01 — R3/RG terminal G-G product-matrix optimization

Approved owner:

```text
Owner module: CartesianFinalBasisRealization
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
```

The first implementation should prefer editing only
`pqs_terminal_residual_gto.jl`. Edits to `pqs_terminal_one_body.jl` are
approved only for a small internal terminal-product workspace or multi-product
helper needed to reuse function-local buffers across consecutive product
assemblies.

Approved product matrices:

- kinetic `K_GG`;
- coordinate moments `x_GG`, `y_GG`, `z_GG`;
- second moments `x2_GG`, `y2_GG`, `z2_GG`.

Approved implementation shapes:

- accumulate the three kinetic-axis product contributions into one destination;
- reuse an already constructed base Hamiltonian kinetic `G-G` block in the
  same-construction path when available and validated equal;
- build coordinate and second-moment `G-G` products one axis at a time and
  transform immediately;
- share function-local scratch/workspace across consecutive terminal product
  assemblies;
- delete or simplify `_r3a_product_matrix(...)` when replaced and no live
  caller remains.

This ID does not approve `G-A`/`A-A` raw-block changes, nuclear raw-block
changes, unit-nuclear `U_A` Gaussian-sum changes, terminal basis realization
changes, residual Gaussian algorithm changes, Qiu-White semantic changes,
IDA/MWG changes, parent construction, persistent caches, metadata,
report/status/payload fields, public API/export, artifact changes, Cr2 facade
support, or Cr2 artifact workflow.

Line budget: at most `100` added `src` lines total. If implementation needs a
broad product-operator framework, persistent workspace/cache object, files
outside the approved source files, or a public/internal payload, stop and
request a new docs-only amendment.

### HP-R3GG-TEST-01 — terminal G-G product validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian usability/performance measurement unchanged except for
  allowed timing/allocation improvement;
- Cr2 q4 `K_GG`, coordinate moment `G-G`, and second-moment `G-G` products
  match the current construction at roundoff as ignored validation;
- augmented exact operators remain finite and symmetric;
- base `G-G` block equality checks in the existing H2 endpoint still pass;
- Cr2 q4 exact-operator allocation is remeasured after parity.

No new committed test file is approved by this ID.

## Approved For R3 Unit-Nuclear U_GG Gaussian Sum

This section approves only the terminal final-basis unit-nuclear `U_GG`
Gaussian-sum optimization recorded in
`r3_unit_nuclear_ugg_gaussian_sum.md`. It is owned by
`CartesianFinalBasisRealization`, not by `CartesianGaussianRawBlocks`.

### HP-R3UN-FN-01 — R3/RG unit-nuclear U_GG Gaussian-sum optimization

Approved owner:

```text
Owner module: CartesianFinalBasisRealization
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

The first implementation should prefer `pqs_terminal_one_body.jl`. Edits to
`pqs_terminal_residual_gto.jl` are approved only for narrow R3 exact-operator
caller wiring needed to use function-local scratch or the optimized helper.

Approved target functions:

```text
_accumulate_terminal_gaussian_sum!
_terminal_gaussian_sum_action
```

Approved implementation shapes:

- reuse function-local scratch/workspace across Gaussian-sum terms and center
  calls;
- accumulate terminal Gaussian-sum contributions in-place into the caller's
  destination;
- reduce avoidable allocation in factor lookup and terminal Gaussian-sum action
  construction;
- add small internal scratch arguments or file-local helpers only if they remain
  inside `CartesianFinalBasisRealization` and create no persistent state;
- simplify or delete allocation-heavy helper code inside the targeted
  Gaussian-sum path after parity.

This ID does not approve neutral raw-block changes, terminal kinetic/moment
`G-G` product changes, residual Gaussian selection/orientation/transform
changes, MWG/IDA changes, Qiu-White semantic changes, route/stage setup
cleanup, raw-block setup cleanup, parent construction, terminal basis
realization changes, persistent caches/workspaces, broad Gaussian-sum
frameworks, metadata/report/status/payload fields, artifacts, public
API/export, Cr2 facade support, or Cr2 artifact workflow.

Line budget: at most `100` added `src` lines total. If implementation needs a
broad Gaussian-sum framework, persistent cache/workspace object, files outside
the approved source files, or source edits outside the terminal unit-nuclear
`U_GG` path, stop and request a new docs-only amendment.

### HP-R3UN-TEST-01 — unit-nuclear U_GG validation

Approved validation:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian facade/readback unchanged except for allowed
  timing/allocation improvement;
- Cr2 q4 exact-operator audit reports before/after unit-nuclear `U_GG`
  allocation and total wrapper allocation;
- Cr2 q4 unit-nuclear `U_GG` block replay parity and final exact augmented
  operator parity at roundoff;
- exact operators remain finite and symmetric.

No new committed test file is approved by this ID.

## Approved For Canonical Cartesian Driver Usability

This section approves only the compact artifact-producing canonical driver
workflow recorded in `cartesian_driver_usability_workflow.md`. It is workflow
authority over approved producer surfaces, not algorithm, kernel, solver,
artifact-schema, or diagnostic authority.

### HP-DRV-FILE-01 — canonical driver file

Approved file:

```text
bin/cartesian_ham_builder.jl
```

No other `bin`, `tools`, `src`, `test`, or committed driver-input fixture file
is approved by this ID.

### HP-DRV-FN-01 — compact functional driver workflow

Approved invocation shape:

```text
julia --project=. bin/cartesian_ham_builder.jl [input.jl] [key=value ...]
```

Approved behavior:

- visible editable defaults near the top of the driver;
- optional trusted local Julia input file for project-specific defaults;
- later command-line `key=value` overrides;
- visible public `system`, `basis`, and optional `supplement` contract
  construction before calling an approved facade;
- compact normalized run summary;
- coarse user-facing phase timing;
- visible physics-level construction stages through the staged producer surface;
- base or supported supplemented Hamiltonian construction through approved
  producer surfaces;
- artifact write;
- optional readback check.

Approved configuration concepts are `basisname`, `system`, base `basis`,
optional `supplement`, `hamfile`, `padding`, `check_file`, `print_contract`,
`print_timing`, and `expected_dimension`.

Compact summary printing and artifact readback checks remain allowed workflow
behavior, but they are not open-ended hooks and must not introduce route,
diagnostic, artifact-schema, or solver controls.

`basisname = nothing` selects base mode. `basisname !== nothing` selects the
supported supplemented diatomic mode and is the visible supplement basis label;
that path must reject `Natom == 1`. Supplemented atoms remain unapproved.

`padding` is a public physical box-padding control. For one-center atoms it
maps to the base facade `radius`. For z-axis diatomics it maps to the existing
facade extents as padding around the two nuclei; under the current origin-based
z-axis contract this means `xmax_parallel = max(abs(z_i)) + padding` and
`xmax_transverse = padding`.

Approved hooks are only `check_file`, `print_contract`, `print_timing`, and
`expected_dimension`. They may support human expert review and
Codex-controlled artifact checks. They must not expose route internals,
stop-after stages, raw-block switches, allocation probes, artifact schema
dumps, solver controls, Cr2-specific workflow, or private helper calls.
`check_file` may contain compact public contract facts, artifact path, final
dimension, expected-dimension result, readback deltas, and coarse timing only.

This ID does not approve private route-stage controls, stop-after internals,
ladder probes, stage markers, fixture hacks, diagnostic knobs, underscored
package helper calls from the driver, raw-block provider switches,
report/status/payload dumps, metadata field clouds, allocation probes,
benchmark harness behavior, solver/RHF/ECP/EGOI/HamV6 workflow, public
API/export changes, artifact schema changes or dumps, committed test files,
committed driver-input fixtures, supplemented atoms, or Cr2-specific workflow
support. Generic explicit homonuclear z-axis Cr2 stress through
`HP-R3U-ZDI-WIRE-01` is separate ignored/user-run validation authority, not
driver-owned Cr2 support.

Line budget: at most `150` added `bin` lines. If implementation needs a parser
framework, source files outside the approved driver and staged producer
surfaces, committed input fixtures, route-stage diagnostics,
status/report/payload expansion, artifact schema changes, or Cr2-specific
workflow support, stop and request a new docs-only amendment.

### HP-DRV-STAGE-FN-01 — visible physics-stage producer surface

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved purpose: expose a small non-exported, non-underscored,
driver-facing staged producer surface so the canonical driver can execute and
time visible physics-level construction stages without calling private
underscored helpers.

Approved visible stages are:

- construct public `system`, `basis`, and optional `supplement`;
- build the base working basis / terminal realization;
- assemble the base Hamiltonian;
- load or build the Gaussian supplement basis when `basisname !== nothing`;
- build residual Gaussian augmentation;
- build exact augmented operators;
- assemble the supplemented Hamiltonian;
- write and check the artifact.

The first and last stages remain driver/writer responsibilities. This ID
approves source factoring needed for the base working-basis/terminal
realization, base-Hamiltonian assembly, Gaussian supplement, residual
augmentation, exact augmented operators, and supplemented-Hamiltonian assembly
stages.

The staged surface may factor the existing `cartesian_base_hamiltonian(...)`
and `cartesian_residual_gto_mwg_hamiltonian(...)` bodies so that those facades
can remain wrappers over the same implementation. It may return existing
domain objects and small fixed-key ephemeral stage products required by the
next approved stage.

The staged surface must be a set of separate named construction-stage
functions. It must not be a single opaque replacement wrapper that hides the
same construction sequence under a new name. The canonical driver must be able
to bind visible local variables for the base realization, base Hamiltonian,
supplement basis, residual augmentation, exact augmented operators, and final
Hamiltonian assembly.

This ID does not approve public exports, public API redesign, route-stage
objects, reports, status/result payloads, metadata field clouds, runtime-keyed
field groups, persistent caches, raw-block switches, solver/ECP workflow,
artifact schema changes, or source files outside
`src/cartesian_base_hamiltonian.jl`.

Line budget: at most `150` added `src` lines. If the staged surface requires a
new module, source files outside `src/cartesian_base_hamiltonian.jl`, a broad
payload object, committed tests, or new artifact keys, stop and request a new
docs-only amendment.

### HP-DRV-STAGE-WIRE-01 — canonical driver staged wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

The canonical driver may call the `HP-DRV-STAGE-FN-01` staged producer surface
as separate top-level stage calls and assign local variables using the approved
physics-stage names. It may print or record coarse user-facing timings for
those stages. Replacing the current facade call with one all-in-one staged
wrapper call is not approved for the canonical driver.

This ID does not approve calls from the driver to underscored package helpers,
old route stages such as `cartesian_parent`, `cartesian_shells`,
`cartesian_units`, `cartesian_pair_terms`, or `cartesian_assembly`, raw-block
provider switches, stop-after controls, route diagnostics, allocation probes,
artifact schema dumps, solver controls, Cr2-specific workflow, or new
committed fixtures/tests.

### HP-DRV-STAGE-TEST-01 — staged driver validation

Approved validation:

- package load;
- H atom or H2 base driver run with visible base-stage timing/summary;
- H2 supplemented driver run with visible supplement/residual/operator/
  Hamiltonian stage timing/summary;
- artifact write/readback still passes for those runs;
- `expected_dimension`, `print_contract`, and `check_file` behavior still
  uses only public contract and coarse stage facts.

No committed test file, committed input fixture, Cr2-specific driver run,
solver run, or diagnostic harness is approved by this ID.

### HP-DRV-TEST-01 — driver workflow validation

Approved validation:

- package load;
- public contract construction and optional `print_contract`/`check_file`
  output for at least one base run when driver construction code changes;
- H2 base driver run writes a `CartesianIDAHamiltonian` artifact and optional
  readback passes;
- H2 supplemented driver run writes a supplemented `CartesianIDAHamiltonian`
  artifact with approved compact `supplement_provenance/` and optional readback
  passes;
- optional ignored Be2 usability run if the implementation touches
  supplemented mode.

Validation input files, if needed, must be ignored `tmp/work` files. No
committed test file, committed driver-input fixture, Cr2-specific driver run,
or solver run is approved by this ID.

## Approved For Canonical Driver Atom Workflow

This section approves only the base atom workflow recorded in
`cartesian_driver_atom_workflow.md`. It is driver authority over the existing
base facade, not new atom physics, Residual Gaussian, artifact-schema, or
solver authority.

### HP-DRV-ATOM-FN-01 — explicit base atom driver workflow

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- accept explicit one-center atom input in `mode = :base`;
- require `atom_symbols`, `nuclear_charges`, `atom_locations`, `nup`, and
  `ndn`;
- require exactly one center at `(0.0, 0.0, 0.0)`;
- require finite positive explicit nuclear charge and neutral all-electron
  count `nup + ndn == round(Int, only(nuclear_charges))`;
- pass explicit one-center base `basis` fields, including required public `d`,
  to the existing base facade;
- use clear `ArgumentError`s for unsupported atom workflow inputs where
  practical.

Current driver validation remains origin-centered H. This driver ID does not
approve changing `src/cartesian_base_hamiltonian.jl`; producer-side
one-center atom support is governed separately by `HP-R1-ATOM-*`.

### HP-DRV-ATOM-WIRE-01 — driver atom-to-base-facade wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- the canonical driver may call
  `cartesian_base_hamiltonian(system; basis, hamfile)` for one-center base
  atom construction;
- artifact write/readback uses the existing base facade and existing
  `producer_provenance/` schema;
- no package-internal route-stage helper, terminal basis object, raw-block
  provider, report/status/payload object, or new artifact field is approved.

Supplemented atom Hamiltonians are not approved. If the requested atom is
outside the existing base facade support, the implementation must stop at a
clear unsupported-input error rather than adding broader atom construction.

Line budget for `HP-DRV-ATOM-FN-01` plus `HP-DRV-ATOM-WIRE-01`: at most `80`
added `bin` lines, with no committed test, tool, or input-fixture file.

### HP-DRV-ATOM-TEST-01 — base atom driver validation

Approved validation:

- package load;
- origin-centered H base driver artifact write/readback with explicit system
  and one-center basis fields;
- optional ignored negative checks for non-origin atom input, nonneutral
  electron count, missing `d`, or unsupported atom input.

No supplemented atom endpoint, translated-atom gate, committed atom fixture,
new committed test file, solver run, artifact-schema validation, or broader
base-atom validation is approved by this ID.

## Approved For Homonuclear Z-Axis Diatomic Supplemented Workflow

This section approves only the molecule-scope relaxation recorded in
`r3_homonuclear_diatomic_supplemented_workflow.md`. It is generic
homonuclear z-axis diatomic authority, not element-specific Cr2 authority.

### HP-R3U-ZDI-FN-01 — homonuclear z-axis diatomic supplemented facade

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- replace hardcoded H/Be supplemented guards with explicit homonuclear z-axis
  diatomic validation;
- require explicit atom symbols, nuclear charges, `nup`, `ndn`, geometry, base
  basis parameters, supplement basis labels, and optional supplement
  `basisfile`;
- support exactly two equal-symbol/equal-charge centers on the Cartesian
  z-axis with distinct finite `z` coordinates;
- require neutral all-electron count
  `nup + ndn == round(Int, sum(nuclear_charges))`;
- throw clear `ArgumentError`s for unsupported systems before expensive
  construction where practical.

This ID does not approve heteronuclear systems, non-z-axis/general orientation,
charged systems, ECP, solver/RHF workflow, public API/export redesign,
artifact schema changes, route diagnostics, metadata/status/report fields, or
Cr2-specific branches/defaults/fixtures.

### HP-R3U-ZDI-WIRE-01 — canonical driver supplemented-mode wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- canonical driver `:supplemented` mode may call the supported
  `cartesian_residual_gto_mwg_hamiltonian(...)` facade;
- driver inputs may carry explicit homonuclear z-axis diatomic system, base
  basis, supplement labels, optional `basisfile`, and `hamfile`;
- Cr2 may be an ignored/user-run stress or usability case through the generic
  path only after H2/Be2 validation.

This ID does not approve package-internal helper composition from the driver,
Cr2-specific workflow, committed Cr2 fixtures, route diagnostics, artifact
schema changes, public exports, solver workflow, or broad driver feature
growth.

Line budget for `HP-R3U-ZDI-FN-01` plus `HP-R3U-ZDI-WIRE-01`: at most `100`
added `src`/`bin` lines total, with net simplification expected where
H/Be-specific checks are removed.

### HP-R3U-ZDI-TEST-01 — homonuclear diatomic validation

Approved validation:

- H2 supplemented facade/driver artifact path remains unchanged;
- Be2 supplemented facade/driver artifact path remains unchanged and acts as
  the non-H correctness/performance gate;
- optional ignored/user-run Cr2 stress or usability run after H2/Be2 pass.

No committed Cr2 fixture, committed Cr2 test, new committed test file,
heteronuclear gate, non-z-axis gate, solver run, or artifact schema validation
is approved by this ID.

## Approved Measurement-Only Authority

These entries authorize ignored measurement/probe work only. They do not
authorize production source edits, committed tests, source files, persistent
objects, metadata/report/status/payload fields, artifacts, public API, or Cr2
workflow support.

### HP-R3REM-AUDIT-01 — remaining exact-operator allocation audit

Approved scope:

- measure the Cr2 q4 R3/RG exact augmented-operator allocation remaining after
  `954c86cd` and the terminal `G-G` product-workspace optimization;
- separate total wrapper allocation from neutral raw-block construction,
  terminal `G-G` product buffers, unit-nuclear `U_GG` Gaussian-sum
  construction, exact augmented nuclear transforms, route/stage setup, and
  audit/replay overhead;
- use ignored `tmp/work` probes only, with H2/Be2 sanity if needed.

Required outcome:

- classify the dominant remaining allocation bucket;
- recommend a future source lane only if the owner, files, functions,
  forbidden surfaces, validation gates, line budget, deletion/simplification
  expectation, and failure rule are specific enough for a separate docs-only
  amendment.

This ID does not extend `HP-R3GG-FN-01` and does not approve unit-nuclear
`U_GG` Gaussian-sum optimization, route/raw-block setup cleanup, final-basis
`G-G` changes, `G-A`/`A-A` raw-block changes, residual Gaussian algorithm
changes, IDA/MWG changes, Qiu-White semantic changes, parent or terminal-basis
changes, persistent caches/workspaces, artifacts, public API/export, Cr2 facade
support, or Cr2 artifact workflow.

## Rejected Or Deferred

### HP-RES-01 — terminal basis build result — rejected

Do not introduce a persistent terminal-basis result wrapper. The realizer
returns `CartesianTerminalBasisRealization` on success.

### HP-CHANGE-01 — return shell overlap from existing shell plan — rejected/deferred

This can be a helper detail under HP-FN-00, but it is not standalone authority.

### HP-OBJ-03 — generic build-result wrapper — rejected

Do not introduce `CartesianHamiltonianBuildResult`, another payload, or a broad
status wrapper around `CartesianIDAHamiltonian`.

### HP-TEST-01 — new committed terminal smoke — rejected

No new committed terminal smoke/probe is approved. Use existing smokes or
ignored `tmp/work` validation unless a later design explicitly approves a test.

### R3 Be2/Cr2 readiness guardrail — deferred

R3-A/B/C are implemented for the narrow H2 residual-GTO/MWG endpoint and
compact supplemented artifact. The following items are closed for that narrow
path and should not be listed as future blockers by default:

- same-construction internal path for the accepted H2 R3 construction;
- deterministic rank-deficient handling in the legacy global-selection
  implementation, now superseded by the approved owner-local selection source
  correction;
- one-shot parent-by-supplement analytic exact-block organization for R3-A
  mixed/self blocks, avoiding repeated CPB-per-terminal-block construction on
  the Be2 proxy;
- independent weight-aware final-basis `V_GM` validation for R3-B;
- compact `supplement_provenance/` artifact group for R3-C.

Do not present Cr2 or broader residual-GTO/MWG supplement support as approved
until a later docs-only amendment chooses and closes the next lane. The
non-exported H2/Be2 R3 usability facade is approved separately by
`HP-R3U-FILE-01`, `HP-R3U-FN-01`, `HP-R3U-WIRE-01`, and `HP-R3U-TEST-01`.
Remaining deferred lanes are:

- implementation and validation of the approved owner-local source correction,
  including the updated H2 scalar and no full Cr2 Hamiltonian;
- public export, public examples, or driver workflow beyond the non-exported
  R3 usability facade;
- a Cr2-readiness lane for measurement-only candidate/rank/memory forecasting,
  with no full Cr2 Hamiltonian yet;
- a basis/supplement-realism lane for validated supplement choices, basis
  labels, and filtering policy beyond the first H2 fixture;
- bounded or streamed residual MWG term storage if higher residual rank makes
  the current dense residual term storage costly;
- allocation-free or bounded-allocation validation reductions for large dense
  matrices.
