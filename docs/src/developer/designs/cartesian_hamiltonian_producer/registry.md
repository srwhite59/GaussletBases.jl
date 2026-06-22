# Registry

Only entries marked approved/implemented authorize source work. Candidate or
rejected entries do not authorize implementation.

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

## Candidate, Not Approved

Candidate entries record proposed design surfaces for review. They do not
authorize source, test, tool, driver, public API, or artifact changes until
explicitly moved into an approved section.

### HP-R3-OBJ-01 — augmented residual-GTO/MWG basis object — candidate

Candidate scope for residual-GTO/MWG augmentation only. A future object may
describe fixed base final basis `G`, Gaussian supplement candidates `A`,
residual block `R`, residual metric facts, and augmented dimension. It must not
rotate or re-Lowdin the base PQS basis, carry route-global matrices in
metadata, or introduce broad status/report payloads.

### HP-R3-FN-01 — residual-GTO/MWG augmentation boundary — candidate

Candidate conceptual boundary:

```julia
augment_cartesian_ida_hamiltonian_with_residual_gto_mwg(
    terminal_basis_realization::CartesianTerminalBasisRealization,
    base_hamiltonian::CartesianIDAHamiltonian{Float64},
    supplement::CartesianGaussianShellSupplementRepresentation3D,
    provider_blocks,
    residual_options,
)::CartesianIDAHamiltonian{Float64}
```

The final name, owner file, and argument types are not approved. The candidate
construction extends the fixed localized final basis with residual GTO/MWG
functions before final supplemented `K`/`U_A`/`V` assembly. It is not a
post-hoc wrapper around an opaque base Hamiltonian.

### HP-R3-FN-02 — augmented one-body assembly — candidate

Candidate scope: assemble augmented kinetic and uncharged by-center nuclear
matrices in basis `[G, R]`. Reuse base `G-G` blocks only when exactly
equivalent, and build mixed/residual blocks from explicit GTO/provider kernels.
No route-global metadata matrices, report clouds, or status-only preflight
graphs are candidate authority.

### HP-R3-FN-03 — residual MWG/IDA assembly — candidate

Candidate scope: assemble supplemented localized IDA terms using the base
`G-G` block plus residual-GTO/MWG descriptors for residual-containing terms.
The base PQS positive final-weight gauge must not be applied blindly to
residual-GTO/MWG integrals or moments; near-zero or negative residual
contributions are allowed where the residual construction or physical operator
permits them.

### HP-R3-ART-01 — supplemented artifact provenance — candidate

Candidate artifact extension only. The in-memory numerical object remains the
existing `CartesianIDAHamiltonian{Float64}` when its matrices satisfy the
current contract. A supplemented artifact may add fixed
`supplement_provenance/` keys exactly as defined in
`r3_residual_gto_mwg_augmentation.md`. This candidate does not approve a
Hamiltonian wrapper, separate manifest, public provenance reader, HamV6 export,
solver export, or consumer API.

### HP-R3-TEST-01 — first supplemented endpoint validation — candidate

Candidate first proxy: z-axis H2 with the existing base H2 public geometry and
a small H-centered Cartesian GTO supplement such as H/cc-pVTZ with `lmax = 1`.
Be2 should follow as the first performance/realism proxy. Cr2 remains a later
stress and consumer-readiness milestone, not the first R3 correctness gate.

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
