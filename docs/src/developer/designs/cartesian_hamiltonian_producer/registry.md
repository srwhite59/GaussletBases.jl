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

Under `HP-HAM-MANIFEST-SRC-FN-01`, a source pass may add one optional compact
terminal source-mode provenance carrier if the implementation chooses to attach
the manifest seam to this realization object. That carrier is artifact
provenance only and must not become a basis, operator, route-report, or
algorithm input.

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

## Approved For White-Lindsey Terminal Basis Implementation

This section approves only the narrow terminal-basis seam recorded in
`white_lindsey_terminal_basis_realization.md`.

### HP-WLTERM-FILE-01 — optional WL terminal realization file

Approved source files:

```text
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

This ID is optional. It approves creating a small WL-specific terminal
realization sibling and adding its include to `CartesianFinalBasisRealization`
only if extending `pqs_terminal_basis_realization.jl` directly would obscure
the distinct WL boundary-stratum construction. No public export, root include,
new module, new basis object, artifact, report, or status/payload object is
approved.

### HP-WLTERM-FN-01 — WL low-order terminal basis realization

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
```

Approved target: realize terminal final-basis blocks for the existing
`:white_lindsey_low_order` route family and return the existing
`CartesianTerminalBasisRealization`.

Allowed behavior:

- support direct identity terminal blocks on their authoritative owned rows;
- realize WL boundary-stratum/product terminal blocks from existing native
  terminal support, retained-rule, and transform records;
- use only support-local coefficients on `support.support_indices` /
  `support.support_states`;
- preserve deterministic terminal support, lowering, retained-record, and
  transform-contract order;
- validate disjoint owned supports and block-local identity overlaps under the
  same structural terminal-basis policy as PQS.

This ID does not approve old WL H1/H1+J materialization adaptation, new
Hamiltonian objects, new route-stage objects, route reports, status/result
payloads, diagnostics, artifact/schema changes, public API/export changes,
raw-block changes, Residual Gaussian/MWG/IDA changes, terminal
shellification-policy changes, retained-selection-policy changes, route
skeleton construction changes, or source files outside the approved surfaces.

Failure rule: if WL boundary-stratum final basis cannot be materialized from
existing terminal lowering, retained-unit, and transform records without
broader route redesign, make no source commit and report the exact missing
native fact.

### HP-WLTERM-WIRE-01 — WL route helper terminal-basis wiring

Approved source file:

```text
src/pqs_source_box_route_driver_helpers.jl
```

Approved behavior:

- remove or narrow the PQS-only terminal-basis guard so
  `route_family = :white_lindsey_low_order` can call the approved WL terminal
  realizer when native terminal records are available;
- keep `route_family = :pqs_source_box` behavior unchanged;
- keep route skeleton construction semantics, route recipe semantics,
  shellification behavior, terminal lowering order, retained-rule order,
  public driver contract, and artifact schema unchanged;
- reject or return a clear unsupported route error when the WL route lacks the
  native facts required for terminal-basis realization.

This ID does not approve adapting old WL materialization, broad route redesign,
supplemented WL behavior, driver input changes, diagnostics/status/report
expansion, raw-block switches, stop-after controls, or source files outside the
approved route helper file.

### HP-WLTERM-TEST-01 — WL terminal-basis validation

Approved validation:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` atom or H2 base artifact/readback remains
  unchanged;
- `nesting = :wl` base atom artifact/readback;
- `nesting = :wl` base H2 artifact/readback if the existing WL diatomic route
  has sufficient native terminal records;
- H2 residual-GTO/MWG PQS endpoint remains unchanged if terminal realization
  code is touched;
- clear unsupported-input/blocker report if WL H2 cannot be realized from
  current native records.

No Cr2 run, supplemented WL run, committed fixture, committed test file,
solver/RHF/ECP/EGOI workflow, artifact schema validation, or broad WL workflow
validation is approved.

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
- `core_spacing` is the authoritative physical near-nucleus spacing for each
  center;
- public `d` is deprecated and is not part of the durable producer contract;
- one-center H maps `reference_spacing` and `core_spacing` separately:
  `spacing_inputs.reference_spacing = basis.reference_spacing`,
  `parent_inputs.parent_mapping_rule = :white_lindsey_atomic_mapping`, and
  `parent_inputs.parent_mapping_d = basis.core_spacing`;
- the reviewed H baseline uses explicit public `core_spacing = 0.3` and
  `reference_spacing = 1.0`;
- if a temporary compatibility path accepts public `d`, it must require
  `d == resolved core_spacing`; z-axis H2 and durable public examples must not
  use `d`;
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

### HP-R1-CORE-FN-01 — unified core-spacing producer contract

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- `core_spacing` is the single public physical near-nucleus spacing for each
  center after explicit input or preset resolution;
- public `d` is deprecated as a producer field;
- a temporary compatibility acceptance of `d` may exist only when
  `d == resolved core_spacing`; mismatches must throw `ArgumentError`;
- one-center White-Lindsey atom wiring sets
  `parent_inputs.parent_mapping_d = resolved core_spacing`;
- the White-Lindsey `Z` dependence is an internal mapping-shape rule:
  `core_range = sqrt(core_spacing / Z)` and
  `mapping_strength = sqrt(core_spacing * Z)`;
- multi-center mappings use the same per-center resolved-core-spacing model
  before applying the combined/neighbor mapping effects;
- future automatic presets may derive `core_spacing = core_scale(q or n_s) / Z`
  or an equivalent fixed `core_spacing * Z` family, but once resolved,
  `core_spacing` remains the authoritative scale.
- canonical driver or project-input defaults may choose visible editable values
  such as `core_spacing = 0.3`; these are explicit resolved inputs, not hidden
  universal producer defaults, and normal overrides such as quick-test
  `core_spacing = 0.5` remain allowed.
- routine correctness tests may override driver physics defaults, but any
  asserted scalar must be tied to the exact test input and not described as a
  physics-default result.

This ID does not approve a public `d`, public `parent_mapping_d`, public
mapping-strength/range knobs, element-table defaults, ECP behavior, solver
workflow, artifact schema changes, or source files outside
`src/cartesian_base_hamiltonian.jl`. `reference_spacing`, `tail_spacing`, and
physical box padding remain separate concepts.

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

### HP-ROUTE-RECIPE-FN-01 — family-selective route recipe cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- `cartesian_recipe(route_inputs)` may construct only the subrecipe selected by
  `route_inputs.route_family`;
- for `route_family = :pqs_source_box`, route inputs must not require inactive
  `white_lindsey_*` fields; the produced recipe may set the inactive
  `white_lindsey` subrecipe to `nothing` while retaining the existing field
  name for caller compatibility;
- for `route_family = :white_lindsey_low_order`, explicit White-Lindsey route
  support must be preserved and the selected `white_lindsey` subrecipe must
  continue to be built from the existing WL route fields; the inactive
  `source_box` subrecipe may be `nothing` if no live WL caller requires it;
- `_cartesian_base_route(kind)` in `src/cartesian_base_hamiltonian.jl` may
  remove unused `white_lindsey_*` fields because the live base producer route
  uses `route_family = :pqs_source_box`;
- existing precomposed recipes that already provide `source_box` and
  `white_lindsey` fields may remain accepted if that compatibility path is
  still live, but it must not force new PQS-only route inputs to carry inactive
  WL vocabulary.

This ID preserves real WL/PQS algorithm differences while removing inactive WL
route-family fields from the current PQS base producer contract. It does not
approve canonical-driver changes, numerical kernel changes, terminal lowering
policy changes, shellification behavior changes, materialization or artifact
schema changes, route-stage diagnostics, status/report expansion, deletion of
WL materialization, or source files outside the two approved files.

Line budget: at most `80` added `src` lines, with net simplification expected.

Failure rule: if `cartesian_recipe(...)` cannot be made family-selective
without broader route-driver, report, materialization, or stage-object changes,
make no source commit and report the blocker.

### HP-ROUTE-RECIPE-TEST-01 — route recipe cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H atom/base artifact readback;
- H2 base artifact readback;
- compact H2 supplemented facade or driver path;
- focused route recipe smoke for explicit `:white_lindsey_low_order` if still
  practical, or a report of the exact live test/tool callers that block further
  WL route-input cleanup.

Existing committed tests may be adjusted only where they directly construct
route inputs that now no longer need inactive family fields. Known direct
route-recipe tests such as
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` may drop inactive
WL fields if required by the source cleanup. No new committed test file, Cr2
run, driver workflow change, route diagnostic, or physics-reference scalar is
approved by this ID.

### HP-ROUTE-INV-FN-01 — retained-unit route inventory type cleanup

Approved source file:

```text
src/pqs_source_box_route_driver_helpers.jl
```

Approved cleanup targets:

- `_pqs_source_box_route_driver_named_tuple_from_units(...)`;
- runtime-keyed retained-unit inventory fields derived from unit labels,
  including `source_boxes`, `source_dimensions`, `retained_counts`, and
  `ranges`;
- runtime-keyed `pair_family_counts = NamedTuple{families}(...)`;
- same-file internal consumers that currently expect those runtime-keyed
  `NamedTuple` shapes.

Approved replacements:

- vector-backed records or tables with stable field names;
- stable dictionaries keyed by unit or pair-family labels only where lookup by
  label is genuinely needed;
- helper accessors that hide the storage shape from same-file callers;
- compact summaries that expose counts/order without encoding route size in the
  concrete type.

The retained-unit vector remains the ordered inventory authority. Unit labels
and pair-family labels may remain data values, but they must not become type
parameters.

This ID does not approve edits to `RawProductBoxPlan.source_mode_indices`,
`source_mode_column_indices`, `source_mode_indices(...)`,
`TerminalLoweringPlan.available_contracts`, `TerminalLoweringPlan.contracts`,
or `RetainedUnitTransformContractPlan.contracts`. It also does not approve
public input `NamedTuple` changes, fixed `NTuple{3,Int}` coordinate/dimension
changes, artifact sidecar table changes, numerical kernels, route recipe
behavior changes, shellification, terminal lowering, terminal basis, Residual
Gaussian, raw-block changes, canonical driver changes, Hamiltonian object
changes, matrix-key changes, reader changes, public API/export changes,
report/status/payload expansion, compatibility adapters, new committed tests,
Cr2 runs, or Cr2-specific workflow.

Line budget: at most `120` added `src` lines, with net simplification expected.
Failure rule: if the cleanup requires source files outside the approved file,
broader route/stage rewiring, public API changes, artifact changes, or an
adapter that preserves the old runtime-keyed type shape, make no source commit
and report the blocker.

### HP-ROUTE-INV-TEST-01 — retained-unit route inventory cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader or
  canonical driver path;
- focused search confirming no `NamedTuple{unit_keys}` or
  `NamedTuple{families}` route inventory remains in
  `src/pqs_source_box_route_driver_helpers.jl`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
runtime-keyed inventory shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, or route-diagnostic test is approved by this
ID.

### HP-RAW-SRCMODE-FN-01 — raw product source-mode inventory cleanup

Approved source files:

```text
src/cartesian_raw_product_sources/records.jl
src/cartesian_raw_product_sources/source_mode_indices.jl
src/cartesian_raw_product_sources/summaries.jl
```

Approved narrow consumer files, only as required by the storage change:

```text
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_base_hamiltonian.jl
```

Approved cleanup targets:

- `RawProductBoxPlan.source_mode_indices::Tuple{Vararg{NTuple{3,Int}}}`;
- `RawProductBoxPlan.source_mode_column_indices::Tuple{Vararg{Int}}`;
- `source_mode_indices(...)` / source-mode summary accessors only to the extent
  required to hide vector-backed storage from approved callers;
- same-file and listed narrow consumers that currently depend on the
  tuple-backed storage shape.

Approved replacement:

- vector-backed source-mode coordinate storage;
- vector-backed source-mode column storage, or no stored column vector when the
  column sequence is exactly `1:count` and accessors provide the same ordered
  column numbers;
- stable accessors preserving deterministic source-mode order, mode values,
  length, indexing/iteration where currently used, and retained-rule parity.

The fixed `NTuple{3,Int}` coordinate and dimension values remain valid. The
variable-length source-mode inventory must not be encoded in `RawProductBoxPlan`
field types or accessor return types. Accessor compatibility means same facts
and order, not the old concrete `Tuple{Vararg{...}}` return shape.

This ID does not approve terminal-lowering `contracts` /
`available_contracts` tuple cleanup, retained-unit transform-contract tuple
cleanup outside the listed narrow consumer wiring, broad pair-block/source-box
rewrites, numerical kernel changes, route semantic changes, public API/export
changes, canonical driver changes, Hamiltonian object changes, matrix-key
changes, reader changes, artifact schema changes, route-stage objects,
report/status/payload expansion, persistent caches, compatibility layers that
preserve the old tuple-backed shape, new committed tests, Cr2 runs, or
Cr2-specific workflow.

Line budget: at most `150` added `src` lines, with net simplification expected.
Failure rule: if vectorizing the raw product plan requires source files outside
the approved surfaces, broad pair-block/source-box rewrites, public API or
artifact changes, numerical changes, or compatibility layers preserving the old
tuple-backed shape, make no source commit and report the exact caller/blocker.

### HP-RAW-SRCMODE-TEST-01 — raw product source-mode inventory validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint;
- focused raw-product source order and retained-rule parity;
- manifest source-mode and final-basis source-relation inspection;
- focused search confirming `RawProductBoxPlan` no longer stores source-mode
  inventories as `Tuple{Vararg{...}}`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
tuple-backed source-mode inventory shape. No new committed test file, Cr2
fixture, driver-input fixture, benchmark, or route-diagnostic test is approved
by this ID.

### HP-CONTRACT-VEC-FN-01 — contract-plan vector cleanup

Approved source files:

```text
src/cartesian_terminal_lowering/contracts.jl
src/cartesian_terminal_lowering/selection.jl
src/cartesian_terminal_lowering/summaries.jl
src/cartesian_retained_unit_transform_contracts/records.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_retained_unit_transform_contracts/summaries.jl
```

Approved narrow consumer files, only as required by the storage change:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Approved cleanup targets:

- `TerminalLoweringPlan.available_contracts::Tuple{Vararg{TerminalLoweringContract}}`;
- `TerminalLoweringPlan.contracts::Tuple{Vararg{TerminalLoweringContract}}`;
- `RetainedUnitTransformContractPlan.contracts::Tuple{Vararg{RetainedUnitTransformContract}}`;
- same-file and listed narrow consumers that currently depend on those
  tuple-backed plan field shapes.

Approved replacement:

- vector-backed terminal-lowering available-contract storage;
- vector-backed terminal-lowering selected-contract storage;
- vector-backed retained-unit transform-contract storage;
- stable accessors preserving ordered contract facts and current behavior:
  `available_contracts(plan)`, `selected_contracts(plan)`, `contracts(plan)`,
  and `transform_contracts(plan)`.

Accessor compatibility means same ordered facts, iteration order, selected
contract behavior, transform-contract behavior, and summaries. It does not
mean preserving variable-length `Tuple` concrete field types or accessor return
types.

This ID does not approve changing
`source_cpbs::Tuple{Vararg{CoordinateProductBox}}`, raw product source-mode
storage, retained-unit route inventories, public input `NamedTuple`s, fixed
coordinate/product-box value objects, numerical kernels, route semantic
changes, shellification behavior changes, public API/export changes, canonical
driver changes, Hamiltonian object changes, matrix-key changes, reader changes,
artifact/manifest schema changes, route-stage objects,
report/status/payload expansion, persistent caches, compatibility layers that
preserve the old tuple-backed plan field types, new committed tests, Cr2 runs,
or Cr2-specific workflow.

Line budget: at most `150` added `src` lines, with net simplification expected.
Failure rule: if vectorizing the plan inventories requires source files outside
the approved surfaces, broad route/stage rewrites, public API or artifact
changes, numerical changes, or compatibility layers preserving the old
tuple-backed plan field types, make no source commit and report the exact
caller/blocker.

### HP-CONTRACT-VEC-TEST-01 — contract-plan vector cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint;
- focused terminal-lowering contract order parity;
- focused retained-unit transform-contract order parity;
- focused search confirming targeted plan inventories no longer store
  contracts as `Tuple{Vararg{...}}`;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
tuple-backed contract-plan field shape. No new committed test file, Cr2
fixture, driver-input fixture, benchmark, or route-diagnostic test is approved
by this ID.

### HP-ROUTE-STAGE-TYPE-FN-01 — route/stage type-surface cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_terminal_shellification_geometry.jl
```

Approved cleanup targets:

- `_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`;
- `cartesian_units` route/stage return surfaces that carry oversized
  compatibility inventories;
- `_pqs_source_box_route_driver_transform_stage_low_order_summary`;
- `cartesian_transforms` route/stage return surfaces that carry oversized
  compatibility inventories;
- `_cartesian_terminal_shellification_region_unit_inventory`;
- related terminal-region lowering inventory summary surfaces in
  `src/cartesian_terminal_shellification_geometry.jl` only where the same
  runtime-sized type-surface pattern appears.

Approved replacement/deletion shapes:

- delete stale route/stage compatibility inventories with no active approved
  caller;
- replace remaining runtime-sized `NamedTuple` / `Tuple` carriers with
  vector-backed compact internal objects, stable dictionaries, accessors, or
  smaller summaries;
- shrink wide internal stage return signatures only where all live approved
  callers can be updated within the approved source files;
- preserve deterministic terminal shellification/lowering order and existing
  behavior.

Required preservation:

- H2 base artifact/readback behavior;
- H2 supplemented artifact/readback behavior;
- deterministic terminal shellification/lowering order;
- existing public driver contract;
- existing artifact schema and manifest behavior;
- existing numerical matrices.

This ID does not approve source files outside the approved boundary, driver
changes, artifact schema or manifest changes, public API/export changes,
numerical kernel changes, matrix value changes, raw-block changes, Residual
Gaussian/MWG/IDA semantic changes, route semantic changes, shellification
behavior changes, route diagnostic/status/report expansion, broad route-stage
redesign, new public contracts, PackageCompiler/PrecompileTools/sysimage or
precompile workload work, new committed tests, Cr2 runs, or Cr2-specific
workflow. No compatibility adapter may preserve the old runtime-sized type
surface merely under a new name.

Line budget: at most `200` added `src` lines, with net simplification expected.
Failure rule: if cleanup requires source files outside the approved boundary,
broad route-stage redesign, new public contracts, artifact changes, numerical
changes, or a precompile/sysimage mechanism, make no source commit and report
the exact blocker.

### HP-ROUTE-STAGE-TYPE-TEST-01 — route/stage type-surface cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint if the pass touches terminal realization behavior;
- focused terminal shellification/lowering order parity;
- focused scan for newly introduced `NamedTuple{...}`, variable-size
  `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed inventories in the
  approved files;
- optional Be2 q5 compile/timing comparison after correctness passes;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
stale compatibility inventory shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, route-diagnostic test, or precompile workload
is approved by this ID.

### HP-ROUTE-STAGE-CARRIER-FN-01 — route/stage carrier cleanup

Approved source files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Optional source file, only if directly required to slim the approved path:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

Approved cleanup targets:

- `cartesian_shells` stage carrier and return signature;
- `cartesian_units` stage carrier and return signature;
- `cartesian_transforms` stage carrier and return signature;
- terminal topology support-region planning in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- terminal retained-rule planning in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`;
- terminal realization plan carriers in
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl` only
  where directly required to avoid reintroducing a large stage-carried plan
  shape through the approved route/stage path.

Approved replacement/deletion shapes:

- stop carrying giant shellification, route-skeleton, support-plan,
  retained-rule-plan, and terminal-plan `NamedTuple` / tuple shapes across the
  approved stage function signatures;
- replace necessary carriers with compact typed/vector-backed records, stable
  dictionaries, accessors, or smaller summaries;
- recompute small derived summaries from canonical objects inside the approved
  path where that is simpler than carrying wide stage payloads;
- delete stale compatibility carriers with no active approved caller;
- preserve deterministic terminal support, shellification, and lowering order.

Route skeleton construction semantics are not changed by this ID. Edits to
`src/pqs_source_box_route_driver_skeletons.jl` are not approved.

Required preservation:

- H2 base artifact/readback behavior;
- H2 supplemented artifact/readback behavior;
- H2 R3 endpoint if terminal realization is touched;
- deterministic terminal support/shellification/lowering order;
- existing public driver contract;
- existing artifact schema and manifest behavior;
- existing route semantics and numerical matrices.

This ID does not approve source files outside the approved boundary, edits to
`src/pqs_source_box_route_driver_skeletons.jl`, driver changes, artifact schema
or manifest changes, public API/export changes, numerical kernel changes,
matrix value changes, raw-block changes, Residual Gaussian/MWG/IDA semantic
changes, route semantic changes, shellification behavior changes, route
diagnostic/status/report expansion, broad route-stage redesign, new public
contracts, PackageCompiler/PrecompileTools/sysimage or precompile workload
work, new committed tests, Cr2 runs, or Cr2-specific workflow. No compatibility
adapter may preserve the old runtime-sized carrier merely under a new name.

Line budget: at most `250` added `src` lines, with net simplification expected.
Failure rule: if cleanup requires source files outside the approved boundary,
broad route-stage redesign, public API changes, artifact changes, numerical
changes, or precompile/sysimage machinery, make no source commit and report the
exact blocker.

### HP-ROUTE-STAGE-CARRIER-TEST-01 — route/stage carrier cleanup validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- H2 R3 endpoint if terminal realization is touched;
- focused terminal support/shellification/lowering order parity;
- focused scan for newly introduced runtime-sized `NamedTuple{...}`,
  `Tuple(...)`, `Tuple{Vararg{...}}`, and runtime-keyed inventories in the
  approved files;
- optional Be2 q5 post-cleanup compile/timing comparison after correctness
  passes;
- no Cr2 run.

Existing committed tests may be adjusted only if they directly assert the old
stale carrier shape. No new committed test file, Cr2 fixture,
driver-input fixture, benchmark, route-diagnostic test, or precompile workload
is approved by this ID.

### HP-R1-ART-01 — public base producer artifact provenance

Approved artifact extension for R1 public facade writes only. When
`hamfile !== nothing`, the final Hamiltonian file may add the following JLD2
keys under `producer_provenance/`:

```text
provenance_version
producer
nesting
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
`core_spacing = 0.3`, `reference_spacing = 1.0`,
`mapping_kind = :white_lindsey_atomic_mapping`, and
`mapping_d = 0.3` for the reviewed endpoint. The `mapping_d` provenance value
is the resolved internal White-Lindsey mapping parameter and equals
`core_spacing` for one-center atoms; it is not a separate public input.
`nesting` must record the public construction-family input, and `route` must be
the truthful base route label derived from `(input.kind, input.nesting)`, not a
PQS-oriented default string.
Existing
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
reviewed H baseline using explicit public `core_spacing = 0.3` with
`reference_spacing = 1.0` and H2 endpoint facts, validate unknown-key and
malformed input errors including mismatched temporary H `d` and rejected H2
`d`, validate
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
- keep required one-center basis fields `q`, `core_spacing`, and `radius`;
- treat public `d`, if temporarily accepted, as a deprecated compatibility
  alias that must equal resolved `core_spacing`.

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
- map the resolved public `core_spacing` to the private White-Lindsey
  `parent_mapping_d`;
- keep `reference_spacing`, `tail_spacing`, and box/domain controls separate
  from `core_spacing`;
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
  now expressed as `core_spacing = 0.3`, `reference_spacing = 1.0`, and
  internal `parent_mapping_d = core_spacing`;
- optional ignored/user-run Be or Cr one-center base atom artifact
  write/readback using explicit charge, spin sectors, origin geometry, and
  basis controls;
- finite/symmetric `K`, unit `U_A`, and IDA `V` for ignored/user-run non-H
  atom checks;
- clear `ArgumentError` for translated atom input, mismatched temporary `d`,
  noninteger or nonpositive charge, nonneutral electron count, or
  element-table/default requests where practical.

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

## Approved For Compact Hamiltonian Artifact Manifest

This section approves only artifact sidecar groups for existing
`CartesianIDAHamiltonian{Float64}` JLD2 files. It does not approve a new
Hamiltonian object, new matrix keys, public reader API, driver public input
change, solver workflow, CR2-consumer-specific field, Cr2-specific field,
report/status payload, or artifact schema dump in the driver.

### HP-HAM-MANIFEST-FN-01 — compact Hamiltonian artifact manifest

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_ida_hamiltonian.jl
```

`src/cartesian_ida_hamiltonian.jl` is approved only for a small unexported
sidecar writer/helper if needed; it must not change
`write_cartesian_ida_hamiltonian` matrix keys or
`read_cartesian_ida_hamiltonian` behavior.

Approved sidecar groups:

```text
hamiltonian_manifest/
hamiltonian_manifest/final_basis_labels/
hamiltonian_manifest/final_basis_source_relations/
hamiltonian_manifest/source_shells/
hamiltonian_manifest/source_modes/
recipe_provenance/
```

`hamiltonian_manifest/` reuses the source/fixed-column provenance model from
`docs/src/developer/projected_q_shell_policy.md`. Basis identity is a
status-bearing construction label, not `center_xyz`. Centers are
representative metadata only.

Approved `hamiltonian_manifest/` root key:

- `manifest_version = 1`;

Approved required `hamiltonian_manifest/final_basis_labels/` fields:

- `final_basis_col`;
- `sector`;
- `unit_label`;
- `unit_kind`;
- `source_region_label`;
- `source_region_label_status`;
- `source_box_label`;
- `source_box_label_status`;
- `owner_nucleus_index`;
- `owner_label_status`;
- `shell_label_status`;
- `shell_index`;
- `ray_label_status`;
- `ray_id`;
- `ray_family_label`;
- `radial_order_status`;
- `radial_order`;
- `center_x`;
- `center_y`;
- `center_z`;
- `center_definition`;
- `center_status`;
- `lowdin_correction_applied`;
- `supplement_label`;
- `angular_power_x`;
- `angular_power_y`;
- `angular_power_z`;
- `inferred_from_centers`;
- `inferred_from_nearest_grid`;
- `inferred_from_support_order`;
- `inferred_from_support_indices`;
- `inferred_from_raw_to_final_support`.

The final-basis label rows must follow the exact matrix row/column order.
`owner_nucleus_index` uses one-based physical nucleus indices and `0` when no
owner is meaningful; unavailable integer shell/ray/radial labels use `0`;
unavailable angular powers use `-1`; and all `inferred_from_*` flags must be
`false` for production manifest rows. Approved sectors are `:base`,
`:residual`, and `:supplement_derived`.

Approved optional `hamiltonian_manifest/final_basis_source_relations/` fields:

- `final_basis_col`;
- `relation_index`;
- `relation_kind`;
- `source_shell_id`;
- `source_mode_label`;
- `local_axis_x`;
- `local_axis_y`;
- `local_axis_z`;
- `relation_status`;
- `shell_label_status`;
- `ray_label_status`;
- `radial_order_status`;
- `coefficient_status`;
- `weight_status`;
- `span_status`;
- `inferred_from_centers`;
- `inferred_from_nearest_grid`;
- `inferred_from_support_order`;
- `inferred_from_support_indices`;
- `inferred_from_raw_to_final_support`.

Relations, `source_shells/`, and `source_modes/` may be populated only when
the construction natively defines those facts. Source-mode identity is
`(source_shell_id, local_axis_x, local_axis_y, local_axis_z)` in shell-local
coordinates. It is a label, not a coefficient map, support row, parent row, or
operator reconstruction claim. Missing shell, ray, radial, source-box, or
relation facts must be status-bearing `:unavailable` or `:mixed`, not inferred
from centers, nearest-grid snapping, support order, support indices, or
raw-to-final support.

Approved `recipe_provenance/` keys:

- `provenance_version = 1`;
- `producer`;
- `nesting`;
- `route`;
- `q`;
- `core_spacing`;
- `padding`;
- `radius`;
- `xmax_parallel`;
- `xmax_transverse`;
- `extent_source`;
- `parent_axis_counts`;
- `atom_symbols`;
- `nuclear_charges`;
- `atom_locations`;
- `nup`;
- `ndn`;
- `basisname`;
- `basisfile`;
- `lmax`;
- `uncontracted`;
- `width_filtering`;
- `base_dimension`;
- `residual_dimension`;
- `augmented_dimension`.

The recipe group may repeat facts from `producer_provenance/` and
`supplement_provenance/` so consumers have one uniform location. Values must
come from the validated public construction contract and produced dimensions,
not route reports, element tables, solver assumptions, or private diagnostics.
`nesting` records the public construction family (`:pqs` or `:wl`), and
`route` records the truthful base route label derived from `(input.kind,
input.nesting)`.

Center conventions and construction labels must be derived from existing
terminal basis blocks, parent axes, residual metadata, and augmented moment/MWG
descriptors. If the implementation cannot derive a required center convention
or source label from those existing objects without adding algorithmic
metadata, it must stop and report the missing seam.

This ID does not approve `T_G`, `T_A`, dense residual transforms, coefficients,
dense moment matrices, raw inventories, allocation probes, report/status
payloads, public reader APIs, public exports, driver public input changes,
artifact schema dumps in the driver, solver-specific fields,
CR2-consumer-specific fields, Cr2-specific fields, committed Cr2 fixtures, or
Cr2-specific branches.

One-center atom padding is provenance-only for this lane. Source work under
this ID must not change the current one-center atom size policy or parent-axis
counts; diatomic padding-derived extents continue to use the existing facade
contract.

Line budget: at most `150` added `src` lines. Stop for a separate amendment if
the implementation needs source files outside the approved surfaces, new
algorithmic metadata, a public reader, driver contract changes, or changes to
the Hamiltonian matrix writer/reader contract.

### HP-HAM-MANIFEST-TEST-01 — artifact manifest validation

Approved validation:

- `git diff --check`;
- package load;
- H atom or H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- direct JLD2 checks that `hamiltonian_manifest/final_basis_labels/` rows
  match matrix dimension and approved status-bearing fields;
- direct JLD2 checks that unavailable or mixed shell/ray/radial/source labels
  are explicit and that no production row is marked inferred from centers,
  nearest grid, support order, support indices, or raw-to-final support;
- direct JLD2 checks that `recipe_provenance/` records validated public
  system, basis, supplement, route, parent-axis counts, and dimensions;
- optional practical Be2 supplemented artifact manifest inspection;
- no Cr2 run.

No new committed test file, public reader API, artifact schema dump, driver
public input change, Cr2 fixture, or solver/CR2 workflow validation is approved
by this ID.

### HP-HAM-MANIFEST-SRC-FN-01 — source-mode provenance seam

Approved purpose: carry compact construction-native source-mode provenance
from terminal lowering / retained-unit / raw-product source planning to the
base working basis manifest context so artifact writing can populate optional
source provenance groups without route reports or center inference.

Approved source files:

```text
src/cartesian_terminal_lowering/contracts.jl
src/cartesian_terminal_lowering/region_contracts.jl
src/cartesian_raw_product_sources/CartesianRawProductSources.jl
src/cartesian_raw_product_sources/records.jl
src/cartesian_raw_product_sources/source_mode_indices.jl
src/cartesian_retained_units/CartesianRetainedUnits.jl
src/cartesian_retained_units/records.jl
src/cartesian_retained_units/lower_contract_units.jl
src/cartesian_retained_unit_transform_contracts/CartesianRetainedUnitTransformContracts.jl
src/cartesian_retained_unit_transform_contracts/records.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_base_hamiltonian.jl
```

Module wrapper files are approved only for internal exports/includes required
by the compact provenance seam. They are not public API authority.

Approved carrier:

- preferred: one internal `source_mode_provenance` field on the
  `cartesian_base_working_basis(...)` return value;
- allowed if cleaner: one optional compact source-mode provenance field on
  `CartesianTerminalBasisRealization`;
- no other stage object, report, status payload, route summary, artifact
  wrapper, or persistent cache is approved.

Approved source facts:

- source shell IDs, unit links, source-box/source-region labels, construction
  kind, source intervals, source-mode dimensions, source-mode ordering, center
  definition/status, Lowdin-correction status, and shell/ray/radial label
  statuses;
- source mode identities
  `(source_shell_id, local_axis_x, local_axis_y, local_axis_z)`, local
  source-mode ordering, native parent-lattice coordinates only when already
  available, representative center metadata, and status flags;
- final-basis source-relation rows only where the relation is construction
  native: direct identity, boundary source-mode selection, product-axis tuple,
  or explicit `:mixed` / `:unavailable`;
- final-basis label improvements only where the final basis column is directly
  and natively tied to the row's unit/source mode.

Ray, cone, shell, and radial labels may be written only when already natively
defined by the construction. This ID does not approve a repo-chosen ray/cone
grouping policy or inferred labels from centers, nearest-grid snapping, support
order, support indices, or raw-to-final support.

This ID does not approve coefficients, dense transforms, `T_G`, `T_A`, raw
candidate inventories, raw pair inventories, allocation probes, route reports,
diagnostic payloads, metadata/status field clouds beyond the compact approved
provenance object, Hamiltonian object changes, matrix-key changes,
`read_cartesian_ida_hamiltonian` changes, public API/export changes, driver
changes, artifact schema dumps, solver fields, CR2-consumer-specific fields,
Cr2-specific fields, committed Cr2 fixtures, or Cr2-specific branches.

Line budget: at most `180` added `src` lines. Stop for a separate amendment if
the implementation needs new source files, public exports, dense payloads,
driver changes, reader changes, route report plumbing, or a source-mode/ray
policy not already present in construction metadata.

### HP-HAM-MANIFEST-SRC-TEST-01 — source-mode provenance seam validation

Approved validation:

- `git diff --check`;
- package load;
- H2 base artifact write/readback through the existing reader;
- H2 supplemented artifact write/readback through the existing reader;
- direct JLD2 checks that optional `source_shells/`, `source_modes/`, and
  `final_basis_source_relations/` are present only for construction-native
  provenance rows;
- checks that absent shell/ray/radial/source facts are explicitly
  `:unavailable` or `:mixed` and not inferred;
- optional practical Be2 supplemented manifest inspection;
- no Cr2 run.

No new committed test file, public reader API, driver public input change,
artifact schema dump, Cr2 fixture, solver/CR2 workflow validation, or broad
route/report validation is approved by this ID.

### HP-NEST-ART-FN-01 — nesting artifact-truth cleanup

Approved source files:

```text
src/cartesian_base_hamiltonian.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

`src/cartesian_base_hamiltonian.jl` is approved only for artifact provenance
truth for the public `nesting` construction-family input. It may:

- record `nesting` in `producer_provenance/` and `recipe_provenance/`;
- choose truthful base route labels from `(input.kind, input.nesting)`, with
  approved labels `:one_center_pqs_base`, `:one_center_wl_base`, and
  `:z_axis_diatomic_pqs_base`;
- reject supplemented `nesting = :wl` before expensive base-stage
  construction with a clear `ArgumentError`;
- leave unsupported WL H2 without a provenance label until that path succeeds
  under separate authority.

`src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl` is
approved only for a docstring correction so the module description no longer
says it is exclusively PQS-specific now that the WL terminal-basis seam uses the
same final-basis boundary.

This ID does not approve driver public input changes, route skeleton changes,
shellification changes, terminal lowering changes, raw-block changes,
Residual Gaussian/MWG/IDA changes, artifact matrix-key changes,
`read_cartesian_ida_hamiltonian` behavior changes, public API/export changes,
Cr2 workflow, committed tests, diagnostic/report changes, or WL H2 support.

Failure rule: if truthful nesting provenance requires changing reader behavior,
artifact matrix keys, or the broader manifest structure, make no source commit
and report the blocker.

### HP-NEST-ART-TEST-01 — nesting artifact-truth validation

Approved validation:

- `git diff --check`;
- package load;
- small `nesting = :pqs` base artifact write/readback plus direct provenance
  inspection;
- small `nesting = :wl` one-center atom artifact write/readback plus direct
  provenance inspection;
- supplemented `nesting = :wl` rejects before base-stage construction;
- no Cr2 run.

No new committed test file, driver-input fixture, public reader API, artifact
schema dump, WL H2 validation, supplemented WL validation, or Cr2 fixture is
approved by this ID.

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

### HP-RG-ORTHO-FN-01 — residual final-orthogonality robustness

Approved source files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is
approved only for narrow existing internal keyword plumbing if an approved
tolerance/check option must be passed through a compatibility entry point.

Approved target: make final residual orthogonalization and final
`R' S R - I` identity validation robust for small floating-point overshoots
with healthy owner-local selection and final merge spectra. The motivating
strict N2 q5 p10 case at `core_spacing = 0.042857` has `max |G' S R| =
1.776e-14`, `max |R' S R - I| = 1.673e-10`, retained counts `9,9`, and final
merge eigenvalues `7.232e-2 .. 1.928`.

Allowed changes:

- stable symmetric final merge normalization/check;
- explicitly symmetrized final residual-overlap validation;
- combined absolute/relative final identity check
  `err_RR <= 1.0e-10 + 1.0e-10 * max(1, scale_RR)`;
- no public API unless the option is already routed through internal keywords.

This ID does not approve blind broad tolerance relaxation, residual-selection
semantic changes, global residual selection, occupation-cutoff changes,
negative-eigenvalue tolerance changes, final merge eigenvalue flooring, width
filtering as a conditioning repair, MWG/IDA/nuclear/raw-block changes, artifact
schema changes, driver changes, status/report fields, public API/export
changes, new committed tests, Cr2 workflow, or source files outside the two
approved files.

Failure rule: if the strict N2 case requires changing residual selection,
supplement construction, or final-basis construction, make no source commit and
report the blocker.

### HP-RG-ORTHO-TEST-01 — residual final-orthogonality validation

Approved validation:

- existing H2 residual-GTO/MWG endpoint;
- H2 base/supplemented readback if touched through the facade or compatibility
  file;
- ignored strict N2 q5 p10 residual audit or artifact smoke at
  `core_spacing = 0.042857`;
- one passing N2 comparison at `core_spacing = 0.05` or `0.075`;
- report `max |G' S R|`, `max |R' S R - I|`, retained owner counts, and final
  merge eigenvalue range/condition.

No committed fixture/test, Cr2 full Hamiltonian, Cr2 artifact, Cr2 facade
support, driver workflow, artifact schema change, solver/RHF, ECP, or EGOI
work is approved.

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

## Approved For R3 Same-Construction Base K/U Reuse

This section approves only narrow reuse of already-built same-construction
base final-basis kinetic and unit-nuclear blocks in supplemented residual-GTO
/ MWG exact augmented operators. It is an orchestration reuse lane, not a
terminal product, Gaussian-sum, raw-block, residual-basis, or interaction
algorithm lane.

Evidence recorded before approval: a replay that reused same-construction base
`K_GG` and unit `U_GG[A]` blocks had exact operator delta `0.0` and reduced the
exact augmented-operator replay to `0.8620s / 1237.136 MiB`.

### HP-R3BASE-FN-01 — same-construction base K/U reuse

Approved owner:

```text
Owner module: CartesianFinalBasisRealization plus narrow caller wiring
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- `pqs_terminal_residual_gto_augmented_products(...)`, or its approved caller
  wrapper, may accept a trusted same-construction base kinetic matrix and use
  it as the `G-G` kinetic block for `transform_augmented_operator`;
- `pqs_terminal_residual_gto_augmented_unit_nuclear(...)`, or its approved
  caller wrapper, may accept trusted same-construction unit nuclear
  `U_GG[A]` matrices and use them as the `G-G` unit blocks for
  `transform_augmented_operator`;
- `cartesian_residual_gto_mwg_hamiltonian(...)` and staged helpers in
  `src/cartesian_base_hamiltonian.jl` may pass `base_ham.kinetic` and
  `base_ham.nuclear_attraction_unit_by_center` into the augmented operator
  construction;
- current behavior must be preserved when trusted base blocks are not supplied.

Trust condition:

- the base Hamiltonian, terminal basis realization, parent bundles, residual
  basis, and supplement must come from the same
  `cartesian_base_working_basis(...)` construction path;
- implementation must validate matrix dimensions and center count before
  reuse;
- no provenance payload, metadata proof, report field, status object, or
  persistent cache is required or approved for this trust check.

This ID does not approve public API/export changes, canonical-driver changes,
raw-block changes, residual selection/orientation/transform changes, MWG/IDA
convention changes, terminal product or Gaussian-sum kernel rewrites,
persistent cache/workspace objects, metadata/status/report/artifact schema
fields, route/stage setup cleanup, committed tests, Cr2 workflow, or source
files outside the two approved files.

Line budget: target under `100` added `src` lines. If trusted
same-construction provenance cannot be guaranteed by local call shape plus
dimension/center validation, or if implementation needs public payloads,
metadata, or stage objects, make no source commit and report the blocker.

### HP-R3BASE-TEST-01 — base K/U reuse validation

Approved validation:

- `git diff --check`;
- package load;
- existing H2 R3 endpoint unchanged;
- Be2 supplemented facade/readback unchanged except allowed
  timing/allocation improvement;
- Cr2 exact-operator attribution audit or focused ignored replay showing base
  `K_GG` / unit `U_GG[A]` reuse parity and allocation effect;
- final exact operators finite and symmetric.

No new committed test file, Cr2 artifact, Cr2 workflow, public API/export,
driver workflow change, metadata/status/report field, or artifact schema
change is approved by this ID.

### HP-R3BASE-DRV-WIRE-01 — canonical driver K/U reuse call-site wiring

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- in supplemented mode only, pass `base_ham.kinetic` into
  `cartesian_residual_gto_augmented_products(...)` as `base_kinetic`;
- in supplemented mode only, pass
  `base_ham.nuclear_attraction_unit_by_center` into
  `cartesian_residual_gto_augmented_unit_nuclear(...)` as
  `base_unit_nuclear`;
- keep the current public inputs, hooks, timing labels, visible stage sequence,
  artifact schema, and driver contract unchanged.

This ID is only call-site wiring so the canonical driver uses the already
approved same-construction base K/U reuse path. It does not approve source or
kernel changes, diagnostics, new hooks, new timing labels, public API/export
changes, artifact changes, tests/fixtures, Cr2 workflow, or edits outside
`bin/cartesian_ham_builder.jl`.

Failure rule: if the driver call-site update needs any visible driver contract
change, make no source commit and report the blocker.

### HP-R3BASE-DRV-TEST-01 — driver K/U reuse validation

Approved validation:

- `git diff --check`;
- package load;
- H2 supplemented driver artifact/readback;
- Be2 supplemented driver artifact/readback if practical;
- no Cr2 run.

No new committed test file, fixture, diagnostic, hook, timing label, public
input, artifact schema, or Cr2 workflow is approved by this ID.

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
optional `supplement`, `nesting`, `hamfile`, `padding`, `check_file`,
`print_contract`, `print_timing`, and `expected_dimension`.

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

### HP-DRV-NEST-FN-01 — construction-family driver input

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- add a visible public driver input `nesting`;
- accepted values are `:pqs` and `:wl`;
- default is `nesting = :pqs`;
- `nesting = :pqs` means the PQS source-box construction family;
- `nesting = :wl` means the White-Lindsey low-order construction family;
- include `nesting` in public contract construction, optional
  `print_contract`, and optional `check_file` output as a public contract fact.

This is a first-class construction-family choice, not a diagnostic route
switch. The driver must not expose internal route-family names, route
skeletons, retained-rule plans, raw-block switches, stop-after controls,
diagnostic knobs, old route-stage labels, private helper calls, allocation
probes, or route reports.

This ID does not approve public API/export changes, artifact schema changes,
stage-label changes, solver/ECP workflow, Cr2-specific behavior, broad driver
diagnostics, committed fixtures/tests, or source files outside the canonical
driver.

### HP-DRV-NEST-WIRE-01 — construction-family route mapping

Approved source files:

```text
bin/cartesian_ham_builder.jl
src/cartesian_base_hamiltonian.jl
```

Approved behavior:

- map public `nesting = :pqs` to the existing internal `:pqs_source_box`
  route family;
- map public `nesting = :wl` to the existing internal
  `:white_lindsey_low_order` route family;
- keep route skeletons, retained rules, raw-block switches, stop-after
  controls, diagnostics, and internal route-stage vocabulary hidden;
- preserve the existing public stage labels, Hamiltonian object, matrix keys,
  artifact schema, driver hooks, and solver-free workflow;
- reject unsupported combinations with clear `ArgumentError`s.

Supplemented `nesting = :wl` is approved only if it is already valid through
the existing supported supplemented facade/staged path. If it is not already
valid, implementation must reject the combination clearly and report it as a
separate design decision instead of adding broad supplemented White-Lindsey
route behavior.

This ID does not approve new route algorithms, route-skeleton construction
changes, White-Lindsey materialization deletion, terminal lowering policy
changes, shellification behavior changes, numerical kernel changes, raw-block
changes, Residual Gaussian/MWG/IDA changes, artifact/provenance schema changes,
public API/export changes, committed tests, Cr2-specific workflow, or source
files outside the two approved files.

Line budget: at most `80` added source/bin lines, with net simplification
preferred where old hidden assumptions can be removed.

Failure rule: if `nesting = :wl` cannot produce a small base artifact/readback
through the existing White-Lindsey low-order route without broader route or
materialization work, make no source commit and report the exact blocker.

### HP-DRV-NEST-TEST-01 — construction-family validation

Approved validation:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` base driver or facade artifact/readback;
- current default `nesting = :pqs` supplemented H2 driver/facade path if
  supplemented-mode input plumbing is touched;
- one small `nesting = :wl` base artifact/readback using a currently supported
  base geometry;
- explicit negative check or ignored smoke showing unsupported supplemented
  `nesting = :wl` fails clearly if that combination is not already valid;
- no Cr2 run.

No new committed test file, committed input fixture, artifact schema
validation, solver run, Cr2-specific driver run, or broad White-Lindsey
workflow validation is approved.

### HP-DRV-STAGE-FN-01 — visible physics-stage producer surface

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
src/pqs_source_box_low_order_materialization.jl
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

Approved purpose: expose a small non-exported, non-underscored,
driver-facing staged producer surface so the canonical driver can execute and
time visible physics-level construction stages without calling private
underscored helpers.

`src/cartesian_base_hamiltonian.jl` remains the primary driver-facing owner.
The lower-level files listed above are approved only for behavior-preserving
operator-class stage factoring in their existing domains; they are not
approved for new algorithms, raw-block changes, or diagnostics.

Approved visible stages are:

- construct public `system`, `basis`, and optional `supplement`;
- build the base working basis / terminal realization;
- build base product/moment operators;
- build base unit-nuclear attraction operators;
- build base electron-electron / IDA interaction;
- assemble the base Hamiltonian;
- load or build the Gaussian supplement basis when `basisname !== nothing`;
- build residual Gaussian augmentation;
- build augmented product/moment operators;
- build augmented unit-nuclear attraction operators;
- build augmented electron-electron / residual-MWG interaction;
- assemble the supplemented Hamiltonian;
- write and check the artifact.

Approved physical operator classes are:

- product/moment operators: kinetic `K`, Cartesian coordinate moments
  `x`/`y`/`z`, and second moments `x^2`/`y^2`/`z^2` where present;
- unit-nuclear attraction: uncharged by-center `U_A` / `Vnuc` matrices before
  applying physical nuclear charges;
- electron-electron interaction: base localized IDA `Vee` and supplemented
  residual-MWG/IDA `Vee`.

The first and last stages remain driver/writer responsibilities. This ID
approves source factoring needed for the base working-basis/terminal
realization, base product/moment, base unit-nuclear, base IDA, Gaussian
supplement, residual augmentation, augmented product/moment, augmented
unit-nuclear, residual-MWG/IDA, and Hamiltonian assembly stages.

The staged surface may factor the existing `cartesian_base_hamiltonian(...)`
and `cartesian_residual_gto_mwg_hamiltonian(...)` bodies so that those facades
can remain wrappers over the same implementation. It may return existing
domain objects and small fixed-key ephemeral stage products required by the
next approved stage.

The staged surface must be a set of separate named construction-stage
functions. It must not be a single opaque replacement wrapper that hides the
same construction sequence under a new name. The canonical driver must be able
to bind visible local variables for the base realization, base products,
base unit-nuclear operators, base `Vee`, base Hamiltonian, supplement basis,
residual augmentation, augmented products, augmented unit-nuclear operators,
augmented `Vee`, and final Hamiltonian assembly.

This ID does not approve public exports, public API redesign, route-stage
objects, reports, status/result payloads, metadata field clouds, runtime-keyed
field groups, persistent caches, raw-block switches, allocation probes,
per-kernel timing frameworks, solver/ECP workflow, artifact schema changes, or
source files outside the four approved paths listed above.

Line budget: at most `200` added `src` lines across the approved source files.
If the staged surface requires a new module, source files outside the four
approved paths, a broad payload object, committed tests, new artifact keys,
raw-block changes, or kernel rewrites, stop and request a new docs-only
amendment.

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

Driver timing should expose the three physical operator classes separately:
product/moment, unit-nuclear, and electron-electron. These timings are
user-facing stage timings only. They must not become allocation probes,
raw-block timing controls, per-kernel instrumentation, or diagnostic stop
points.

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
- pass explicit one-center base `basis` fields, including `core_spacing`, to
  the existing base facade;
- allow visible, easily edited driver/project defaults such as
  `core_spacing = 0.3` and template `padding`, while treating them as explicit
  resolved input values that may be overridden for quick tests;
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
  electron count, mismatched temporary `d` if accepted, or unsupported atom
  input.

No supplemented atom endpoint, translated-atom gate, committed atom fixture,
new committed test file, solver run, artifact-schema validation, or broader
base-atom validation is approved by this ID.

### HP-DRV-ATOM-CLEAN-01 — remove hidden atom `d` driver residue

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

Approved behavior:

- remove the hidden one-center atom basis field `d = vars[:core_spacing]`;
- keep the visible driver atom basis in terms of `q`, `core_spacing`,
  `radius`, and existing optional public fields only;
- keep public inputs, defaults, overrides, hooks, timing labels, visible stage
  sequence, artifact schema, and driver contract unchanged.

This ID exists because the producer no longer requires public `d` for
one-center atoms. It does not approve source/kernel changes, diagnostics, new
hooks, new timing labels, public input changes, committed tests/fixtures, Cr2
workflow, old `:white_lindsey_low_order` retirement, test/tool route-input
cleanup, or edits outside `bin/cartesian_ham_builder.jl`.

Failure rule: if removing the hidden `d` field requires any visible driver
contract change or producer/source change, make no source commit and report the
blocker.

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
