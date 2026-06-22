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

## Approved For R3-A Implementation

These entries authorize only the R3-A residual-GTO basis plus exact one-body
and moment scope recorded in `r3_residual_gto_mwg_augmentation.md`. R3-A may
implement deterministic residual-basis construction plus exact augmented
kinetic, uncharged by-center nuclear attraction, and moment matrices
`x`/`y`/`z`/`x^2`/`y^2`/`z^2`.

R3-A does not approve MWG/IDA `V`, supplemented
`CartesianIDAHamiltonian` construction, artifact provenance, public API
expansion, driver/bin/tool workflow, broad provider payloads, status/result
objects, report fields, pair/assembly public workflow, Be2 first-gate
validation, Cr2 validation, ECP, or EGOI.

Approved first fixture and spike evidence from manager-log Pass 048:

- public/base z-axis H2;
- contracted H/cc-pVTZ on both physical H centers;
- `lmax = 1`;
- `uncontracted = false`;
- no width filtering;
- 18 supplement candidates total, 9 per center;
- full residual rank `18`;
- residual metric eigenvalue range approximately `3.05e-4` to `1.35e-2`;
- symmetric Lowdin residualization in candidate order is not marginal.

Approved residual thresholds:

```text
tau_abs = 1.0e-10
tau_rel = 1.0e-10
tau_neg_abs = 1.0e-12
tau_neg_rel = 1.0e-12
```

Approved R3-A source owner/path:

```text
Owner module: CartesianFinalBasisRealization
Source file: src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`HP-R3-OBJ-01`, `HP-R3-FN-01`, and `HP-R3-FN-02` may be implemented only in
that file. Existing `CartesianFinalBasisRealization` module include plumbing may
include the file only to expose the approved internal R3-A surfaces. This does
not approve a public API or export, a second source owner, a driver/tool
workflow, or any R3-B/R3-C surface.

### HP-R3-OBJ-01 — residual-GTO augmentation object

Approved R3-A scope. The residual object is a numerical object, not a
status/result payload, and must expose matrices as typed fields rather than
metadata. Fields and invariants are defined in
`r3_residual_gto_mwg_augmentation.md` and include base dimension, candidate
count, residual dimension, deterministic candidate labels/order, derived
candidate owner indices, retained candidate indices, `T_G::Matrix{Float64}`,
`T_A::Matrix{Float64}`, residual metric eigenvalue diagnostics, frozen keep
thresholds, rank rule, orientation rule, and sign rule.

No hidden metadata matrices, route-global field clouds, readiness/status
graphs, or Hamiltonian wrappers are approved.

### HP-R3-FN-01 — deterministic residual-basis construction

Approved R3-A scope. Construction:

```text
X = G' S A
S_AA = A' S A
S_R = S_AA - X'X
A_perp = A - G X
R = (A - G X) T
T_A = T
T_G = -X T
```

`G` is the fixed orthonormal base terminal final basis and `A` is the
deterministically ordered Gaussian supplement sector. `X` must be assembled
from exact mixed overlaps one terminal block at a time; no global parent-space
coefficient matrix is approved.

Residual orientation is selected-candidate-order symmetric Lowdin, not raw
eigenvectors of `S_R`. Rank selection uses
`tau_keep = max(tau_abs, tau_rel * lambda_max)` with negative-eigenvalue error
handling as defined in the R3 note and the frozen thresholds above.

### HP-R3-FN-02 — exact augmented one-body and moment assembly

Approved R3-A scope. For each exact one-body or moment operator `O`, assemble
raw `[G, A]` blocks and transform:

```text
O_GR = O_GG T_G + O_GA T_A
O_RR =
    T_G' O_GG T_G
  + T_G' O_GA T_A
  + T_A' O_AG T_G
  + T_A' O_AA T_A
O_aug = [O_GG  O_GR
         O_GR' O_RR]
```

This applies to kinetic, every uncharged by-center nuclear attraction, `x`,
`y`, `z`, `x^2`, `y^2`, and `z^2`. Moment matrices are not stored in
`CartesianIDAHamiltonian`, so R3-A must produce/consume them in the same
construction boundary and must not recover them from an arbitrary base
Hamiltonian.

`HP-R3-FN-02` does not require CPB providers as the only exact-block
implementation. Inside the approved owner file, R3-A may use the QW analytic
1D-table donor organization to build full parent-by-supplement `G-A` matrices
once for overlap, kinetic, position, second moments, and by-center nuclear
attraction, then project parent rows through terminal blocks. It may reuse the
once-built overlap `G-A` block for `X = G' S A` and use the once-built
`G-A`/`A-A` blocks in `pqs_terminal_residual_gto_augmented_operators`. This is
allowed rectangular cross data, not an approved parent-by-parent global
operator, dense global pair matrix, new shared QW API, persistent provider
bundle, payload/status/report surface, artifact/provenance, public API,
driver/tool workflow, Be2 validation gate, or Cr2 work.

### HP-R3-TEST-01 — first augmented one-body endpoint validation

Approved standalone R3-A/R3-B endpoint gate:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

Invocation:

```text
julia --project=. test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

The gate covers only the H2 residual-GTO endpoint family. It should
validate the frozen H2/H/cc-pVTZ fixture, candidate ownership, `G' S R`,
`R' S R`, base G-G block equality, finite/symmetric augmented `K`, uncharged
`U_A`, and moment matrices, and `E1_aug <= E1_base + epsilon`. It is a
standalone endpoint/integration gate and is not approved for inclusion in
`test/runtests.jl`. Under approved R3-B, the same file may be extended only
with the first in-memory supplemented Hamiltonian checks: finite/symmetric
`V_aug`, unchanged base `V_GG`, returned `CartesianIDAHamiltonian{Float64}`,
augmented dimension `489`, an independent weight-aware `V_GM` comparison, and
lowest-orbital IDA self-Coulomb `0.4574256036192161` within `1.0e-10`. The
independent `V_GM` check must recompute the final-basis density-normalized
mixed block from support weights, final weights, and support-to-M donor values;
it must not only compare the final self-Coulomb scalar against the same
implementation path. It must not assert private pair/assembly/report/status
behavior and must not run Be2 or Cr2.

## Approved For R3-B Implementation

### HP-R3-FN-03 — residual MWG/IDA and in-memory Hamiltonian

Approved source owner/path/function:

```text
Owner module: CartesianFinalBasisRealization
Source file: src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
Function: pqs_terminal_residual_gto_augmented_hamiltonian
```

This ID approves only the narrow R3-B in-memory interaction continuation for
the accepted H2 R3-A path. It may compute residual MWG descriptors from exact
R3-A moments of the actual final residual functions:

```text
c_ralpha = <r | alpha | r>
v_ralpha = <r | alpha^2 | r> - c_ralpha^2
sigma_ralpha = sqrt(2 * v_ralpha)
```

The MWG approximation is separable, uses the repo Gaussian-width convention in
the factor `sqrt(2 * v_ralpha)`, and omits off-diagonal covariance. Nonfinite
moments, nonpositive variances, or nonpositive widths are construction errors.

R3-B uses density-normalized `M-M` pair factors directly. For `G-M`, donor
support-to-M factors are parent-density normalized and must be transformed to
final-basis density normalization for each terminal block before insertion:

```text
support_weights = wx .* wy .* wz
final_weights = C' * support_weights
C_density = C .* support_weights ./ final_weights'
V_GM_block = C_density' * V_support_M
```

Direct blocks use identity/final weights consistently and therefore agree with
the direct insertion formula. PQS shell blocks must use the weight-aware
contraction above. The resulting interaction is inserted into

```text
V_aug = [V_GG_base  V_GM
         V_GM'      V_MM]
```

with no additional downstream division by final weights. `M` denotes
moment-matched effective
Gaussians, not raw supplement candidates and not exact residual-GTO Coulomb
integrals. Term-first pair-factor reuse and bounded workspace are binding
requirements in the R3 note.

The approved function must combine `V_aug` with the accepted R3-A augmented
`K` and uncharged `U_A` blocks, and return the existing
`CartesianIDAHamiltonian{Float64}` directly. It must reuse base
`nup`/`ndn`, nuclear charges, and nuclear positions from the same-construction
base Hamiltonian. An arbitrary dimension-compatible Hamiltonian is not an
approved provenance source.

H2 closure value: for the public/base z-axis H2 plus contracted
two-center H/cc-pVTZ `lmax = 1` fixture, the lowest augmented one-body orbital
must have IDA self-Coulomb `0.4574256036192161` within `1.0e-10` under the
compact R3-A residual basis, `sigma = sqrt(2v)` MWG convention, and
weight-aware `G-M` contraction. The earlier target `0.457435475059184` is
superseded for R3-B because it came from a retired private
`[pre_final_pqs, residual_gto]` density-gauge diagnostic. The intermediate
target `0.4574331709135599` is also superseded because it came from direct
parent-density insertion of `G-M` factors rather than final-basis
density-normalized `G-M` contraction.

Do not add a width scale factor and do not relax tolerance to fit the old
scalar. This ID does not approve artifacts, public API expansion,
driver/bin/tool workflow, broad provider payloads, status/result objects,
report fields, pair/assembly workflow expansion, Be2 validation, Cr2
validation, ECP, EGOI, RHF/solver work, wrappers, or a new test file.

The first R3-B endpoint may keep the approved full R3-A moment matrices
`x`/`y`/`z`/`x^2`/`y^2`/`z^2`; this amendment does not replace them with a
diagonal-only moment contract.

## Candidate, Not Approved

Candidate entries record proposed design surfaces for review. They do not
authorize source, test, tool, driver, public API, or artifact changes until
explicitly moved into an approved section.

### HP-R3-ART-01 — compact supplemented artifact provenance — candidate/deferred

R3-C candidate scope, deferred until after accepted R3-A and R3-B numerical
endpoints. The in-memory numerical object remains
`CartesianIDAHamiltonian{Float64}`. A future artifact may add compact
`supplement_provenance/` keys defined in
`r3_residual_gto_mwg_augmentation.md`, derived from the validated R3
construction specification rather than from an in-memory Hamiltonian.

The first compact schema must not store full residual eigenvalue vectors, MWG
center matrices, MWG width matrices, `T_G`, `T_A`, or candidate labels. If
those arrays become consumer-critical, a later explicit residual-basis artifact
group must promote them together. This candidate does not approve a Hamiltonian
wrapper, separate manifest, public provenance reader, HamV6 export, solver
export, or consumer API.

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

### R3 Be2/Cr2 hardening guardrail — deferred

Do not present Be2, Cr2, or broader residual-GTO/MWG supplement support as
approved until a later docs-only amendment closes these items:

- an exact same-construction guard for independently supplied basis, bundle,
  supplement, residual, and base-Hamiltonian objects. Acceptable candidates are
  a single internal construction orchestrator or a numerical consistency check
  such as recomputing `X = G' S A` and validating `T_G + X*T_A`;
- deterministic rank-deficient residual selection under the already-approved
  rank/orientation policy;
- owned-support mixed overlap/operator providers that avoid full-parent
  identity allocation and avoid bounding CPB interior materialization for
  hollow shells;
- bounded or streamed residual MWG term storage for higher residual rank;
- allocation-free or bounded-allocation validation reductions for large dense
  matrices.
