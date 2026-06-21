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
global coefficients, or global self-overlap matrices.

### HP-FILE-01 — terminal realization file

Approved file:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
```

### HP-FN-00 — projected terminal shell realization

Approved as a file-local/internal helper under HP-FILE-01. It realizes one PQS
terminal shell by projecting raw boundary seed columns against all previously
accepted blocks, forming the projected Gram matrix, applying shell-local Lowdin,
and returning effective support plus coefficients.

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
    projection_atol = 1.0e-12,
    weight_atol = 1.0e-14,
)
```

Returns `CartesianTerminalBasisRealization` on success. It owns direct-sector
checks, PQS shell realization, positive final-integral sign canonicalization,
and construction of `CartesianTerminalBasisBlock`.

### HP-FN-02 — cross-block overlap audit

Approved as internal Slice A logic. It computes `C_i' * S_ij * C_j` one block
pair at a time and returns the largest infinity norm. It must not form global
self-overlap or coefficient matrices.

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
