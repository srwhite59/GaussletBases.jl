# Experimental High-Order Doside Nesting Plan

Date: 2026-04-25

Status: active experimental implementation plan, not a production merge plan

## Purpose

Record the bounded repo plan for the experimental high-order doside nesting
lane described in:

- `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/higher-order-doside-sandbox/reports/high_order_doside_nesting_handoff_2026-04-25.md`

This note is the durable repo-side planning memo for that lane.

## Main Judgment

The cleanest repo insertion point is a new quarantined Cartesian experimental
lane parallel to the current ordinary Cartesian IDA layer, not an extension of
the existing face-first nested fixed-block / QW line.

Why:

- the current nested fixed-block line in `src/cartesian_nested_faces.jl` is a
  face-first, support-local, packet-carrying shell contract
- the new lane is instead:
  - full tensor shell, not face-only shell
  - full-block-supported
  - inner-to-outer orthogonalized
  - intended first as an experimental basis-design / He-style validation lane
  - not a QW consumer route
- the repo already has the right parent-grid one-electron and IDA substrate in
  `src/ordinary_cartesian_ida.jl`

The missing seam is therefore not another `_NestedFixedBlock3D` route. It is a
small experimental contracted-pure-Cartesian basis lane.

## Scope Boundaries

This plan is intentionally bounded.

Phase-1 implementation must enforce:

- uniform undistorted `MappedUniformBasis`
- 3D only
- `doside = 5` only
- odd side lengths only
- full tensor shell, not face-only shell
- inner-to-outer orthogonalization
- existing repo one-electron machinery
- IDA only for electron-electron terms

Phase 1 must not broaden to:

- distorted grids
- rectangular boxes
- even side lengths
- arbitrary `doside`
- public API commitments
- QW route integration
- bundle export
- production docs/manual exposure

## Where Current Repo Contracts Fit Cleanly

Good fits:

- `src/ordinary_mapped_backends.jl`
  - `mapped_ordinary_one_body_operators(...)`
- `src/ordinary_cartesian_ida.jl`
  - `_mapped_cartesian_one_body_matrix(...)`
  - parent-grid pair-factor / interaction machinery
  - `_normalized_density_transfer(...)` style density-transfer logic as a
    reference for the two-electron action shape
- `src/ordinary_pgdg.jl`
  - `_cleanup_comx_transform(...)` for local sign/localization cleanup if
    needed
- `src/atomic_ida_two_electron.jl`
  - a model for a tiny explicitly experimental two-electron consumer / solver
    layer

Poor fits:

- `src/cartesian_nested_faces.jl`
  - wrong top-level shell contract
- `src/ordinary_qw_operator_assembly.jl`
  - wrong consumer line and too much route baggage
- `_NestedFixedBlock3D`
  - this experiment is not another fixed-block/QW route

## Proposed File Placement

Add two new internal files:

- `src/cartesian_high_order_doside_experimental.jl`
- `src/cartesian_high_order_doside_ida_experimental.jl`

They should be included in `src/GaussletBases.jl` immediately after
`src/ordinary_cartesian_ida.jl`.

These files should remain internal-only in phase 1:

- no public exports
- no manual/reference entries
- no public route promises

Test placement:

- `test/ordinary/high_order_doside_experimental_runtests.jl`
- included from `test/ordinary/runtests.jl`

Scratch validation drivers:

- `tmp/work/high_order_doside_heplus_validation.jl`
- `tmp/work/high_order_doside_he_smoketest.jl`

## Proposed Internal Data Structures

The experimental lane should use its own small explicit data structures.

Suggested shape:

```julia
struct _ExperimentalHighOrderAxisData1D
    centers::Vector{Float64}
    weights::Vector{Float64}
    overlap::Matrix{Float64}
    position::Matrix{Float64}
    kinetic::Matrix{Float64}
    gaussian_factors::Vector{Matrix{Float64}}
    pair_factors_1d::Vector{Matrix{Float64}}
end

struct _ExperimentalHighOrderBlock1D
    side::Int
    interval::UnitRange{Int}
    coefficients::Matrix{Float64}
    local_overlap::Matrix{Float64}
    local_position::Matrix{Float64}
    local_weights::Vector{Float64}
    local_centers::Vector{Float64}
    localized_centers::Vector{Float64}
end

struct _ExperimentalHighOrderTensorShell3D
    side::Int
    interval::NTuple{3,UnitRange{Int}}
    full_block_coefficients::SparseMatrixCSC{Float64,Int}
    shell_coefficients::SparseMatrixCSC{Float64,Int}
    shell_labels::Vector{NTuple{3,Int}}
    shell_kind_counts::NamedTuple
    column_range::UnitRange{Int}
end

struct ExperimentalHighOrderDosideStack3D
    parent_basis::MappedUniformBasis
    backend::Symbol
    doside::Int
    sides::Vector{Int}
    parent_side::Int
    coefficient_matrix::SparseMatrixCSC{Float64,Int}
    block_column_ranges::Vector{UnitRange{Int}}
    block_labels::Vector{Symbol}
    shell_layers::Vector{_ExperimentalHighOrderTensorShell3D}
    contracted_weights::Vector{Float64}
    diagnostics::NamedTuple
end

struct _ExperimentalHighOrderDosideIDAData3D
    stack::ExperimentalHighOrderDosideStack3D
    expansion::CoulombGaussianExpansion
    parent_overlap::Matrix{Float64}
    contracted_overlap::Matrix{Float64}
    one_body_hamiltonian::Matrix{Float64}
    parent_interaction::Matrix{Float64}
end
```

Design intent:

- `ExperimentalHighOrderDosideStack3D` is the primary basis object for the
  experimental lane
- `_ExperimentalHighOrderDosideIDAData3D` is the narrow He-style validation
  carrier
- neither object should try to pretend to be `_NestedFixedBlock3D`

## Constructor and Driver Shape

Primary stack constructor:

```julia
_experimental_high_order_doside_stack_3d(
    basis::MappedUniformBasis;
    backend::Symbol = :numerical_reference,
    doside::Int = 5,
    sides::AbstractVector{<:Integer} = [5, 7, 9, 11],
)
```

Phase-1 constructor rules:

- `basis` must be uniform and undistorted
- cubic parent only
- `doside == 5`
- `sides[1] == 5`
- all `sides` odd
- `sides` strictly increasing by `2`
- `maximum(sides) == length(basis)` in the first pass

Validation drivers:

```julia
_experimental_high_order_doside_ida_data_3d(stack; expansion, Z)
_experimental_high_order_doside_heplus_energy(stack; expansion, Z = 2.0)
_experimental_high_order_doside_he_singlet_lanczos(
    stack;
    expansion,
    Z = 2.0,
    krylovdim = ...,
    maxiter = ...,
)
```

## Core Construction Plan

For side ladder `[5, 7, 9, 11, ...]` with fixed `doside = 5`:

1. Build one normalized 1D axis data bundle from the chosen backend.
2. For each side, build one centered 1D `5`-function local contraction.
3. Tensor the 1D contraction to a full `5^3 = 125` local block.
4. Label tensor columns by `(ix, iy, iz)` with each coordinate in `1:5`.
5. Define shell columns as all tensor labels with at least one coordinate in
   `{1, 5}`.
6. Define residual columns as all labels with every coordinate in `2:4`.
7. Use the side-`5` full block as the initial stack.
8. For each larger side:
   - take only the `98` shell columns
   - orthogonalize them against the accumulated inner stack in the parent
     overlap metric
   - Lowdin-orthonormalize the residual shell space
   - sign-fix / center-order consistently
   - append to the stack

The implementation should reuse:

- `_nested_interval_data(...)`
- `_nested_retained_span(...)`
- `_cleanup_comx_transform(...)`

as local algebra helpers where useful, but it should not reuse the face-first
shell orchestration types as the host contract.

## One-Electron Reuse Plan

Use the repo one-electron machinery directly on the parent grid:

1. Build:

   ```julia
   one_body = mapped_ordinary_one_body_operators(
       basis;
       exponents = expansion.exponents,
       center = 0.0,
       backend = backend,
   )
   ```

2. Build parent matrices:

   ```julia
   S_parent, h_parent = _mapped_cartesian_one_body_matrix(
       one_body,
       expansion;
       Z = Z,
   )
   ```

3. Project into the stack basis:

   ```text
   S_stack = C' * S_parent * C
   h_stack = C' * h_parent * C
   ```

This is the intended contract:

- one-electron is handled by the mature repo path
- the experiment changes only the basis
- phase 1 should require `S_stack ≈ I` as a validation invariant

## Electron-Electron IDA Path

Phase 1 should use IDA only for the electron-electron term.

Use a dense parent-grid two-index kernel:

```text
W = parent-grid ordinary Cartesian IDA interaction matrix
```

and apply the two-electron action on a symmetric coefficient matrix `A` by:

```text
P_parent = C * A * C'
V(A) = C' * (W .* P_parent) * C
H(A) = h*A + A*h + V(A)
```

Important constraints:

- no four-index tensor
- no attempt to present this as a general many-electron framework
- matrix-free action is acceptable for He validation
- dense parent `W` is acceptable in phase 1

This keeps the electron-electron lane explicitly experimental while using the
existing ordinary Cartesian IDA substrate.

## Invariants and Tests

Structural invariants:

- `C' * S_parent * C ≈ I`
- `rank(C) == expected dimension`
- stack span matches full-block union span
- shell counts are correct
- contracted weights are finite and sign-fixed
- center drift remains reasonable
- moment diagnostics through degree `4` are finite

Expected dimensions:

- for 3D `doside = 5`, one full block has `125`
- one shell contributes `98`
- for sides `[5, 7, 9, 11]`, the expected dimension is `419`

Required structural tests:

- side `[5]`: dimension `125`
- sides `[5, 7]`: dimension `223`
- sides `[5, 7, 9, 11]`: dimension `419`
- shell-kind counts per outer layer:
  - faces `54`
  - edges `36`
  - corners `8`
- stack-vs-union span agreement on the bounded `3D [5, 7, 9, 11]` ladder

Two-electron action tests:

- `V(A)` preserves symmetry for symmetric `A`
- Frobenius Hermiticity:
  - `<A, V(B)> == <V(A), B>` on random symmetric inputs
- no NaN/Inf under representative ladders

## Validation Sequence

### Gate 1: structural basis validation

Use `backend = :numerical_reference`, parent side `11`, sides `[5, 7, 9, 11]`.

Required:

- exact expected dimensions
- near-identity overlap defect
- shell counts correct
- stack-vs-union span recovers

### Gate 2: He+

Use projected repo one-electron Hamiltonian only.

Required:

- finite energies
- shell additions variationally nonincreasing
- no obvious pathologies in overlap / conditioning

### Gate 3: He singlet

Use:

- projected one-electron repo Hamiltonian
- parent-grid IDA electron-electron action only

Required:

- energies improve smoothly with shell additions
- moderate spacing / parent-box perturbations do not cause obvious instability
- no severe conditioning or outer-shell domination by a tiny set of bad moment
  directions

## Failure Criteria / Stop Criteria

Stop immediately if:

- structural dimension or rank is wrong
- `C' S C` is not near identity
- shell-kind counts are wrong
- span-vs-union validation fails
- moments or contracted weights become non-finite

Stop the He+ lane if:

- energies are not finite
- shell additions are not variationally nonincreasing
- stack and full-union reference disagree materially on the same parent box

Stop the He lane if:

- the two-electron action fails symmetry / Hermiticity checks
- shell additions do not improve smoothly
- modest spacing/box perturbations cause unstable electron-electron behavior
- severe conditioning or near-null modes dominate

## What Must Remain Explicitly Experimental

Keep clearly quarantined:

- all new types and drivers should use `_experimental_...` naming
- no public export
- no bundle export
- no `basis_representation(...)` integration unless needed later
- no attempt to merge into the current nested fixed-block / QW public surface
- no attempt to claim distorted-grid or general-`doside` support

## Repo-Doer Execution Plan

### Chunk 1: structural stack core

Implement:

- `src/cartesian_high_order_doside_experimental.jl`
- 1D block build
- full 3D block / shell split
- inner-to-outer orthogonalized stack
- structural diagnostics

Tests:

- dimensions
- shell counts
- orthogonality
- span-vs-union on `[5, 7, 9, 11]`

Stop if any structural gate fails.

Suggested commit message:

- `Add experimental high-order doside stack core`

### Chunk 2: projected one-electron lane

Implement:

- projection of repo one-electron parent matrices into the stack basis
- He+ validation helper / scratch driver

Tests:

- `S_stack ≈ I`
- He+ finite and variational over shell additions

Suggested commit message:

- `Project one-body operators into experimental doside stack`

### Chunk 3: ee-only IDA He lane

Implement:

- `src/cartesian_high_order_doside_ida_experimental.jl`
- dense parent `W`
- matrix-free `A -> H(A)` action
- tiny He singlet Lanczos driver

Tests:

- action symmetry
- He shell-addition smoothness
- basic spacing / box stability smoke

Suggested commit message:

- `Add experimental He IDA action for high-order doside stack`

## Current Bottom Line

This experiment should be treated as:

- a bounded Cartesian research lane
- built on top of existing repo one-electron machinery
- using IDA only for electron-electron terms
- validated first on He+ and He

It should not be merged into the current public nested Cartesian contract
until the electron-electron IDA behavior is shown to be stable enough on the
bounded phase-1 cases.
