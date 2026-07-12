# Terminal Basis And Base Assembly

Status: implemented internal producer infrastructure under `HP-OBJ-01`,
`HP-OBJ-02`, `HP-FILE-01`, `HP-FN-00`, `HP-FN-01`, `HP-FN-02`,
`HP-WIRE-01`, `HP-FN-03`, `HP-FN-04`, and `HP-FN-05`.

This page is the canonical contract for the terminal-basis objects, PQS
terminal realization, blockwise exact one-body assembly, localized IDA matrix,
and final `CartesianIDAHamiltonian` construction boundary. The registry owns
permission, lifecycle, and exact source/test surfaces. Route geometry,
shellification, artifacts, residual augmentation, and public driver behavior
remain in their own contracts.

## Numerical Boundary

The foundational producer sequence is:

```text
typed terminal support, retained, and transform records
    -> support-local terminal basis blocks
    -> block-pair kinetic and unit-nuclear matrices
    -> block-pair localized IDA interaction
    -> CartesianIDAHamiltonian
```

The terminal basis is block sparse as a representation. Its exact one-body and
IDA matrices are generally dense. Structural support orthogonality must not be
misread as block-diagonal operator physics.

## Terminal Basis Objects

The implemented objects are exactly:

```julia
struct CartesianTerminalBasisBlock
    unit_key::Symbol
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    coefficients::Union{Nothing,Matrix{Float64}}
    column_range::UnitRange{Int}
end

struct CartesianTerminalBasisRealization
    blocks::Vector{CartesianTerminalBasisBlock}
    final_dimension::Int
    max_cross_overlap::Float64
end
```

For each block:

- `unit_key` preserves the terminal unit identity and deterministic order.
- `support_indices` and `support_states` are the authoritative owned parent
  rows, in matching order.
- `coefficients === nothing` means direct identity on those support rows.
- Otherwise `coefficients` has one row per owned support row and one column per
  retained final function.
- `column_range` is the block's contiguous native range in the final basis.

`final_dimension` is the sum of the realized block column counts. The object
does not own a global coefficient matrix, global overlap, parent bundle,
route-stage report, or artifact metadata.

`max_cross_overlap` remains an exact live field for compatibility with the
implemented object shape, but structural support correction made it legacy
debt. PQS construction currently returns `0.0`. It is not a physical residual,
a quality score, or permission to project one block into another. Removing the
field requires separate source authority.

## Owned Supports And Structural Overlap

Terminal support records partition authoritative parent rows. Every realized
block must match its corresponding support record exactly, contain no duplicate
row, and be disjoint from every earlier block.

Because parent gausslet rows are orthonormal and terminal supports are
disjoint, cross-block overlap is zero by construction. A nonzero structural
overlap means one of:

- duplicated support rows;
- incorrect row restriction;
- wrong support ownership;
- inconsistent state/index ordering;
- an indexing error.

It is not a physical residual to minimize. Previous-block projection,
recursive projection, an accumulated projection basis, and effective-support
growth onto earlier terminal regions are rejected.

Cross-block kinetic, nuclear-attraction, and IDA matrix elements can still be
nonzero. Those operators must be assembled over all terminal block pairs.

## PQS Terminal Realization

The implemented entry point is:

```julia
pqs_terminal_basis_realization(
    support_records,
    retained_records,
    transform_contracts,
    bundles;
    identity_atol = 1.0e-8,
    weight_atol = 1.0e-14,
)
```

There is no current `cross_atol` keyword. Cross-block overlap is structural,
not a computed tolerance gate.

The realizer consumes support, retained-rule, transform-contract, and mapped
axis facts already owned by earlier stages. It does not infer shell geometry or
retained policy.

### Direct Identity Blocks

A direct record keeps `coefficients = nothing`. Before accepting it, the
realizer checks that the support-local overlap is identity and that all
product IDA weights are finite and positive. The block remains implicit
identity on exactly its owned rows.

### PQS Shell Blocks

For a PQS shell, the durable sequence is:

```text
construct full source-box product coefficients
    -> generate boundary COMX product-mode columns
    -> validate retained/source ordering and dimensions
    -> restrict rows to support_indices/support_states
    -> form the shell-local Gram matrix
    -> apply inv(sqrt(Symmetric(Gram)))
    -> canonicalize column signs from final product weights
    -> append the block on unchanged owned support
```

Full source-box coefficients may come from carried materialized axis facts or
the established projected-shell product construction. In either case, they
only generate source columns before restriction to owned rows. They do not
enlarge terminal support.

The shell-local Gram must produce an identity overlap within `identity_atol`.
Each final column must have a finite nonzero product weight so its sign can be
canonicalized deterministically.

### Compact Thin Slabs

The same PQS realizer also accepts the implemented compact thin-slab transform
contract. Coefficient construction remains owned by the thin-slab lowering
contract; this boundary validates its support-local identity, preserves the
owned rows, and appends it in terminal order. This page does not redefine slab
geometry or lowering policy.

Unsupported transform kinds, missing contracts, retained/source mismatches,
support mismatches, duplicate rows, nonidentity local overlaps, and
ungaugeable final weights are construction failures.

`HP-CHANGE-01` is rejected as standalone authority. Returning a shell overlap
may exist only as a private implementation detail of `HP-FN-00` inside its
already approved source surface. It creates no independent helper, result
field, object, module, or source permission.

## Route Wiring

`cartesian_transforms` reaches terminal realization through the current helper
in `src/pqs_source_box_route_driver_helpers.jl`.

For `:pqs_source_box`, it passes the typed terminal support plan, retained
records, transform contracts, and parent axis bundle directly to
`pqs_terminal_basis_realization(...)`. It does not reconstruct terminal facts
from reports or summaries.

White-Lindsey uses the same `CartesianTerminalBasisRealization` object through
its separate terminal realizer. Its boundary-stratum semantics are owned by
[White-Lindsey terminal basis realization](white_lindsey_terminal_basis_realization.md),
not by this PQS algorithm.

## Blockwise One-Body Assembly

The general product-matrix entry point is:

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

It evaluates every upper-triangular terminal block pair and mirrors the result.
Direct blocks act as implicit selectors. Compact blocks apply their
support-local coefficient matrices on the left and/or right. No global parent
operator or global final-basis coefficient matrix is formed.

Kinetic energy is the sum of the three product terms:

```text
T_x (x) S_y (x) S_z
S_x (x) T_y (x) S_z
S_x (x) S_y (x) T_z.
```

Unit nuclear attraction uses the file-local term-first Gaussian-sum
accumulator. The Coulomb expansion term is the inner reduction over reusable
one-dimensional factor matrices. Each center produces one uncharged matrix
`U_A = -1/r_A`; physical charges are applied only when assembling H1.

Support-pair actions are tiled under the established 64 MiB local workspace
cap. Dense direct identity matrices and persistent one-body caches are not part
of this contract.

The destination and factors must have compatible dimensions. Factors and
coefficients must be finite; one-dimensional factor matrices must be symmetric
within their active numerical checks.

## Localized IDA Assembly

The implemented IDA entry point is:

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

For each terminal block, the final IDA weights are:

```text
direct block:   product support weights
compact block:  coefficients' * product support weights.
```

Every final weight must be finite and greater than `weight_atol`. The raw
Gaussian-expanded pair action is assembled blockwise and divided elementwise
by the left and right final weights. Upper-triangular block pairs are mirrored,
and both diagonal blocks and the complete destination must satisfy the symmetry
gate.

The result is the final-basis two-index `electron_electron_ida` matrix. This is
the localized IDA density-density convention. It is not a four-index tensor,
not a cache packet, and not a matrix to transform later with `C' V C`.

## Cartesian IDA Hamiltonian Boundary

The existing Hamiltonian object is:

```julia
struct CartesianIDAHamiltonian{T}
    kinetic::Matrix{T}
    nuclear_attraction_unit_by_center::Vector{Matrix{T}}
    electron_electron_ida::Matrix{T}
    nup::Int
    ndn::Int
    nuclear_charges::Vector{T}
    nuclear_positions::Matrix{T}
    nuclear_repulsion::T
end
```

Construction uses the existing `CartesianIDAHamiltonian(...)` constructor.
Kinetic, every unit-nuclear matrix, and `electron_electron_ida` must be finite,
symmetric, square, and dimensionally identical. Electron counts must be valid
for the basis. Charges and positions must be finite and center-aligned.
Nuclear repulsion is derived from the stored charges and positions.

The physical one-body matrix is assembled on demand as:

```text
H1 = K + sum_A Z_A U_A.
```

The current staged base producer obtains terminal products, unit-nuclear
matrices, and IDA through `cartesian_base_products`,
`cartesian_base_unit_nuclear`, and `cartesian_base_vee`, then calls
`cartesian_base_hamiltonian_assembly(...)`. Optional artifact writing delegates
to the existing writer and artifact contracts. No wrapper result or route
report is part of the Hamiltonian construction boundary.

## Historical Slice D Warning

The former route-driver wrapper:

```text
cartesian_materialization(report, terminal_basis_realization,
    materialization_inputs)
```

was removed by `e2e164e9b` under the completed materialization-retirement
authority. It is not a compatibility interface. Do not restore the wrapper,
its report/save choreography, or adapters around its old shape. Current code
uses staged producer functions, direct `CartesianIDAHamiltonian` construction,
and the existing artifact path. See
[route-driver materialization retirement](route_driver_materialization_retirement.md).

## Ownership

| Contract | Source owner |
| --- | --- |
| Terminal objects and PQS realization | `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl` |
| Exact product and Gaussian-sum one-body assembly | `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl` |
| Localized IDA assembly | `src/cartesian_final_basis_realization/pqs_terminal_ida.jl` |
| Terminal-stage route wiring | `src/pqs_source_box_route_driver_helpers.jl` |
| Base product/unit-nuclear/IDA composition | `src/pqs_source_box_low_order_materialization.jl`, `src/cartesian_base_hamiltonian.jl` |
| Hamiltonian object and one-body/nuclear accounting | `src/cartesian_ida_hamiltonian.jl` |

## Validation

Current bounded validation is owned by:

- `test/driver_public/cartesian_base_hamiltonian_runtests.jl` for base H/H2,
  compact/high Coulomb construction, finite/symmetric matrices, endpoint
  values, and artifact/readback behavior;
- `test/ida/cartesian_ida_hamiltonian_runtests.jl` for object construction,
  one-body and nuclear-repulsion accounting, ownership, and readback;
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` as a downstream
  consumer of product, unit-nuclear, and IDA assembly.

The current public base gate includes:

```text
H lowest H1        -0.49877574806444014 Ha
H2 final dimension 487
H2 lowest H1       -0.79460371733658908 Ha
H2 self-Coulomb     0.4569012290840094 Ha
```

Validation also requires exact compact omitted/explicit parity, finite and
symmetric high-accuracy matrices, package load, and `git diff --check` when
these surfaces change. Endpoint interpretation follows the terminal
due-diligence contract.

## Explicit Non-Goals

These foundational IDs do not approve:

- terminal shell geometry, retained-selection, or White-Lindsey policy changes;
- recursive projection or global Lowdin repair;
- a global coefficient matrix or global parent/final overlap;
- residual-Gaussian, MWG augmentation, protected, EGOI, or screened-Hartree
  behavior;
- Coulomb-policy changes, new caches, reports, status objects, or stage fields;
- artifact schema, public driver, solver, ECP, or Cr2 workflow changes;
- restoration of retired Slice D wrappers;
- interaction rotation such as `C' V C`.
