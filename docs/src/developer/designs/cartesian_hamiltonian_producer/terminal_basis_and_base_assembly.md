# Terminal Basis And Base Assembly

Status: implemented internal producer infrastructure under `HP-OBJ-01`,
`HP-OBJ-02`, `HP-FILE-01`, `HP-FN-00`, `HP-FN-01`, `HP-FN-02`,
`HP-WIRE-01`, `HP-FN-03`, `HP-FN-04`, and `HP-FN-05`. Pass 526 completed the
direct support-local shell-seed replacement; `HP-FN-00` is in maintenance.
Pass 532 completed the bounded call-local terminal-buffer optimization below;
`HP-FN-03` is again in maintenance.

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
resolve ordinary or carried source-axis coefficient matrices
    -> generate and validate retained boundary COMX modes
    -> assemble support_states x retained_modes directly
    -> form the shell-local Gram matrix
    -> apply inv(sqrt(Symmetric(Gram)))
    -> canonicalize column signs from final product weights
    -> append the block on unchanged owned support
```

The three axis coefficient matrices may come from carried materialized axis
facts or the established ordinary projected-shell side construction. They are
consumed by one local loop: retained modes remain in
`modes.mode_indices`/`modes.column_indices` order, rows remain in
`support.support_states` order, and each selected axis entry is converted to
`Float64` at the existing route boundary. If any selected axis entry is zero,
the coefficient is literal `0.0`; otherwise it is evaluated exactly as
`vx * vy * vz`. The state order must continue to match
`support.support_indices`.

No parent-row by complete-source-mode coefficient matrix is part of terminal
realization. Pass 526 replaced that global sparse tensor-product construction
and restriction with direct allocation of only the final support-row by
retained-mode seed inside `_shell_seed`. It deleted
`_shell_seed_full_coefficients_from_axis_facts` and
`_shell_seed_full_coefficients` without changing `_nested_product_coefficients`,
`_realize_shell`, `_support_action`, the retained-rule validator, Lowdin, sign
canonicalization, or the independent post-Lowdin metric check.

Accepted clean comparisons covered ordinary and mapped source paths for both
one-center and bond-aligned diatomic block sets. Every coefficient matrix was
byte-identical. Supports, retained modes, unit keys, column ranges,
fingerprints, topology, dimensions, captures, energies, warnings, and column
accounting were unchanged.

Under the order-controlled Julia `1.12.6` protocol, eight shell seeds improved
from `134--216 ms` and `969--1,028 MB` to `2.56--2.61 ms` and
`13.26--13.51 MB`. Warm PQS construction improved from `0.622--0.729 s` and
`1.756--1.764 GiB` to `0.280--0.384 s` and `0.874 GiB`. Fresh-process and
post-White-Lindsey allocations also fell substantially, while the smaller
fresh-process wall-time reduction confirms that cold-compilation specialization
is an independent optimization question.

The accepted validation used the unchanged public Cartesian base, focused H2+
release, mapped-COMX, and source-q owners; the bounded
`core,ida,cartesian,examples` groups; authority check/self-test; documentation
tests and build; and all three remote numerical gates. No permanent baseline
helper, committed performance fixture, test, API, metadata, cache, or file was
added. Cold-compilation specialization, scalar-loop optimization, Gram-policy
changes, and compatibility cleanup require separate authority.

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
in `src/cartesian/pqs_source_box_route_driver_helpers.jl`.

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

### Implemented Call-Local Terminal Buffer Reuse

Pass 532 implemented one storage-only optimization in
`src/cartesian/cartesian_final_basis_realization/pqs_terminal_one_body.jl` and its base
caller in `src/cartesian/pqs_source_box_low_order_materialization.jl`. The existing
action, tile, and block buffers may be sized once from the live terminal blocks
and reused within one base-product or unit-nuclear construction. The three
kinetic product terms may share one lexical buffer set, and all physical
nuclear centers may share a separate lexical buffer set. Each buffer remains a
plain call-local matrix held through the current internal buffer mechanism; it
must be fully overwritten before use and must not alias a live input,
destination, accumulator, or another simultaneously live scratch view.

The existing 64 MiB tile bound, terminal-block order, three-term kinetic order,
upper-triangle/mirroring behavior, `mul!` calls, scaling, and destination
accumulation order are unchanged. In particular, the exact 135-term scalar loop
in `_fill_terminal_gaussian_sum_action!` remains byte-for-byte outside this
authority: term order, multiplication association, signs, and `Float64`
conversion points may not change.

Buffer ownership is lexical and nonescaping. No global or module-mutable
storage, task-local state, pool, lock, retained stage/result field, metadata,
cache, public workspace, or new workspace type is permitted. Independent
constructions must own disjoint buffers. If safe reuse requires a public API,
another source owner, persistent state, a new file, arithmetic reordering, or
disproportionate machinery, implementation stops without a source commit.

The accepted Julia `1.12.6` implementation preserved every matrix bit. PQS
kinetic construction changed from `0.637 s / 2.741 GB` to
`0.435 s / 0.140 GB`; PQS unit-nuclear construction changed from
`10.553 s / 1.926 GB` to `10.171 s / 0.267 GB`. Combined stage allocation fell
by `3.968 GiB`. The warmed complete comparison changed from
`32.774 s / 6.816 GB` to `32.028 s / 2.536 GB`, a `2.27%` wall-time
improvement. These measurements are acceptance evidence, not cross-machine
performance thresholds.

Before performance interpretation, the old and candidate kinetic matrices and
every by-center unit-nuclear matrix must be bitwise identical for the matched
PQS and White-Lindsey constructions. Parent and terminal fingerprints,
dimensions, `275 + 960 + 50` accounting, energies, captures, eigen-residuals,
symmetry, topology, due-diligence warnings, and all release tolerances must
remain unchanged. Validation uses transient comparison code only, package
load, the existing public Cartesian and matched-H2+ owners, normal bounded CI,
authority/docs gates, and explicit terminal due-diligence inspection.

The implementation changed the two authorized files by `+25/-10` lines, with
no test, file, API, type, metadata, or cache addition. Repeated zero-sized
buffer initialization was removed from the production paths. At that closeout,
scalar-loop optimization, cold-compilation/provenance cleanup, Gram policy,
compatibility deletion, route construction, release work, and the
duplicate-example question remained separate.

### Implemented Four-Element Gaussian-Sum Reduction

Pass 534 authorized, and commit `94ec277d954b5435a04b0ad68ae352c95b0434c7`
implemented, replacement of only the exact 135-term scalar hot loop in
`_fill_terminal_gaussian_sum_action!`. The implementation precomputes one
call-local `Float64` table containing
`coefficients[term] * fx[term, ix, jx]` per complete Gaussian-sum action or
accumulation call and reuses it across the existing block pairs. That table is
plain lexical scratch: it does not escape, enter a stage result, or become a
cache, field, metadata item, public workspace, or shared mutable object.

Each complete batch processes exactly four independent support-pair output
elements with four explicit scalar accumulators and ordinary three-index array
access. A final zero-to-three-element remainder stays scalar. For each output
element, terms remain in `1:nterms` order and each contribution retains the
same left-associated `((coefficients * fx) * fy) * fz` `Float64` operations,
conversion points, and addition order. Batching changes only when independent
outputs are visited. It does not permit tuples, eight-lane batching, explicit
full unrolling, term-major traversal, layout-dependent offsets, or another
algorithm. The old scalar inner loop was deleted; no parallel fallback remains.

The exact implementation delta is `+50/-14` lines in
`src/cartesian/cartesian_final_basis_realization/pqs_terminal_one_body.jl`, reaching but
not exceeding the authorized `50`-line hard limit. It added no helper, tuple
machinery, type, file, API, test, cache, metadata, dependency, or other owner
edit. The design evidence is the transient report with SHA-256
`6c33a5d264c92a45cd2d66e419ceacca901259463a6515dc6cc0a7c8fc8b703a`.
Its controlled scratch measurements were bitwise exact and changed isolated
PQS time from `4.9875` to `2.7444 s` and White-Lindsey time from `4.1457` to
`2.1969 s`; preweighting added `491,600` call-local bytes and fresh compilation
about `0.15 s`. Production validation then matched the frozen SHA-256 hashes
of both nuclei's PQS and White-Lindsey unit-nuclear matrices bitwise. Dimensions,
fingerprints, energies, captures, eigen-residuals, symmetry, topology,
due-diligence warnings, column accounting, and release tolerances were
unchanged.

Paired Julia `1.12.6` production measurements improved warmed isolated PQS
from `5.021` to `2.687 s` (`46.50%`) and White-Lindsey from `4.125` to
`2.182 s` (`47.09%`). The warmed complete comparison improved from `32.286`
to `24.020--24.102 s` (`25.35--25.60%`), and fresh comparison improved from
`60.770` to `51.911--52.961 s` (`12.85--14.58%`). Maximum additional fresh
compilation was `0.134 s`. Scalar accumulation remained allocation-free;
call-local preweighting used `491,600` bytes, below the `1 MiB` gate.

Validation passed the augmented owner `464/464`, matched-H2+ release owner
`18/18`, public Cartesian owner `232/232`, residual-GTO/MWG owner `80/80`,
docs `110/110 + 10/10`, package load, authority check/self-test, Documenter,
and diff checks. CI run `33023091521` and Docs run `33023091519` passed. The
matched fixture retained `21 x 21 x 29` parent axes, dimension `12789`,
terminal dimension `1285`, `10`-bohr padding, `275 + 960 + 50` columns, eight
shells, and two slabs. The inspected augmented endpoint retained `9 x 9 x 15`
axes, `4`-bohr padding, and terminal dimension `487`.

Eight-lane batching and the measured no-more-than-`2%` remaining loop gain are
not a follow-on target. Further loop restructuring, cold aggregate-wrapper
compilation, compatibility cleanup, Gram changes, route work, release work,
and duplicate-example policy remain separate. The next performance
investigation should examine the cold reporting boundary under new authority.

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

This existing object is a compatibility bundle: `nup` and `ndn` select the
intended electron sector but do not participate in terminal-basis or operator
construction. At fixed nuclei, basis input, and numerical policy, changing
only those counts must leave `kinetic`, every unit-nuclear matrix,
`electron_electron_ida`, the assembled physical one-body matrix, and nuclear
repulsion exactly unchanged. A separate operator/problem type split is future
architecture, not part of this boundary.

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
| Terminal objects and PQS realization | `src/cartesian/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl` |
| Exact product and Gaussian-sum one-body assembly | `src/cartesian/cartesian_final_basis_realization/pqs_terminal_one_body.jl` |
| Localized IDA assembly | `src/cartesian/cartesian_final_basis_realization/pqs_terminal_ida.jl` |
| Terminal-stage route wiring | `src/cartesian/pqs_source_box_route_driver_helpers.jl` |
| Base product/unit-nuclear/IDA composition | `src/cartesian/pqs_source_box_low_order_materialization.jl`, `src/cartesian/cartesian_base_hamiltonian.jl` |
| Hamiltonian object and one-body/nuclear accounting | `src/cartesian/cartesian_ida_hamiltonian.jl` |

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
