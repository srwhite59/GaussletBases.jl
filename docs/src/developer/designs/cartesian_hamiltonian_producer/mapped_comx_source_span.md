# Mapped-COMX Source Span

Status: implemented opt-in source-span, terminal-consumption, and canonical
driver facility under the six `HP-MCOMX` source IDs. The three corresponding
test IDs are completed evidence with no continuing permission.

This page is the canonical numerical and workflow contract. The registry owns
individual ID lifecycle and source/test surfaces.

## Purpose

Mapped-COMX is one alternative raw one-dimensional source span inside the
existing nested doside/COMX path:

```text
pqs source-axis transform facts
-> _nested_doside_1d
-> _nested_retained_span
-> _cleanup_comx_transform
-> carried AxisSourceTransformFact
-> PQS terminal shell seed
```

It changes only the columns entering existing physical-coordinate COMX
cleanup. It is not a second COMX wrapper, a high-order workflow branch, a
Hamiltonian correction, or a numerical builder owned by
`CartesianRawProductSources`.

Ordinary source spans remain the default. High-order/physics workflows are
consumers of the installed mainline option, not duplicate owners.

## Numerical Contract

The implemented mapped span is:

```text
protected physical P2
+ mapped Chebyshev enrichment T_k(s_lambda(u))
lambda = 0.5
no sqrt(J)
physical-coordinate COMX localization
```

For source-interval centers `x`, define dimensionless local coordinates:

```text
u = (x - x_mid) / x_half
s_lambda(u) = (1 + lambda) * u / (1 + lambda*u^2)
```

Construction:

1. build the protected physical polynomial block `1,u,u^2`;
2. add mapped Chebyshev columns `T_k(s_lambda(u))` to the requested count;
3. project enrichment columns out of protected `P2` in the local parent
   metric;
4. orthonormalize the combined source span;
5. run existing `_cleanup_comx_transform(...)` with the physical position
   matrix.

Applying `s_lambda` to raw physical centers is forbidden. The nonlinear map
uses normalized `u`; localization remains physical. The implemented first
facility fixes `protected_degree = 2`. Other degrees require separate
parity/order policy.

Compact provenance records source family, protected degree, lambda, mapped
family/orders, no-sqrt-J convention, physical localization coordinate,
requested/resolved count, rank, and overlap/orthogonality diagnostics.
Metadata is descriptive; it is not an algorithmic data bus.

## Terminal Consumption

PQS terminal `_shell_seed(...)` prefers carried materialized
`AxisSourceTransformFact` data when present. It validates:

- exactly three axis facts with materialized coefficients;
- interval agreement with shell support;
- source dimensions and coefficient matrix shapes;
- consistency with the contract source shape.

It then forms full source coefficients from the carried axis matrices and
continues through the existing boundary-mode selection, owned-support
restriction, shell-local Lowdin, sign canonicalization, and support checks.
When carried facts are absent, the ordinary projected-q fallback remains.

The Hamiltonian, one-body, IDA/MWG, residual, and artifact layers consume the
same terminal basis and do not branch on source-span family.

## Driver Convention

The canonical driver and staged base input accept:

```julia
source_span = :ordinary      # default
source_span = :mapped_comx
```

The selector is a construction choice, not a route diagnostic. It is visible
in editable defaults, trusted input keys, overrides, compact contract output,
and construction provenance. `:mapped_comx` is currently PQS-only and
rejects with `nesting = :wl`.

The driver forwards only the normalized symbol. It does not build transforms,
inspect mapped orders, expose high-order controls, or add another route.

## Implemented Ownership

Source-span construction:

- `src/cartesian_nested_faces.jl`;
- narrow source-axis plumbing in
  `src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl`.

Terminal consumption:

- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`;
- module import/include support in
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`.

Driver/staged selection:

- `bin/cartesian_ham_builder.jl`;
- `src/cartesian_base_hamiltonian.jl`;
- `src/pqs_source_box_route_driver_helpers.jl`.

No mapped-COMX production file under `cartesian_raw_product_sources/` is
owned or needed.

## Validation And Current Evidence

Accepted source validation covered:

- local `ns = 5,6,7` rank, protected-P2, normalized-coordinate, overlap, and
  physical-COMX checks;
- ordinary fallback parity;
- mapped terminal coefficients matching carried axis facts and differing from
  ordinary when expected;
- canonical driver artifact/readback and provenance;
- bounded H/He/H2 and supplemented endpoint smokes;
- package load and `git diff --check`.

Post-installation He evidence found `ns = 5`, `P2+T1+T2`,
`lambda=0.5` did not improve the bounded all-electron H1/IDA checks and was
worse than ordinary PQS for the harder `core_spacing=0.2` case. Therefore:

- mapped-COMX remains opt-in;
- ordinary remains the default;
- no `ns=5` physics promotion is justified;
- the next consumer evidence, if pursued, is bounded He `ns=6/7` with
  shell-restricted scalar-capture diagnostics before Cr or molecular claims.

This limitation is a physics result, not a wiring failure.

## Failure Behavior

Stop if the option would require:

- a new source file or second COMX wrapper;
- numerical source-span ownership under `CartesianRawProductSources`;
- changed shell ownership, boundary selection, support restriction, Lowdin,
  or terminal sign conventions;
- Hamiltonian, raw Gaussian, RG/MWG/IDA, artifact, solver, or ECP changes;
- route records, diagnostic payloads, stop-after switches, or high-order
  workflow controls.

If a faithful mapped span fails a physics gate, report it. Do not tune defaults,
add injection, change protected degree, add `sqrtJ`, or change localization in
the same lane.

## Explicit Non-Goals

This contract does not approve:

- default promotion or automatic source-span tuning;
- public exports or API beyond the existing driver/facade construction symbol;
- White-Lindsey mapped source spans;
- explicit spherical-harmonic injection;
- artifacts or reader schema changes;
- EGOI, screened-Hartree, residual, solver, DMRG, ECP, or Cr2 workflow;
- committed Cr/Cr2 fixtures or production claims.
