# Mapped-COMX Source Span

Status: approved narrow mainline source-span authority under
`HP-MCOMX-FILE-01`, `HP-MCOMX-OBJ-01`, `HP-MCOMX-FN-01`,
`HP-MCOMX-WIRE-01`, and `HP-MCOMX-TEST-01`.

## Purpose

The high-order lane is an experimental proving ground. It produced the
mapped-COMX source-span idea and scratch evidence that the idea is worth
installing once in mainline so real atom gates can test it through the current
carried-space / PQS / PGDG machinery.

Mainline owns the facility and contract. High-order remains a consumer and
benchmark lane, like CR2: it should call the installed mainline option rather
than maintain a duplicate scratch implementation.

Correction after source review: this is not a new
`CartesianRawProductSources` numerical facility. It is a new source-span option
inside the existing nested doside / COMX path:

```text
pqs_source_axis_transform_facts_from_pgdg_axes(...)
-> _nested_doside_1d(...)
-> _nested_retained_span(...)
-> _cleanup_comx_transform(...)
```

Mapped-COMX changes only the raw one-dimensional span passed into the existing
physical-coordinate COMX cleanup. It must not add a second COMX wrapper, a
parallel axis-transform constructor, or a numerical transform builder under
`CartesianRawProductSources`.

## Numerical Rule

The approved first source-span option is:

```text
protected physical P2
+ mapped Chebyshev enrichment T_k(s_lambda(u))
+ lambda = 0.5
+ no sqrtJ
+ physical-u COMX localization
```

For a one-dimensional interval with physical centers `x`, the nonlinear map
uses a dimensionless local coordinate:

```text
u = (x - x_mid) / x_half
```

where `x_mid` and `x_half` come from the source interval. The map is:

```text
s_lambda(u) = (1 + lambda) * u / (1 + lambda * u^2)
```

The final COMX localization still uses the physical position matrix. Applying
`s_lambda` directly to raw physical centers is not the tested construction and
is forbidden because it makes the source span unit- and location-dependent.

Default source-span specification:

```text
protected_degree = 2
lambda = 0.5
mapped_family = :chebyshev_s
include_sqrt_jacobian = false
localization_coordinate = :physical_u
```

Construction semantics:

1. Build the protected physical polynomial block `1, u, u^2`.
2. Add mapped Chebyshev columns `T_k(s_lambda(u))` until the requested source
   mode count is reached.
3. Project mapped columns against the protected physical block in the local
   parent metric.
4. Orthonormalize the combined source span.
5. Continue through the existing `_cleanup_comx_transform(...)` using the
   physical position matrix.

The first implementation is restricted to `protected_degree = 2`. General
protected degrees require a later parity-balanced mapped-order fill rule; they
must not be implied by blindly adding `T_1`, `T_2`, and so on for every
protected degree.

## Evidence Summary

High-order-manager scratch scans on 2026-06-25 found that pure mapped columns
improve angular center spacing on cube faces, while protected `P2 + T` gives a
better first occupied-mode source span than stronger `P3`/`P4` protection at
small retained counts.

The proposed default `lambda = 0.5` was a balanced setting across the reported
low angular proxy rows:

```text
n_s = 5: L1 = 1.431e-3, L2 = 7.529e-3, Eg = 2.238e-3, T2g = 9.547e-3
n_s = 6: L1 = 3.391e-4, L2 = 1.446e-3, Eg = 2.238e-3, T2g = 3.863e-4
n_s = 7: L1 = 9.702e-5, L2 = 2.818e-4, Eg = 4.105e-4, T2g = 1.417e-4
```

No Hamiltonian, IDA, HF, DMRG, ECP, EGOI, MWG, or Cr2 calculation was used as
approval evidence. Those are later consumer gates after the source-span option
is installed in mainline.

Post-installation high-order-manager He/PQS evidence on 2026-06-26 found that
the installed `n_s = 5` mapped-COMX recipe is not robust enough for
all-electron scalar capture. In the tested He framework,
`P2 + T1(s_lambda), T2(s_lambda)` with `lambda = 0.5` was mechanically wired
through the real driver, but it did not improve the bounded all-electron He
H1/IDA checks and was worse than ordinary PQS on the harder
`core_spacing = 0.2` stress case.

That evidence does not retract the mapped-COMX source-span option. It does
block promoting mapped-COMX, or `n_s = 5` mapped-COMX in particular, as a
default construction rule. Ordinary PQS remains the baseline/default. The next
mapped-COMX evidence gate is bounded He `n_s = 6` and `n_s = 7` H1/IDA testing,
with shell-restricted scalar-capture diagnostics, before any Cr or
molecule-facing promotion.

## Approved IDs

- `HP-MCOMX-FILE-01` - source files for the mainline mapped-COMX source-span
  option.
- `HP-MCOMX-OBJ-01` - compact `MappedCOMXSourceSpec` or equivalent source
  specification object.
- `HP-MCOMX-FN-01` - source-span and axis-transform construction.
- `HP-MCOMX-WIRE-01` - narrow wiring into existing raw product source / PQS
  axis-transform construction.
- `HP-MCOMX-TEST-01` - validation gates.
- `HP-MCOMX-TERM-FN-01` - terminal-basis consumption of carried materialized
  source-axis transform facts.
- `HP-MCOMX-TERM-TEST-01` - terminal seam validation gates.
- `HP-MCOMX-DRV-FN-01` - canonical driver and staged-facade selection of the
  source-span family.
- `HP-MCOMX-DRV-TEST-01` - driver-level validation gates.

## Approved Source Surface

Primary owner:

```text
existing nested doside / COMX source-span seam
```

Approved files:

```text
src/cartesian_nested_faces.jl
src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl
src/cartesian_raw_product_sources/axis_transform_facts.jl
src/cartesian_raw_product_sources/records.jl
```

`src/cartesian_nested_faces.jl` owns the existing low-level doside span and
physical COMX cleanup. `pqs_source_axis_transforms.jl` is approved only for
narrow keyword/spec plumbing into that existing seam and for reporting the
returned source-span facts. `CartesianRawProductSources` files are approved
only for compact provenance/accessors on existing `AxisSourceTransformFact`
records if needed.

No new source file is approved. In particular,
`src/cartesian_raw_product_sources/mapped_comx_source_span.jl` is not an
approved production surface.

## Approved Behavior

The implementation may:

- add a compact internal `MappedCOMXSourceSpec` or equivalent source-span
  selector at the existing doside seam;
- extend `_nested_doside_1d(...)` / `_nested_retained_span(...)` with an
  internal keyword or spec that changes only the raw source-span columns;
- keep the current ordinary span as the default behavior;
- construct mapped-COMX columns from normalized local `u`, not raw physical
  centers;
- project mapped columns against the protected physical block in the local
  parent metric;
- use the existing physical-coordinate `_cleanup_comx_transform(...)` after the
  mapped source span is built;
- return the usual materialized `AxisSourceTransformFact`s compatible with
  existing `RawProductBoxPlan` and PQS boundary product-mode retained rules;
- record compact metadata:
  - `source_span_family = :mapped_comx`;
  - `protected_degree = 2`;
  - `lambda`;
  - `mapped_family`;
  - `mapped_orders`;
  - `include_sqrt_jacobian = false`;
  - `localization_coordinate = :physical_u`;
  - requested/resolved source mode count;
  - rank and overlap/orthogonality diagnostics;
- leave ordinary polynomial source spans available and unchanged;
- expose the option only through internal construction controls needed by the
  approved validation gates.

Internal consumers must not read provenance metadata as a data bus. If a later
source pass needs mapped-order data after construction, the doside/span helper
should return it as a real result field or accessor, and metadata should remain
reporting/provenance.

The Hamiltonian, operator, and artifact layers should consume the same
carried-space / raw product source facts as before. They should not branch on
whether a source span was ordinary polynomial or mapped-COMX, except for
descriptive provenance already carried by source facts.

## Terminal-Basis Wiring

`HP-MCOMX-TERM-FN-01` approves a narrow terminal-basis wiring pass so carried
mapped-COMX `AxisSourceTransformFact`s become basis-defining for PQS shell
realization.

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

The second file is approved only for import/include cleanup if directly
required.

Approved behavior in `_shell_seed(...)`:

- prefer `contract.metadata.raw_product_source_axis_transform_facts` when it
  is present;
- validate exactly three axis facts;
- validate each fact is an `AxisSourceTransformFact`;
- validate `coefficient_status === :materialized`;
- validate source intervals match `support.outer_box`;
- validate source mode dimensions match `source_shape`;
- validate coefficient matrix row/column sizes match interval length and
  source mode dimension;
- build `full_coefficients` from the carried axis coefficient matrices;
- keep the existing boundary mode selection, support restriction,
  shell-local Lowdin, canonicalization, and support validation;
- preserve the ordinary fallback through
  `_nested_projected_q_shell_full_sides(...)` when materialized facts are
  absent.

This lane makes mapped-COMX source-axis facts terminal-basis inputs. It does
not change source-span construction, retained-rule semantics, shell ownership,
or artifact/Hamiltonian behavior.

## Driver Selection

`HP-MCOMX-DRV-FN-01` approves a narrow public construction choice in the
canonical driver and staged base/facade path:

```julia
source_span = :ordinary      # default
source_span = :mapped_comx
```

Approved source files:

```text
bin/cartesian_ham_builder.jl
src/cartesian_base_hamiltonian.jl
src/pqs_source_box_route_driver_helpers.jl
```

`src/pqs_source_box_route_driver_helpers.jl` is approved only for narrow
propagation of the normalized selector to the already-approved source-axis
transform fact path. It must not add route records, route-stage diagnostics,
or new terminal-lowering contracts.

Approved behavior:

- add `source_span` to the canonical driver's editable defaults, trusted input
  file keys, command-line overrides, and compact `print_contract` output;
- normalize and validate `source_span` in the staged base/facade input path;
- accept only `:ordinary` and `:mapped_comx` as driver-facing values;
- treat omitted `source_span` as `:ordinary`;
- pass `:mapped_comx` through the existing PQS source-box path so the mapped
  doside source span produces materialized axis facts consumed by terminal
  realization;
- reject `source_span = :mapped_comx` clearly when `nesting = :wl`, unless a
  later WL-specific source-span amendment approves otherwise;
- preserve ordinary driver behavior and artifact/readback by default.

This is a construction choice, not a diagnostic route switch. It exposes the
same kind of visible, copyable driver control as `nesting`; it must not expose
route skeletons, retained-rule dumps, raw-block switches, stop-after controls,
allocation probes, or high-order benchmark controls.

## Forbidden

This amendment does not approve:

- changing default source spans;
- public API or export changes;
- canonical driver input changes;
- artifact schema, manifest, or reader changes;
- Hamiltonian, one-body, IDA, MWG, Residual Gaussian, raw Gaussian block, or
  solver changes;
- ECP, EGOI, RHF, ED, DMRG, or Cr2 workflow;
- explicit `Y_lm` / angular injection;
- `sqrtJ` weighting;
- mapped-`s` COMX as the production localization gauge;
- applying `s_lambda` to raw physical centers;
- `protected_degree != 2` in the first implementation;
- a parallel mapped-COMX axis-transform route outside the existing
  `_nested_doside_1d(...)` / `_cleanup_comx_transform(...)` path;
- numerical source-span builders under `CartesianRawProductSources`;
- `src/cartesian_raw_product_sources/mapped_comx_source_span.jl`;
- importing high-order branch scaffolding, scripts, route wrappers, status
  objects, diagnostics, or reports;
- duplicate high-order-maintained implementation of the same mainline option;
- changing terminal shell ownership, retained-rule semantics, support
  restriction, shell-local Lowdin, sign canonicalization, or support
  validation under `HP-MCOMX-TERM-FN-01`;
- route records, terminal-lowering changes, route-stage diagnostics, driver
  hooks beyond `source_span`, artifact/schema/manifest/reader changes,
  another COMX path, or high-order workflow controls under
  `HP-MCOMX-DRV-FN-01`;
- committed Cr fixtures or broad high-order benchmark fixtures.

## Validation

`HP-MCOMX-TEST-01` approves only:

- `git diff --check`;
- package load;
- local source-span validation for `n_s = 5`, `6`, and `7`:
  - full retained rank;
  - protected `P2` span preserved;
  - mapped columns use normalized local `u in [-1, 1]`;
  - source columns/centers match the high-order scratch convention within a
    reviewed tolerance;
  - source-axis overlap approximately identity after construction;
  - physical-`u` COMX off-diagonal residual reported;
  - metadata records the approved source-span rule;
- a bounded cubic H and He+ one-electron gate comparing ordinary polynomial
  and mapped-COMX source spans with fixed support and retained count;
- a bounded He `1s^2` fixed-orbital IDA gate if the already-supported
  analytic path can run without new solver or artifact workflow;
- high-order-manager consumer benchmarks on the installed mainline option for
  Cr occupied capture, reported back as evidence rather than committed mainline
  fixtures;
- no Cr2 run.

`HP-MCOMX-TERM-TEST-01` approves later source-pass validation:

- `git diff --check`;
- package load;
- ordinary PQS H2 endpoint/regression unchanged;
- mapped source-span probe still passes;
- focused He or H terminal seam check showing mapped terminal shell
  coefficients differ from ordinary and match the carried materialized axis
  facts;
- H2 supplemented RG endpoint if the touched path crosses it;
- no Cr2 run.

`HP-MCOMX-DRV-TEST-01` approves later source-pass validation:

- `git diff --check`;
- package load;
- default ordinary driver artifact/readback still passes;
- mapped-COMX He or H PQS driver smoke proves carried facts are
  basis-defining;
- ordinary versus mapped He supplemented/MWG/IDA comparison through the real
  driver if bounded;
- H2 RG endpoint still passes;
- no Cr2 run.

No committed test file is approved by default. A later implementation blurb may
name a small standalone script or ignored probe when needed for the approved
gates. Committed fixtures, public driver tests, solver tests, and Cr/Cr2
benchmark fixtures require a separate amendment.

## Failure Rule

If the source option cannot be installed as a small branch in the existing
doside source-span seam before `_cleanup_comx_transform(...)`, make no source
commit and report the exact missing mainline seam.

If implementation requires a new source file, a second COMX wrapper, a
`CartesianRawProductSources` numerical builder, Hamiltonian assembly changes,
artifact schemas, public driver inputs, or high-order-specific workflow, make
no source commit and request a separate design amendment.

If terminal realization cannot consume the carried axis facts without changing
shell ownership, retained-rule semantics, Lowdin realization, artifact schema,
or driver inputs, make no source commit and report the exact blocker.

If making `source_span` driver-selectable requires new route records,
terminal-lowering changes, artifact schema changes, or another COMX path, make
no source commit and report the exact blocker.

If the real atom gates fail after a faithful implementation of the approved
source-span rule, do not tune defaults or add injection in the same pass.
Report the measured failure and request a separate design amendment.
