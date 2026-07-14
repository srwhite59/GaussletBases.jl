# Historical High-Order Doside Experimental Plan

Status: historical evidence for the QW/high-order cluster approved for
retirement under `HP-RETIRE-QW-DONOR-*`. This page grants no source authority
and is not a current implementation plan.

## Original Question

The experiment tested whether a higher-order Cartesian nested construction
could improve core resolution by combining:

- a doside refinement ladder;
- full tensor-shell basis (FSB) increments rather than only the full basis
  union (FBU);
- inner-to-outer overlap-metric orthogonalization;
- endcap and panel treatments around elongated molecular support;
- projected one-body operators and an ee-only IDA control.

The implementation lived primarily in:

- `cartesian_high_order_doside_experimental.jl`;
- `cartesian_high_order_doside_ida_experimental.jl`;
- `ordinary_qw_experimental_paths.jl`;
- `cartesian_qw_operator_carried_spaces.jl`.

## Evidence Retained

The bounded experimental sequence established:

- a working doside ladder and tensor-shell construction;
- a meaningful FSB/FBU distinction;
- successful undistorted He+ and He controls;
- consistent projected one-body and ee-only IDA behavior in that control
  regime;
- bounded Cr capture measurements;
- useful endcap/panel concepts that later migrated into active adjacent
  geometry and retained-unit code.

It did not establish that the large implementation improved the distorted
mapped-parent regime where the approach would need to be useful. The
distorted-parent benefit remained unresolved, and no current producer or test
consumer adopted the experimental exported APIs.

## Retirement Meaning

The four-file implementation is retired rather than kept as an oracle. Detailed
numerical history remains in this page's follow-up and validation notes and in
Git history. New source must not depend on the retiring files while deletion is
pending.

The retirement preserves:

- endcap/panel and doside helper concepts already owned by
  `cartesian_nested_faces.jl`, `cartesian_nested_owned_units.jl`, and
  `cartesian_nested_diatomic.jl`;
- chain/square basis types and constructors in
  `ordinary_qw_types_and_bases.jl`;
- active QW residual, raw-block, operator, and hybrid-representation kernels.

Most importantly, implementation retirement does **not** abandon homonuclear
linear atomic chains as a scientific target. A future chain route requires a
new design grounded in the preserved basis/geometry contracts and current
producer architecture; it must not restore this experimental cluster by
default.

See
[QW and high-order experimental cluster retirement](designs/cartesian_hamiltonian_producer/qw_high_order_experimental_retirement.md)
for exact source surfaces, caller evidence, guardrails, and validation.
