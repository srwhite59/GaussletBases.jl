# White-Lindsey Terminal Basis Realization

Status: approved narrow source authority under `HP-WLTERM-FILE-01`,
`HP-WLTERM-FN-01`, `HP-WLTERM-WIRE-01`, and `HP-WLTERM-TEST-01`.

## Reason

The canonical driver now exposes `nesting = :pqs` and `nesting = :wl` as
public construction-family choices. The driver/facade plumbing can route
`nesting = :wl` to `route_family = :white_lindsey_low_order`, but the staged
Hamiltonian path still requires a `CartesianTerminalBasisRealization`.

The current terminal-basis seam deliberately drops non-PQS routes:

```julia
recipe.route_family === :pqs_source_box || return nothing
```

This is the real missing seam. White-Lindsey complete shells lower to native
White-Lindsey boundary-stratum concepts, including
`:white_lindsey_boundary_strata` and
`:white_lindsey_boundary_stratum_product`, while the existing terminal
realizer supports direct identity blocks and PQS source-mode shell
realization. The next source pass should make the existing White-Lindsey
low-order route produce the same terminal-basis object consumed by the staged
Hamiltonian producer.

## Approved IDs

- `HP-WLTERM-FILE-01` - optional White-Lindsey terminal-realization file and
  include.
- `HP-WLTERM-FN-01` - White-Lindsey low-order terminal basis realization.
- `HP-WLTERM-WIRE-01` - route helper wiring from the WL route to the terminal
  basis realization.
- `HP-WLTERM-TEST-01` - validation gates for the WL terminal-basis lane.

## Approved Source Surface

Approved files:

```text
src/pqs_source_box_route_driver_helpers.jl
src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl
src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
```

The route helper file is approved only for removing the PQS-only terminal-basis
guard where the WL route has sufficient native terminal lowering and retained
records, and for calling the approved terminal realizer.

`pqs_terminal_basis_realization.jl` may be extended only if the WL lowering
kinds fit the existing terminal-realizer structure cleanly. If the
implementation is clearer as a separate sibling, `HP-WLTERM-FILE-01` approves
creating `white_lindsey_terminal_basis_realization.jl` and adding its include
to `CartesianFinalBasisRealization.jl`. The sibling must return the same
`CartesianTerminalBasisRealization` type and must not introduce a new basis
object, route-stage object, report payload, or public API.

No other `src`, `bin`, `test`, `tools`, or artifact files are approved.

## Approved Behavior

`HP-WLTERM-FN-01` approves only terminal-basis realization for the existing
`:white_lindsey_low_order` route family. The output must be the existing:

```julia
CartesianTerminalBasisRealization
```

with `CartesianTerminalBasisBlock` entries supported on authoritative owned
terminal rows.

Allowed construction:

- preserve direct identity terminal blocks as identity on their owned support;
- realize White-Lindsey boundary-stratum/product terminal blocks only from
  existing construction-native terminal support, retained-rule, and transform
  records;
- produce support-local coefficients on `support.support_indices` /
  `support.support_states`;
- preserve deterministic terminal support, lowering, retained-record, and
  transform-contract order;
- validate disjoint owned supports and shell/block-local identity overlaps
  under the same structural terminal-basis policy as the PQS path;
- use existing module-local helper structure or a small WL-specific sibling
  helper, whichever gives the clearer implementation without new persistent
  vocabulary.

The lane is deliberately not approval to adapt the old White-Lindsey
H1/H1+J materialization path. The common downstream boundary is terminal basis
realization followed by the existing final-basis product, unit-nuclear, IDA,
Hamiltonian construction, artifact, and driver machinery.

## Forbidden

This amendment does not approve:

- changes to route skeleton construction semantics;
- changes to terminal shellification behavior or retained-selection policy;
- old WL H1/H1+J materialization revival or adaptation;
- new Hamiltonian object, artifact schema, manifest/provenance keys, reader
  behavior, or public export/API;
- driver public input changes beyond the already-approved `nesting` field;
- raw-block, Residual Gaussian, MWG/IDA, Qiu-White, supplement, ECP, solver, or
  Cr2 workflow changes;
- report/status/payload expansion, diagnostic route switches, stop-after
  controls, retained-rule dumps, raw-block switches, or route-stage labels;
- committed test files or committed driver input fixtures.

If White-Lindsey boundary-stratum final basis cannot be materialized from
existing terminal lowering, retained-unit, and transform records without
broader route redesign, implementation must make no source commit and report
the exact missing native fact.

## Validation

`HP-WLTERM-TEST-01` approves only:

- `git diff --check`;
- package load;
- current default `nesting = :pqs` atom or H2 base artifact/readback remains
  unchanged;
- `nesting = :wl` base atom artifact/readback;
- `nesting = :wl` base H2 artifact/readback if the existing WL diatomic route
  has sufficient native terminal records;
- H2 residual-GTO/MWG PQS endpoint remains unchanged if terminal realization
  code is touched;
- clear unsupported-input/blocker reporting if WL H2 cannot be realized from
  current native records.

No Cr2 run, supplemented WL run, committed fixture, committed test file,
solver/RHF/ECP/EGOI workflow, artifact schema validation, or broad WL workflow
validation is approved.
