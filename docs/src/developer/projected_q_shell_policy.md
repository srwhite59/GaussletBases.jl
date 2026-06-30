# Projected q-Shell Policy

This file is retained as a compatibility pointer for historical links. The
current PQS source-box framework lives in
[`pqs_source_box_operator_framework.md`](pqs_source_box_operator_framework.md),
with implementation authority under
`docs/src/developer/designs/cartesian_hamiltonian_producer/`.

Durable facts preserved from the old policy:

- PQS raw operator authority starts from a filled source CPB and boundary
  COMX/product-mode selection.
- Shell support is owned support, not a CPB; shell realization is a later
  projection plus shell-local Lowdin step.
- The useful count picture is
  `dim PQS(q, L) = q^2 L - (q - 2)^2 * (L - 2)`, so cubic `PQS(5, 5)` has
  retained count `98`.
- PGDG analytic pair-factor paths remain required unless an explicit
  debug/reference path is selected.
- Artifact provenance that historically cited this file should treat basis
  identity as a status-bearing construction/source-column label, not
  `center_xyz`; centers are representative metadata only.

The removed long-form history about first-gate evidence, private smoke tests,
sidecar prototypes, and route-shadow migration was transitional narrative, not
current authority.
