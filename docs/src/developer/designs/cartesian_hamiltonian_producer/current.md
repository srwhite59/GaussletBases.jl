# Current Authority

Status: Slice A, Slice B, Slice C1, Slice C2, and Slice D base handoff are
implemented for the internal base PQS Hamiltonian producer.

This authority covers the base all-electron PQS path:

```text
terminal support and retained contracts
-> terminal localized final basis
-> final-basis kinetic matrix
-> final-basis unit nuclear-attraction matrices by center
-> localized IDA electron-electron matrix
-> CartesianIDAHamiltonian
-> existing minimal artifact writer
```

This is internal base-Hamiltonian authority, not public API polish. The visible
driver shape may call the implemented path, but this design does not approve a
new public workflow, new artifact format, solver integration, or broad driver
redesign.

Current implementation boundary:

- One-center atomic and bond-aligned diatomic terminal plans share the same
  terminal-basis realization entry point once typed terminal support, retained,
  and transform records exist.
- `cartesian_transforms` owns terminal basis realization for supported PQS
  terminal plans.
- `cartesian_materialization(report, terminal_basis_realization,
  materialization_inputs)` receives `transforms.terminal_basis_realization`
  directly.
- Requested base PQS materialization returns `CartesianIDAHamiltonian{Float64}`
  directly. No-request materialization returns `nothing`.
- Optional base-Hamiltonian artifact writing uses the existing
  `write_cartesian_ida_hamiltonian` shape.

Deferred lanes:

- public-driver polish and examples;
- Cr2-scale stress and performance validation;
- residual-GTO/MWG supplements and other non-base Hamiltonians;
- solver integration;
- White-Lindsey pair-framework completion;
- distorted-product COMX realization;
- EGOI or other Hamiltonian corrections.

Normal startup reading for this lane is this file, `registry.md`,
`invariants.md`, `implementation_slices.md`, and
`docs/src/developer/algorithm_implementation_index.md`.
