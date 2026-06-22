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
-> existing minimal artifact writer plus R1 producer provenance when requested
```

This is internal base-Hamiltonian authority plus the narrow approved R1 public
base producer surface recorded in `r1_public_base_producer.md` and
`registry.md`. The visible driver shape may call the implemented path, but this
design does not approve a new artifact format except the `HP-R1-ART-01`
`producer_provenance/` keys in the final Hamiltonian file, solver integration,
broad driver redesign, or public workflow outside the R1 H/H2 scope.

Current implementation boundary:

- One-center atomic and bond-aligned diatomic terminal plans share the same
  terminal-basis realization entry point once typed terminal support, retained,
  and transform records exist.
- Terminal basis realization is block-local. A PQS shell uses the full source
  box only to generate boundary product-mode columns, then restricts rows to
  the shell-owned `support_indices` / `support_states` before shell-local Gram
  construction, symmetric Lowdin, final sign canonicalization, and appending
  the block with unchanged owned support.
- Previous-block projection, recursive projection, projection-basis repair, and
  effective-support growth onto previous terminal regions are forbidden.
  Cross-block overlap remains an audit; large overlap after shell-owned
  realization is a parent metric or shell construction problem, not a reason to
  mix coefficients into previous supports.
- `cartesian_transforms` owns terminal basis realization for supported PQS
  terminal plans.
- `cartesian_materialization(report, terminal_basis_realization,
  materialization_inputs)` receives `transforms.terminal_basis_realization`
  directly.
- Requested base PQS materialization returns `CartesianIDAHamiltonian{Float64}`
  directly. No-request materialization returns `nothing`.
- Optional base-Hamiltonian artifact writing uses the existing
  `write_cartesian_ida_hamiltonian` shape.

Base pair/assembly role decision:

- The future base public workflow should be:

  ```text
  system / specification
  -> parent and route geometry
  -> terminal basis realization
  -> Hamiltonian production
  -> artifact
  ```

- `cartesian_pair_terms` and `cartesian_assembly` are not required
  base-public concepts. The current base Hamiltonian construction path already
  uses terminal basis realization, blockwise `K` and unit `U_A`, term-first
  localized IDA `V`, and direct `CartesianIDAHamiltonian` construction.
- No new base-route consumer should be added to `cartesian_pair_terms` or
  `cartesian_assembly`.
- The existing stages may remain temporarily for legacy script and report
  compatibility until R1 rewires the public facade.
- Their direct report dependency is narrow: `cartesian_assembly` currently
  exists chiefly so `cartesian_report` can recover `route_skeleton` and a
  low-order shellification summary. That is not numerical assembly authority.
- Pair modules remain donor/oracle inventory pending R2/R3 file-level
  classification. Useful local product-box and 1D factor kernels should move
  to the module that owns their scientific consumer rather than justify empty
  public stages.
- Future pair authority requires an explicitly approved, factorized, local,
  consumer-owned contract with a scale/workspace model and immediate numerical
  consumption. Metadata-only all-pairs inventories, status frameworks, and
  payload graphs are not future pair authority.
- Quantitative R0 baselines should be recorded before deleting or rewiring
  these stages.

Deferred lanes:

- public-driver polish and examples outside the approved R1 origin-centered H
  and z-axis H2 base producer scope;
- Cr2-scale stress and performance validation;
- residual-GTO/MWG supplements and other non-base Hamiltonians;
- solver integration;
- White-Lindsey pair-framework completion;
- distorted-product COMX realization;
- EGOI or other Hamiltonian corrections.

Normal startup reading for this lane is this file, `registry.md`,
`invariants.md`, `implementation_slices.md`, and
`docs/src/developer/algorithm_implementation_index.md`.
