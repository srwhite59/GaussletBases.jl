# Implementation Slices

## Slice A — Terminal Basis Realization

Status: implemented.

Implemented boundary:

- `CartesianTerminalBasisBlock`
- `CartesianTerminalBasisRealization`
- `pqs_terminal_basis_realization(...)`
- terminal-basis wiring from `cartesian_transforms`

Corrected construction contract:

```text
full source-box product modes
-> boundary product-mode column selection
-> restrict rows to support.support_indices / support.support_states
-> shell-local Gram on that owned support
-> symmetric shell-local Lowdin
-> final sign canonicalization
-> append block with unchanged owned support
```

Previous-block projection, recursive projection, and effective-support growth
onto previous terminal regions are rejected. Cross-block overlap is structurally
zero because parent gausslet rows are orthonormal and terminal regions own
disjoint parent rows. A nonzero structural overlap means duplicated support
rows, incorrect row restriction, wrong support ownership, or an indexing error;
it is not a physical residual to compute or repair.

Cross-block kinetic, nuclear-attraction, and IDA interactions may still be
nonzero and remain assembled over terminal block pairs.

Validation gates used:

- one-center atomic terminal basis through the shared entry point;
- H2 terminal basis dimension `471`;
- Cr2 terminal basis dimension `4291` during design/implementation validation;
- positive final integrals;
- structurally disjoint terminal supports and shell-local identity overlaps;
- no global Lowdin.

Deletion obligations completed:

- terminal source-realization preflight path removed;
- source-plan summary mirrors removed;
- old blocked-route smoke assertions reduced.

Deferred gates:

- larger performance/stress work is not part of Slice A.

Implementation correction required:

- remove `_subtract_previous` and recursive projection;
- make `_shell_seed` use `support.support_indices` /
  `support.support_states`, not all `outer_box` rows;
- remove `projection_atol` plumbing when it has no remaining construction use;
- replace production cross-overlap audit / `max_cross_overlap` with structural
  checks: every block support equals its authoritative terminal support,
  terminal support sets are pairwise disjoint, and each shell-local overlap is
  identity;
- H2 realized block supports should return to local counts such as `(275, 362,
  578)`, not cumulative supports like `(275, 637, 1215)`;
- rerun Slice A/B/C and R1 validation after source correction.

The recursive-projection Be2 measurements are superseded. Post-`d2bf139c` Be2
measurements are the valid optimization baseline.

## Slice B — Final-Basis One-Body Operators

Status: implemented.

Implemented boundary:

- `assemble_terminal_product_operator!(...)`
- file-local term-first Gaussian-sum nuclear attraction helper

Validation gates used:

- one-center H one-body baseline around `-0.49855234726272035`;
- H2 one-body lowest energy `-0.79460371733658908`;
- light separated-diatomic N2 validation for topology/performance;
- finite/symmetric `K` and `U_A` matrices;
- 64 MiB local workspace cap.

Deletion obligations completed:

- dense direct identity allocation removed from terminal overlap helper;
- no production global support one-body operator was introduced.

Deferred gates:

- Cr2 one-body stress/performance validation.

## Slice C1 — Localized IDA Matrix

Status: implemented.

Implemented boundary:

- `assemble_terminal_ida_interaction!(...)`

Validation gates used:

- H2 one-body energy remains at the Slice B value;
- final IDA weights positive and finite;
- `electron_electron_ida` finite and symmetric;
- reviewed H2 self-Coulomb `0.4569117646737212`.

Deletion obligations completed:

- no CPBM route revival;
- no IDA payload/cache/status framework.

Deferred gates:

- non-base/supplement IDA variants.

## Slice C2 — CartesianIDAHamiltonian Construction

Status: implemented.

Implemented boundary:

- existing `CartesianIDAHamiltonian` construction from `K`, uncharged
  by-center `U_A`, localized IDA `V`, charges, positions, and electron counts.

Validation gates used:

- H2 final dimension `471`;
- constructed type `CartesianIDAHamiltonian{Float64}`;
- `one_body_hamiltonian(ham)` reproduces direct Slice B H1;
- H2 one-body lowest `-0.79460371733658908`;
- self-Coulomb from `ham.electron_electron_ida`
  `0.4569117646737212`.

Deferred gates:

- public workflow polish and larger molecule stress.

## Slice D — Base Materialization Handoff

Status: HP-WIRE-02 approved and implemented.

Implemented boundary:

```julia
cartesian_materialization(report, terminal_basis_realization, materialization_inputs)
```

Validation gates used:

- no-request path returns `nothing`;
- requested H2 path returns `CartesianIDAHamiltonian{Float64}`;
- final dimension `471`;
- H2 one-body lowest `-0.79460371733658908`;
- H2 self-Coulomb `0.4569117646737212`;
- optional artifact write/readback with one-body delta `0.0`.

Deletion obligations completed:

- materialization wrapper/result-kind path removed for the PQS base handoff;
- blocked physical-gausslet source-plan/report mirrors removed;
- physical-gausslet target/supplement payload chain removed;
- H2 and Cr2 stage probes no longer treat blocked source-plan summaries as
  route authority.

Deferred gates:

- public-driver polish;
- public export or driver workflow beyond the approved non-exported R3
  usability facade;
- Cr2 stress/performance and any full Cr2 supplemented Hamiltonian run;
- measurement-only Cr2-readiness forecasting, consumer workflow beyond the
  approved H2/Be2 R3 usability facade, and basis/supplement-realism decisions
  until separately approved;
- non-base Hamiltonian variants.
