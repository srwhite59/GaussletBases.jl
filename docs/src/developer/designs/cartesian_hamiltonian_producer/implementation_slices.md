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
- Residual Gaussian migration cleanup: residual basis construction, exact
  augmented operators, and residual MWG interaction now live under
  `src/cartesian_residual_gaussians/`; keep deleting old R3 wrappers when live
  callers move;
- compact supplemented artifact writing and facade parsing remain outside
  `CartesianResidualGaussians` as terminal/facade workflow glue unless a later
  amendment identifies a real duplication or consumer need;
- Cr2 stress/performance and any full Cr2 supplemented Hamiltonian run;
- measurement-only Cr2-readiness forecasting, consumer workflow beyond the
  approved H2/Be2 R3 usability facade, and basis/supplement-realism decisions
  until separately approved;
- non-base Hamiltonian variants.

## Residual Gaussian Module Migration

Status: residual-basis construction, exact augmented operator transformation,
and matched-width Gaussian residual interaction are migrated to the
`CartesianResidualGaussians` domain module.

Canonical algorithm authority:

- `residual_gaussian_domain_module.md`

Current compatibility boundary:

- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` remains
  the terminal/facade compatibility file for existing R3 entry points, artifact
  writing, and facade wiring;
- moved physics helpers should not be reintroduced there;
- RG does not own `supplement_provenance/`, JLD2 artifact workflow, facade
  input parsing, basis loading, parent lattice construction, or public exports.

Completed cleanup:

- residual-basis helpers moved to `residual_basis.jl`;
- exact `[G,A] -> [G,R]` operator transform moved to
  `augmented_operators.jl`;
- residual MWG descriptor and interaction assembly moved to
  `mwg_interaction.jl`;
- standalone R3 test no longer depends on old test-only R3-B wrapper names.

Remaining cleanup:

- delete any future compatibility wrapper once the exact live caller moves;
- keep artifact and facade hooks outside RG unless a separate design amendment
  approves moving them.

## Cartesian Gaussian Raw Blocks - Nuclear Slice

Status: approved for implementation; not part of R3/RG public workflow.

Approved boundary:

- neutral owner files under `src/cartesian_gaussian_raw_blocks/`;
- exact uncharged by-center Cartesian Gaussian nuclear raw blocks:
  parent-supplement `G-A` and supplement-supplement `A-A`;
- analytic one-dimensional nuclear factor construction, unique nuclear
  coordinate reuse, upper-triangular `A-A` assembly with mirroring,
  function-local scratch reuse, and term-first contraction.

Implementation sequence:

1. Extract the current Residual Gaussian/Qiu-White nuclear `G-A`/`A-A`
   behavior into the neutral owner without changing numerical conventions.
2. Rewire Residual Gaussian and Qiu-White callers to consume the neutral
   kernel.
3. Delete the duplicate route-local nuclear loops once parity is established.
4. Optimize allocation inside the neutral owner only after extraction parity.

Validation gates:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged as an ignored performance/usability
  measurement if needed;
- Cr2 q4 exact nuclear blocks match the current implementation at roundoff as
  an ignored measurement only;
- a small Qiu-White nuclear parity fixture passes.

Forbidden in this slice:

- overlap, kinetic, coordinate moments, second moments;
- pair factors or matched-width Gaussian interaction;
- terminal projection or residual-basis transformation;
- Qiu-White route objects or parent construction;
- persistent caches, metadata, reports, status fields, artifacts, public API,
  Cr2 facade support, or Cr2 artifact workflow.

Follow-on low-level optimization approval:

- `HP-CGAI-FN-01` approves only
  `src/cartesian_gaussian_axis_integrals.jl` for an in-place
  `_cartesian_gaussian_axis_integral_table!(...)` helper and
  `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl` for consuming that
  helper in the neutral nuclear raw-block kernel;
- the bottleneck addressed is one-dimensional analytic axis-integral table
  allocation measured in Cr2 q4 nuclear raw blocks, not Residual Gaussian
  selection, Qiu-White route objects, or wrapper/result allocation;
- the later source pass must preserve table parity, H2/Be2/Qiu-White/Cr2
  nuclear parity, and materially reduce Cr2 q4 nuclear raw-block allocation
  from the `44552.840 MiB` baseline without adding caches, metadata, route
  objects, reports, artifacts, public API, or semantic changes.
