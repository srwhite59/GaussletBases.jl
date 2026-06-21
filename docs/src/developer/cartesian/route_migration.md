# Cartesian Route Migration

This is the current developer note for the Cartesian/PQS route on `main`.

## Current State

- The thin Cartesian/PQS route has been promoted to `main`.
- The old `demolition/pqs-thin-route` branch name is retired.
- Old Cartesian helper/schema/status test groups were deleted.
- Validation now uses driver ladders, not private helper payload tests.
- Old Cartesian code remains feature-donor inventory until its features are
  migrated, explicitly abandoned, or proven unused.
- The base H2 PQS Hamiltonian endpoint now materializes
  `CartesianIDAHamiltonian{Float64}` and validates artifact readback through
  `tools/h2_pqs_base_hamiltonian_smoke.jl`.
- Residual-GTO/MWG supplements and other non-base Hamiltonians are no longer
  the immediate PQS route status; they are future roadmap work.
- `cartesian_pair_terms` and `cartesian_assembly` remain in the public stage
  sequence, but their role in the base public workflow is an open R0 decision.

The target architecture is a thin, staged driver:

```text
input file
-> system
-> recipe
-> parent
-> shells
-> units
-> transforms
-> pair terms
-> assembly
-> report
-> materialization
-> print/save
```

Required construction objects should be built and passed as objects. Missing
objects should fail plainly at the construction point rather than flowing
through `status`, `readiness`, `available`, `blocker`, or `probe` payloads.

## Validation Policy

Run the full Cartesian smoke matrix with:

```text
julia --project=. tools/run_cartesian_driver_ladder.jl
```

Run current line-specific ladders with:

```text
julia --project=. tools/run_cartesian_line_ladder.jl --line=wl_atomic
julia --project=. tools/run_cartesian_line_ladder.jl --line=wl_diatomic
julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_atomic
julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic
```

These are route smoke validators. They are not `Test.jl` helper/schema suites.
They should answer whether a surviving Cartesian sub-line can still execute the
driver far enough to prove the line exists and to expose the first real missing
construction object.

Use `--dry-run` or `--list` for quick command-discovery checks when a pass does
not need numerical validation.

Run the current H2 base PQS Hamiltonian endpoint smoke with:

```text
julia --project=. tools/h2_pqs_base_hamiltonian_smoke.jl
```

That smoke validates the materialized `CartesianIDAHamiltonian`, H2 final
dimension, H1 lowest value, localized IDA self-Coulomb value, and artifact
readback. It is the current endpoint evidence for the recovered base H2 PQS
path.

## R0 Status Audit

This audit classifies stale pre-recovery claims for the roadmap R0 closure lane.

Update docs now:

- this route migration/status page, especially the `pqs_diatomic` row;
- [Feature donor inventory](feature_donor_inventory.md), especially
  residual-GTO/MWG and artifact rows that previously implied no base PQS
  producer exists;
- [Cartesian route retirement ledger](../cartesian_route_retirement_ledger.md),
  with a current-status note separating old replacement-pressure bullets from
  the recovered base Hamiltonian endpoint.

Leave as historical:

- `docs/src/developer/pqs_manager_running_log.md`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/history/`;
- `docs/src/developer/designs/cartesian_hamiltonian_producer/reviews/`;
- archived demolition notes, handoffs, blurb logs, and old pass responses.

Active pressure for later R0/repo-manager work:

- record the quantitative R0 baseline: line counts, H2 smoke timing/allocation,
  light separated-diatomic one-body/IDA timing if practical, and current Cr2
  stress status;
- decide whether `cartesian_pair_terms` and `cartesian_assembly` gain real
  downstream authority or collapse out of the base public workflow;
- keep Cr2 as stress/performance and consumer-readiness evidence, not the next
  base-PQS correctness gate.

## Current Driver Lines

| line | current driver inputs | purpose | expected coarse facts | deletion/merge condition |
|---|---|---|---|---|
| `wl_atomic` | `test/driver_inputs/he_wl_q5_pure_gausslet_h1.jl`; `test/driver_inputs/he_wl_q5_gto_h1.jl` | Preserve WL atomic pure-gausslet and supplement-capable driver entry. | Driver recognizes WL; parent/basis construction runs; H1 materialization is finite when requested; GTO path reaches its current intended stage. | Delete or fold into main matrix once WL atomic is fully represented by the common route and no separate line ladder adds information. |
| `wl_diatomic` | `test/driver_inputs/h2_wl_q5_pure_gausslet_h1.jl`; `test/driver_inputs/h2_wl_q5_gto_h1.jl` | Preserve WL diatomic route capability while WL and PQS are merged into common staged construction. | Driver recognizes WL diatomic; parent axes and route stages execute; H1/GTO stages fail only at real missing construction objects. | Delete or fold into main matrix once diatomic WL uses the common route without old donor wrappers. |
| `pqs_atomic` | `test/driver_inputs/he_pqs_q5_wlmap.jl`; `test/driver_inputs/he_pqs_q5_gto.jl` | Keep an atomic PQS/source-box route smoke while PQS support/source construction is generalized. | Driver recognizes PQS; route executes through the current atomic stage; GTO case reaches supplement staging or a real missing object. | Delete or fold into main matrix once atomic PQS shares the same route stages as the diatomic PQS path. |
| `pqs_diatomic` | `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl`; `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_base_hamiltonian.jl` | Protect the independent H2 PQS source-box line and the recovered base Hamiltonian endpoint. | H2 reaches terminal basis realization, materializes `CartesianIDAHamiltonian{Float64}`, and validates H1/self-Coulomb/artifact readback through the endpoint smoke. Cr2 remains later stress/performance evidence, not the immediate correctness blocker. Supplement materialization is deferred roadmap work. | Delete or fold into main matrix once the base public producer covers H2 PQS without private route-stage vocabulary and pair/assembly role is settled. |

## Feature Donor Inventory

The current feature-donor migration table lives in:

- [Feature donor inventory](feature_donor_inventory.md)

The highest-priority donor features are:

1. P1 residual-GTO / MWG supplement materialization: H2 residual-GTO
   materialization is no longer the current base-PQS route blocker. Remaining
   work is a generic final-basis supplement design, arbitrary electron/spin and
   non-H2 source dimensions, external consumer agreement/smoke, performance
   review, and donor-wrapper deletion after shared kernels and the public
   producer absorb the capability.
2. P1 Ham/JLD2 artifact contract and basis transfer/roundtrip: the Ham side has
   a public Cartesian IDA writer/reader and the H2 base PQS endpoint writes and
   reads that Hamiltonian artifact. Any future basis/provenance artifact should
   be consumer-driven rather than reviving the deleted private H2 sidecar.
   Private H1-J/RHF solver diagnostics are not part of the producer contract.
3. P2 hydrogenic-core / ESOI corrections.
4. P2 EGOI / density-density correction.
5. P3 branch / fragment Hamiltonian workflow.
6. P3 high-order slab/endcap/panel geometry variants.
7. P4 legacy nested fixed-source oracles, only if still needed as references.

## Guardrails

- Do not resurrect deleted status/readiness/probe/helper-schema tests.
- Do not add exact print-string tests.
- Do not add new helper-payload schema tests.
- Do not preserve old donor lines after their feature is migrated or
  intentionally abandoned.
- Old Cartesian code is feature-donor inventory, not target architecture.
- Each donor feature should get a tiny driver ladder case before or while it is
  migrated.
- Once a feature is migrated, delete the donor line or mark it intentionally
  retained with a clear reason.

## Historical Note

The thin route was promoted to `main` after a non-adabatic demolition branch
removed old route-state scaffolding, flat payload/report field clouds, and
Cartesian helper/schema tests. The driver ladder became the validation
authority for this private route. Historical details are archived in:

- [PQS thin route demolition history](../archive/pqs_thin_route_demolition_history.md)
