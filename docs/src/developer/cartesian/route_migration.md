# Cartesian Route Migration

This is the current developer note for the Cartesian/PQS route on `main`.

## Current State

- The thin Cartesian/PQS route has been promoted to `main`.
- The old `demolition/pqs-thin-route` branch name is retired.
- Old Cartesian helper/schema/status test groups were deleted.
- Validation now uses driver ladders, not private helper payload tests.
- Old Cartesian code remains feature-donor inventory until its features are
  migrated, explicitly abandoned, or proven unused.

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

## Current Driver Lines

| line | current driver inputs | purpose | expected coarse facts | deletion/merge condition |
|---|---|---|---|---|
| `wl_atomic` | `test/driver_inputs/he_wl_q5_pure_gausslet_h1.jl`; `test/driver_inputs/he_wl_q5_gto_h1.jl` | Preserve WL atomic pure-gausslet and supplement-capable driver entry. | Driver recognizes WL; parent/basis construction runs; H1 materialization is finite when requested; GTO path reaches its current intended stage. | Delete or fold into main matrix once WL atomic is fully represented by the common route and no separate line ladder adds information. |
| `wl_diatomic` | `test/driver_inputs/h2_wl_q5_pure_gausslet_h1.jl`; `test/driver_inputs/h2_wl_q5_gto_h1.jl` | Preserve WL diatomic route capability while WL and PQS are merged into common staged construction. | Driver recognizes WL diatomic; parent axes and route stages execute; H1/GTO stages fail only at real missing construction objects. | Delete or fold into main matrix once diatomic WL uses the common route without old donor wrappers. |
| `pqs_atomic` | `test/driver_inputs/he_pqs_q5_wlmap.jl`; `test/driver_inputs/he_pqs_q5_gto.jl` | Keep an atomic PQS/source-box route smoke while PQS support/source construction is generalized. | Driver recognizes PQS; route executes through the current atomic stage; GTO case reaches supplement staging or a real missing object. | Delete or fold into main matrix once atomic PQS shares the same route stages as the diatomic PQS path. |
| `pqs_diatomic` | `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl`; `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl`; `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_materialized.jl` | Protect the independent H2 PQS source-box line, supplement preflight, and H2 PQS Ham/Basis plus residual-GTO sidecar materialization. | Independent route remains fake-free; final dimension is 471 for the materialized H2 PQS sidecar case; overlap error is finite/small; H1 is finite/symmetric with finite lowest energy; H1-J self-Coulomb is finite/positive; basis/ham sidecar artifacts are written and reloadable; the materialized artifact includes narrow one-body-only provider blocks with `provider_blocks_included = :one_body_only`. | Delete or fold into main matrix once independent PQS and residual-GTO sidecar materialization are ordinary route functionality with compact consumer checks. |

## Feature Donor Inventory

The current feature-donor migration table lives in:

- [Feature donor inventory](feature_donor_inventory.md)

The highest-priority donor features are:

1. P1 residual-GTO / MWG supplement materialization: first sidecar slice and
   narrow one-body-only provider-block artifact completed; density/pair/H1-J
   provider blocks and supplemented values remain.
2. P1 Ham/JLD2 artifact contract and basis transfer/roundtrip: H2 PQS
   Ham/Basis plus residual-GTO sidecar artifacts, including the one-body-only
   provider-block fields, write and roundtrip; broader downstream consumer
   coverage and density/pair/H1-J provider-block artifact coverage remain.
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
