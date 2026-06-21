# Cartesian Long-Range Roadmap

Status: strategic roadmap, not implementation authority.

This roadmap is based on `main` after the design-governed Cartesian
Hamiltonian producer implementation merged back to main. It records the next
program after the internal base PQS Hamiltonian path was restored.

The compact implementation authority remains:

- `docs/src/developer/designs/cartesian_hamiltonian_producer/README.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/current.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md`
- `docs/src/developer/designs/cartesian_hamiltonian_producer/invariants.md`

This roadmap sequences future work. It does not approve source implementation,
new production surfaces, tests, report fields, metadata keys, public APIs, or
artifact shapes.

## North Star

Geometry and lowering may differ by method, but everything downstream of
localized final-basis construction should converge:

```text
method-specific geometry/lowering
-> localized final basis
-> common operators and corrections
-> common Hamiltonian object and artifact contract
-> common consumer handoff
```

The current base PQS path already reaches:

```text
terminal basis
-> K and by-center unit U_A
-> localized IDA V
-> CartesianIDAHamiltonian
-> existing artifact writer/readback
```

The next era is not recovery of a PQS Hamiltonian algorithm. It is public
producer hardening, method unification, larger-system performance, supplement
and correction migration, and donor retirement.

## Full Cartesian Functionality

Full Cartesian functionality means all of these are available through a small,
documented producer surface:

- base atomic and molecular Hamiltonians;
- independent PQS and supported WL/QW routes;
- residual-GTO/MWG augmentation;
- hydrogenic-core/ESOI and EGOI corrections;
- fragment and counterpoise workflows;
- high-order slab/endcap/panel geometry;
- Cr2-scale consumer readiness;
- one documented public producer with durable artifacts and stable downstream
  contracts.

## Strategic Rules

- The compact Hamiltonian producer design remains the implementation gate.
- New production surfaces still need approved design IDs or a later compact
  design amendment.
- Each migrated feature must name the old authority it replaces.
- Old and new production paths may coexist for at most one milestone unless a
  manager explicitly extends the migration window.
- Line-positive numerical work must be followed by donor retirement.
- Tests and docs count as carrying cost.
- Final Cartesian/WL/QW/PQS source plus tests should be smaller than the R0
  baseline despite broader functionality.
- The `pair terms -> assembly` stages are not required base-public concepts.
  Future pair authority needs a named downstream consumer and an approved local,
  factorized, immediately consumed contract.

## Current Planning Integration

Slice D base materialization handoff is closed under the compact Hamiltonian
producer authority. Do not keep extending "Stage D" as a catch-all for later
work.

Short-term follow-up planning should classify each task by roadmap lane:

- stale authority, route-status, fixture, baseline, and pair/assembly-role
  reconciliation belongs to R0;
- public driver shape, examples, artifact naming, and user-facing producer
  polish belongs to R1;
- WL/QW/PQS downstream sharing belongs to R2;
- residual-GTO/MWG and other non-base Hamiltonians belong to R3;
- Cr2-scale stress, performance, and consumer-readiness gates belong to R6.

Volatile pass sequencing belongs in repo-manager blurbs and the manager running
log, not in this roadmap. Any task that needs a new production surface, report
field, metadata key, public API, artifact shape, committed test, or source file
still needs a compact design amendment before implementation.

## Roadmap

### R0 - Reconcile Authority, Tests, And Baselines

Status: partially complete.

Done:

- Slice D authority ratified in the compact Hamiltonian producer design.
- Stale H2 supplement-preflight fixture/smoke replaced with a base H2 PQS
  Hamiltonian endpoint smoke.
- Design authority compacted into
  `docs/src/developer/designs/cartesian_hamiltonian_producer/`.
- Full June 2026 design and review rounds preserved as history.
- Route migration/status docs, feature donor inventory, and route retirement
  ledger refreshed so they no longer describe base PQS as blocked at
  source-plan or supplement-preflight recovery.
- Pair/assembly role decision recorded in the compact Hamiltonian producer
  authority: `cartesian_pair_terms` and `cartesian_assembly` are not required
  base-public concepts, though useful local-box kernels remain donor/oracle
  inventory for R2/R3 classification.
- Quantitative R0 baseline recorded in
  `docs/src/developer/roadmaps/cartesian_r0_baseline_2026-06.md`.

Still open:

- decide whether any remaining stale roadmap/design references should be
  archived or left as historical breadcrumbs.

Exit condition:

- Current docs, smoke tests, and design authority describe the real base
  Hamiltonian path without preserving stale supplement-preflight or blocked
  source-plan stories.

### R1 - Public Base Producer

Status: candidate design drafted; implementation not approved.

Goal:
Expose a real public producer for base atomic and molecular Cartesian
Hamiltonians.

Scope:

- document the intended public call shape;
- use `r1_public_base_producer.md` as the candidate design for first H/H2
  public base producer review;
- route atomic and molecular base cases through the common final-basis
  Hamiltonian path;
- keep `CartesianIDAHamiltonian` as the output boundary unless a later design
  explicitly approves a different public object;
- define expected artifact names and minimal readback validation.

Exit condition:

- A user-facing H or H2 example constructs a base Hamiltonian through the public
  producer without relying on private route-stage vocabulary.

### R2 - WL/QW And PQS Downstream Unification

Goal:
Unify WL/QW and PQS downstream of final-basis construction.

Scope:

- identify where WL/QW already owns equivalent final-basis data;
- route WL/QW base operators and Hamiltonian construction through the common
  final-basis operator/IDA/Hamiltonian seams where mathematically compatible;
- retire duplicate operator/Hamiltonian paths as parity lands.

Exit condition:

- Base WL/QW and PQS Hamiltonians share the same downstream Hamiltonian
  construction boundary, with method-specific differences confined to
  geometry/lowering/final-basis realization.

### R3 - Generic Residual-GTO/MWG Augmentation

Goal:
Restore supplement augmentation as a generic final-basis augmentation, not an
H2-specific residual-GTO sidecar.

Scope:

- design the supplement representation against the current terminal/final-basis
  boundary;
- keep residual-GTO/MWG near-zero weight behavior separate from the base PQS
  positive-weight gauge;
- avoid report metadata matrices, status payloads, and source-plan wrappers.

Exit condition:

- A reviewed supplement endpoint augments a base Hamiltonian through a generic
  final-basis boundary and deletes its migration donor.

### R4 - Corrections, Branches, Fragments, And Counterpoise

Goal:
Move corrections and chemistry workflow variants onto the common Hamiltonian
producer boundary.

Scope:

- hydrogenic-core/ESOI and EGOI corrections;
- branch-specific Hamiltonian variants;
- fragment and counterpoise workflows;
- clear provenance for charges, centers, fragment labels, and correction
  terms.

Exit condition:

- Corrections and fragment workflows consume the same base Hamiltonian
  structures rather than private pair/materialization frameworks.

### R5 - High-Order Geometry Integration

Goal:
Integrate high-order slab/endcap/panel geometry through shellification and
lowering while preserving the common downstream Hamiltonian path.

Scope:

- keep high-order geometry-specific code upstream of final-basis construction;
- validate q-ladders and shape variants;
- avoid importing old high-order operator duplication when a common final-basis
  operator path exists.

Exit condition:

- High-order geometry can feed the common final-basis Hamiltonian producer for
  reviewed base cases.

### R6 - Cr2 Performance And Consumer Readiness

Goal:
Meet Cr2-scale performance and downstream-consumer gates.

Scope:

- Cr2 stress/performance validation;
- memory and allocation budgets for K/U/V and artifact IO;
- downstream consumer handoff, such as RHF/DMRG setup, without forcing a public
  solver API in this milestone.

Exit condition:

- Cr2-scale Hamiltonian construction is measured, bounded, and usable by the
  next consumer without private manual staging.

### R7 - Public API Stabilization And Final Retirement

Goal:
Stabilize the public Cartesian producer and finish donor retirement.

Scope:

- public docs and examples;
- stable argument names and artifact conventions;
- final removal of migration-only donors and compatibility wrappers;
- final test stratification around scientific endpoints and compact module
  contracts.

Exit condition:

- The public Cartesian producer is documented and stable, and old donor paths
  that no longer serve live validation are deleted or archived as history.

## Authority Map

- Roadmap: strategic sequencing and milestone boundaries.
- Compact Hamiltonian producer design: current implementation authority.
- Algorithm implementation index: reusable kernels and donor navigation.
- Feature-donor inventory and retirement ledger: migration/deletion evidence.
- Manager running log: accepted history and strategic interpretation.
- Manual/reference: public truth once the public producer is stable.

If these documents disagree, do not infer implementation permission from this
roadmap. Update the compact design authority first.
