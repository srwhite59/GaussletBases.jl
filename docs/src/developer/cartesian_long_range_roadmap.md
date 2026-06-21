# Cartesian Long-Range Roadmap

Status: **current strategic roadmap**

Baseline commit: `02c97685` (`Delete physical-gausslet payload chain`)

This page is the durable project-level plan for the Cartesian line. It turns the
long-term goals in the Cartesian/PQS manager log into a sequenced program for a
smaller, unified, fully functional, public Cartesian producer.

This roadmap is **not implementation authority**. New production types, files,
functions, stage fields, artifact fields, or status vocabulary still require the
applicable frozen design document and the `AGENTS.md` approval process. The
roadmap says what outcome the project is pursuing and in what order.

Use this page as the current strategic plan. The initiation-era medium-term goals
in `pqs_manager_running_log.md` remain useful history, but independent H2 PQS
recovery is no longer the active frontier.

See also:

- [Architecture and current direction](architecture.md)
- [Algorithm implementation index](algorithm_implementation_index.md)
- [Cartesian Hamiltonian producer design](cartesian_hamiltonian_producer_design.md)
- [Cartesian route migration](cartesian/route_migration.md)
- [Cartesian feature donor inventory](cartesian/feature_donor_inventory.md)
- [Cartesian route retirement ledger](cartesian_route_retirement_ledger.md)
- [Cartesian source/build unification plan](cartesian_source_build_unification_plan.md)
- [Performance review contracts](performance_review_contracts.md)
- [Test-suite reorganization plan](test_suite_reorganization_plan.md)

## 1. North star

The Cartesian line should provide a scientifically transparent, efficient,
public workflow from physical system specification to a stable Hamiltonian
artifact:

```text
physical system + route/basis specification
-> geometry and parent basis
-> disjoint terminal support
-> route-specific lowering and retained rules
-> localized final Cartesian basis
-> K + separated unit U_A + localized IDA V
-> optional supplements and corrections
-> CartesianIDAHamiltonian
-> versioned artifact
-> CR2 / HFDMRG / DMRG or other consumer
```

Atomic, diatomic, chain, lattice, White--Lindsey/Qiu--White, PQS, and high-order
routes may differ where the mathematics really differs: geometry,
shellification, lowering, and retained-space construction. Once a localized
final basis exists, one-body construction, IDA construction, Hamiltonian
assembly, artifact semantics, and consumer handoff should be common unless a
measured algorithmic reason requires a different kernel.

The project is done only when this workflow is both usable and smaller than the
collection of donor, compatibility, and migration paths it replaces.

## 2. Definition of full Cartesian functionality

“Full functionality” means all of the following are available through one
coherent architecture, with unsupported combinations documented rather than
represented by blocker payloads.

1. **Base all-electron Hamiltonians**
   - one-center atomic;
   - contact-core and separated bond-aligned diatomic;
   - independent PQS and supported White--Lindsey/Qiu--White construction;
   - blockwise kinetic, separated unit nuclear attraction, and localized IDA.

2. **Supplemented bases**
   - residual-GTO / MWG augmentation for atomic and molecular systems;
   - mixed gausslet/GTO and GTO/GTO operator blocks;
   - the same final Hamiltonian and artifact boundary as the base route.

3. **Corrections and branches**
   - hydrogenic-core / ESOI corrections;
   - EGOI or stationary-Fock correction where scientifically justified;
   - fragment, branch, and counterpoise-style Hamiltonians using separated
     center terms.

4. **High-order geometry choices**
   - slab/endcap/panel and related high-order shellification/lowering choices;
   - q-ladder convergence and route-owned provenance;
   - no separate downstream Hamiltonian implementation.

5. **Consumer readiness**
   - stable `CartesianIDAHamiltonian` semantics;
   - versioned read/write contract;
   - enough basis/provenance information for a real downstream consumer when
     that consumer demonstrates the need;
   - Cr2-scale performance and memory evidence.

6. **Public use**
   - one recommended high-level producer;
   - a documented advanced staged workflow;
   - user-facing route support matrix and examples;
   - no requirement to understand private payloads, report fields, or historical
     route names.

## 3. Baseline at `02c97685`

### Working numerical spine

- Generic terminal-basis realization is live for one-center, H2, and Cr2
  terminal topologies.
- Blockwise final-basis one-body construction is live.
- H and H2 one-body values have reviewed parity.
- A separated N2 topology has exercised the one-body path and term-first
  Gaussian contraction.
- Localized terminal-basis IDA is live and reproduces the reviewed H2
  self-Coulomb value.
- The existing `CartesianIDAHamiltonian` constructor is the in-memory C2
  boundary.
- H2 materialization returns the Hamiltonian directly and can write/read the
  existing artifact format.
- The old physical-gausslet target/source/supplement payload chain has been
  deleted.

### Still incomplete

- Slice D implementation is live, but its design/`AGENTS.md` authority and the
  actual accepted signature still need formal reconciliation.
- The committed H2 smoke and fixtures still encode stale supplement-preflight,
  blocker, route-role, and source-plan-era vocabulary instead of the real
  Hamiltonian endpoint.
- The canonical user documentation does not yet present the new terminal PQS
  Hamiltonian producer as a supported public workflow.
- One-center atomic materialization is not yet fully normalized onto the same
  public base-Hamiltonian route as H2.
- Cr2 has a real terminal basis, but not yet a reviewed full Hamiltonian
  stress/performance result.
- White--Lindsey/Qiu--White materialization still has separate legacy/donor
  orchestration.
- Residual-GTO/MWG supplements, corrections, branch workflows, and high-order
  policies remain donor or experimental capabilities rather than unified
  terminal-producer features.
- Pair-term and assembly stages still carry historical structure; their future
  role is unresolved because the base terminal producer no longer needs them to
  construct K/U/V.
- Several developer status pages, migration tables, fixtures, and public
  ordinary-branch pages describe an older state.

## 4. Long-term goals and exit criteria

The manager-log goals remain valid, but their current interpretation is:

| Goal | Current status | Exit criterion |
|---|---|---|
| LT1 Reliable Hamiltonian construction | Partial: complete for the base H2 PQS endpoint | Atomic, representative molecular, supplemented, corrected, and supported high-order routes produce the same stable Hamiltonian contract |
| LT2 Reduce repo bloat and complexity | Active: major payload deletion completed | Cartesian/QW/PQS source plus tests are net smaller than the roadmap baseline, with no unclassified donor or compatibility authority |
| LT3 Public Cartesian driver line | Early | One recommended public producer, one documented advanced staged workflow, route support table, and runnable examples |
| LT4 High efficiency | Partial | H/H2 correctness, N2 topology/allocation gate, and Cr2 peak-memory/time gate pass without silent fallback or global support matrices |
| LT5 Scientific correctness and provenance | Strong for base PQS | Every supported route and artifact records real construction authority; oracle and diagnostic routes cannot masquerade as production |
| LT6 Stable consumer contracts | Partial | Versioned Hamiltonian artifact is consumed by at least one downstream workflow; any basis/provenance artifact exists only because a consumer needs it |
| LT7 Stratified testing | Partial | Kernel, public workflow, physics, artifact, slow stress, and benchmark tiers are distinct; internal-vocabulary smokes are gone |
| LT8 Extensible architecture | Partial | Geometry/lowering differences remain local; final-basis operators, supplements, corrections, Hamiltonian, and artifacts use common contracts |

## 5. Target architecture

### 5.1 Public facade

The public surface should have two levels:

1. A recommended high-level producer that accepts physical system, basis/route,
   and output options and returns `CartesianIDAHamiltonian`.
2. An advanced staged workflow for scientists who need to inspect geometry,
   terminal support, retained spaces, or the final basis.

Exact names require a design freeze. Conceptually:

```julia
ham = build_cartesian_hamiltonian(system, basis_or_spacing, route; options...)
write_cartesian_ida_hamiltonian(path, ham)
```

The public producer must not return route payloads or readiness wrappers. Expected
unsupported input should throw a clear argument/capability error before expensive
construction.

### 5.2 Internal dependency direction

```text
geometry policy
-> shellification / owned support
-> lowering / source CPBs
-> retained rules
-> final-basis realization
-> blockwise K/U/V kernels
-> optional basis augmentation
-> optional Hamiltonian corrections / branches
-> CartesianIDAHamiltonian
-> artifact / consumer
```

Rules:

- geometry code does not own operator assembly;
- retained/source plans do not carry final operators in metadata;
- final-basis operator code does not branch on atom count or route names;
- supplements extend a final basis rather than creating a parallel Hamiltonian
  family;
- corrections operate on stable matrices/bases rather than route payloads;
- artifacts contain consumer data, not the construction-stage object graph.

### 5.3 Pair-stage decision

The current `pair terms -> assembly` stages must earn their place.

Before broad public stabilization, make one explicit decision:

- **retain and make real** if supplements, high-order mixed blocks, or another
  active consumer needs typed pair planning/materialization; or
- **collapse/delete** from the base public workflow if they remain pass-through
  inventory while K/U/V are built directly from the terminal basis.

Do not preserve empty stages merely because an older architecture diagram listed
them.

## 6. Sequenced roadmap

### R0 — Truth alignment and branch closure

Purpose: make docs, authority, tests, and the live base endpoint tell the same
story.

Work:

- formally approve/record the implemented `HP-WIRE-02` Slice D boundary;
- add `HP-WIRE-02` to `AGENTS.md` or otherwise reconcile the authority record;
- ratify the accepted materialization signature and its narrowly allowed use of
  direct report fields, or move computation to parent/system/recipe inputs;
- replace the stale H2 supplement-preflight fixture and internal-vocabulary smoke
  with a base H2 Hamiltonian endpoint check;
- update current route migration/status pages that still claim missing terminal
  shell projection or missing source-plan materialization;
- record a measurable roadmap baseline: relevant source/test/docs line counts,
  exported Cartesian names, active route wrappers, test times, H2/N2 metrics,
  and Cr2 basis memory facts.

Done when:

- branch authority is internally consistent;
- committed H2 validation checks a real Hamiltonian and artifact semantics;
- no active fixture claims the base H2 PQS endpoint is unavailable;
- the design branch is ready for broad review/merge.

Deletion target:

- stale fixtures/blocker variables;
- obsolete smoke assertions;
- superseded D-candidate wording.

### R1 — Public base producer

Purpose: make the working A/B/C/D base lane a real package workflow.

Work:

- put one-center atomic and bond-aligned diatomic base Hamiltonians through the
  same final-basis/operator/Hamiltonian/materialization boundary;
- choose and freeze the recommended public high-level producer;
- document the advanced staged path without exposing private route vocabulary;
- publish H and H2 examples plus an artifact roundtrip example;
- document supported/experimental/reference route status;
- establish a light separated-diatomic public integration case;
- run Cr2 base construction as a later stress/performance gate.

Done when:

- a user can build H and H2 Hamiltonians from documented public input without a
  private harness;
- one-center and diatomic code share final-basis K/U/V and artifact machinery;
- the artifact is readable and useful outside the producing process;
- current user docs no longer describe the Cartesian producer only as hidden or
  route-internal research machinery.

Deletion target:

- migration-only atomic common-H1/final-basis adapters after parity;
- route-specific base-Hamiltonian materializers;
- stale private fixtures and report aliases;
- duplicate public overloads superseded by the chosen producer.

### R2 — White--Lindsey/Qiu--White and PQS final-basis unification

Purpose: keep algorithmic differences below the final-basis boundary and remove
parallel operator/Hamiltonian implementations.

Work:

- define the honest common final-basis consumption contract for WL/QW and PQS;
- route supported WL atomic and diatomic bases through common K/U/V,
  `CartesianIDAHamiltonian`, and artifact code;
- preserve WL/QW-specific lowering and retained rules where mathematically real;
- decide the pair-stage future based on actual WL/supplement consumers;
- classify every old fixed-source/fixed-block path as production, oracle,
  migration, or deletion.

Done when:

- route family changes basis construction, not the downstream Hamiltonian API;
- line ladders validate physical/public endpoints rather than route internals;
- no old fixed-source path remains active authority accidentally.

Deletion target:

- separate White--Lindsey materialization wrapper where no live consumer needs it;
- duplicate one-body/IDA/Hamiltonian builders;
- route-global safe-term pilots and pair scaffolding that have no chosen role;
- fixed-source oracles whose equivalence duty is complete.

### R3 — Generic residual-GTO / MWG augmentation

Purpose: recover supplement capability without reviving the deleted payload
chain.

Architecture:

```text
base terminal basis
+ explicit supplement specification
-> project supplement against base basis
-> residual rank selection and local cleanup
-> augmented final basis
-> common K/U/V and Hamiltonian construction
```

Work:

- define one typed supplement input and one final-basis augmentation contract;
- reuse existing analytic Gaussian, mixed-block, and MWG donor kernels;
- keep center count and species as data;
- support atomic and molecular cases through common algebra where real;
- add basis/provenance artifact data only after a downstream consumer states the
  exact need.

Done when:

- an H/He atomic and H2 molecular supplement case produce a real Hamiltonian;
- mixed and supplement/supplement blocks are finite, symmetric, and validated;
- residual ranks and weight conventions are explicit;
- the public artifact roundtrips and is consumable.

Deletion target:

- legacy residual-GTO/QW orchestration wrappers;
- supplement request/representation/preflight vocabulary;
- duplicate GA/AA and residual final-space mixers after shared kernels take over.

### R4 — Corrections, branches, and fragments

Purpose: move existing useful correction science onto the common Hamiltonian
boundary.

Work:

- expose hydrogenic-core/ESOI correction as a post-base-Hamiltonian operation;
- expose EGOI/stationary-Fock correction with explicit target/provenance limits;
- express full/fragment/counterpoise branches using separated unit `U_A` and
  shared basis/operator data;
- keep dense exact Gaussian Coulomb machinery reference-only and size-guarded.

Done when:

- corrections are basis- and route-neutral where their mathematics is neutral;
- branch results preserve center/charge provenance;
- public docs state approximation and scaling limits;
- correction paths do not manufacture new operator payload families.

Deletion target:

- ordinary-QW-only correction adapters after generic matrix-level entry points
  cover their real consumers;
- duplicate branch charge/orchestration code;
- stale correction result wrappers that only mirror matrices.

### R5 — High-order geometry and q-ladder integration

Purpose: make slab/endcap/panel and related high-order choices first-class
geometry/lowering policies, not separate Hamiltonian routes.

Work:

- move selected high-order donor algorithms into shellification/lowering-owned
  choices;
- run q-family convergence rather than treating q4 as an accuracy label;
- preserve support ownership, no-fallback PGDG behavior, capture, H1, and
  chemistry diagnostics;
- validate H2 chemistry and Cr2 occupied-space capture before broader claims;
- decide whether chain/square routes have a real consumer or remain explicitly
  experimental.

Done when:

- high-order choices feed the same final-basis and Hamiltonian producer;
- q-specific provenance and convergence are visible;
- no high-order-specific operator or artifact implementation exists.

Deletion target:

- experimental high-order orchestration duplicated by route-owned policies;
- one-off q4 wrappers and probes after durable q-ladder validation exists;
- unsupported geometry families with no scientific consumer, after explicit
  retirement review.

### R6 — Cr2 and downstream consumer readiness

Purpose: demonstrate that the architecture is useful at the scale that motivated
it.

Work:

- produce a base Cr2 Hamiltonian with measured wall time, allocations, and peak
  RSS;
- add required supplements/corrections only after their preceding milestones;
- validate basis capture, H1, IDA conventions, and consumer input semantics;
- run an external CR2/HFDMRG/DMRG smoke that reads the artifact and performs a
  meaningful operation;
- version the artifact if consumer-driven provenance additions are required.

Performance gates:

- no global support operator or normal-working global dense coefficient matrix;
- local simultaneous workspace at or below the reviewed 64 MiB cap unless a
  later design changes it;
- Cr2 peak RSS target at or below 1.2 GiB and hard cap 1.5 GiB for the reviewed
  base two-center fixture, unless a measured design amendment revises it;
- no silent numerical-reference fallback on a PGDG production route;
- compilation/specialization pressure is measured as well as runtime memory.

Done when:

- a downstream consumer can use a documented, versioned artifact;
- Cr2 performance matches the scale model;
- producer and consumer agree on charges, centers, electron counts, basis order,
  and localized IDA convention.

### R7 — Public stabilization and final retirement wave

Purpose: finish the transition from a research collection to a maintainable
public Cartesian subsystem.

Work:

- classify exported Cartesian/QW names as stable, experimental, deprecated, or
  internal;
- deprecate/remove duplicate front doors and diagnostic types from the public
  export list;
- publish a concise public route support matrix, tutorials, examples, and scale
  limits;
- archive stale developer notes while preserving the algorithm implementation
  index and feature-donor history;
- complete the route retirement ledger;
- reduce tests to durable scientific, public-contract, kernel, and performance
  gates.

Done when:

- there is one recommended producer and one stable Hamiltonian artifact;
- no active code depends on blocker/status payload vocabulary;
- every retained donor/oracle has a named live purpose and owner;
- Cartesian/QW/PQS source plus tests are net smaller than the R0 baseline;
- public docs and exported names match actual supported functionality.

## 7. Continuous repository-shrink program

Shrinking is not a final cleanup phase. It is a condition on every milestone.

### Rules

1. Every feature migration names the old authority it makes unnecessary.
2. A replacement switches live callers and deletes the replaced implementation
   in the same merge whenever possible.
3. Parallel old/new production paths may survive for at most one milestone and
   require a named external caller and deletion condition.
4. New numerical capability may be line-positive within its approved budget;
   the enclosing milestone should still be net-negative after donor retirement.
5. No compatibility wrapper is kept merely because it is small.
6. Variable-size basis/unit/pair/center inventories use vectors or indexed/lazy
   views, not specialization-heavy tuples or runtime-keyed named tuples.
7. Tests and docs count as carrying cost too.

### Baseline metrics to record at R0

- source lines and files under `src/cartesian*`, `src/pqs*`,
  `src/ordinary_qw*`, and related route-driver files;
- corresponding test/tool lines and files;
- number of exported Cartesian/QW types and functions;
- number of route-specific materialization/front-door wrappers;
- number of active compatibility/oracle paths;
- package load latency and representative first-call compilation latency;
- H2 smoke time, N2 one-body time/allocations, and Cr2 basis time/peak workspace.

### Final shrink criterion

The program is not complete if new unified code has merely been added beside all
old routes. At R7, relevant source plus tests must be net smaller than the R0
baseline, notwithstanding supplements, corrections, high-order support, and the
public producer.

## 8. Public API and artifact policy

### Public API

- Prefer one high-level producer over a family of route-specific `gethams`
  equivalents.
- Keep an advanced staged API only for scientifically meaningful inspection.
- Public inputs should be typed physical/basis/route specifications, not report
  objects or arbitrary named tuples.
- Public outputs should be basis objects, Hamiltonians, corrections, or stable
  diagnostics—not construction payload graphs.
- Unsupported combinations belong in a documented support table and clear
  exceptions.

### Artifact

- `CartesianIDAHamiltonian` remains the base in-memory boundary.
- The existing Hamiltonian writer/reader remains the default artifact contract.
- Add a basis/provenance artifact only for a demonstrated consumer requirement.
- Version format changes and provide migration/readback tests.
- Do not serialize internal stages, raw pair tensors, report graphs, or caches.

## 9. Validation matrix

| Tier | Purpose | Representative cases | Cadence |
|---|---|---|---|
| Kernel | Local algebra and shape contracts | Lowdin, projection, factor contraction, IDA normalization | Per relevant change |
| Base physics | Public numerical correctness | H, H2 | Per numerical milestone |
| Topology/performance | Separated mixed terminal topology | N2 or similarly light diatomic | Per operator/augmentation milestone |
| Route parity | Algorithm family comparison without conflating authority | supported WL/QW vs independent PQS | At route-unification gates |
| Artifact | Stable producer/consumer contract | H2 write/read and downstream operation | Per artifact change |
| Stress | Scale and memory | Cr2 | Milestone/nightly/manual, not every pass |
| High-order | Convergence and geometry policy | H2 q-ladder, Cr2 capture | High-order milestones |
| Supplements/corrections | Added science | atomic plus H2 supplemented/corrected fixtures | Feature milestones |

Durable tests should assert physical quantities, stable public types, artifact
semantics, symmetry/finiteness, and scale limits. Terminal-role tuples, blocker
names, helper names, and report-field presence belong only in temporary local
probes unless they are themselves a stable module contract.

## 10. Performance invariants

- Term-first Gaussian/Coulomb contraction is the default organization.
- Reuse parent 1D factors; do not rebuild them inside block-pair loops.
- Direct identity sectors remain implicit.
- Production does not form global support-space K/U/V or a normal-working global
  dense parent-to-final coefficient matrix.
- At most one bounded support-pair workspace is live.
- Numerical-reference and diagnostic dense paths are explicit and cannot be
  silent fallbacks.
- Public routes document expected scaling and representative measured behavior.
- Compilation time, method specialization, and invalidation pressure are part of
  the performance review, not only peak numerical memory.

## 11. Developer-document policy

The old developer documents are valuable algorithm and donor indexes. Preserve
that value without allowing old status text to override current architecture.

- `algorithm_implementation_index.md` is required navigation before numerical
  work.
- `cartesian/feature_donor_inventory.md` owns donor status and deletion
  conditions.
- `cartesian_route_retirement_ledger.md` owns retirement state.
- This roadmap owns project sequence and current strategic goals.
- Frozen design documents own implementation authority.
- The manager running log owns accepted-decision history.
- User-facing manual/reference pages own supported public behavior.

When a note becomes historical:

1. mark its status clearly;
2. link to the current authority;
3. move it under an archive index if it no longer guides active work;
4. keep source anchors and lessons discoverable through the algorithm index.

Do not create another disconnected long-range plan. Update this page when the
project direction changes.

## 12. Governance and update cadence

At every roadmap milestone:

- update the baseline/status table here;
- update the feature-donor and retirement ledgers;
- update public support documentation if user-visible capability changed;
- report source/test/docs line changes and retired surfaces;
- report scientific and performance evidence;
- state the next milestone and the explicit non-goals.

Every five accepted Cartesian passes, the manager should check whether work still
advances the current roadmap milestone. Every major milestone should end with a
broad architecture, numerical, performance, deletion, documentation, and public
API review.

A roadmap entry is not permission to code an unapproved surface. If execution
needs a new type, persistent object, public function, stage field, artifact
field, or status vocabulary, first amend and approve the relevant design.

## 13. Immediate next sequence

1. Complete R0 authority reconciliation for Slice D.
2. Replace the stale H2 supplement-preflight fixture/smoke with a physical base
   Hamiltonian endpoint test.
3. Refresh the route-migration and current-ordinary status pages.
4. Record the R0 line/export/performance baseline.
5. Perform a broad review and merge decision for the current producer branch.
6. Begin R1 with a docs-only public-producer/API freeze before further source
   expansion.
