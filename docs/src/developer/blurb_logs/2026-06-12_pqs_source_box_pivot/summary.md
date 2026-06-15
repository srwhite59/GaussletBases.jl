# PQS Source-Box Pivot Baton Log

This log records the curated manager/doer blurbs, responses, and reviews for
the PQS source-box-first pivot.

The live polling files for this loop are not stored here. They live directly in
`.agent_handoffs/` under the root live baton convention described in
`docs/src/developer/root_live_baton_loop.md`.

## Current Theme

The pivot should develop PQS as a source-box-first operator route:

```text
raw product source boxes
-> source-space operator factors
-> retained source-mode rules/transforms
-> retained source-space operator blocks
-> optional shell realization only as oracle/adapter
```

Shell realization must not become the algorithm. It may compare, diagnose, or
adapt, but the route authority should remain the source-box retained transform.

## First Target

The first live target is a read-only PQS restart audit. It should identify the
smallest first implementation step for:

```text
one-center raw source box
one retained PQS boundary-mode rule
one self pair
safe one-body terms first: overlap and kinetic
```

No Coulomb/IDA/RHF work should start until the retained source-mode transform
and safe one-body block contraction are clear.

## Test And Artifact Policy

Use `tmp/work` probes for exploratory audits. Do not add tests by default.
Long-lived tests should be compact module-contract tests or physics/workflow
acceptance gates once the contract is stable.

Every implementation response must include deletion/shrinkage reporting.

## Consolidated Pass Summary: 001-100

Passes 1-100 were mainly the first PQS source-box restart and maturation phase.
The work began by re-establishing module ownership:
`CartesianRawProductSources` was identified as the owner of raw product
source-box facts, while older raw-product and CCPM/CCP fixture surfaces were
explicitly treated as oracle/reference or migration material, not new route
authority. The line then added the key source-box operator pieces: PQS/PQS safe
one-body source-space blocks, source-mode Gaussian-factor projection, centered
electron-nuclear block conventions, and final-basis/shell-realization bridges.
By pass 020, the new final-basis object matched the old shell-realization
oracle exactly on the cubic q=5 projected-shell fixture, with final overlap
identity error about `2.2e-14`, while still correctly blocking H1/IDA/RHF/driver
adoption.

The middle of the range shifted from source-box kernel plumbing toward physical
interpretation and complete/multilayer PQS. Boundary-shell-only H1 was shown to
be mathematically consistent with its shell-support oracle but not a physical H
acceptance fixture: the H1 value was only about `-0.0817`, far from `-0.5`,
because the basis was just the retained boundary shell. A complete core/shell
density/J attempt exposed a serious IDA-weight convention problem: signed
final-weight division gave an unphysical `J ~= 18.13`, so density-density/RHF
was correctly blocked rather than accepted.

The big constructive advance was the multilayer PQS shell source plan. By pass
060 it could build repeated one-cell projected-q-shell layers around a core,
cover a side-13 parent with 343 core rows plus 1854 shell rows, produce a
1549-dimensional final basis with overlap error about `5.5e-13`, and give
plausible H/He+ H1 values (`-0.4942`, `-1.9756`) in a route-owned planning
seam. Later passes around 70-80 clarified authority: support kinetic could be
promoted, electron-nuclear needed by-center convention review, and explicit-box
multilayer entry points were downgraded to compatibility/probe bridges while
shellification/lowering-backed region plans became the route-authority
direction.

The detailed per-pass files for passes 001-100 were removed after this
consolidation to keep the curated log from accumulating obsolete scaffolding.

## Consolidated Pass Summary: 101-150

Passes 101-150 pivoted from one-center private H1/J and RHF diagnostics into a
Hamiltonian-payload direction, then into Be2/PQS diatomic route construction.
The one-center line gained private RHF input, one-step and initial-density
checks, SCF/residual/Fock-DIIS controls, and finally a private object-carrying
Hamiltonian payload. The medium-term target then shifted toward comparable
Be2 WL/PQS payloads for downstream CR2 inspection, using route-owned readiness,
parent-axis, source-plan, support-window/order, raw-box, source-realization,
final-basis, and H1 payloads.

Durable decisions:

- Compact `7^3 / 5^3 / dim 223` one-center fixtures are route-smoke and
  convention diagnostics, not physics endpoints.
- RHF is private diagnostic/prototype only, distinct from serious HF.
- RHF route wiring was rejected in favor of private Hamiltonian-constructor
  payloads.
- One-center PQS Hamiltonian payload authority is object-carrying route data,
  not report aliases.
- Be2/PQS must not fake the old `:pqs_multilayer_shell_source_plan`; it gets
  its own `:pqs_diatomic_complete_core_shell_source_plan`.
- After private Be2/PQS H1, H1-J is not the production-facing seam; Ham-input
  and electron-electron payloads should come first.

Validation landmarks:

- Pass 101 one-center H1/J dry run: final dimension 223, H1
  `-5.6629907690725245`, self-Coulomb `1.8691288063594704`, cold about 120s
  and warm in-process about 0.21s.
- Passes 105-109 added focused synthetic RHF tests for input contract, one-step,
  initial density, and SCF.
- Passes 113-128 used local ignored compact RHF probes. Fixed-point and simple
  damping did not converge; Fock-DIIS with history 8 converged at iteration 34,
  energy `-10.032119189804888`, commutator residual about `2.95e-9`.
- Pass 127 committed default private Fock-DIIS history 8.
- Pass 131 committed the one-center private PQS Hamiltonian payload; focused
  test passed 21/21.
- Passes 135-149 advanced Be2/PQS readiness/final-basis/H1 fingerprints:
  probe-enabled final retained count reached 221 and private H1 lowest energy
  was `-0.27746109235228694`.

Cleanup and guardrails:

- Pass 142 removed direct/ad hoc raw-box producer assertions from focused tests
  in favor of route-owned `diatomic_raw_box_route_payload`.
- Most passes in this window built new private seams rather than deleting large
  surfaces, but report aliases, H1/J scalar diagnostics, RHF validators, and
  raw-box probes were kept out of route authority.
- Do not promote compact route-smoke fixtures to physics endpoints.
- Do not infer RHF electron count from nuclear metadata; keep it explicit.
- Do not route-wire RHF or treat it as serious HF.
- Do not use report aliases as Hamiltonian-construction authority.
- Keep public API, exports/artifacts, CR2, HFDMRG, WL promotion, RHF, and
  IDA/MWG promotion out of these private seams unless explicitly assigned.

Raw `blurb`, `response`, and `review` files for passes 101-150 were removed
after this consolidation.

## Consolidated Pass Summary: 151-200

Passes 151-200 completed much of the Be2/PQS private handoff line, then pivoted
to atom-first He driver-owned physics and finally to H2 target classification.
Be2/PQS gained private Ham-input, Hamiltonian handoff, CR2 read-only
inspection, and JLD2/TSV artifact support, but CR2 solver/export/HFDMRG
readiness stayed false. The Be2/WL/PQS artifact was correctly labeled
read-only and not same-basis comparison-ready because PQS and WL dimensions
differed. After pass 178, work pivoted away from Be2/CR2 toward atom-first He.

Durable decisions:

- CR2 decides downstream handoff shape; GaussletBases should not become the
  CR2/HF runner.
- The correct first atom target is fixed-q multi-shell He, not compact 223 and
  not q-ladder growth.
- He q=5/n_s=5 PQS is aligned to
  `white_lindsey_atomic_mapping(Z=2,d=0.3,tail_spacing=10.0)`.
- He route authority moved from hand-built fixture construction to visible
  driver input, `bin/cartesian_ham_builder.jl`, saved JLD2 artifact, and thin
  readback endpoint.
- Private RHF became driver-owned, default-off, and explicit-endpoint only.
- Old H2 WL/QW HF totals include H/cc-pVTZ S/P residual supplements, so
  gausslet-only PQS must not compare directly to those totals.
- The H2 221 route is diagnostic-only, not a physics endpoint.
- The physical H2 target needs a new route kind with atom-contact core,
  retained atom-core interiors, and shared shell layers, not mutation of the
  221 boundary diagnostic.

Validation landmarks:

- Be2 Ham payload and handoff tests repeatedly passed in about 40-50s.
- The CR2 artifact generator produced JLD2/TSV, shrank from about 242 MB to
  about 82 MB, and recorded clean producer metadata after regeneration.
- He fixed-q inventory was `125 + 3*98 = 419`.
- He WL-mapped PQS fingerprints: H1 lowest `-1.991334820314074`, H1-J
  self-Coulomb `1.2420423900074902`, PQS/WL H1 delta about `+9.65e-6`, and
  self-Coulomb delta about `-5.0e-6`.
- He private RHF endpoint: PQS total `-2.850817886618113`, WL 419
  gausslet-only total `-2.85080350301779`, delta about `-1.44e-5`.
- H2 221 diagnostic: parent axis counts `(9,9,15)`, final dimension 221, H1
  lowest `+0.14582426982296057`, treated as a warning.
- H2 physical audit recovered old gausslet-only fixed block `(1215,463)` and
  supplemented final dimension 481.

Cleanup and guardrails:

- Major deletions removed old RHF/SCF, seed/WL, Be2 artifact, compact 223 Ham,
  direct He fixture, materializer/timing probes, CPB overlap-placement
  facts/pilot tests, and many mixed one-body scaffold tests.
- Multi-minute He endpoints were removed from the default runner and retained
  only as explicit/manual endpoints.
- Do not compare gausslet-only PQS H2 to supplemented WL/QW HF/ED totals.
- Do not continue the H2 221 diagnostic route into H1-J/RHF as if physical.
- Keep private RHF default-off and artifact-scalar-only; do not write matrices,
  densities, orbitals, Fock matrices, or histories.
- Do not save final-basis self-overlap as downstream working data; use scalar
  identity diagnostics only.
- Prefer compact summary objects over scalar alias clouds.
- New physical H2 work should add a route kind and compact inventory payload
  first; no H1-J/RHF/supplement/comparison in that first pass.

Raw `blurb`, `response`, and `review` files for passes 151-200 were removed
after this consolidation.

## Consolidated Pass Summary: 201-229

Passes 201-229 cover the H2 463 source-backed route leadup, the fake-PQS
correction, and the first quarantine work. H2 R=4 q5 physical target inventory
was initially defined as 463-dimensional, no supplement, support counts
`(275, 578, 362)`, retained counts `(251, 98, 114)`, and retained order
`(:atom_contact_core, :shared_shell_1, :shared_shell_2)`. The old H2 221 route
was downgraded to source-box diagnostic, not a physics endpoint. Source-backed
WL/QW fixed-source data was allowed only as a checked private adapter/candidate,
not independent PQS authority.

Fake-PQS leadup and correction:

- Passes 210-211 promoted a source-backed fixed-source candidate into
  `:private_source_backed_adapter_authority`.
- Passes 212-217 built a working H2 463 route through final basis, H1, H1-J,
  RHF input, and RHF execution.
- Passes 218-221 compared it to no-supplement WL/QW values and made it appear
  endpoint-ready.
- Pass 229 corrected the interpretation: the H2 463 "PQS" artifact was a
  WL/QW fixed-source retained-transform reproduction.
- The route was renamed fake-PQS and marked `fake_pqs/enabled = true`,
  `independent_pqs_transform = false`, `delete_after_independent_pqs = true`,
  and `physics/endpoint_ready = false`.
- Its blocker became
  `:fake_pqs_source_backed_wl_reproduction_not_independent_pqs`.

Authority facts:

- Independent PQS authority must come from PQS source-box retained transform
  construction.
- Source-backed WL/QW fixed-source rows/coefs may be a regression oracle or
  private adapter, but not scientific PQS route authority.
- Stable artifact groups must carry fake markers when fake/source-backed paths
  are reported; hiding fake status only under a side group is unsafe.
- Future independent work must start from a separate target with
  `fake_pqs/enabled = false`.

Validation landmarks:

- Pass 212 final basis: dimension 463, final overlap identity error about
  `1.6e-13`.
- Pass 213 H1 lowest: about `-0.7946609179724462`.
- Pass 214 H1-J/self-Coulomb: about `0.45696639804337114`.
- Pass 217 private RHF: electronic energy about `-1.1589518556683855`,
  converged in 8 iterations.
- Pass 220 WL no-supplement reference: electronic RHF
  `-1.1589518556651142`, total with nuclear repulsion
  `-0.9089518556651142`.
- Pass 221 fake route vs WL deltas were tiny, later explained by source-backed
  WL reproduction rather than independent PQS agreement.

Cleanup and guardrails:

- Major cleanup included deleting private global-overlap driver hooks, splitting
  `CartesianContractedParentMetrics.jl`, splitting/deleting mixed
  `source_box_route_shadow.jl`, extracting low-order materialization, shrinking
  stale route-shadow assertions, retiring the H2 221 diagnostic input/test,
  retiring component-smoke sidecars, and deleting obsolete probe/smoke tests.
- Do not use the fake-PQS artifact as a PQS/WL scientific comparison.
- Do not infer private-RHF electron count from nuclear charges; require explicit
  route input.
- Keep RHF and supplement paths private/diagnostic until route authority is
  real.
- Preserve the H2 fake-PQS artifact only as a golden driver regression until a
  real independent PQS target replaces it.
- Supplement provider-block work remained blocked on concrete H2 support/source
  ordering, support-to-retained transform placement, mixed gausslet-GTO
  accumulation, GTO/GTO rules, and raw moment matrices.

Raw `blurb`, `response`, and `review` files for passes 201-229 were removed
after this consolidation.
