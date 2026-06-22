# Current Authority

Status: Slice A, Slice B, Slice C1, Slice C2, and Slice D base handoff are
implemented for the internal base PQS Hamiltonian producer. R1 public base
producer implementation is approved for the narrow H/H2 scope. R3-A, R3-B,
and R3-C are implemented for the narrow internal H2 residual-GTO/MWG endpoint:
deterministic residual-GTO basis construction, exact augmented one-body and
moment matrices, same-construction in-memory MWG/IDA Hamiltonian construction,
and compact supplemented artifact provenance in the existing Hamiltonian file.
The R3 owner-local residual-selection source correction is approved under the
existing R3 residual/Hamiltonian IDs. The R3 usability lane is approved only as
a non-exported supported facade for H2 and internal/performance-supported Be2
supplemented artifacts; it must use the corrected owner-local residual
selection and the updated H2 scalar below.

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
`registry.md`, plus the narrow implemented R3-A/R3-B/R3-C residual-GTO/MWG
surfaces recorded in `r3_residual_gto_mwg_augmentation.md` and `registry.md`,
and the R3 usability facade recorded in
`r3_usability_supplemented_workflow.md`. The visible driver shape may call the
implemented base path, but this design does not approve a new artifact format
except the `HP-R1-ART-01` `producer_provenance/` keys and the `HP-R3-ART-01`
`supplement_provenance/` keys in the final Hamiltonian file, solver
integration, broad driver redesign, or public workflow outside the R1 H/H2
scope.

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
  Parent gausslet rows are orthonormal to machine precision and terminal
  regions own disjoint parent rows, so block-local terminal basis supports are
  structurally orthogonal across blocks. Cross-block overlap is zero by
  construction, not a physical residual to compute or repair. A nonzero
  structural overlap means duplicated support rows, incorrect row restriction,
  wrong support ownership, or an indexing error.
- Cross-block kinetic, nuclear-attraction, and IDA interactions may still be
  nonzero and remain assembled over terminal block pairs. Structural
  cross-overlap zero does not imply block-diagonal operators.
- `cartesian_transforms` owns terminal basis realization for supported PQS
  terminal plans.
- `cartesian_materialization(report, terminal_basis_realization,
  materialization_inputs)` receives `transforms.terminal_basis_realization`
  directly.
- Requested base PQS materialization returns `CartesianIDAHamiltonian{Float64}`
  directly. No-request materialization returns `nothing`.
- Optional base-Hamiltonian artifact writing uses the existing
  `write_cartesian_ida_hamiltonian` shape.
- Once the Cartesian parent lattice and axis bundle are realized, that
  construction is the authority for reusable parent-only one-dimensional
  numerical data: overlap, kinetic, coordinate, second moment, integral
  weights, Gaussian factor terms, raw pair-factor terms, and exponent ordering.
  Supplement-dependent parent-by-supplement cross tables are derived
  construction-local work data, built from that parent-axis source for a
  validated supplement, expansion, and physical centers. They are not metadata,
  report fields, route-stage fields, artifacts, public API, or global mutable
  caches.

Approved R3-A residual-GTO exact one-body/moment scope:

- approved source owner/path: `CartesianFinalBasisRealization` owns
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` for
  `HP-R3-OBJ-01`, `HP-R3-FN-01`, and `HP-R3-FN-02`;
- first fixture: public/base z-axis H2 plus contracted two-center H/cc-pVTZ,
  `lmax = 1`, `uncontracted = false`, no width filtering;
- frozen residual thresholds: `tau_abs = 1.0e-10`,
  `tau_rel = 1.0e-10`, `tau_neg_abs = 1.0e-12`,
  `tau_neg_rel = 1.0e-12`, residual-occupation cutoff
  `eta_RG = 1.0e-8`, final-merge thresholds
  `tau_merge_abs = 1.0e-12` and `tau_merge_rel = 1.0e-12`, and final
  orthogonality tolerance `1.0e-10`;
- allowed numerical work: deterministic residual-basis construction plus exact
  augmented `K`, uncharged `U_A`, and moment matrices `x`/`y`/`z`/`x^2`/`y^2`/
  `z^2`;
- implementation organization: R3-A may use the QW analytic 1D-table donor
  pattern inside the approved owner file to build full parent-by-supplement
  `G-A` blocks once and project them through terminal blocks. This does not
  approve parent-by-parent global operators, a new shared QW API, persistent
  provider bundles, payloads, or edits outside the approved file;
- corrected residual-selection invariant: residual-GTO candidates are
  partitioned by physical owner center before residual-content selection.
  Owner-local residual Gram eigenvalues are residual occupations and control
  retention. If every owner-local mode is retained, preserve donor-style
  full-rank owner orientation by symmetric Lowdin in owner candidate order. If
  rank is lost, use deterministic owner-local natural modes. Retained
  owner-local sectors are concatenated and merged by one final symmetric
  Lowdin over inter-owner overlap. Global raw-candidate symmetric Lowdin and
  global raw-column pivoted-Cholesky selection are superseded and are not the
  R3 residual algorithm;
- first validation gate: H2 augmented one-body/moment endpoint only, checking
  `G' S R`, `R' S R`, base G-G block equality, finite/symmetric augmented
  operators and moments, and `E1_aug <= E1_base + epsilon`.

R3-A does not approve MWG/IDA `V`, supplemented
`CartesianIDAHamiltonian` construction, artifact provenance, public API
expansion, driver/bin/tool workflow, broad provider payloads, status/result
objects, report fields, Be2 first-gate validation, or Cr2 validation.

Approved R3-B residual-MWG/IDA in-memory Hamiltonian scope:

- approved source owner/path/function:
  `CartesianFinalBasisRealization` owns
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, function
  `pqs_terminal_residual_gto_augmented_hamiltonian`, for `HP-R3-FN-03`;
- `HP-R3-FN-03` is extended as the same-construction internal path for the H2
  endpoint. The approved function may accept the same-construction
  `base_hamiltonian`, `CartesianTerminalBasisRealization`, bundles,
  supplement, atom locations, and nuclear charges, then construct the residual
  augmentation object, exact augmented `K`/`U_A`/moments, residual MWG
  descriptors, weight-aware `V_GM`, direct `V_MM`, and the existing
  `CartesianIDAHamiltonian{Float64}` inside one call;
- residual MWG centers and widths are computed from exact R3-A moment matrices
  using `c_ralpha = <r | alpha | r>`,
  `v_ralpha = <r | alpha^2 | r> - c_ralpha^2`, and
  `sigma_ralpha = sqrt(2 * v_ralpha)`;
- `V_aug = [V_GG_base V_GM; V_GM' V_MM]`, with `V_GG_base` unchanged from the
  base Hamiltonian, `V_MM` using density-normalized MWG/IDA factors directly,
  and `V_GM` transformed from parent density normalization to final-basis
  density normalization block by block:
  `support_weights = wx .* wy .* wz`,
  `final_weights = C' * support_weights`,
  `C_density = C .* support_weights ./ final_weights'`, and
  `V_GM_block = C_density' * V_support_M`;
- the returned object is the existing in-memory
  `CartesianIDAHamiltonian{Float64}`;
- corrected owner-local H2 closure value: lowest augmented one-body orbital
  IDA self-Coulomb `0.4574265214362075` within `1.0e-10`;
- historical global-selection compact-path H2 closure value
  `0.4574256036192161` belongs to the superseded global candidate-order
  residual basis and is not a future target;
- superseded R3-B targets: `0.457435475059184`, from the retired private
  `[pre_final_pqs, residual_gto]` density gauge, and
  `0.4574331709135599`, from direct parent-density `G-M` insertion, are not
  future acceptance targets.

R3-B does not approve artifact provenance, public API expansion,
driver/bin/tool workflow, broad provider payloads, status/result objects,
report fields, parent-stage fields, Be2 validation, Cr2 validation,
RHF/solver work, rank-loss implementation, width scaling, or tolerance
relaxation.

Approved R3-C compact supplemented artifact provenance scope:

- approved source owner/path: `CartesianFinalBasisRealization` owns
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` for
  `HP-R3-ART-01`;
- R3-C writes the existing `CartesianIDAHamiltonian{Float64}` artifact shape
  with `write_cartesian_ida_hamiltonian`, then adds only the compact
  `supplement_provenance/` group defined in `registry.md`;
- `read_cartesian_ida_hamiltonian` is used for validation/readback only and
  remains a Hamiltonian reader, not a public provenance API;
- R3-C does not approve a Hamiltonian wrapper, payload/status/report object,
  public API/export, driver/bin/tool workflow, solver/RHF/Cr2, broad
  residual-basis serialization, or any artifact format beyond the existing
  Hamiltonian file plus provenance group.

R3 closeout status:

- R3-A/B/C provide the first narrow supplemented Hamiltonian path for the H2
  fixture. The owner-local residual-selection correction is now the approved
  R3 source target: augmented dimension `489` and lowest-orbital IDA
  self-Coulomb `0.4574265214362075` within `1.0e-10`.
- The Be2 R3-A donor-kernel measurement closed the first exact-operator
  scaling blocker. Repeated CPB-per-terminal-block construction measured about
  `43.2 s` and `35.4 GiB`; the one-shot parent-by-supplement analytic block
  organization measured about `1.94 s` and `2.1 GiB` with roundoff agreement
  for tested `G-A` and `A-A` blocks.
- Be2 R3-B residual rank `26` showed modest MWG/IDA storage and runtime at
  that proxy size. Bounded/streamed MWG storage is therefore not urgent before
  the next planning lane, but remains a high-rank/Cr2 guardrail.
- Cr2 remains deferred. It is a stress/consumer-readiness milestone, not the
  next correctness gate.
- Owner-local residual-selection measurement for H2, Be2, Cr2 q4, and Cr2 q5
  found final orthogonality below `1.0e-10` and no rank loss under trial
  residual-occupation cutoffs `1.0e-8` or `1.0e-7`. The approved cutoff is
  `eta_RG = 1.0e-8`, preserving the less aggressive retained content while the
  evidence shows no first-scope rank difference.

Approved owner-local residual-selection correction:

- For each owner center, compute `M_a = S_AaAa - X_a' X_a` after projecting
  that owner's candidates against the fixed orthonormal terminal gausslet
  basis. The eigenvalues of `M_a` are residual occupations, not just numerical
  rank diagnostics.
- Retain owner-local modes using residual-occupation cutoff
  `eta_RG = 1.0e-8`; numerical negative-eigenvalue tolerance and physical
  residual occupation cutoff are separate policies.
- If all owner-local modes are retained, preserve donor-style full-rank owner
  orientation with owner-local symmetric Lowdin in candidate order. If rank is
  lost, use deterministic owner-local natural residual modes with sign
  canonicalization.
- Concatenate retained owner-local sectors and perform one final symmetric
  Lowdin over the inter-owner merge overlap.
- Final merge failure rule: with
  `tau_merge = max(1.0e-12, 1.0e-12 * max(lambda_max(S_merge), 1.0))`,
  any merge eigenvalue below `-tau_merge` is a construction error and any
  eigenvalue `<= tau_merge` is a near-singular merge error. Do not floor merge
  eigenvalues to preserve directions. Final `G' S R` and `R' S R - I` errors
  must be below `1.0e-10`.
- MWG centers and widths are computed from the final merged residual
  functions. Residual integral weights may be near zero or sign-changing and
  must not be treated as base-PQS IDA weights.
- Do not implement a stabilized global raw-candidate Lowdin pass, do not use
  global raw-column pivoted-Cholesky as residual-content selection, do not use
  eigenvalue flooring to retain tiny residual modes, and do not use width
  filtering as conditioning repair.
- The active R3-B/R3U H2 scalar is `0.4574265214362075`. The old R3-B scalar
  `0.4574256036192161` is a global-selection baseline, not a target to
  preserve.

Approved R3 usability workflow scope:

- approved IDs: `HP-R3U-FILE-01`, `HP-R3U-FN-01`, `HP-R3U-WIRE-01`, and
  `HP-R3U-TEST-01`;
- approved facade: non-exported
  `cartesian_residual_gto_mwg_hamiltonian(system; basis, supplement,
  hamfile = nothing)::CartesianIDAHamiltonian{Float64}`;
- approved primary owner file: `src/cartesian_base_hamiltonian.jl`;
- existing R3 owner file
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` may be
  adjusted only to reuse the R3 same-construction path and R3-C writer without
  recomputing residual objects or adding new payload/status/artifact shapes;
- approved validation path:
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`, as a
  standalone gate only;
- first supported systems: z-axis H2 as the committed endpoint and z-axis Be2
  as an internal/performance-supported proxy. Cr2 remains unsupported;
- supplement input is a compact `NamedTuple` with `basis_by_center`, `lmax`,
  optional `uncontracted`, and optional `width_filtering`;
- optional artifact output writes the existing Hamiltonian artifact plus the
  approved `supplement_provenance/` group and returns the Hamiltonian.

R3 usability does not approve a public export, driver/bin/tool workflow, Cr2
artifact or full run, ECP/EGOI/RHF/solver work, a Hamiltonian wrapper,
new artifact format, broad report/status/payload object, or exposed internal
stage objects.

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
- public export or driver workflow for supplemented Hamiltonians beyond the
  approved non-exported R3 usability facade;
- full Cr2 residual-GTO/MWG Hamiltonian or artifact construction;
- Cr2-readiness lane: measurement-only candidate/rank/memory forecast, with no
  full Cr2 Hamiltonian yet;
- high-rank R3 performance guardrails: bounded/streamed residual MWG storage
  if residual rank grows, nonallocating large-matrix validation checks, and
  consumer-scale timing;
- basis/supplement-realism lane: validated supplement choices, basis labels,
  and filtering policy beyond the first H2 fixture;
- Cr2-scale stress and performance validation;
- R3 cleanup beyond compact artifact provenance and other non-base
  Hamiltonians;
- solver integration;
- White-Lindsey pair-framework completion;
- distorted-product COMX realization;
- EGOI or other Hamiltonian corrections.

Normal startup reading for this lane is this file, `registry.md`,
`invariants.md`, `implementation_slices.md`, and
`docs/src/developer/algorithm_implementation_index.md`.
