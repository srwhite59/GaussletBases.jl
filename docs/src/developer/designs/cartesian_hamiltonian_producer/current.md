# Current Authority

Status: Slice A, Slice B, Slice C1, Slice C2, and Slice D base handoff are
implemented for the internal base PQS Hamiltonian producer. R1 public base
producer implementation is approved for the narrow H/H2 scope. R3-A
residual-GTO basis plus exact augmented one-body/moment implementation is
approved for the first H2 endpoint. R3-B in-memory residual MWG/IDA
Hamiltonian construction is reapproved for the first H2 endpoint with the
corrected weight-aware compact-path baseline.

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
`registry.md`, plus the narrow R3-A residual-GTO exact one-body/moment and
R3-B in-memory residual MWG/IDA Hamiltonian surface
recorded in `r3_residual_gto_mwg_augmentation.md` and `registry.md`. The
visible driver shape may call the implemented base path, but this design does
not approve a new artifact format except the `HP-R1-ART-01`
`producer_provenance/` keys in the final Hamiltonian file, solver integration,
broad driver redesign, public workflow outside the R1 H/H2 scope, or R3-C
artifact/cleanup work.

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
  `tau_rel = 1.0e-10`, `tau_neg_abs = 1.0e-12`, and
  `tau_neg_rel = 1.0e-12`;
- allowed numerical work: deterministic residual-basis construction plus exact
  augmented `K`, uncharged `U_A`, and moment matrices `x`/`y`/`z`/`x^2`/`y^2`/
  `z^2`;
- implementation organization: R3-A may use the QW analytic 1D-table donor
  pattern inside the approved owner file to build full parent-by-supplement
  `G-A` blocks once and project them through terminal blocks. This does not
  approve parent-by-parent global operators, a new shared QW API, persistent
  provider bundles, payloads, or edits outside the approved file;
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
- corrected weight-aware compact-path H2 closure value: lowest augmented
  one-body orbital IDA self-Coulomb `0.4574256036192161` within `1.0e-10`;
- superseded R3-B targets: `0.457435475059184`, from the retired private
  `[pre_final_pqs, residual_gto]` density gauge, and
  `0.4574331709135599`, from direct parent-density `G-M` insertion, are not
  acceptance targets for the compact R3-A residual basis.

R3-B does not approve artifact provenance, public API expansion,
driver/bin/tool workflow, broad provider payloads, status/result objects,
report fields, parent-stage fields, Be2 validation, Cr2 validation,
RHF/solver work, rank-loss implementation, width scaling, or tolerance
relaxation.

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
- R3 hardening before Be2/Cr2: independent weight-aware `V_GM` validation,
  same-construction consistency, deterministic rank-loss implementation,
  owned-support mixed providers, bounded MWG term storage, and nonallocating
  large-matrix validation checks;
- Cr2-scale stress and performance validation;
- R3-C artifact provenance/cleanup, Be2/Cr2 validation, and other non-base
  Hamiltonians;
- solver integration;
- White-Lindsey pair-framework completion;
- distorted-product COMX realization;
- EGOI or other Hamiltonian corrections.

Normal startup reading for this lane is this file, `registry.md`,
`invariants.md`, `implementation_slices.md`, and
`docs/src/developer/algorithm_implementation_index.md`.
