# Current Cartesian Hamiltonian Producer Authority

This page owns live implementation status, active lanes, blockers, and next
work. It is not the algorithm manual and does not independently grant source
authority. It remains a transitional long-form ledger pending the approved
topology compression; consult only the task-relevant status after startup.

Normal startup reading:

- `README.md` for orientation;
- `current.md` for live status;
- the assigned ID entry in `registry.md` for permission, lifecycle, ownership,
  and approved surfaces;
- `invariants.md` for architecture-wide guardrails;
- the one subsystem contract linked by that registry entry.

Before numerical implementation, also consult
`docs/src/developer/algorithm_implementation_index.md`. Endpoint/probe work has
the additional operational and due-diligence reading required by `AGENTS.md`.
The longer link catalog in `README.md` is task-specific navigation, not startup
reading.

Historical design and review material remains under `history/`, `reviews/`,
completed retirement records, and planning/implementation ledgers. Those files
are evidence and rationale, not normal startup authority when they conflict
with the compact current files.

## Live Status

The internal base PQS Hamiltonian lane is implemented for origin-centered H and
Cartesian z-axis H2. `HP-R1-ATOM-*` approves a narrow source relaxation for
explicit origin-centered all-electron one-center atoms beyond H through the
same base facade. This is internal numerical producer authority, not broad
public API polish.

Implemented base path:

- Slice A terminal basis realization over owned terminal supports;
- Slice B final-basis one-body assembly;
- Slice C localized IDA matrix assembly and existing
  `CartesianIDAHamiltonian` construction;
- direct staged base Hamiltonian construction and optional artifact writing;
- historical Slice D route-driver materialization wrapper retired;
- R1 public base facade for the approved H/H2 scope and fixed
  `producer_provenance/` artifact group.
- R1 one-center base atom relaxation for explicit origin-centered all-electron
  atoms through the existing base facade; current committed endpoint remains H.
- `HP-COMP-NS-*` approves public size-parameter normalization: `ns` is the
  durable user-facing cube/source/nesting size, while route-local `q` is
  derived as `q = ns` for `nesting = :pqs` or `q = ns - 2` for
  `nesting = :wl`.
- `HP-COMP-WLNS-*` approves the narrow WL z-axis diatomic `ns` contract
  cleanup: normalized `nesting = :wl`, `Natom = 2`, `ns < 4` should reject
  early in `src/cartesian_base_hamiltonian.jl`, and WL diatomic retained
  support may saturate across `ns` ranges when physical parent extent
  dominates.
- `HP-PQS-MAP-SFACTOR-FN-01` and `HP-PQS-MAP-SFACTOR-TEST-01` approve a
  narrow expert mapping-strength scalar. Optional positive `s_factor` defaults
  to `1.0`; omitted/default behavior must remain unchanged. For one-center
  White-Lindsey/atom mapping,
  `effective_s = s_factor * sqrt(Z * core_spacing)` while `core_spacing`
  remains the physical near-core scale. Multicenter PQS mapping may apply the
  analogous per-center factor only if the combined-invsqrt semantics are
  unambiguous and provenance records the actual mapping. This does not revive
  public `d`, `parent_mapping_d`, element defaults, automatic tuning, solver
  workflow, EGOI, rho0/P0, or protected-localized interaction changes.
- `HP-PQS-COULOMB-ACCURACY-FN-01` and
  `HP-PQS-COULOMB-ACCURACY-TEST-01` approve an expert
  `coulomb_accuracy = :compact | :high` producer choice with default
  `:compact`. The producer resolves one existing
  `CoulombGaussianExpansion` before parent/PGDG construction and carries it
  through base unit-nuclear/IDA, residual-GTO exact Coulomb-expanded blocks,
  and MWG. New artifacts record one Hamiltonian-wide expansion summary;
  protected/ladder readback exposes it, while missing legacy provenance is
  never inferred as high accuracy. Atomic reference packets remain a separate
  role-qualified exception: packet RHF is high accuracy, while current
  density/self-energy and fitted-potential scaffold evaluations are compact.
  A narrow amendment also stabilizes cancellation-prone analytic Gaussian
  factor formulas in `GaussianAnalyticIntegrals.jl` and the neutral nuclear
  raw-block owner. High-exponent failure is not authority for scaled/log PGDG
  carriers or terminal-contraction redesign.
  The canonical driver may now expose the same expert symbol, default it to
  `:compact`, validate `:compact | :high`, forward it in `common_basis`, and
  print it. This does not approve other CLI controls, custom parameters,
  shellification, residual/injection/EGOI/screened-Hartree policy, solver
  workflow, or Cr2-specific behavior.
- `HP-COMP-NSCORE-*` approves direct nucleus-centered core side parity:
  route-local `q` derivation remains unchanged, but direct core side comes
  from public `ns` as `isodd(ns) ? ns : ns + 1`. This oddization rule is only
  for direct core identity blocks; boundary retained sizes keep route-local
  construction.
- `HP-COMP-SHELLGEOM-*` approves only common terminal shell decomposition as a
  shared route-family-free first step. PQS and White-Lindsey may diverge in
  retained-construction geometry only after common direct core and shell-owned
  support regions exist.
- `HP-COMP-SHELLGEOM-DIAT-*` extends that rule to z-axis diatomics: PQS and
  White-Lindsey must enter the common shellifier with the same first-step
  arguments, including public `ns`, before central-gap/contact and shared-shell
  regions are planned.
- `HP-COMP-THINSLAB-*` supersedes the outer-mismatch-only lane. It approves a
  narrow unified z-axis diatomic thin-slab stack lowering repair. For both
  PQS and White-Lindsey, `:direct_midpoint_slab` and `:outer_mismatch_slab`
  must call the same compact slab lowering function with the same region,
  public `ns`, slab normal/thickness, and native source/support facts. The unit
  slice scale is `ns x ns x 1`; a thickness-`t` outer-mismatch region with
  `t <= ns` should scale about `t * ns * ns`, not become one identity block.
  Direct/core sectors remain identity, and real shell regions still diverge
  after common shellification. If slab thickness exceeds `ns`, source work
  must stop and request a policy decision rather than silently deleting slabs
  or using identity blocks.
- `HP-COMP-THINSLAB-META-*` approves only live metadata/scaffold inventory
  cleanup in `cartesian_terminal_shellification_geometry.jl`. The route helper
  still consumes `_cartesian_terminal_shellification_region_unit_inventory(...)`,
  so `_cartesian_terminal_region_unit_mapping(region)` must describe midpoint,
  outer-mismatch, and angular z-extension slabs as compact thin-slab lowering
  regions, not direct identity CPBs. This does not approve coefficient
  materialization, shellification geometry changes, route skeleton redesign, or
  a new reporting framework.
- `HP-COMP-FACEPROD-*` approves a neutral internal
  `CartesianFinalBasisRealization` face-product helper seam in
  `terminal_face_product_blocks.jl`. The helper owns compact face-like
  terminal blocks over one fixed normal-axis index or an ordered stack of
  fixed normal indices, using `_nested_doside_1d(...)` and
  `_nested_face_product(...)`. WL facet realization and future
  `HP-COMP-THINSLAB-*` slab realization should reuse this helper; it must not
  live in a WL- or PQS-specific terminal file.
- `HP-COMP-ANGBOX-*` approves the shellification side of the angular-balanced
  z-axis diatomic repair in `terminal_geometry.jl`. Each shared-shell step
  computes a physical outer-nucleus angular-balanced target box. The ordinary
  index-layer shell body plus planned z-extension thin-slab stacks, not the
  ordinary body alone, realizes that target coverage. When ordinary
  shared-shell growth stops with transverse axes saturated and bond-axis parent
  support remaining, shellification emits the bond-axis leftovers as planned
  `:angular_z_extension_slab` stack regions. Lowering those slabs remains
  deferred to `HP-COMP-THINSLAB-*`.
- `HP-PQS-SCREEN-HARTREE-CORR-FN-01` and
  `HP-PQS-SCREEN-HARTREE-CORR-TEST-01` implement the internal in-memory
  correction described by
  [the durable residual-density formalism](screened_hartree_residual_density.md)
  and [the assembly contract](screened_hartree_correction_assembly.md).
  `Vnuc_G` remains Galerkin; represented occupied determinants define
  `P0/q0`; the density fit defines the compressed cloud and `E0_G`; and the
  potential fit is only a fast `J0_G` evaluator. The H/Be/Be2, Ne
  fitted-cloud, and potential-fit audit IDs are completed historical evidence;
  the original direct Ne lane is historical and operationally superseded by
  the packet path. The facility remains internal and in memory. Its current
  boundary is downstream workflow and endpoint validation: no public default,
  corrected artifact, solver integration, exchange, interaction transform,
  `C' V C`, or Cr2 production claim is implemented.
- `HP-PQS-ATOMREF-PACKET-FN-01` and
  `HP-PQS-ATOMREF-PACKET-TEST-01` implement reusable one-center atomic HF
  reference packets under
  [the canonical packet contract](atomic_hf_reference_packets.md). Current
  bounded references are Be core `2e` and Ne all-electron `10e`, cc-pV5Z,
  `lmax = 1`. Unconverged packets are rejected; fitted density/potential terms
  are not protected orbitals; and density-fit `J0_G` receives the explicit
  role-qualified compact Coulomb expansion. Packet overlap fingerprints remain
  exact self-integrity checks, while translated owner-local molecular overlap
  blocks use exact identity/order mapping plus the unchanged numerical
  `norm(..., Inf) <= 1e-10` equivalence gate.
- `HP-PQS-ATOMREF-POTMOM-FN-01` and
  `HP-PQS-ATOMREF-POTMOM-TEST-01` implement the fixed-policy determinant-
  moment polish for those packet potential fits. It preserves all `33`
  exponents and the first `5` broad-tail coefficients, adjusts only
  coefficients `6:33` against the fixed `13`-distance Coulomb-moment grid, and
  keeps old packets readable without inferring polish. It does not weaken the
  screened-Hartree anchor or add a public/tunable fitting policy.
- `HP-REP-XGTO-IMPORT-FN-01` and `HP-REP-XGTO-IMPORT-TEST-01` approve only a
  representation-transfer facility for importing explicit external Gaussian AO
  orbitals into an orthonormal GaussletBases final working basis. The rule is
  `S_FG = <F|G_external>` and `C_F = S_FG * C_G`. External AO self-overlap
  `S_GG` is validation-only for checks such as `C_G' * S_GG * C_G ~= I`; it is
  not a generalized final-basis transfer metric. The lane may add a compact
  explicit packet reader/writer, import result, and capture diagnostics around
  existing `gto_overlap_matrix(...)`. It does not approve Hamiltonian
  transforms, `C' V C`, `Vee`/source transforms, solver workflow,
  screened-Hartree/EGOI changes, residual/injection policy changes, PySCF
  dependency in repo tests, or Cr2 production claims.
- `HP-MCOMX-*` approves a protected-`P2` plus mapped Chebyshev source-span
  option at the existing nested doside / COMX seam. The nonlinear map uses
  normalized local `u`, while `_cleanup_comx_transform(...)` still uses the
  physical position matrix. High-order may benchmark the installed option as a
  consumer, but the high-order branch does not own the mainline implementation
  shape.
- `HP-MCOMX-TERM-*` approves only terminal-basis wiring in
  `pqs_terminal_basis_realization.jl`: `_shell_seed(...)` may consume carried
  materialized `AxisSourceTransformFact`s as the basis-defining shell seed when
  they are present and validated, while ordinary
  `_nested_projected_q_shell_full_sides(...)` fallback remains unchanged.
- `HP-MCOMX-DRV-*` approves only a compact `source_span` construction choice in
  `bin/cartesian_ham_builder.jl` and the staged base/facade path. Public values
  are `:ordinary` and `:mapped_comx`; ordinary remains the default, and
  `:mapped_comx` is currently PQS-only.
- `HP-DRV-INV-*` implements a compact terminal-region / shellification inventory
  summary in the canonical driver output. This is bounded human-facing output:
  region label, kind, lowering or realization kind, support rows, final
  columns, compression ratio, shell index, index ranges for `x`/`y`/`z`,
  physical coordinate ranges for `x`/`y`/`z`, identity-vs-compact/product
  realization, and slab axis/side/thickness/stack facts when applicable. The
  physical `x`/`y` ranges are required because angular-balance failures compare
  transverse physical scale against the bond-axis margin. This does not
  approve new driver inputs, route diagnostics, source-mode or pair dumps,
  artifact schema changes, numerical construction changes, or Cr2 workflow.
- `HP-DRV-SHELLDD-*` implements a standard terminal due-diligence report for
  Cartesian/PQS terminal bases. The report is the live replacement for the old
  "driver printed enough for a human to inspect" practice: consumers are
  expected to inspect the derived system, parent axes/box, gausslet/IDA weight
  statistics, dimensions, and shell-by-shell geometry before interpreting
  energies, residual occupations, or injection behavior. It must expose
  normalized atom/charge/electron/geometry facts, validated atom locations,
  bond axis/length, snapped nuclear indices, parent physical bounds, axis
  counts, 1D center summaries or bounded tables, spacing summaries,
  gausslet/IDA weight statistics, dimension/compression accounting, terminal
  order/key, role, region kind, shell index, owner/contact/shared
  classification, index and physical boxes, physical side lengths and aspect
  ratios, actual and expected aspect-balanced source-mode shapes, source-mode
  count, retained count, final column range, lowering/retained/realization
  rules, slab metadata, and advisory warning flags. The implemented report
  extends
  `_cartesian_terminal_inventory_rows(...)` in
  `src/cartesian_base_hamiltonian.jl` and joining existing terminal inventory
  rows with terminal retained-rule plan/support records; the canonical driver
  prints the resulting bounded report. This does not approve
  artifact schema changes, public semantics changes, shellification policy
  changes, aspect-balanced source-mode implementation, dense coefficient/
  transform/support dumps, or Cr2 workflow.
- `HP-PQS-ASPECTSHELL-*` approves the separate future source-policy lane for
  PQS complete-shell source modes. Current complete shells hard-code cubic
  `(q,q,q)` source modes in `_pqs_complete_shell_contract(...)`; the restored
  policy should choose explicit aspect-aware `(q,q,L)` dimensions for z-axis
  diatomic complete shells, with `L` derived from the older angular-resolution
  rule or a documented validated equivalent. The approved seam is in the
  terminal low-order route-driver path after shellification has produced
  complete-shell regions and parent/bundle facts, but before lowering,
  retained-unit, transform-contract, and terminal retained-rule records are
  frozen. `region_contracts.jl` is too early to own the `L` choice by itself,
  and `pqs_multilayer_shell_source_plan.jl` is too late to be the only fix.
  This lane may change retained counts, final dimensions, Hamiltonian matrices,
  and energies. It does not approve artifact schema changes, public driver
  inputs, WL policy changes, residual/MWG/IDA changes, old route-global
  materialization, or Cr2 production claims.
- 2026-06-26 He/PQS evidence found `n_s = 5` mapped-COMX not robust enough for
  all-electron scalar capture. This keeps mapped-COMX opt-in only and blocks
  default promotion until bounded He `n_s = 6`/`7` H1/IDA evidence is reviewed.

Implemented Residual Gaussian path:

- Residual Gaussian basis construction, exact augmented one-body/moment
  operators, and residual MWG/IDA interaction live in
  `src/cartesian_residual_gaussians/`;
- compatibility entry points and artifact/facade hooks remain in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- the H2 owner-local residual-GTO/MWG endpoint has augmented dimension `489`
  and lowest-orbital IDA self-Coulomb `0.4574265214362075` within `1.0e-10`;
- Be2 remains an internal/performance-supported usability proxy, not a
  committed public gate;
- `HP-R3U-ZDI-*` relaxes the supported supplemented molecule scope to explicit
  homonuclear z-axis diatomics with no element-specific defaults or
  Cr2-specific branch;
- Cr2 may be used only as an explicit generic homonuclear z-axis ignored/user
  stress or usability run after H2/Be2 validation, not as a committed gate or
  special workflow.

Approved Residual Gaussian robustness lane:

- `HP-RG-ORTHO-FN-01` approves only final residual normalization/validation
  robustness in `src/cartesian_residual_gaussians/residual_basis.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` allowed
  only for narrow internal keyword plumbing if needed;
- the target is strict N2 q5 p10 at `core_spacing = 0.042857`, where
  `G' S R` passes, owner metrics are full rank, retained counts are `9,9`, the
  merge metric is healthy, and only `R' S R - I` slightly exceeds the old
  absolute tolerance;
- the lane may use a symmetric final-overlap validation and a combined
  absolute/relative final identity check, but it must not change residual
  selection semantics, occupation cutoff, merge failure rule, MWG/IDA, raw
  blocks, artifacts, driver workflow, or public API.
- `HP-RG-IDTOL-FN-01` approves only changing the default final residual
  `R' S R` identity validation tolerance to `1.0e-8` in the same RG owner,
  with optional narrow compatibility keyword plumbing in
  `pqs_terminal_residual_gto.jl`. This older Be tolerance policy is superseded
  for production defaults by `HP-RG-CUTOFF-FN-01` and then
  `HP-RG-CUTOFF-FN-02`.
- The Be atom cc-pV5Z `lmax = 1` evidence for `HP-RG-IDTOL-*` is a tiny final
  identity overshoot: `21` retained residual directions, minimum occupation
  `6.151e-6`, final merge condition `1.0`, `max |G' S R| = 1.776e-14`, and
  `max |R' S R - I| = 2.183e-10` against an old allowed error of about
  `2.000e-10`.
- `HP-RG-CUTOFF-FN-01` set the prior production defaults:
  `residual_occupation_cutoff = 5.0e-8` and `identity_atol = 5.0e-8`. That
  pass addressed the Cr atom `basis_ns = 9`, `map_ns = 11`, `lmax = 1`
  marginal direction at occupation `3.637e-8`.
- `HP-RG-CUTOFF-FN-02` supersedes the residual occupation default for
  production: `residual_occupation_cutoff = 1.0e-6`, while `identity_atol`
  remains `5.0e-8`. The Cr2 residual audit found low-H1 modes built from
  marginal directions at occupations about `1.27e-7` to `8.98e-7`; the first
  production move is to discard those directions by default before considering
  kinetic/`H1_RR` spectral guards.
- `HP-RG-CUTOFF-FN-02` does not change owner-local grouping, merge checks,
  `G' S R` validation, width/zeta filtering, MWG/IDA, artifacts, driver
  workflow, public API, or source files outside the approved RG owner/plumbing
  surface.
- `HP-RG-SPECTRAL-AUDIT-01` is approved as measurement-only follow-up after
  the cutoff cleanup. The cited Cr2 residual-only replay now has owner counts
  `62 + 62`, but still shows a low two-owner residual mode:
  `min eig(K_RR) = 0.3700413519`,
  `min eig(H1_RR) = -7.1647854052`, with owner weights about `0.5 / 0.5`.
  The audit may classify low residual-sector modes by owner weights and
  residual-occupation composition using ignored probes and durable text/TSV
  output. It does not approve pruning, kinetic/`H1_RR` guards, cutoff or
  tolerance changes, full HF, Vee/solver work, artifacts, source
  instrumentation, or driver changes.
- `HP-RG-INJECT-AUDIT-01` approves the first measurement-only lane for the
  optional injection-plus-RG scheme recorded in
  `residual_gaussian_injection_hybrid.md`. Ignored probes may classify
  owner-local principal modes, merge the trial injected subspace, report
  `B = G' S Y_inj` rank/condition, report true RG counts and `K_RR`/`H1_RR`
  spectra after injection, and measure injected-sector one-body projection
  errors for `K`, each unit `U_A`, and `H1`. This does not approve
  injected-sector source behavior, defaults, artifacts, driver input, public
  API, automatic pruning, MWG/IDA convention changes, full HF, or Cr2 workflow.
- `HP-RG-INJECT-FN-01` approves a narrow default-off in-memory implementation
  of the injection-plus-RG hybrid. The audit did not remove the current Cr2
  low-H1 residual sector, but injection remains the better general
  construction because RG alone has a singular-complement limit. The source
  lane may implement owner-local principal-mode classification, global
  injected-subspace merge, replacement base sector `F`, exact one-body
  transformation into `[F, R]`, and true-residual MWG/IDA only for directions
  not injected. It must preserve current behavior when
  `residual_injection_cutoff <= 0` and must not change driver input, public
  API, artifact schema, production defaults, MWG channels for injected
  directions, spectral pruning, full HF, or Cr2 workflow.
- `HP-RG-OCC-FIRST-INJECT-FN-01` and
  `HP-RG-OCC-FIRST-INJECT-TEST-01` implement the internal
  [occupied-first injection geometry](occupied_first_injection.md). It makes
  identified `Y_occ` mandatory, reports pre-inclusion base capture separately
  from roundoff post-inclusion recovery, validates physical complement/capture
  spectra, and rejects weak optional directions without creating MWG residual
  channels. The helper is tested but not wired into the protected-localized
  builder; composition over `M = [G, R_compact]` belongs to
  `HP-RG-PROTECT-ADDREF-*`.
- `HP-RG-PROTECT-ADDREF-FN-01` and
  `HP-RG-PROTECT-ADDREF-TEST-01` implement that first real internal consumer for
  homonuclear protected-localized molecules. It is explicit internal opt-in;
  protected members without packet references remain unchanged. The staged
  geometry consumes the already-built compact residual basis
  `M = [G,R_compact]`, makes the
  full-rank union of all placed packet occupied spaces mandatory before current
  optional staged selection, and preserves the original packet blocks
  separately for additive `P0`. Native-order packet blocks define
  `P0 = sum_a P_a`; placed fitted potentials define `J0 = sum_a J_a`; packet
  density fits define
  `E0 = sum_a E_aa + 2*sum_{a<b} E_ab`. Existing raw `GG/GA/AA`, protected
  fixed-sector transforms, localized `W`, inherited-site `Vee_L`, and
  `ScreenedHartreeCorrection` algebra must be reused. The first acceptance
  gate is a physically padded Be2 construction with two Be `2e` packets,
  separate packet recovery/trace checks, total charge `4`, explicit cross
  energy, anchors, unchanged unscreened matrices, and terminal due diligence.
  Stored packet fingerprints and structural mapping remain hard failures;
  mapped overlap hash equality is diagnostic when the numerical equivalence
  gate passes.
  No public input, corrected artifact, solver, endpoint claim, protected atom,
  counterpoise, compact/high transfer, EGOI, exchange, or Cr2 production
  workflow is approved. The strict padded Be2 gate passed; CR2 may consume the
  same in-memory path for measurement only.
- `HP-RG-PROTECT-INJECT-DESIGN-01`,
  `HP-RG-PROTECT-INJECT-FN-01` / `TEST-01`, and
  `HP-RG-PROTECT-ONEBODY-FN-01` / `TEST-01` govern the internal,
  default-off [protected-localized basis convention](protected_localized_basis.md)
  implemented in `residual_basis.jl` and `augmented_operators.jl`. The
  completed one-body audit is historical evidence. The completed
  `HP-RG-PROTECT-VEE-AUDIT-01`
  rejects direct `C' V C` interaction rotation and fixes localized `L`, exact
  `H1_L`, and inherited pre-injection site-order `Vee_M` as the viable
  convention. Public/default workflow, artifacts, EGOI, ladder behavior,
  screening, solver, and production claims remain separately governed.
- `HP-RG-PROTECT-ART-FN-01` / `HP-RG-PROTECT-ART-TEST-01` and
  `HP-RG-PROTECT-ARTLOC-FN-01` / `HP-RG-PROTECT-ARTLOC-TEST-01` are
  implemented. The opt-in artifact identity, native-order sector and matrix
  law, row-locality metadata, compatibility, and rejection behavior are
  canonical in
  [protected-localized artifact contract](protected_localized_artifact.md).
  Ordinary artifact semantics, public/default workflow, solver behavior,
  selection policy, screening, and Cr2 production claims remain excluded.
- `HP-RG-PROTECT-EGOI-AUDIT-01` is a completed historical measurement,
  summarized in manager running-log Passes 302-308.
- `HP-RG-PROTECT-EGOI-FN-01` and `HP-RG-PROTECT-EGOI-TEST-01` are approved
  pending internal source and validation authority. The retained-original
  `s1+s2`, local-product, `M2`, exactly local `DeltaV` contract is canonical
  in [retained-GTO local-product EGOI](retained_gto_egoi.md). No protected
  retained-GTO helper or focused test is implemented in committed source;
  current uncommitted `hamiltonian_corrections.jl` work is not authoritative.
- `HP-RG-PROTECT-LADDER-XFER-AUDIT-01` is completed historical
  same-parent transfer evidence.
- `HP-RG-PROTECT-LADDER-BUNDLE-FN-01` and
  `HP-RG-PROTECT-LADDER-BUNDLE-TEST-01` implement the internal opt-in
  [protected-localized ladder bundle](protected_localized_ladder.md). It owns
  versioned manifests, protected member references, exact adjacent
  final-basis cross overlaps, optional native-order restarts, target-member
  fixed-density evaluation, and bounded summaries. Generalized overlap,
  source-Hamiltonian/`Vee` transforms, interaction rotation, new
  representation sidecars, solver workflow, and Cr2 production remain
  excluded.
- `HP-RG-RHO0-GAL-AUDIT-01` approves only a measurement audit for rho0/Galerkin
  IDA correction over the protected-localized inherited-site baseline. It may
  use ignored probes, `/Users/srw/dmrgtmp` outputs, analytic IDA/Coulomb
  sanity checks, small H/He/H2 checks, Cr2 fixed-density diagnostics, and one
  bounded Cr2 HF replay only if static rho0/Galerkin diagnostics are sane. It
  does not approve source edits, public wiring, artifacts/provenance,
  production Hamiltonian workflow, `C' V C` revival, rejected broad directions
  as MWG residual channels, Vee scaling as the fix, screened-reference
  production claims, Cr2 production energy claims, or publication-scale
  validation sweeps.
  Subsequent row-gauge audits showed this line was algebraically
  under-specified: `(J*w)/w`, `diag(J)`, and the IDA/MWG density-proxy
  potential are different objects. Treat this ID as historical measurement
  evidence, not the final correction formulation.
- `HP-RHO0-REFDENS-AUDIT-01` is the current rho0 successor lane. It approves
  only a measurement audit for a fixed reference density matrix `P0` in the
  protected-localized final basis. The correction target is
  `Delta_F0_sigma = F_exact0_sigma[P0] - F_app0_sigma[P0]` plus
  `C0 = E_exact0[P0] - E_app0[P0] - sum_sigma Tr(P0_sigma * Delta_F0_sigma)`,
  so the corrected model matches both reference energy and first derivative
  at `P0`. First audit is Hartree-only. It must report `P0`
  construction/normalization, representability in the protected-localized
  basis, exact and approximate finite-difference Fock checks, corrected
  energy anchoring, and sector safety diagnostics. It does not approve source
  edits, public driver/API/export changes, artifact/provenance changes,
  production Hamiltonian or solver workflow, direct `C' V C`, residual/MWG
  default changes, basis-fate policy changes, broad rejected directions as MWG
  residuals, Cr2 production claims, or publication-scale sweeps.
  `HP-RHO0-REFDENS-FN-01` and `HP-RHO0-REFDENS-ERI-01` are candidate future
  IDs only and are not approved.
- `HP-RHO0-REFDENS-MIXH-AUDIT-01` approves only the measurement-first exact
  mixed Hartree seam audit needed before Be/Be2 protected-localized `P0`
  work. It may use ignored probes and `/Users/srw/dmrgtmp` outputs to find or
  prototype `(final final | reference reference)` for H, Be, and Be2 only. It
  must identify existing kernels/helpers if available, run Be/Be2 fixed-`P0`
  only if the exact mixed Hartree seam is feasible in ignored code, or report
  the smallest future source owner if not. It does not approve source edits,
  public workflow, artifacts/provenance, solver workflow, Cr/Cr2 diagnostics,
  HF exchange, direct `C' V C`, row action, `diag(J)`, `q0`, center metadata,
  IDA proxy shortcuts, residual/MWG default changes, basis-fate changes, broad
  rejected directions as MWG residuals, or committed tests.
- `HP-RHO0-MIXH-GG-FN-01` approves the first source-backed exact mixed Hartree
  seam only for terminal/base `GG` blocks from one-center atomic `P_A`.
  Authority is limited to a neutral `CartesianGaussianRawBlocks` helper file,
  internal module include wiring, and narrow Gaussian Coulomb oracle/pair-term
  reuse. The pass may build same-center reference pair-density terms, separable
  Coulomb one-body packets, and exact `GG` output plus compact diagnostics. It
  does not approve `GA`/`AA`, protected-localized transforms, `F_app[P0]`,
  `Delta_F0`, `C0`, artifacts, public workflow, Cr/Cr2, HF exchange, dense
  final ERIs, or residual/MWG/basis-fate changes. `HP-RHO0-MIXH-GG-TEST-01`
  requires H/Be/Be2-scale `GG` validation, dense-oracle spot checks, and an
  angular/off-diagonal same-center pair check.
- `HP-RHO0-MIXH-GAAA-FN-01` approves extending the same neutral mixed Hartree
  owner from `GG` to exact `GA = <G|v_P_A|A>` and
  `AA = <A|v_P_A|A>` blocks. It may reuse the source-backed `P_A` validation,
  same-center pair-density term stream, and Coulomb-expanded factor packets,
  but still only for one-center atomic reference contributions. It does not
  approve protected/final transforms, `F_app[P0]`, `Delta_F0`, `C0`,
  artifacts, public workflow, Cr/Cr2, exchange, dense final ERIs, cross-atom
  reference density products, or residual/MWG/basis-fate changes.
  `HP-RHO0-MIXH-GAAA-TEST-01` requires bounded H/Be `GA`/`AA` validation and
  dense-oracle spot checks including angular reference pairs and angular
  supplement rows.
- `HP-RHO0-MIXH-FEXACT-FN-01` approves only the exact-side transform from raw
  mixed Hartree `GG`/`GA`/`AA` blocks into the current final/protected-localized
  fixed sector using existing protected one-body transform helpers. It may
  return in-memory dense `F_exact_Hartree[P0]` plus compact diagnostics. It
  does not approve new raw kernels, geometry selection changes, `F_app[P0]`,
  `Delta_F0`, `C0`, artifacts, public workflow, Cr/Cr2, exchange, IDA/MWG
  interaction transforms, approximate Fock construction, or residual/MWG/
  basis-fate changes. `HP-RHO0-MIXH-FEXACT-TEST-01` requires H/Be/Be2-only
  final/protected transform validation and dense-oracle spot checks.
- `HP-RHO0-FAPP-AUDIT-01` approves only a measurement audit for the
  approximate-side fixed-`P0` Fock seam. It must identify `F_app[P0]` as the
  derivative of the actual current IDA/MWG approximate energy convention for
  represented `P0_final`, and finite-difference validate that derivative on
  H/Be/Be2. It may use ignored probes and `/Users/srw/dmrgtmp` output only. It
  does not approve source edits, `Delta_F0`, `C0`, corrected Hamiltonian
  assembly, artifacts/public workflow, Cr/Cr2, residual/MWG default changes, or
  row-action/`diag(J)`/`q0`/center-metadata/direct-`C' V C` shortcuts.
- `HP-RHO0-FAPP-FN-01` approves only a narrow source seam in
  `src/cartesian_ida_hamiltonian.jl` for paired approximate energy/Fock
  helper(s) on `CartesianIDAHamiltonian`: fixed represented `P_alpha` and
  `P_beta` in the Hamiltonian basis, the actual two-index IDA/MWG
  density-density convention, and alpha/beta Fock matrices
  finite-difference checked against the matching energy helper. It does not
  approve public API/export/default changes, artifacts, `Delta_F0`, `C0`,
  corrected Hamiltonian assembly, Cr/Cr2, exact exchange extensions, or
  row-action/`diag(J)`/`q0`/center-metadata/direct-`C' V C` substitutes.
  `HP-RHO0-FAPP-TEST-01` requires package load plus compact alpha/beta
  finite-difference validation, with H/Be/Be2-only ignored endpoint replay if
  the helper is consumed by the rho0 audit.
- `HP-RHO0-ANCHOR-FN-01` is implemented, superseded, and closed to new source
  work. It subtracted the full approximate interaction Fock, including the
  current same-spin exchange-like term, from `F_exact_Hartree[P0]`; its
  `Delta_F0` must not be used as the Hartree reference-density correction.
- `HP-RHO0-CORR-AUDIT-01` is historical measurement evidence. Its original
  full-interaction-anchor run is invalid. `HP-RHO0-JANCHOR-*` is implemented,
  and the later direct-Hartree H/Be/Be2 behavior remains a stop-signal pending
  `HP-RHO0-XPAIR-AUDIT-01`.
- `HP-RHO0-JANCHOR-FN-01` / `HP-RHO0-JANCHOR-TEST-01` approve the direct-
  Hartree replacement. `src/cartesian_ida_hamiltonian.jl` may add private
  direct-only helpers with `q = diag(P_alpha) + diag(P_beta)`,
  `E_app_direct = 1/2 * q' * V * q`, and
  `F_app_direct = Diagonal(V * q)`. `src/cartesian_residual_gaussians/
  augmented_operators.jl` may replace or supplement the anchor helper to form
  `Delta_J0 = F_exact_Hartree[P0] - F_app_direct[P0]` and
  `C0_J = E_exact_Hartree[P0] - E_app_direct[P0] -
  Tr((P0_alpha + P0_beta) * Delta_J0)`. The corrected full interaction keeps
  the existing approximate exchange-like contribution:
  `E_corr_int = E_app_full_int + Tr((P_alpha + P_beta) * Delta_J0) + C0_J`.
  Validation must prove the direct anchor, shared spin-independent `Delta_J0`,
  the full corrected finite-difference derivative
  `F_exact_Hartree[P0] - K_app_sigma[P0]`, and rerun H/Be/Be2 corrected
  audit behavior with `Delta_J0`/`C0_J`. It does not approve artifacts,
  public workflow, production Hamiltonian integration, solver workflow,
  Cr/Cr2, exact exchange correction, or changes to the current approximate
  exchange convention.
- Implementation evidence now shows supplement-space atomic `P0` is
  represented at roundoff in H/Be/Be2 and the direct-Hartree anchor algebra is
  viable. The corrected Hamiltonian with inherited approximate exchange-like
  terms remains blocked: H/Be/Be2 replay still shows sizable negative low-mode
  shifts. `HP-RHO0-XPAIR-AUDIT-01` is the next approved rho0 lane. It is
  measurement-only, H/Be/Be2 only, and should compare direct-Hartree
  correction, inherited approximate exchange-like contributions, exact/
  supplement-space exchange diagnostics where feasible, and an explicitly
  non-production direct-only corrected-operator diagnostic. It does not
  approve source edits, artifacts, public workflow, solver integration,
  Cr/Cr2, exact exchange implementation, changes to the current approximate
  exchange convention, or direct-only physics as a Hamiltonian.
- `rho0_reference_density_implementation_plan.md` is a review memo for the
  likely fast separable atomic-reference Hartree source shape. The current
  approved source target is raw exact mixed Hartree `GG` plus `GA`/`AA` from
  one-center atomic `P_A`, followed by the exact-side final/protected transform;
  the paired Cartesian IDA `F_app[P0]` seam is source-backed. The old
  full-interaction anchor is superseded by the approved direct-Hartree anchor
  replacement in `HP-RHO0-JANCHOR-*`; corrected-Hamiltonian behavior must be
  rerun with `Delta_J0`/`C0_J`, and its current H/Be/Be2 result is a
  stop-signal pending exchange/direct pairing diagnosis. Cr/Cr2, artifacts,
  public workflow, solver integration, exact exchange correction, and
  production Hamiltonian integration remain later lanes.

Completed stale complete-core-shell RHF retirement:

- `28e9b2c84` removed the include and deleted
  `src/pqs_multilayer_complete_core_shell_rhf.jl`;
- `HP-RETIRE-CCS-RHF-*` is now historical deletion evidence, not active
  source/test authority;
- do not restore the payload stack or add a replacement/compatibility layer
  without a new docs-only amendment.

Completed route-driver materialization/report/save retirement:

- `e2e164e9b` removed the old materialization/report/save wrappers and their
  stale tool/test pressure;
- `77fa2700b` removed the two dangling ladder runners;
- `HP-RETIRE-DRV-MAT-*` and `HP-RETIRE-LADDER-RUNNERS-*` are historical
  deletion evidence, not active source/test/tool/docs authority;
- the canonical staged driver and current producer/artifact path remain the
  only live workflow. Do not restore wrappers, runners, or adapters without a
  new docs-only amendment.

Approved neutral Cartesian Gaussian raw-block owner:

- `src/cartesian_gaussian_raw_blocks/` is approved for exact uncharged
  by-center Cartesian Gaussian nuclear `G-A` and `A-A` raw blocks;
- `HP-CGRB-FN-02` approves only a neutral-kernel reorganization around unique
  one-dimensional supplement axis families and term-first table reuse;
- `HP-CGRB-NN-*` approves only non-nuclear overlap, kinetic, coordinate
  moment, and second-moment `G-A`/`A-A` raw blocks;
- `src/cartesian_gaussian_axis_integrals.jl` remains optional support only for
  `HP-CGAI-FN-01` helper work needed by that family-reuse kernel;
- the owner may be consumed by Residual Gaussian and Qiu-White code after
  behavior-preserving parity;
- the already-crossed cleanup applies to the main diatomic Qiu-White
  non-nuclear path. Remaining QW-local non-nuclear helpers used by atomic QW
  reference paths, factor-term outputs, hybrid sidecars, dense-parent probes,
  and CPB/provider surfaces are retained reference/sidecar/provider surfaces,
  not dead duplicates under `HP-CGRB-NN-*`. Rewiring them requires later
  authority for factor-block ownership or for those specific sidecar/provider
  callers;
- it does not own terminal projection, residual Gaussian transforms,
  Qiu-White route objects, final-basis `G-G` product-matrix optimization,
  caches, reports, artifacts, or public API.

Approved R3/RG terminal `G-G` product-matrix lane:

- `HP-R3GG-FN-01` approves only `K_GG`, coordinate moment `G-G`, and
  second-moment `G-G` product matrices used by
  `pqs_terminal_residual_gto_augmented_operators(...)`;
- owner module is `CartesianFinalBasisRealization`;
- approved files are
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` and,
  only if needed for small internal terminal-product workspace/helper reuse,
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`;
- it does not approve `G-A`/`A-A` raw-block changes, unit-nuclear Gaussian-sum
  work, IDA/MWG changes, route setup, public API, artifacts, metadata/status
  fields, persistent caches, or Cr2 workflow.

Approved remaining exact-operator allocation decision:

- after `954c86cd`, terminal `G-G` product workspace is considered crossed for
  the current Cr2 q4 proxy;
- `HP-R3REM-AUDIT-01` attributed the largest in-wrapper remaining allocation to
  unit-nuclear `U_GG` factor lookup plus Gaussian-sum construction;
- `HP-R3UN-FN-01` now approves only that narrow terminal final-basis
  unit-nuclear `U_GG` Gaussian-sum allocation lane under
  `CartesianFinalBasisRealization`;
- `HP-R3BASE-FN-01` approves same-construction reuse of already-built base
  final-basis `K_GG` and unit `U_GG[A]` blocks in supplemented exact augmented
  operators, with dimension and center-count validation and no metadata or
  provenance payload;
- `HP-R3BASE-DRV-WIRE-01` approves only the canonical driver supplemented-mode
  call-site wiring that passes `base_ham.kinetic` as `base_kinetic` and
  `base_ham.nuclear_attraction_unit_by_center` as `base_unit_nuclear`;
- the driver wiring must not change public inputs, hooks, timing labels,
  visible stage sequence, artifact schema, or driver contract;
- route/stage setup, raw-block setup, neutral raw-block kernels,
  residual/MWG/IDA changes, public workflow, and Cr2 facade/artifact work
  remain unapproved.

Approved compact Hamiltonian artifact manifest lane:

- `HP-HAM-MANIFEST-FN-01` approves only sidecar JLD2 groups for existing
  `CartesianIDAHamiltonian{Float64}` artifact files;
- approved source files are `src/cartesian_base_hamiltonian.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, and
  `src/cartesian_ida_hamiltonian.jl` only for a small unexported sidecar writer
  helper if needed;
- `hamiltonian_manifest/` follows the earlier PQS fixed-column/source-mode
  provenance model: basis identity is a status-bearing construction label, not
  a representative center;
- `hamiltonian_manifest/final_basis_labels/` records one row per matrix-order
  final basis column with sector, unit/source labels, shell/ray/radial status,
  representative center metadata, owner nucleus index where meaningful,
  locality/freezing labels, and supplement angular labels where available;
- optional `final_basis_source_relations/`, `source_shells/`, and
  `source_modes/` subgroups may be written only for native construction facts;
- `recipe_provenance/` records the validated public construction recipe and
  base/residual/augmented dimensions, including the public `nesting`
  construction family and a truthful route label derived from `(input.kind,
  input.nesting)`;
- `HP-NEST-ART-FN-01` approves only the narrow source cleanup that records
  `nesting` in base and recipe provenance and writes `:one_center_wl_base`
  rather than a PQS-oriented label for WL one-center artifacts. Its original
  supplemented-WL early-rejection boundary is historical and is superseded for
  supported z-axis diatomic supplemented WL construction by `HP-COMP-SUPPWL-*`;
- existing matrix keys and `read_cartesian_ida_hamiltonian` behavior must not
  change;
- no public reader API, driver public input change, route report/status
  payload, dense transform, raw inventory, inferred shell/ray labels, solver
  field, Cr2-specific field, or committed Cr2 fixture is approved.
- `HP-HAM-MANIFEST-SRC-FN-01` separately approves a compact construction-native
  source-mode provenance seam from terminal lowering / retained-unit /
  raw-product source plans to the base working basis manifest context, only for
  optional `source_shells/`, `source_modes/`, native source relations, and
  native final-label improvements.
- the source-mode seam must not serialize coefficients, dense transforms,
  `T_G`, `T_A`, raw inventories, route reports, allocation probes, diagnostic
  payloads, or repo-chosen ray/radial labels.

Approved canonical driver usability lane:

- `HP-DRV-FILE-01` approves only `bin/cartesian_ham_builder.jl`;
- `HP-DRV-FN-01` approves a compact functional driver workflow with visible
  defaults, one optional trusted project input file, command-line `key=value`
  overrides, visible public `system`/`basis`/optional `supplement` contract
  construction, coarse user-facing timing/summary printing, artifact write,
  optional readback check, and only the compact hooks `check_file`,
  `print_contract`, `print_timing`, and `expected_dimension`;
- `HP-DRV-STAGE-FN-01` approves a narrow non-exported, non-underscored staged
  producer surface so the driver can execute visible physics-level stages and
  operator-class timings: base working basis / terminal realization, product
  / moment operators, unit-nuclear attraction, electron-electron / IDA or
  residual-MWG interaction, base Hamiltonian assembly, Gaussian supplement,
  residual augmentation, and supplemented Hamiltonian assembly; this must be
  separate named stage functions, not one all-in-one replacement wrapper;
- approved source files for `HP-DRV-STAGE-FN-01` are
  `src/cartesian_base_hamiltonian.jl`,
  `src/pqs_source_box_low_order_materialization.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, only
  for behavior-preserving operator-class stage factoring;
- `HP-DRV-STAGE-WIRE-01` allows the canonical driver to call that staged
  surface as visible top-level stage calls and time/print those workflow stages
  without calling underscored package helpers;
- `HP-DRV-NEST-FN-01` and `HP-DRV-NEST-WIRE-01` approve one visible
  construction-family input, `nesting = :pqs` or `nesting = :wl`, in the
  canonical driver and narrow base-facade plumbing; `:pqs` maps to the existing
  PQS source-box route and remains the default, while `:wl` maps only to the
  existing White-Lindsey low-order route;
- `nesting` is not a diagnostic route switch and must not expose route
  skeletons, retained rules, raw-block switches, stop-after controls,
  diagnostics, route reports, or route-stage labels;
- supplemented `nesting = :wl` is now governed by `HP-COMP-SUPPWL-*` for the
  supported homonuclear z-axis diatomic composition cell; unsupported geometry
  or supplement combinations must still reject clearly;
- artifact provenance must record the public `nesting` input truthfully under
  `HP-NEST-ART-FN-01`; a WL artifact must not be labeled as PQS merely because
  helper code remains PQS-named;
- the target producer contract is the three-choice composition recorded in
  `nesting_supplement_composition_plan.md`: geometry (`atom` or z-axis
  diatomic), nesting (`:pqs` or `:wl`), and supplement state (`off` or `on`);
  this is planning authority except for explicitly promoted cells, and
  unsupported cells must keep clear rejection until separately approved;
- the driver may call only approved base, staged, and supported supplemented
  producer surfaces and the approved artifact writer/readback;
- `basisname = nothing` selects base mode; `basisname !== nothing` selects a
  supported supplemented mode; `HP-COMP-SUPPATOM-FN-01` separately approves
  relaxing the old supplemented `Natom == 2` guard so `Natom == 1` and
  `Natom == 2` can use the same supplemented staged path;
- `padding` is physical box padding around the atom or two nuclei and maps
  internally to the existing public facade fields;
- route diagnostics, stop-after internals, ladder probes, raw-block switches,
  underscored package helper calls, status/report/payload fields, allocation
  probes, artifact schema dumps, solver work, public API/export changes,
  artifact schema changes, old route-stage choreography, and Cr2-specific
  workflow remain unapproved in the canonical driver.

Approved White-Lindsey terminal-basis seam:

- `HP-WLTERM-FN-01` approves only terminal-basis realization for the existing
  `:white_lindsey_low_order` route family, returning the same
  `CartesianTerminalBasisRealization` consumed by the staged Hamiltonian path;
- `HP-WLTERM-WIRE-01` approves only the route-helper seam that lets a WL route
  with sufficient native terminal lowering/retained records produce that
  terminal basis instead of being discarded by the PQS-only guard;
- approved files are `src/pqs_source_box_route_driver_helpers.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  optional
  `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`,
  and the include in
  `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`;
- this is not approval to adapt the old WL H1/H1+J materialization path,
  change route skeleton semantics, change shellification/retained-selection
  behavior, add diagnostics, change artifacts, or add supplemented WL behavior.

Approved first composition lane:

- `HP-COMP-WLDIAT-FN-01` approves the native WL z-axis diatomic base path for
  `Natom = 2`, `nesting = :wl`, `basisname = nothing`;
- the implementation must produce native WL diatomic terminal records and then
  use the same `CartesianTerminalBasisRealization`, staged base Hamiltonian
  construction, writer, and reader path as the PQS producer;
- `HP-WLDIAT-COMPACT-FN-01` separately approves the compact retained-basis
  correction for WL z-axis diatomics: the current elongated shared-shell
  boundary-stratum identity realization is classified as a mechanical endpoint,
  not the intended production compact WL retained basis;
- each WL diatomic unit must carry or realize compact retained columns from
  products of one-dimensional contractions, not full-support identity retention;
- identity realization remains valid only for true direct/core identity units,
  not for WL boundary-stratum retained units;
- deleted WL coefficient helpers are historical donor/reference material only:
  source work may re-express the essential CPB-local product-of-1D coefficient
  primitive behind the current terminal-basis boundary, but must not revive the
  route-global WL stack, reports, adapters, or H1/H1+J materialization;
- `HP-WLDIAT-PARITY-FN-01` separately approves the boundary-stratum retained
  count parity cleanup: odd-side enforcement belongs to direct
  nucleus-centered core blocks, not to boundary shell strata; public WL
  `ns = 4`
  therefore targets `4^3 - 2^3 = 56` boundary columns rather than the inherited
  symmetric-odd donor result `26`;
- approved source files are the narrow diatomic complete-core-shell,
  terminal-shellification, terminal-lowering, route-helper, terminal-basis,
  final-basis include, and base-facade surfaces named in `registry.md`;
- `:z_axis_diatomic_wl_base` is approved only as a truthful route-provenance
  value under existing artifact keys, not as an artifact schema change;
- no driver special cases, fake compactness by dropping rows, old WL H1/H1+J
  materialization, RG/MWG/supplement work, route diagnostics, public
  API/export changes, committed tests, or Cr2 workflow is approved.
- `HP-COMP-BASEDIAT-FN-01` separately approves only
  `src/cartesian_base_hamiltonian.jl` to relax base two-center validation from
  H2-only to explicit homonuclear z-axis all-electron diatomics;
- the base diatomic relaxation requires equal symbols, equal finite positive
  integer-valued charges, two finite distinct z-axis centers, and neutral
  `nup + ndn == sum(charges)`;
- it preserves both `nesting = :pqs` and `nesting = :wl`, but it does not
  approve route/shellification/terminal-lowering changes, supplement/RG/MWG
  work, driver changes, element tables, heteronuclear support, translated or
  non-z-axis geometry, committed tests, or Cr2 workflow.
- `HP-COMP-SUPPWL-FN-01` approves the supplemented White-Lindsey z-axis
  diatomic composition lane for `Natom = 2`, `basisname !== nothing`, and
  `nesting = :wl`;
- approved source surface is `src/cartesian_base_hamiltonian.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  optional only for a direct RG/MWG compatibility genericity blocker;
- the implementation may remove the two early supplemented-WL blockers only if
  the existing Residual Gaussian/MWG path works with the WL
  `CartesianTerminalBasisRealization`;
- it must preserve the supplement contract, residual selection, exact
  augmented operators, MWG/IDA convention, base K/U reuse, artifact keys,
  manifest/provenance, driver inputs, and stage labels;
- it does not approve driver changes, supplemented atoms, route
  skeleton/shellification/terminal lowering changes, raw-block changes,
  residual-selection changes, MWG/IDA convention changes, artifact schema or
  reader changes, public API/export changes, old WL H1/H1+J materialization,
  committed tests, or Cr2 workflow.
- `HP-COMP-SUPPATOM-FN-01` approves the supplemented one-center atom
  composition lane for `Natom = 1`, `basisname !== nothing`, and
  `nesting = :pqs` or `nesting = :wl`;
- approved source surface is `src/cartesian_base_hamiltonian.jl` and
  `bin/cartesian_ham_builder.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`
  optional only for a direct one-owner RG/MWG genericity blocker;
- the implementation must use the existing base atom validation, terminal
  basis construction, residual Gaussian augmentation, exact augmented
  operators, residual MWG/IDA interaction, base K/U reuse, assembly, writer,
  readback, manifest, and provenance;
- it may select `legacy_atomic_gaussian_supplement(...)` for one-center inputs
  and keep the existing diatomic supplement loader for two-center inputs;
- it may relax only the canonical driver's supplemented `Natom == 2` guard and
  must otherwise preserve driver public inputs, ordering, hooks, spacing,
  stage labels, and artifact contract;
- it does not approve a separate atom-only Hamiltonian builder, new driver
  inputs, route switches, diagnostics, stop-after controls, route
  skeleton/shellification/terminal lowering changes, raw-block changes,
  residual-selection changes, MWG/IDA convention changes, artifact schema or
  reader changes, public API/export changes, solver/ECP work,
  heteronuclear/general geometry, translated atoms, committed tests, or Cr2
  workflow.

Composition lane status:

- the explicit `atom | z-axis diatomic`, `:pqs | :wl`, and
  `supplement = off | on` composition lanes now all have approved
  implementation authority under the current origin-centered atom and
  homonuclear z-axis diatomic geometry constraints;
- `HP-COMP-ATOMBOX-FN-01` separately approves one-center atom parent sizing in
  `src/cartesian_base_hamiltonian.jl`: public `basis.radius` is the physical
  box extent authority, parent axis counts must depend on radius plus
  `core_spacing` / existing spacing policy, and the fixed `2*q + 1` atom
  parent side count artifact must be removed;
- `HP-COMP-NS-FN-01` separately approves public size naming: `ns` is the
  requested source/nesting size, and route-local `q` is derived from `ns` and
  `nesting`;
- this atom box lane preserves origin-centered atom validation, explicit
  charge/electron-count validation, `nesting = :pqs` and `nesting = :wl`,
  supplemented atoms, artifact keys, manifest/provenance, and canonical driver
  inputs;
- it does not approve driver changes, route-family switches, raw-block
  changes, residual-selection changes, MWG/IDA convention changes, artifact
  schema or reader changes, public API/export changes, solver/ECP work,
  diagnostics/status/report payloads, committed tests, Cr2-specific workflow,
  translated atoms, non-origin atom support, element lookup/default tables,
  broad parent-construction rewrites, or diatomic sizing changes;
- remaining geometry, solver, ECP, public export, and Cr2-specific work still
  need later docs-only amendments before implementation may begin.

Approved R1 one-center base atom relaxation:

- `HP-R1-ATOM-FN-01` relaxes `cartesian_base_hamiltonian(...)` from
  one-center H-only validation to explicit origin-centered all-electron
  one-center atoms;
- `HP-R1-ATOM-WIRE-01` maps the public nuclear charge to the existing
  White-Lindsey atomic mapping `Z` and maps resolved public `core_spacing` to
  the private `parent_mapping_d`; public `d` is deprecated and, if
  temporarily accepted, must equal resolved `core_spacing`;
- `HP-R1-CORE-FN-01` freezes `core_spacing` as the single public
  near-nucleus physical scale while keeping `reference_spacing`, tail spacing,
  and box padding separate; White-Lindsey `Z` behavior is a mapping-shape rule
  or preset rule, not a second public knob except for the later narrow expert
  `s_factor` override under `HP-PQS-MAP-SFACTOR-*`; visible driver/project defaults
  such as `core_spacing = 0.3` are explicit resolved inputs and may be
  overridden for quick tests; routine correctness-test scalars must be tied to
  their exact override inputs, not described as physics-default results;
- atoms and diatomics must share the same producer workflow after the narrow
  geometry/shellification differences; no atom-only Hamiltonian builder,
  materialization path, route-stage object, report/status payload, or artifact
  shape is approved;
- `HP-R1-ATOM-TEST-01` keeps the committed H endpoint as the regression gate
  and allows ignored/user-run Be or Cr atom artifact checks;
- the R1 atom IDs do not approve translated atoms, ECP, solver workflow,
  artifact schema change, element lookup/default table, public API redesign,
  or source files outside `src/cartesian_base_hamiltonian.jl`;
- supplemented atoms are separately governed by `HP-COMP-SUPPATOM-*`.

Approved route-recipe cleanup:

- `HP-ROUTE-RECIPE-FN-01` allows `cartesian_recipe(...)` in
  `src/pqs_source_box_route_driver_helpers.jl` to construct only the
  `route_family`-selected subrecipe and to set the inactive subrecipe to
  `nothing` where needed for compatibility;
- `:pqs_source_box` route inputs must no longer require inactive
  `white_lindsey_*` fields, and `_cartesian_base_route(kind)` in
  `src/cartesian_base_hamiltonian.jl` may delete those unused fields;
- explicit `:white_lindsey_low_order` route support remains live and must keep
  using its WL route fields where selected;
- no canonical-driver changes, numerical kernel changes, shellification or
  terminal-lowering changes, materialization/artifact schema changes,
  report/status expansion, WL materialization deletion, new tests, or Cr2 run
  are approved by this cleanup.

Approved route-inventory type-surface cleanup:

- `HP-ROUTE-INV-FN-01` approves only
  `src/pqs_source_box_route_driver_helpers.jl`;
- the target is removal of runtime-keyed retained-unit inventory
  `NamedTuple{unit_keys}` shapes and runtime-keyed `pair_family_counts =
  NamedTuple{families}(...)`;
- approved replacements are vector-backed records/tables, stable dictionaries,
  or helper accessors with stable concrete types;
- public input `NamedTuple`s, fixed `NTuple{3,Int}` coordinates/dimensions,
  artifact sidecar tables, `RawProductBoxPlan.source_mode_indices`, terminal
  lowering contract tuples, and retained-unit transform-contract tuples remain
  out of scope.

Approved raw product source-mode inventory cleanup:

- `HP-RAW-SRCMODE-FN-01` approves vector-backed source-mode inventory cleanup
  in `src/cartesian_raw_product_sources/records.jl`,
  `src/cartesian_raw_product_sources/source_mode_indices.jl`, and
  `src/cartesian_raw_product_sources/summaries.jl`;
- narrow consumer wiring is approved only in
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`,
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`,
  and `src/cartesian_base_hamiltonian.jl`;
- the target is removal of `RawProductBoxPlan.source_mode_indices` and
  `source_mode_column_indices` variable-length tuple storage, while preserving
  deterministic source-mode ordering, fixed `NTuple{3,Int}` mode coordinates,
  retained-rule behavior, and manifest source-mode/relation facts;
- accessor semantics mean the same ordered facts and column associations, not
  preserving a variable-length `Tuple` concrete return type;
- terminal-lowering contract tuple cleanup and broader source-box/pair-block
  rewrites remain out of scope.

Approved contract-plan vector cleanup:

- `HP-CONTRACT-VEC-FN-01` approves vector-backed plan-level contract
  inventories in `src/cartesian_terminal_lowering/contracts.jl`,
  `src/cartesian_terminal_lowering/selection.jl`,
  `src/cartesian_terminal_lowering/summaries.jl`,
  `src/cartesian_retained_unit_transform_contracts/records.jl`,
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`, and
  `src/cartesian_retained_unit_transform_contracts/summaries.jl`;
- narrow consumer wiring is approved only if needed in
  `src/pqs_source_box_route_driver_helpers.jl`,
  `src/cartesian_base_hamiltonian.jl`, and
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`;
- the target is removal of variable-length tuple storage for
  `TerminalLoweringPlan.available_contracts`, `TerminalLoweringPlan.contracts`,
  and `RetainedUnitTransformContractPlan.contracts`;
- accessors `available_contracts(plan)`, `selected_contracts(plan)`,
  `contracts(plan)`, and `transform_contracts(plan)` must preserve ordered
  behavior and semantics without preserving variable-length `Tuple` concrete
  field types;
- `source_cpbs::Tuple{Vararg{CoordinateProductBox}}` remains out of scope.

Approved route/stage type-surface cleanup:

- `HP-ROUTE-STAGE-TYPE-FN-01` approves only
  `src/pqs_source_box_route_driver_helpers.jl` and
  `src/cartesian_terminal_shellification_geometry.jl`;
- approved targets are Be2 q5 compile-attributed compatibility inventories and
  wide route/stage return surfaces around
  `_pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan`,
  `cartesian_units`,
  `_pqs_source_box_route_driver_transform_stage_low_order_summary`,
  `cartesian_transforms`, and
  `_cartesian_terminal_shellification_region_unit_inventory`;
- stale compatibility inventories may be deleted, and remaining runtime-sized
  `NamedTuple` / `Tuple` carriers may be replaced with vector-backed compact
  internal objects, stable dictionaries, accessors, or smaller summaries;
- H2 base/supplemented artifact behavior, terminal shellification/lowering
  order, artifact/manifest schema, public driver contract, route semantics, and
  numerical matrices must remain unchanged;
- driver changes, artifact/manifest changes, numerical/raw-block/RG/MWG/IDA
  changes, route diagnostic/status/report expansion, committed tests,
  PackageCompiler/PrecompileTools/sysimage work, and Cr2 workflow remain
  unapproved.

Approved route/stage carrier cleanup:

- `HP-ROUTE-STAGE-CARRIER-FN-01` approves only
  `src/pqs_source_box_route_driver_helpers.jl` and
  `src/pqs_source_box_diatomic_complete_core_shell.jl`, with
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  optional only where directly required to slim terminal realization plan
  carriers in the approved path;
- approved targets are broad stage signatures and plan carriers around
  `cartesian_shells`, `cartesian_units`, `cartesian_transforms`, terminal
  topology support-region planning, terminal retained-rule planning, and
  directly required terminal realization plan carriers;
- giant shellification, route-skeleton, support-plan, retained-rule-plan, and
  terminal-plan `NamedTuple` / tuple shapes may stop crossing stage function
  signatures; necessary facts should move to compact typed/vector-backed
  carriers, smaller summaries, accessors, or local recomputation from canonical
  objects;
- route skeleton construction semantics and
  `src/pqs_source_box_route_driver_skeletons.jl` remain out of scope;
- H2 base/supplemented artifact behavior, H2 R3 endpoint if touched,
  terminal support/shellification/lowering order, public driver contract,
  artifact/manifest schema, route semantics, and numerical matrices must remain
  unchanged;
- driver changes, artifact/manifest changes, numerical/raw-block/RG/MWG/IDA
  changes, route diagnostic/status/report expansion, committed tests,
  PackageCompiler/PrecompileTools/sysimage work, and Cr2 workflow remain
  unapproved.

Approved canonical driver atom workflow:

- `HP-DRV-ATOM-FN-01` approves explicit one-center atom input normalization in
  `bin/cartesian_ham_builder.jl` for `mode = :base`;
- `HP-DRV-ATOM-WIRE-01` lets the driver call the existing
  `cartesian_base_hamiltonian(system; basis, hamfile)` facade for
  origin-centered base atom construction where the base facade already
  supports the requested atom;
- current validation remains origin-centered H; these driver atom IDs do not
  approve broader base atoms, translated atoms, element tables, ECP, solver
  workflow, public API/export changes, or artifact schema changes;
- supplemented atom driver wiring is separately governed by
  `HP-COMP-SUPPATOM-*`.
- `HP-DRV-ATOM-CLEAN-01` approves only removing the stale hidden
  `d = core_spacing` atom-basis field from the canonical driver; public
  inputs, defaults, overrides, hooks, timing labels, stage sequence, artifact
  schema, and driver contract must remain unchanged.

Approved homonuclear z-axis diatomic supplemented workflow:

- `HP-R3U-ZDI-FN-01` relaxes the non-exported supplemented facade from
  hardcoded H2/Be2 checks to explicit homonuclear two-center z-axis diatomic
  validation;
- `HP-R3U-ZDI-WIRE-01` lets canonical driver supplemented mode call the
  supported facade instead of package-internal helper paths;
- required inputs remain explicit: symbols, charges, spin sectors, geometry,
  base basis controls, supplement labels, and optional trusted `basisfile`;
- no heteronuclear, non-z-axis, ECP, solver, public export, artifact schema,
  route diagnostic, or Cr2-specific branch is approved.

## Residual Gaussian Authority

`residual_gaussian_domain_module.md` is the canonical current RG algorithm
contract. Do not duplicate the full residual-selection, exact-transform, or
MWG-interaction algorithm in this file.

Essential live guardrails:

- residual Gaussian directions are selected separately on each atom and then
  merged once;
- residual occupation is a physical retained-direction measure, not numerical
  rank;
- exact augmented one-body/moment transforms are not the MWG approximation;
- MWG descriptors are rotation-dependent and must be computed from the final
  merged residual basis;
- RG does not own artifact writing, artifact provenance, basis loading, facade
  parsing, parent lattice construction, terminal topology, raw Gaussian block
  formula ownership, or public exports.
- final residual identity validation may use the approved
  `HP-RG-ORTHO-FN-01` absolute/relative check only for small floating-point
  overshoots after owner-local selection and a healthy final merge; it is not
  permission to relax residual selection or retain low-occupation directions.
- `HP-RG-IDTOL-FN-01` sets the default final residual identity validation
  tolerance to `1.0e-8` in the older Be tolerance lane. It is superseded for
  production defaults by `HP-RG-CUTOFF-FN-01` and then
  `HP-RG-CUTOFF-FN-02`.
- `HP-RG-CUTOFF-FN-02` supersedes the residual occupation default:
  `residual_occupation_cutoff = 1.0e-6`; `identity_atol = 5.0e-8` remains
  unchanged. Owner grouping, merge metric failure rules, `G' S R` validation,
  width/zeta filtering, MWG/IDA, artifacts, driver workflow, and public API
  remain unchanged.
- `HP-RG-SPECTRAL-AUDIT-01` is measurement-only authority for residual-sector
  spectra after that cutoff: count residuals by owner, compute low `K_RR` and
  `H1_RR`, report owner weights and residual-occupation composition for
  low/flagged modes, and compare with one-center atom baselines when
  available. It is not permission to implement a spectral guard or change
  selection policy.

Exact Cartesian Gaussian raw blocks are separate neutral kernel authority:
uncharged by-center nuclear blocks under
`cartesian_gaussian_raw_blocks_nuclear.md`, and non-nuclear
overlap/kinetic/moment blocks under
`cartesian_gaussian_raw_blocks_non_nuclear.md`. RG consumes those blocks as
exact operator inputs; RG does not own their raw construction.

The old R3-A/R3-B/R3-C labels are implementation-history labels. The current
domain names are:

- `build_residual_gaussian_basis`;
- `transform_augmented_operator`;
- `moment_matched_gaussians`;
- `assemble_residual_ida_interaction`.

## Current Boundaries

Base pair/assembly role:

- Cartesian pair terms remain a structural precursor for the base route;
- `cartesian_assembly` remains topology/report context only;
- production Hamiltonian construction must use realized terminal basis and
  final-basis matrix kernels, not source-plan/report reconstruction.

Artifact ownership:

- base Hamiltonian files use the existing `CartesianIDAHamiltonian` artifact
  shape plus the fixed `producer_provenance/` group;
- supplemented R3/RG writes may add only the approved compact
  `supplement_provenance/` group through the existing helper outside the RG
  module;
- `read_cartesian_ida_hamiltonian` ignores provenance and reads matrices
  normally.

Deferred lanes:

- broad driver diagnostics and public-driver polish beyond the approved
  compact artifact-producing workflow;
- translated atoms; supported one-center supplemented atoms are governed by
  `HP-COMP-SUPPATOM-*`;
- Cr2-specific workflow, committed Cr2 gates, and Cr2 support decisions beyond
  the generic explicit homonuclear z-axis path;
- non-base/supplement public workflow;
- ECP, broad EGOI workflow beyond the approved retained-GTO helper, public RHF
  workflow, solver, and HamV6 work;
- artifact/public API decisions beyond the approved compact provenance groups.

Before implementation, confirm the approved ID and source surface in
`registry.md`; if the needed surface is absent, do a docs-only amendment first.
