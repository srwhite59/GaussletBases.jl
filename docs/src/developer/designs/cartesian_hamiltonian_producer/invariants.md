# Cartesian Hamiltonian Producer Invariants

These are architecture-wide guardrails. Detailed algorithms live in the module
contract pages, especially `residual_gaussian_domain_module.md` for Residual
Gaussian work.

## Authority

- `AGENTS.md` supplies the deny-by-default execution whitelist and repo-wide
  operating rules. Presence there is necessary but not sufficient for source
  work.
- `registry.md` owns each ID's permission, lifecycle, owner, approved source or
  test surface, and canonical contract link. Registry section nesting is
  navigation only; an entry's own status governs.
- `current.md` owns concise live status: implemented facilities, active lanes,
  blockers, and next work. It does not restate full contracts.
- subsystem design pages own normative mathematics, object semantics,
  invariants, failure behavior, and validation contracts. In particular,
  `residual_gaussian_domain_module.md` owns the canonical Residual Gaussian
  algorithm contract.
- `README.md` and `algorithm_implementation_index.md` are navigation. They do
  not independently authorize source work.
- `implementation_slices.md`, planning memos, `history/`, `reviews/`, reports,
  and the manager running log are decision/evidence records unless explicitly
  promoted by current authority.

An ID authorizes source work only when the `AGENTS.md` whitelist, registry
permission/lifecycle, and canonical subsystem contract agree and `current.md`
does not explicitly contradict them. Silence in the concise current-status
page is neutral; an explicit conflict is not. On conflict, fail closed: make no
source edit and request a docs-only authority reconciliation. This transition
rule remains in force until a later amendment establishes one machine-readable
authority source.

## Terminal Basis

The detailed object, realization, one-body, IDA, and Hamiltonian boundary is
[terminal basis and base assembly](terminal_basis_and_base_assembly.md).

- Terminal basis blocks are supported only on their owned terminal support
  regions.
- Previous-block projection, recursive projection, and effective-support growth
  onto earlier terminal regions are forbidden.
- Cross-block overlap is zero by structure because parent rows are orthonormal
  and terminal supports are disjoint. A nonzero structural overlap means
  duplicated support rows, wrong row restriction, wrong ownership, or indexing
  error.
- Cross-block kinetic, nuclear, one-body, and IDA terms may be nonzero and must
  be assembled blockwise.
- Do not form a global parent overlap, global support operator, global dense
  coefficient matrix, global final overlap, or global Lowdin path.

## Lowdin

- Shell-local and residual merge Lowdin operations use
  `inv(sqrt(Symmetric(overlap)))`.
- Do not hand-roll eigendecomposition transforms as the production Lowdin
  implementation.
- Do not floor small eigenvalues to preserve directions that the approved
  residual or merge policy rejects.

## One-Body And IDA Assembly

- Final-basis one-body operators are assembled from terminal block pairs using
  realized block supports and coefficients.
- Nuclear attraction uses term-first Gaussian contraction and reusable 1D factor
  data where possible.
- Neutral Cartesian Gaussian raw-block nuclear kernels may own only exact
  uncharged by-center `G-A` and `A-A` raw matrices, `U_A = -1/r_A`. They do not
  apply physical nuclear charges, perform terminal projection, transform into
  residual bases, or own overlap/kinetic/moment blocks.
- Neutral Cartesian Gaussian raw-block non-nuclear kernels may own only exact
  overlap, kinetic, coordinate moment, and second-moment `G-A` and `A-A` raw
  matrices. They do not perform terminal projection, optimize final-basis
  `G-G` product matrices, transform into residual bases, or own nuclear/MWG
  interaction blocks.
- R3/RG terminal `G-G` product-matrix optimization owns only final-basis
  product matrices for kinetic, coordinate moments, and second moments used by
  residual-Gaussian exact augmented operators. It must not change terminal
  basis realization, `G-A`/`A-A` raw blocks, unit-nuclear Gaussian sums, IDA,
  MWG, residual selection, or route/public/artifact surfaces.
- Remaining R3/RG exact-operator allocation after terminal `G-G` product
  workspace reuse was attributed by completed `HP-R3REM-AUDIT-01`.
  `HP-R3UN-FN-01` owns only terminal final-basis unit-nuclear `U_GG`
  Gaussian-sum construction. Trusted base-block reuse is separately limited
  to one construction and changes no operator. Route/raw-block setup, neutral
  raw blocks, residual/MWG/IDA work, public workflow, and Cr2 facade/artifact
  work remain outside these authorities.
- Slice C/base IDA assembly produces a real final-basis matrix and then uses the
  existing `CartesianIDAHamiltonian`; no Hamiltonian wrapper or result payload is
  approved.
## Residual Gaussian Guardrails

The canonical RG algorithm contract is
`residual_gaussian_domain_module.md`. Keep only these non-negotiable rules here:

- residual basis directions are selected separately on each physical owner atom
  and then merged once;
- the ordinary production residual-occupation cutoff is not numerical rank, an
  integral weight, or a tolerance repair knob. The explicit numerical-complete
  opt-in is the sole current exception and uses the owner-local residual-metric
  spectrum only to discard numerical nulls at its fixed `1e-10` threshold;
- exact augmented one-body/moment transformation is not the MWG approximation;
- raw Cartesian Gaussian block construction is not Residual Gaussian basis
  selection or augmented-operator transformation;
- MWG centers and widths are computed from the final merged residual basis and
  are not invariant under arbitrary residual rotations;
- weight-aware final-basis density normalization is required for PQS-shell
  `V_GM` blocks;
- final residual identity validation may use the narrow
  `HP-RG-ORTHO-FN-01` combined absolute/relative check only for small
  floating-point overshoots after owner-local selection and a healthy final
  merge; it must not become an occupation-cutoff change, eigenvalue flooring,
  or global residual-selection repair;
- `HP-RG-IDTOL-FN-01` sets the default final residual identity validation
  tolerance to `1.0e-8` in the older Be tolerance lane. It is superseded for
  production defaults by `HP-RG-CUTOFF-FN-01` and then
  `HP-RG-CUTOFF-FN-02`;
- `HP-RG-CUTOFF-FN-02` supersedes the RG residual occupation default:
  `residual_occupation_cutoff = 1.0e-6`; `identity_atol = 5.0e-8` remains
  unchanged. Occupations below the new default cutoff, including the Cr atom
  `3.637e-8` marginal direction and the cited Cr2 `1.27e-7` to `8.98e-7`
  directions, are not production-retained residual directions by default.
  Owner-local grouping, merge checks,
  `G' S R` validation, width/zeta filtering, MWG/IDA, artifacts, driver
  workflow, and public API remain unchanged;
- RG does not own artifact writing, artifact provenance, basis loading, facade
  parsing, public exports, driver workflow, or route-stage/report fields.

## Protected Additive References

- Protected molecular geometry uses the exact compact residual basis consumed
  by augmented operators. Do not rebuild compact selection inside the geometry
  helper.
- The full-rank orthonormal union of placed packet occupied spaces is used only
  to make their span mandatory in the protected basis.
- Each packet's original occupied columns and occupations remain separate and
  define its additive contribution to `P0`. Do not globally orthogonalize the
  packet blocks when constructing the reference density.
- Mandatory packet occupied directions are never optional capture candidates.
  Rank or representability failure means insufficient compact-main-basis
  support and must stop the construction.
- Molecular reference fields and energies are additive:
  `P0 = sum_a P_a`, `J0 = sum_a J_a`, and
  `E0 = sum_a E_aa + 2*sum_{a<b}E_ab`.
- Packet density fits define self/cross Coulomb energies. Fitted potentials are
  only fast evaluators of the same placed reference field.
- Transform `J0` through existing raw `GG/GA/AA`, protected fixed-sector, and
  localized `W` owners. Keep native `L` ordering and inherited-site `Vee_L`;
  do not rotate or transform `Vee`.
- `Delta_J0` and `C` remain an in-memory screened direct-Hartree correction.
  They do not mutate unscreened `H1_L/Vee_L` and are not a corrected artifact
  without later authority.

## Provenance And Metadata

- Metadata may contain provenance only. It must not carry transforms,
  coefficients, matrices, source plans, runtime inventories, or numerical data
  needed by algorithms.
- The base `producer_provenance/` group and supplemented
  `supplement_provenance/` group are compact artifact provenance, not algorithm
  inputs.
- Artifact provenance must record public construction-family truth. `nesting`
  is a public recipe fact, and route labels must be derived from
  `(input.kind, input.nesting)`, not from PQS-oriented helper names or default
  route-stage vocabulary.
- Existing Hamiltonian readers must continue to read matrices normally and
  ignore provenance unless a later approved reader contract says otherwise.
- One Cartesian/PQS Hamiltonian construction has one Coulomb Gaussian
  expansion. The producer must resolve `coulomb_accuracy` before parent/PGDG
  construction and carry the same `CoulombGaussianExpansion` through base
  unit-nuclear/IDA, residual-GTO exact Coulomb-expanded blocks, and MWG.
  Construction stages must not independently choose compact or high accuracy.
- A new base or supplemented artifact records one Hamiltonian-wide
  `coulomb_expansion/` summary. Protected-localized members and ladder
  manifests preserve that summary on readback. Missing legacy provenance is
  unavailable, not evidence of high accuracy. Atomic reference packets are
  different: their RHF, density/self-energy, and fitted-potential scaffold
  expansion provenance is role-qualified because those are separate evaluated
  reference objects.
- High Coulomb exponents must be handled by algebraically stable analytic
  Gaussian identities. Do not infer a carrier limitation from cancellation in
  determinant or weighted-variance formulas, and do not repair such failures
  by clamping, exponent truncation, scaled/log PGDG carriers, or terminal
  contraction changes without separate authority.

## Public Boundaries

- The implemented R1 facade, exact inputs, and base provenance are canonical in
  [R1 public base producer](r1_public_base_producer.md). Its scope is
  origin-centered all-electron one-center atoms and homonuclear Cartesian
  z-axis diatomics.
- Atoms and diatomics must share the same base producer workflow after
  geometry/shellification normalization. Do not introduce atom-only
  Hamiltonian builders, materialization paths, route-stage/report/status
  objects, or artifact shapes when the shared terminal-basis, one-body, IDA,
  Hamiltonian-construction, and writer machinery applies.
- Base z-axis homonuclear diatomic relaxation uses explicit public symbols,
  nuclear charges, centers, and electron counts. Symbols are provenance labels;
  charges and explicit `nup`/`ndn` are authority. Do not add element lookup
  tables or inferred electron counts.
- One-center atom parent sizing uses physical extent. Public `basis.radius`
  / driver `padding` is the atom box-size authority.
- Public size input uses `ns` as the requested cube/source/nesting size.
  Route-local `q` is derived only after selecting `nesting`: `q = ns` for PQS
  and `q = ns - 2` for White-Lindsey. Do not treat `q` as the common public
  cube-size field, direct parent side-count rule, or physical box extent.
- White-Lindsey z-axis diatomics with normalized `ns < 4` are unsupported and
  should reject before route construction. Working WL diatomic `ns` ranges may
  saturate the final retained support when physical parent extent dominates;
  this is not an input-ignored bug and must not be "fixed" by changing
  shellification, terminal lowering, or driver semantics without separate
  authority.
- `nesting = :wl` must converge to the same terminal-basis downstream boundary
  as `nesting = :pqs`: a `CartesianTerminalBasisRealization` with disjoint
  owned terminal supports. White-Lindsey boundary-stratum realization must not
  revive old WL H1/H1+J materialization, route reports, or route-stage
  diagnostics.
- WL z-axis diatomic boundary-stratum/product units must not treat
  full-support identity retention as the production compact retained basis.
  The approved compact WL diatomic basis is generated from products of
  one-dimensional contractions on the owned unit support. Do not fake
  compactness by dropping support rows, relabeling identity units, or changing
  the driver comparison.
- White-Lindsey odd-side centering is required only for direct
  nucleus-centered core blocks, where odd side length keeps the nucleus
  centered on a grid point. Boundary shell strata, WL boundary-stratum retained
  products, and non-direct support regions must not inherit this rule. A
  boundary retained count must preserve the requested shell contraction count;
  for WL public `ns = 4`, route-local `q = 2` gives `4^3 - 2^3 = 56`
  boundary columns, not `26`.
- Direct nucleus-centered core side must be derived from public `ns`, not
  route-local `q`: `direct_core_side = isodd(ns) ? ns : ns + 1`. Route-local
  `q` remains `q = ns` for PQS and `q = ns - 2` for White-Lindsey, but it does
  not define the shared direct core side.
- Terminal shell decomposition is common geometry. Direct core regions,
  shell regions, owned support rows, ordering, and coverage must be computed
  by one route-family-free operation before PQS or White-Lindsey lowering.
  For z-axis diatomics, PQS and White-Lindsey must call that common shellifier
  with the same first-step arguments when the public system, parent axes,
  public `ns`, direct core side, nuclear centers, and bond axis match.
  Central-gap/contact, shared-shell, and outer-mismatch ownership are common
  shell geometry, not retained-construction geometry.
  PQS may then use shell support plus a full source CPB for retained-mode
  realization; White-Lindsey may then split shell boundaries into faces, edges,
  corners, and strata for product-of-1D contractions. Those are
  retained-construction geometries, not separate first-step shellifiers.
- A private semantic shell `source_q` is a source-span refinement, not a second
  public/route `q`. It may be applied only after common shellification and must
  not change public `ns`, PQS route `q`, direct-core side, parent axes, shell
  boxes, support ownership, or region order. Shared-shell longitudinal `L`
  remains owned by the existing angular-band selector; callers may not supply
  `L` or a full shape. The first override lane is ordinary-span only and does
  not establish mapped-COMX or artifact semantics.
- Z-axis diatomic thin slab stacks are not real shell regions and are not
  direct identity sectors. For both PQS and White-Lindsey,
  `:direct_midpoint_slab` and `:outer_mismatch_slab` must lower through the
  same compact slab function with the same terminal region, public `ns`, slab
  normal/thickness, and native source/support facts. The compact unit-slice
  retained scale is `ns x ns x 1` after one-dimensional COMX/product
  compression; a thickness-`t <= ns` outer-mismatch region should scale about
  `t * ns * ns`. Direct/core sectors remain identity; real shell regions
  remain route-specific. If slab thickness exceeds `ns`, source work must stop
  for a policy decision; it must not silently delete slabs, retain them as
  identity rows, infer slab geometry from labels, or invent
  route-family-specific lowering.
- Terminal-shellification metadata inventories must agree with this thin-slab
  contract. Metadata/scaffold summaries may describe planned compact
  thin-slab lowering, but they must not preserve stale direct identity mapping
  for midpoint slabs, outer-mismatch fallback slabs, or
  `:angular_z_extension_slab`. Those inventories are not numerical authority:
  they must not materialize coefficients, carry Hamiltonian data, or create a
  parallel report/payload workflow.
- Compact face-product terminal coefficient assembly is neutral
  `CartesianFinalBasisRealization` machinery. White-Lindsey facets and
  thin-slab stacks should reuse the same internal face-like product helper
  over fixed normal-axis indices. Shared face-product assembly must not live
  in `white_lindsey_terminal_basis_realization.jl`, must not be duplicated as
  a PQS-specific slab projection path, and must not relabel thin slabs as WL
  boundary strata merely to reuse names.
- Z-axis diatomic shellification, not lowering, owns the angular-balanced
  shared molecular box rule. Each shared-shell step should compute a target
  box in physical parent-axis coordinates so the longitudinal margin from each
  outer nucleus is balanced against the selected transverse scale. The
  ordinary index-layer shell body may be underextended along the bond axis;
  the ordinary body plus planned z-extension thin-slab stacks, not the
  ordinary body alone, realizes the angular-balanced target coverage. The
  thin-slab category applies uniformly to midpoint slabs, planned
  non-boundary angular z-extension slabs, planned boundary angular z-extension
  slabs, and any unexpected outer-mismatch fallback slab. Planned z-extension
  stack slices lower through the same compact thin-slab machinery as midpoint
  slabs; real shells remain route-specific only after common shellification.
- Geometry (`atom` or z-axis diatomic), `nesting` (`:pqs` or `:wl`), and
  supplement state (`off` or `on`) are intended to compose through shared
  producer boundaries. Do not implement the 2 x 2 x 2 matrix as eight
  driver-level special cases or parallel Hamiltonian builders.
- Source-span options are mainline carried-space/raw-source facts, not
  high-order workflow branches. A mapped-COMX source span changes only the
  raw one-dimensional span passed into the existing physical-coordinate COMX
  cleanup. It may change source-axis transform facts and provenance, but
  Hamiltonian, operator, artifact, and driver layers must not branch on mapped
  versus ordinary source except through descriptive source provenance.
- Materialized source-axis transform facts become basis-defining only at the
  terminal shell seed seam after validation against the shell `outer_box`,
  source shape, and coefficient matrix sizes. Consuming carried facts must not
  change boundary mode selection, owned-support row restriction, shell-local
  Lowdin, sign canonicalization, support validation, artifacts, or driver
  inputs.
- The canonical `source_span` selector is a construction choice over source
  spans, not a diagnostic route switch. It must not expose route records,
  retained-rule dumps, raw-block switches, stop-after controls, allocation
  probes, high-order workflow controls, artifact schema changes, or another
  COMX path.
- R3/RG supplemented usability remains internal unless a later public/export
  amendment approves it.
- The canonical Cartesian driver is an artifact-producing workflow over
  approved producer surfaces. It may own compact editable defaults, trusted
  local input-file loading, command-line overrides, coarse timing, compact
  printing, artifact writing, and readback checks. It must not own route
  diagnostics, private stage controls, raw provider switches, report/status
  payloads, solver workflow, artifact schema, public exports, or
  Cr2-specific support.
- Canonical driver terminal-region inventories are bounded human-facing
  summaries, not route reports. They may print region kind, support/final
  counts, compression, identity-vs-compact realization, and native slab facts
  needed to catch accidental identity sectors. They must not dump source modes,
  pair inventories, raw blocks, all rows, full metadata, artifact schema, or
  private route-stage objects.
- The canonical driver may expose an explicit one-center base atom workflow
  only through the existing base facade. Current atom validation remains
  origin-centered H, while `HP-R1-ATOM-*` permits explicit origin-centered
  all-electron base atoms in the facade. `HP-COMP-SUPPATOM-*` separately
  permits supplemented one-center atoms through the common RG/MWG boundary.
  Translated atoms, element lookup/default tables, ECP, solver workflow, and
  artifact-schema changes require separate authority.
- The supplemented R3 usability facade may support explicit homonuclear
  two-center z-axis diatomics under `HP-R3U-ZDI-*`. This is not
  heteronuclear, non-z-axis, ECP, solver, public export, artifact-schema, or
  Cr2-specific workflow authority.
- Be2 is an internal performance/usability proxy, not a committed public gate.
- Cr2 may be used only as an explicit homonuclear z-axis ignored/user-run
  stress or usability case after H2/Be2 validation. Cr2 diagnostics do not
  authorize Cr2-specific branches, committed Cr2 gates, full production Cr2
  artifacts, ECP, solver workflow, or broad molecule support.

## Carrying Cost

- New source surfaces require an approved ID and owner in `registry.md`.
- Replacement work must switch live callers and delete obsolete helper pressure
  in the same lane unless an exact live caller blocks deletion.
- Tests should validate physics endpoints or small active module contracts, not
  stale helper names, blocker vocabulary, route metadata, or transition-only
  wrappers.
