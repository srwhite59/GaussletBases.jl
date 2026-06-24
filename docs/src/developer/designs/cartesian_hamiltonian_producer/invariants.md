# Cartesian Hamiltonian Producer Invariants

These are architecture-wide guardrails. Detailed algorithms live in the module
contract pages, especially `residual_gaussian_domain_module.md` for Residual
Gaussian work.

## Authority

- `registry.md` owns approved IDs, files, functions, and source ownership.
- `current.md` owns live status and startup direction.
- `residual_gaussian_domain_module.md` owns the canonical Residual Gaussian
  math/algorithm contract.
- Historical design and review files are rationale unless explicitly promoted
  by the compact current authority.

## Terminal Basis

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
- R3/RG terminal `G-G` product-matrix optimization may own only final-basis
  product matrices for kinetic, coordinate moments, and second moments used by
  residual-Gaussian exact augmented operators. It must not change terminal
  basis realization, `G-A`/`A-A` raw blocks, unit-nuclear Gaussian sums, IDA,
  MWG, residual selection, or route/public/artifact surfaces.
- Remaining R3/RG exact-operator allocation after terminal `G-G` product
  workspace reuse was attributed by `HP-R3REM-AUDIT-01`. `HP-R3UN-FN-01` may
  own only terminal final-basis unit-nuclear `U_GG` Gaussian-sum allocation
  reduction. Route/raw-block setup, neutral raw blocks, residual/MWG/IDA work,
  public workflow, and Cr2 facade/artifact work remain outside this authority.
- Slice C/base IDA assembly produces a real final-basis matrix and then uses the
  existing `CartesianIDAHamiltonian`; no Hamiltonian wrapper or result payload is
  approved.

## Residual Gaussian Guardrails

The canonical RG algorithm contract is
`residual_gaussian_domain_module.md`. Keep only these non-negotiable rules here:

- residual basis directions are selected separately on each physical owner atom
  and then merged once;
- residual occupation is not numerical rank, an integral weight, or a tolerance
  repair knob;
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
- RG does not own artifact writing, artifact provenance, basis loading, facade
  parsing, public exports, driver workflow, or route-stage/report fields.

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

## Public Boundaries

- R1 public base scope is origin-centered H, Cartesian z-axis H2, and explicit
  origin-centered all-electron one-center atoms under `HP-R1-ATOM-*`.
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
- Geometry (`atom` or z-axis diatomic), `nesting` (`:pqs` or `:wl`), and
  supplement state (`off` or `on`) are intended to compose through shared
  producer boundaries. Do not implement the 2 x 2 x 2 matrix as eight
  driver-level special cases or parallel Hamiltonian builders.
- R3/RG supplemented usability remains internal unless a later public/export
  amendment approves it.
- The canonical Cartesian driver is an artifact-producing workflow over
  approved producer surfaces. It may own compact editable defaults, trusted
  local input-file loading, command-line overrides, coarse timing, compact
  printing, artifact writing, and readback checks. It must not own route
  diagnostics, private stage controls, raw provider switches, report/status
  payloads, solver workflow, artifact schema, public exports, or
  Cr2-specific support.
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
