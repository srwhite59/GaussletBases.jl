# Residual Gaussian Domain Module

Status: canonical algorithm authority for the internal
`CartesianResidualGaussians` domain module. This is the current home for RG
residual selection, exact augmented operators, and MWG interaction convention.
It approves no public export, no artifact schema change, no Cr2 facade support,
and no new production behavior beyond the approved owner-local
residual-GTO/MWG construction and domain migration.

## Decision

Residual Gaussian construction is a domain concept, not a terminal-basis
detail. The current implementation in
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` works for
the narrow H2 path, but the file name and helper names flatten three different
physics operations:

- residual Gaussian directions are selected by atom-local residual occupation
  and then merged;
- one-body and moment operators are transformed exactly into the augmented
  basis;
- matched-width Gaussians are only the residual-containing IDA interaction
  approximation.

Production code should make those meanings explicit in the internal module:

```text
src/cartesian_residual_gaussians/
  CartesianResidualGaussians.jl
  residual_basis.jl
  augmented_operators.jl
  mwg_interaction.jl
```

`R3-A`, `R3-B`, and `R3-C` are implementation-history labels. They are not
permanent source concepts and should not be used as production function,
module, or object names in the new owner.

## Approved IDs

- `HP-RG-FILE-01` - Residual Gaussian domain module files and root include
  plumbing.
- `HP-RG-OBJ-01` - residual Gaussian basis object and invariant fields.
- `HP-RG-FN-01` - owner-local residual Gaussian basis construction.
- `HP-RG-FN-02` - exact augmented one-body and moment transformation.
- `HP-RG-FN-03` - moment-matched Gaussian descriptor construction.
- `HP-RG-FN-04` - residual-containing MWG/IDA interaction assembly.
- `HP-RG-WIRE-01` - migration from `pqs_terminal_residual_gto.jl` to the
  domain module.
- `HP-RG-TEST-01` - validation gates for the migration.
- `HP-RG-ORTHO-FN-01` - narrow robust final residual
  orthogonalization/identity validation.
- `HP-RG-ORTHO-TEST-01` - validation gates for the robustness pass.
- `HP-RG-IDTOL-FN-01` - final residual identity tolerance default.
- `HP-RG-IDTOL-TEST-01` - validation gates for the tolerance-default pass.
- `HP-RG-CUTOFF-FN-01` - residual occupation cutoff and identity tolerance
  defaults.
- `HP-RG-CUTOFF-TEST-01` - validation gates for the cutoff/tolerance policy.
- `HP-RG-CUTOFF-FN-02` - production residual occupation cutoff tightening.
- `HP-RG-CUTOFF-TEST-02` - residual-only validation for the tightened cutoff.
- `HP-RG-SPECTRAL-AUDIT-01` - measurement-only residual-sector spectral audit.
- `HP-RG-INJECT-AUDIT-01` - measurement-only optional injection audit.
- `HP-RG-INJECT-FN-01` - default-off in-memory injection-plus-RG
  implementation.
- `HP-RG-PROTECT-INJECT-DESIGN-01` - completed compact-main rationale for the
  implemented protected-localized basis convention.
- `HP-RG-OCC-FIRST-INJECT-AUDIT-01` - completed historical measurement for
  occupied-first global injection.
- `HP-RG-OCC-FIRST-INJECT-FN-01` / `HP-RG-OCC-FIRST-INJECT-TEST-01` -
  implemented occupied-first geometry governed by
  `occupied_first_injection.md`.
- `HP-RG-PROTECT-INJECT-FN-01` / `HP-RG-PROTECT-INJECT-TEST-01` and
  `HP-RG-PROTECT-ONEBODY-FN-01` / `HP-RG-PROTECT-ONEBODY-TEST-01` -
  implemented internal geometry and one-body convention governed by
  `protected_localized_basis.md`.
- `HP-RG-PROTECT-ONEBODY-AUDIT-01` and
  `HP-RG-PROTECT-VEE-AUDIT-01` - completed historical evidence for that
  convention.
- `HP-RG-PROTECT-ART-FN-01` - opt-in protected-localized Hamiltonian artifact
  variant.
- `HP-RG-PROTECT-ART-TEST-01` - validation gates for the protected artifact
  variant.
- `HP-RG-PROTECT-ARTLOC-FN-01` - row-locality metadata for protected
  artifacts.
- `HP-RG-PROTECT-ARTLOC-TEST-01` - validation gates for protected artifact
  row-locality metadata.
- `HP-RG-PROTECT-EGOI-AUDIT-01` - measurement-only EGOI audit over
  protected-localized artifacts.
- `HP-RG-PROTECT-EGOI-FN-01` - retained-GTO local-product EGOI helper.
- `HP-RG-PROTECT-EGOI-TEST-01` - validation gates for the retained-GTO EGOI
  helper.
- `HP-RG-PROTECT-LADDER-XFER-AUDIT-01` - measurement-only same-parent
  protected-localized ladder transfer audit.
- `HP-RG-PROTECT-LADDER-BUNDLE-FN-01` - opt-in protected-localized ladder
  bundle writer/reader facility.
- `HP-RG-PROTECT-LADDER-BUNDLE-TEST-01` - validation gates for the ladder
  bundle facility.
- `HP-RG-RHO0-GAL-AUDIT-01` - measurement-only rho0/Galerkin IDA correction
  audit over the protected-localized baseline.
- `HP-RHO0-REFDENS-AUDIT-01` - measurement-only reference-density-matrix IDA
  correction audit; successor to the row-gauge rho0/Galerkin probes.
- `HP-RHO0-REFDENS-MIXH-AUDIT-01` - measurement-only exact mixed Hartree seam
  audit for `(final final | reference reference)`.
- `HP-RHO0-MIXH-GG-FN-01` - first source-backed exact Hartree `GG` block
  builder for one-center atomic `P_A`.
- `HP-RHO0-MIXH-GG-TEST-01` - validation gates for the first `GG` mixed
  Hartree source seam.
- `HP-RHO0-MIXH-GAAA-FN-01` - exact mixed Hartree `GA`/`AA` extension for
  one-center atomic `P_A`.
- `HP-RHO0-MIXH-GAAA-TEST-01` - validation gates for the `GA`/`AA` mixed
  Hartree extension.
- `HP-RHO0-MIXH-FEXACT-FN-01` - exact-side final/protected Hartree transform
  from raw mixed Hartree `GG`/`GA`/`AA` blocks.
- `HP-RHO0-MIXH-FEXACT-TEST-01` - validation gates for the exact-side
  final/protected Hartree transform.
- `HP-RHO0-FAPP-AUDIT-01` - measurement-only approximate-side `F_app[P0]`
  seam audit; this is not source authority.
- `HP-RHO0-FAPP-FN-01` - narrow Cartesian IDA approximate energy/Fock source
  seam for fixed represented spin densities.
- `HP-RHO0-FAPP-TEST-01` - validation gates for the Cartesian IDA approximate
  energy/Fock seam.
- `HP-RHO0-ANCHOR-FN-01` - superseded in-memory Hartree reference
  correction-anchor lane.
- `HP-RHO0-ANCHOR-TEST-01` - validation gates for the superseded Hartree
  correction anchor.
- `HP-RHO0-CORR-AUDIT-01` - suspended measurement-only corrected-Hamiltonian
  small-system audit.
- `HP-RHO0-JANCHOR-FN-01` - direct-Hartree reference anchor replacement.
- `HP-RHO0-JANCHOR-TEST-01` - validation gates for the direct-Hartree anchor
  replacement.
- `HP-RHO0-XPAIR-AUDIT-01` - measurement-only exchange/direct pairing audit
  after direct-Hartree anchor validation.

Implementation IDs in this list are approved only within the surfaces below.
Design-only IDs record authority for future source blurbs but do not approve
implementation by themselves.

## Ownership

### Owned By CartesianResidualGaussians

The module owns:

- owner-local residual Gaussian basis selection;
- residual occupation spectra and cutoff policy;
- atom-local orthonormalization or localization;
- final inter-owner merge;
- final `T_G` and `T_A` transforms;
- exact augmented `K`, uncharged by-center `U_A`, `x`/`y`/`z`, and
  `x^2`/`y^2`/`z^2` transformation;
- moment-matched Gaussian descriptors;
- residual-containing IDA interaction blocks.

### Not Owned By CartesianResidualGaussians

The module must not own:

- basis-set file loading or named-basis lookup;
- parent lattice construction;
- terminal shell topology or terminal support ownership;
- raw analytic Gaussian integral formulas;
- exact uncharged by-center Cartesian Gaussian nuclear `G-A`/`A-A` raw-block
  construction;
- exact Cartesian Gaussian overlap, kinetic, coordinate-moment, and
  second-moment `G-A`/`A-A` raw-block construction;
- facade input parsing;
- terminal or base stage orchestration;
- Hamiltonian artifact writing;
- `supplement_provenance/` schema ownership;
- route reports, route stages, status symbols, driver parsing, public API, or
  public exports.

The module may consume already-validated base Hamiltonian facts, terminal
basis realization, parent axis/bundle data, supplement representations, raw
exact blocks, and donor factor tables. Exact uncharged nuclear raw blocks are
owned by the neutral Cartesian Gaussian raw-block nuclear owner when that
kernel is present. Exact non-nuclear overlap/kinetic/moment raw blocks are
owned by the neutral non-nuclear raw-block owner when that kernel is present.
RG consumes these blocks as exact operator inputs to
`transform_augmented_operator(...)`; it must not turn those inputs into a new
report or provider payload framework.

## Approved Source Surfaces

### HP-RG-FILE-01

Approved source files:

```text
src/cartesian_residual_gaussians/CartesianResidualGaussians.jl
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_residual_gaussians/augmented_operators.jl
src/cartesian_residual_gaussians/mwg_interaction.jl
```

`src/GaussletBases.jl` may add only the internal include needed to load
`cartesian_residual_gaussians/CartesianResidualGaussians.jl`. No public export
is approved. Include-order adjustments are allowed only as needed to place the
new module after its existing internal dependencies and before approved
callers; they must not reorder unrelated subsystems for style.

### HP-RG-OBJ-01

The residual Gaussian basis object is a numerical domain object, not a
status/result payload. It must carry the final residual transformation
authority and enough compact policy facts to make the basis reproducible:

- base dimension;
- candidate count;
- residual dimension;
- candidate owner indices;
- residual source owner indices;
- owner retained counts;
- retained residual occupations;
- `occupation_cutoff = 1.0e-6`;
- `tau_neg_abs = 1.0e-12`;
- `tau_neg_rel = 1.0e-12`;
- `tau_merge_abs = 1.0e-12`;
- `tau_merge_rel = 1.0e-12`;
- selection rule;
- orientation rule;
- sign rule;
- `T_G::Matrix{Float64}`;
- `T_A::Matrix{Float64}`.

It must not carry route metadata, report fields, status flags, raw candidate
payload clouds, full discarded spectra, MWG center/width matrices, dense
moment matrices, artifacts, or public API state.

### HP-RG-FN-01

Approved production name:

```julia
build_residual_gaussian_basis(...)
```

The exact Julia signature may be adjusted to local types, but the constructor
must require candidate owner indices. Parsing owner identity from labels is
not approved.

The algorithm is the approved owner-local residual-selection rule:

```text
for each owner a:
    X_a = G' S A_a
    M_a = S_AaAa - X_a' X_a
    select by residual occupation lambda > eta_RG
    form an owner-local orthonormal residual sector
concatenate owner sectors
merge by final symmetric Lowdin over inter-owner overlap
return final T_G and T_A
```

If all owner-local modes are retained, preserve donor-style full-rank owner
orientation by owner-local symmetric Lowdin in candidate order. If rank is
lost, use deterministic owner-local natural residual modes ordered by
decreasing residual occupation with deterministic sign canonicalization.

Final merge failure rule: with
`tau_merge = max(1.0e-12, 1.0e-12 * max(lambda_max(S_merge), 1.0))`, any merge
eigenvalue below `-tau_merge` is a construction error and any merge eigenvalue
`<= tau_merge` is a near-singular merge error. Eigenvalue flooring must not be
used to preserve directions.

`HP-RG-ORTHO-FN-01` additionally permits robust final residual identity
validation after this merge when the owner metrics are positive/full rank and
the merge spectrum is healthy. The final `R' S R` check may symmetrize the
computed overlap and use
`err_RR <= 1.0e-10 + 1.0e-10 * max(1, scale_RR)`, where `err_RR` is the maximum
absolute identity error and `scale_RR` is the maximum absolute entry or
equivalent infinity-norm scale of the symmetrized residual overlap. This is
not an occupation-cutoff change, a residual-selection change, or permission to
floor merge eigenvalues.

`HP-RG-IDTOL-FN-01` sets the default final residual `R' S R` identity
validation tolerance to `1.0e-8`. This updates only the final identity
acceptance threshold in the older Be tolerance lane. That production default
is superseded by `HP-RG-CUTOFF-FN-01`.

`HP-RG-CUTOFF-FN-01` set the prior production defaults:
`residual_occupation_cutoff = 5.0e-8` and `identity_atol = 5.0e-8`, addressing
the Cr atom `basis_ns = 9`, `map_ns = 11`, `lmax = 1` direction at occupation
`3.637e-8`.

`HP-RG-CUTOFF-FN-02` supersedes the residual occupation default:
`residual_occupation_cutoff = 1.0e-6`; `identity_atol = 5.0e-8` remains
unchanged. The selection cutoff is an owner-local residual occupation policy:
directions below the default cutoff, including the cited Cr2
`1.27e-7` to `8.98e-7` marginal directions, are discarded by default. This
does not change owner grouping, negative eigenvalue tolerances, final merge
metric checks, `G' S R` validation, width/zeta filtering, MWG/IDA, artifacts,
driver workflow, public API, or the approved source owner.

`HP-RG-SPECTRAL-AUDIT-01` is measurement-only authority after the cutoff
cleanup. It records that the Cr2 retained count drops to `62 + 62`, but
residual-only spectra still show a low two-owner mode with
`min eig(K_RR) = 0.3700413519`, `min eig(H1_RR) = -7.1647854052`, and owner
weights about `0.5 / 0.5`. Ignored probes may classify low residual-sector
modes by owner weights, residual-occupation composition, and one-center atom
baselines when available. This audit does not approve a kinetic/`H1_RR` guard,
automatic pruning, cutoff/tolerance changes, source instrumentation, full HF,
artifacts, driver work, or MWG/IDA changes.

The optional injection-plus-RG idea is recorded separately in
`residual_gaussian_injection_hybrid.md`. That memo proposes classifying
near-gausslet supplement modes by owner-local residual norm, globally merging
the injected subspace, replacing the corresponding gausslet-sector directions,
and then returning to owner-local residual selection for true RGs.
`HP-RG-INJECT-AUDIT-01` approves only ignored measurement probes for that
proposal. It is not production source authority and does not change the
approved RG defaults, MWG/IDA convention, artifacts, driver workflow, or
public API.

`HP-RG-INJECT-FN-01` approves the historical default-off in-memory direct
`G`-injection implementation in the RG owner. It may add the numerical
authority needed for an injected replacement base sector and exact one-body
transformation into `[F, R]`, but injected directions must not become
residual-GTO/MWG channels. When `residual_injection_cutoff <= 0`, the current
production RG behavior must remain unchanged within roundoff. This ID does not
approve driver input, public API, artifact schema/provenance changes,
production-default changes, full HF, Cr2 workflow, or spectral pruning.

`HP-RG-PROTECT-INJECT-DESIGN-01` supplied the completed compact-first
rationale now governed by
[Protected-localized basis convention](protected_localized_basis.md). It is not
source authority by itself; the implemented geometry and one-body IDs below
own the source surface.

`HP-RG-OCC-FIRST-INJECT-AUDIT-01` is completed historical measurement
evidence. The implemented geometry and validation contract is owned by
[Occupied-first injection geometry](occupied_first_injection.md) under
`HP-RG-OCC-FIRST-INJECT-FN-01` and
`HP-RG-OCC-FIRST-INJECT-TEST-01`.

The helper makes identified `Y_occ` mandatory before optional capture
selection, distinguishes pre-inclusion capture from post-inclusion recovery,
rejects malformed complement/capture geometry, and does not turn rejected weak
directions into MWG residual channels. It is not wired into staged
protected-original geometry; composition over `M = [G, R_compact]` belongs to
`HP-RG-PROTECT-ADDREF-*`.

`HP-RG-PROTECT-INJECT-FN-01` / `HP-RG-PROTECT-INJECT-TEST-01`
and `HP-RG-PROTECT-ONEBODY-FN-01` / `HP-RG-PROTECT-ONEBODY-TEST-01`
implement the internal/default-off geometry and exact one-body parts of
[Protected-localized basis convention](protected_localized_basis.md).
`residual_basis.jl` owns staged geometry; `augmented_operators.jl` owns exact
fixed-sector and localized one-body transforms.

The one-body audit is completed historical evidence. The completed
`HP-RG-PROTECT-VEE-AUDIT-01` rejected direct `C' V C` interaction rotation
and established localized `L`, exact `H1_L`, and inherited pre-injection
site-order `Vee_M` as the viable convention. Rejected broad directions remain
basis-insufficiency diagnostics and never become MWG residual channels.

`HP-RG-PROTECT-ART-FN-01` approves a narrow opt-in artifact variant for that
protected-localized baseline. The artifact may persist `H1_L`, inherited-site
`Vee_L`, electron counts, a recognized convention/version ID, recipe and
commit provenance, basis controls, compact-RG and injection counts, localized
ordering, sector maps, representability diagnostics, and
orthogonality/localization diagnostics. Primary ownership is
`augmented_operators.jl` and `cartesian_ida_hamiltonian.jl`; `residual_basis.jl`
and `cartesian_base_hamiltonian.jl` are optional only for already-computed
geometry, diagnostics, or producer plumbing. The lane does not approve driver
or public API work, default producer changes, rho0, solver methods, artifact
semantic changes for ordinary PQS/WL/RG files, selection-policy changes,
`C' V C`, or Cr2 production claims. Readback must reject missing or
unrecognized convention/version fields before any consumer uses the file.

`HP-RG-PROTECT-ARTLOC-FN-01` approves row-locality metadata for the same
artifact. Centers must be diagonal position expectations in the actual
protected-localized `L` basis, computed from inherited main-space position
operators and the native `ML` transform. The artifact may add native-order
`center_x/y/z`, deterministic `z_order_to_native` and `native_to_z_order`
permutations, per-row sector labels or native-sector indices, and optional
spread diagnostics when second moments already exist. `H1_L` and `Vee_L`
remain native-order canonical matrices; z sorting is metadata for consumers
and must not mutate matrix order or native sector ranges.

`HP-RG-PROTECT-EGOI-AUDIT-01` approves only ignored measurement probes for
existing matrix-level EGOI routines over the protected-localized artifact
objects. It may reconstruct `Qtarget` from current protected/injection
geometry, build exact Gaussian target Coulomb for selected target orbitals,
and report H/Be/Be2 diagnostics before any optional bounded Cr2 diagnostic.
It does not approve source edits, artifact variants, public workflow, solver
workflow, RG/injection selection changes, protected-artifact convention
changes, rho0/P0 revival, or Cr2 production claims.

`HP-RG-PROTECT-EGOI-FN-01` approves the first source-backed retained-GTO EGOI
helper. The target is owner-balanced retained original supplement `s1+s2`
GTOs mapped from compact retained source indices, not broad protected-`Z`,
atom-HF orbitals, final rows, or residualized RG functions. The source lane is
internal and in-memory: local atom products are first-class targets, the
AA-BB local-product Coulomb block is included in the acceptance metric, AB
overlap products are not default targets, and disallowed/long-range `DeltaV`
must remain zero. Primary ownership is `hamiltonian_corrections.jl`, with
`augmented_operators.jl` and `residual_basis.jl` optional only for retained
source mapping or transform-ready `Qtarget`. The lane does not approve
artifacts, public workflow, solver integration, selection changes, rho0/P0,
or `s3`/`p`/`d` target promotion.

`HP-RG-PROTECT-LADDER-XFER-AUDIT-01` approves only a measurement-first
same-parent protected-localized ladder transfer audit. It may use ignored
probes and `/Users/srw/dmrgtmp` outputs to build same-parent `ns = 7 -> 9`
and optional `9 -> 11` protected-localized inherited-site Hamiltonians,
compute exact final-basis cross overlaps, transfer occupied orbitals with
`C_B = <B|A> C_A`, evaluate the transferred density with target `H1_L` and
`Vee_L`, and run short bounded UHF continuation only after trace and
orthonormality checks pass. It does not approve source edits, public workflow,
durable artifact schema, transformed source Hamiltonians, transformed source
`Vee`, `C' V C`, rho0/P0 revival, EGOI expansion, or Cr2 production claims.
Final working bases are treated as orthonormal, so the audit must not use a
generalized self-overlap transfer.

`HP-RG-PROTECT-LADDER-BUNDLE-FN-01` approves a reusable opt-in
protected-localized ladder bundle facility. Preferred ownership is a small
`cartesian_protected_ladder_bundle.jl` orchestrator that reuses
`cartesian_ida_hamiltonian.jl` for protected-localized Hamiltonian artifacts
and `cartesian_representation_transfer.jl` for exact final-basis cross
overlaps and orbital transfer. The bundle may write member Hamiltonian
artifacts, exact `S_BA = <L_B|L_A>` transfer sidecars, optional transferred
occupied-orbital restart sidecars, a versioned manifest with shared-parent
proof and diagnostics, and bounded summary TSVs. It does not approve package
exports, driver defaults, solver workflow, transformed source Hamiltonians,
transformed source `Vee`, `C' V C`, rho0/P0, EGOI expansion, protected-`Vee`
convention changes, or Cr2 production claims. `HP-RG-PROTECT-LADDER-BUNDLE-TEST-01`
approves only package load, `git diff --check`, small bundle/readback smoke if
feasible, and ignored Cr2 `ns = 7 -> 9` bundle parity/readback validation.

`HP-RG-RHO0-GAL-AUDIT-01` approves only an ignored measurement audit for a
rho0/Galerkin IDA correction on top of that protected-localized baseline. It
may use existing protected-localized geometry and one-body/Vee data, analytic
IDA/Coulomb sanity checks, small H/He/H2 checks, Cr2 fixed-density
diagnostics, and one bounded Cr2 HF replay only if static rho0/Galerkin
diagnostics are sane. It does not approve source edits, source-backed IDA/MWG
implementation, artifact support, public wiring, production Hamiltonian
workflow, `C' V C` revival, broad rejected directions as MWG channels, Vee
scaling as the fix, screened-reference production claims, Cr2 production
energy claims, or publication-scale validation sweeps.

Later row-gauge measurements retired that lane as the target formulation.
`HP-RHO0-REFDENS-AUDIT-01` is the active screened-correction measurement
target: use a fixed reference density matrix `P0`, compare exact and
approximate reference Fock/energy for the same represented `P0`, and include
the energy anchoring constant
`C0 = E_exact0 - E_app0 - sum_sigma Tr(P0_sigma * Delta_F0_sigma)`.
The first audit is Hartree-only and must validate `P0` representability in the
protected-localized basis plus exact/approximate finite-difference Fock checks.
This ID does not approve source edits, artifact support, public wiring,
production Hamiltonian or solver workflow, Cr2 production claims, direct
`C' V C`, residual/MWG default changes, basis-fate changes, or broad rejected
directions as MWG residual channels. Candidate future IDs
`HP-RHO0-REFDENS-FN-01` and `HP-RHO0-REFDENS-ERI-01` are named only for
planning and are not approved.

`HP-RHO0-REFDENS-MIXH-AUDIT-01` approves only an ignored measurement audit for
the exact mixed Hartree seam needed by the reference-density lane:
`(final final | reference reference)`. The scope is H/Be/Be2 only, exact
Hartree only, with a helper/kernel inventory and a Be/Be2 fixed-`P0` audit
only if the seam is feasible without source edits. If no existing seam can
compute the operator, the audit must report the smallest neutral source owner
needed. It does not approve source edits, artifacts, public workflow,
Cr/Cr2 diagnostics, HF exchange, direct `C' V C`, row action, `diag(J)`, `q0`,
center metadata, IDA proxy shortcuts, residual/MWG defaults, basis-fate
changes, or broad rejected directions as MWG residuals.

`HP-RHO0-MIXH-GG-FN-01` is source authority, but not RG ownership. It approves
only a neutral `CartesianGaussianRawBlocks` exact Hartree `GG` helper for
one-center atomic `P_A`: same-center reference pair-density terms,
Coulomb-expanded separable one-body packets, and terminal/base `GG` output.
`GA`/`AA`, protected transforms, `F_app[P0]`, correction constants, artifacts,
public workflow, Cr/Cr2, exchange, residual/MWG default changes, and
basis-fate policy remain outside this source lane.

`HP-RHO0-MIXH-GAAA-FN-01` remains in the neutral raw-block owner and extends
that exact Hartree source seam only to `GA` and `AA`. It is not authority for
Residual Gaussian transforms, protected-localized correction assembly,
approximate Fock construction, artifacts, public workflow, Cr/Cr2, exchange,
or residual/MWG/basis-fate policy.

`HP-RHO0-MIXH-FEXACT-FN-01` is Residual Gaussian transform authority, but only
for the exact-side reference Hartree one-body operator. It may consume raw
mixed Hartree `GG`/`GA`/`AA` blocks and existing protected-original geometry to
produce in-memory `F_exact_Hartree[P0]`. It is not authority for approximate
Fock construction, correction constants, artifacts/public workflow, IDA/MWG
interactions, geometry selection changes, Cr/Cr2, or solver work.

`HP-RHO0-FAPP-AUDIT-01` is measurement-only authority to identify the
approximate-side derivative seam. It may use ignored probes to compute
`E_app[P0]` and validate `F_app[P0]` as the derivative of the actual current
IDA/MWG approximate energy convention for represented `P0_final`. It is not
source authority and does not approve correction constants, artifacts/public
workflow, Cr/Cr2, row-action/`diag(J)`/`q0`/center-metadata/direct-`C' V C`
shortcuts, or residual/MWG/basis-fate changes.

`HP-RHO0-FAPP-FN-01` is not Residual Gaussian ownership. It approves only a
paired approximate energy/Fock seam on `CartesianIDAHamiltonian` in
`src/cartesian_ida_hamiltonian.jl`, for fixed represented spin densities and
the Hamiltonian's stored `electron_electron_ida` convention. The helper may be
used later by rho0 audits, but it does not approve correction assembly,
artifacts/public workflow, Cr/Cr2, residual selection changes, or broad
rejected directions as MWG residuals.

`HP-RHO0-ANCHOR-FN-01` is Residual Gaussian adjacent only because the
protected/final exact-side helper currently lives in
`augmented_operators.jl`. It is implemented but superseded for Hartree-
correction physics/stability interpretation: the old anchor subtracted the
full approximate interaction Fock, including the current same-spin
exchange-like term, from an exact Hartree operator. Do not use its
`Delta_F0_alpha/beta` as the Hartree reference-density correction.

`HP-RHO0-CORR-AUDIT-01` is measurement-only, but suspended until
`HP-RHO0-JANCHOR-*` replaces the old full-interaction anchor. A corrected-
Hamiltonian audit using the old `Delta_F0_alpha/beta` is invalid as
Hartree-correction physics/stability evidence.

`HP-RHO0-JANCHOR-FN-01` is not broad Residual Gaussian ownership. It approves
only a direct-Hartree anchor replacement across
`src/cartesian_ida_hamiltonian.jl` and
`src/cartesian_residual_gaussians/augmented_operators.jl`. The direct
approximate side is `E_app_direct = 1/2 * q' * V * q` and
`F_app_direct = Diagonal(V * q)`, with
`q = diag(P_alpha) + diag(P_beta)`. The anchor forms
`Delta_J0 = F_exact_Hartree[P0] - F_app_direct[P0]` and `C0_J`, then reruns
H/Be/Be2 corrected behavior with `Delta_J0`/`C0_J`. It does not approve
artifacts, public workflow, solver integration, Cr/Cr2, exact exchange
correction, changes to the current approximate exchange convention, residual
selection changes, or production correction use.

`HP-RHO0-XPAIR-AUDIT-01` is measurement-only. It records the current rho0
interpretation: supplement-space atomic `P0` is viable, direct-Hartree anchor
algebra is viable, but corrected Hamiltonian behavior with inherited
approximate exchange-like terms is not ready. The audit may use ignored
H/Be/Be2 probes to compare direct-Hartree correction, inherited approximate
exchange-like expectations, exact/supplement-space exchange diagnostics where
feasible, and an explicitly non-production direct-only corrected-operator
diagnostic. It does not approve source edits, artifacts/public workflow,
solver integration, Cr/Cr2, exact exchange implementation, changes to the
current approximate exchange convention, or direct-only physics as a
Hamiltonian.

Do not approve a vague global entry point such as
`stabilize_residual_metric(...)`. Global raw-candidate symmetric Lowdin and
global raw-column pivoted-Cholesky selection are not the Residual Gaussian
basis algorithm.

### HP-RG-FN-02

Approved production name:

```julia
transform_augmented_operator(...)
```

This function family owns exact transformation of raw `[G, A]` blocks into
the augmented `[G, R]` basis. For any exact one-body or moment operator `O`:

```text
O_GR = O_GG T_G + O_GA T_A
O_RR =
    T_G' O_GG T_G
  + T_G' O_GA T_A
  + T_A' O_AG T_G
  + T_A' O_AA T_A
O_aug = [O_GG  O_GR
         O_GR' O_RR]
```

This applies to kinetic, every uncharged by-center nuclear attraction,
`x`, `y`, `z`, `x^2`, `y^2`, and `z^2`. It is exact operator transformation,
not an MWG approximation.

### HP-RG-FN-03

Approved production name:

```julia
moment_matched_gaussians(...)
```

This builds residual matched-width Gaussian descriptors from the final merged
residual functions and exact transformed moment matrices:

```text
c_ralpha = <r | alpha | r>
v_ralpha = <r | alpha^2 | r> - c_ralpha^2
sigma_ralpha = sqrt(2 * v_ralpha)
```

Nonfinite moments, nonpositive variances, or nonpositive widths are
construction errors. The descriptors are for residual-containing interaction
approximation only. They must not be treated as exact residual-GTO Coulomb
integrals or base-PQS IDA weights.

### HP-RG-FN-04

Approved production name:

```julia
assemble_residual_ida_interaction(...)
```

This owns residual-containing MWG/IDA interaction blocks:

```text
V_aug = [V_GG_base  V_GM
         V_GM'      V_MM]
```

`V_GG_base` is unchanged from the base Hamiltonian. `V_GM` must use the
approved weight-aware final-basis density-normalized contraction for PQS shell
blocks, and `V_MM` uses the direct density-normalized matched-Gaussian
interaction. Term-first pair-factor reuse and bounded workspace remain
binding requirements.

## Do Not Confuse

- Residual occupation is not numerical rank. It is the atom-local amount of
  candidate Gaussian content left outside the gausslet span.
- Residual occupation is not residual integral weight. Residual integral
  weights may be near zero or sign-changing.
- Atom-local orthonormalization or localization is not the final inter-owner
  merge.
- Exact one-body transformation is not the MWG approximation.
- MWG interaction is not invariant under arbitrary residual rotations; the
  residual orientation policy is part of the interaction convention.

## Migration And Cleanup Sequence

1. Create `CartesianResidualGaussians` and move residual-basis construction
   first. The existing H2 endpoint must keep the approved owner-local scalar
   `0.4574265214362075`.
2. Move exact augmented operator transformation next. The H2 one-body/moment
   checks must remain unchanged.
3. Move moment-matched Gaussian descriptors and residual IDA interaction
   assembly. The independent weight-aware `V_GM` check remains required.
4. Rewire the same-construction residual-GTO/MWG Hamiltonian path and the R3U
   facade to call the new domain functions.
5. Delete old R3-named helpers from
   `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` once
   no live callers remain.

The migration may use small compatibility wrappers in
`pqs_terminal_residual_gto.jl` only for existing callers while the callers are
rewired. Wrappers must delegate to the new domain module, must not preserve
the old global-selection algorithm, must not add new behavior, and must be
deleted when the last live caller moves. If a source pass can switch all
callers in one commit, prefer deleting the wrappers immediately.

## Validation

Future source migration must validate:

- existing standalone H2 residual-GTO/MWG endpoint, including augmented
  dimension `489` and lowest-orbital IDA self-Coulomb
  `0.4574265214362075` within `1.0e-10`;
- exact one-body/moment checks: `G' S R`, `R' S R`, base G-G block equality,
  finite/symmetric `K`, uncharged `U_A`, and moment matrices;
- independent weight-aware `V_GM` check;
- R3U facade section, if touched, against the same scalar and artifact
  readback deltas;
- ignored Be2 owner-local usability/performance measurement under `tmp/work`
  when the source pass changes the interaction path or facade wiring.
- for `HP-RG-ORTHO-FN-01`, ignored strict N2 q5 p10 residual audit/artifact
  smoke at `core_spacing = 0.042857`, plus one passing N2 comparison at
  `core_spacing = 0.05` or `0.075`, reporting `G' S R`, `R' S R - I`, retained
  counts, and merge spectrum.
- for `HP-RG-IDTOL-FN-01`, Be atom cc-pV5Z `lmax = 1` residual audit/artifact
  validation with the same `21` retained residual directions, Be atom cc-pVDZ
  `lmax = 1` comparison, unchanged H2 residual-GTO/MWG endpoint, and reporting
  of `R' S R - I`, allowed tolerance, retained count, minimum retained
  occupation, final merge condition, and `G' S R`.
- for `HP-RG-CUTOFF-FN-01`, Cr atom
  `basis_ns = 9`, `map_ns = 11`, `lmax = 1` residual construction passes or
  cleanly drops the marginal `s4` direction at occupation `3.637e-8` as
  intended, Be atom cc-pV5Z still passes, H2 residual-GTO/MWG endpoint remains
  unchanged, and the handoff reports retained counts, minimum retained
  occupation, `G' S R`, `R' S R - I`, allowed tolerance, and final merge
  condition.
- for `HP-RG-CUTOFF-TEST-01`, the existing H2 endpoint test
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` may update only
  its two cutoff assertions from `1.0e-8` to `5.0e-8`:
  `residual.occupation_cutoff` and `values[:occupation_cutoff]`.
- for `HP-RG-CUTOFF-FN-02`, the default residual occupation cutoff changes to
  `1.0e-6` while `identity_atol` remains `5.0e-8`.
- for `HP-RG-CUTOFF-TEST-02`, run a Cr2 residual-only audit, not full HF or a
  Cr2 artifact/workflow run. Owner retained counts should drop from
  `68 + 68` to `62 + 62`; recompute and report `min eig(K_RR)`,
  `min eig(H1_RR)`, and low-mode candidate composition. If low-H1 ghost modes
  remain, stop and request separate kinetic/`H1_RR` spectral-guard authority.
  Be high-zeta and H2 residual-GTO/MWG endpoints must still pass. The existing
  H2 endpoint test may update only its two cutoff assertions from `5.0e-8` to
  `1.0e-6`: `residual.occupation_cutoff` and
  `values[:occupation_cutoff]`.
- for `HP-RG-SPECTRAL-AUDIT-01`, run only measurement probes under ignored
  `tmp/work` scripts, with durable text/TSV output under `/Users/srw/dmrgtmp`
  or CR2 run directories. Report residual count by owner, low `K_RR`, low
  `H1_RR`, owner weights, residual-occupation composition, and one-center atom
  baselines when available. Do not run full HF, write a new Hamiltonian
  artifact, add source instrumentation, or implement pruning/guards.

No Cr2 full Hamiltonian, Cr2 artifact, Cr2 facade support, public export,
driver/bin/tool workflow, artifact schema expansion, report/status/payload
object, solver/RHF, ECP, or EGOI work is approved by this amendment.

## Artifact Writing Boundary

Compact supplemented artifact writing remains outside
`CartesianResidualGaussians`. Artifact writing is workflow/provenance glue
attached to the supported R3 usability path, not residual Gaussian
mathematics.

`CartesianResidualGaussians` may produce residual basis objects, exact
augmented operators, moment-matched Gaussian descriptors, residual interaction
blocks, and in-memory Hamiltonian ingredients consumed by a workflow writer.
It must not own the artifact schema, JLD2 file workflow, facade input parsing,
or `supplement_provenance/` policy.

The current writer location remains acceptable:

```julia
write_pqs_terminal_residual_gto_augmented_hamiltonian(...)
```

That writer stays with the terminal/facade workflow under the existing R3-C
artifact authority unless a later docs-only amendment approves moving or
splitting it. Future movement requires a named duplication or consumer reason;
"closer to RG" is not sufficient justification.
