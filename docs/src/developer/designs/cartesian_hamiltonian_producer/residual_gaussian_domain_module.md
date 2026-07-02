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
- `HP-RG-PROTECT-INJECT-DESIGN-01` - design-only protected-original
  injection over compact main space.
- `HP-RG-PROTECT-INJECT-FN-01` - narrow internal staged
  protected-original geometry prototype.
- `HP-RG-PROTECT-INJECT-TEST-01` - validation gates for the staged geometry
  prototype.
- `HP-RG-PROTECT-ONEBODY-AUDIT-01` - measurement-only audit for exact
  one-body transformation into the protected fixed sector.
- `HP-RG-PROTECT-ONEBODY-FN-01` - narrow internal exact one-body transform for
  the protected fixed sector.
- `HP-RG-PROTECT-ONEBODY-TEST-01` - validation gates for the protected
  one-body transform replay.
- `HP-RG-PROTECT-VEE-AUDIT-01` - measurement-only protected fixed-sector Vee
  interaction audit.

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

`HP-RG-PROTECT-INJECT-DESIGN-01` records the current compact-first direction:
build compact/narrow RGs first, define `M = [G, R_compact]`, then inject
original supplement Gaussians by replacement `F = [Z, M Q_perp]`. Protected
narrow originals are orthonormalized first in original GTO overlap, remaining
originals are Gram-cleaned separately from injection representability, and
`B = M' S Z` must be full rank and well conditioned. If a good-norm broad
original is not represented by `M`, that is an insufficient-main-basis
diagnostic, not permission to make the candidate a MWG residual Gaussian. This
ID is design-only and requires a later source amendment before implementation.

`HP-RG-PROTECT-INJECT-FN-01` approves only an internal geometry prototype for
the staged protected-original construction in `residual_basis.jl`. It may
source-back the measured sequence: construct `M`, protected originals, and the
broad original subspace `W`; filter by representability singular values of
`B = M' S W`; optionally localize and classify shape for diagnostics; filter
by fake-RDM eigenspace occupancy; and report geometry diagnostics for
`Z = [Z_protected, Z_broad]` and `F = [Z, M Q_perp]`. It does not approve
operator/Hamiltonian transformation, inherited IDA/MWG for protected-original
injection, artifact support, public wiring, Cr2 HF, or default behavior
changes.

`HP-RG-PROTECT-ONEBODY-AUDIT-01` approves only an ignored measurement audit
that consumes the source-backed staged geometry and existing exact one-body
data to test in-memory fixed-sector transforms such as `F' K F`, `F' U_A F`,
and `F' H1 F`. It does not approve source helpers, artifact/provenance
changes, public wiring, IDA/MWG interaction transforms, Cr2 HF, residual
default changes, or production Hamiltonian claims. If a later source lane is
needed, exact one-body ownership likely belongs in `augmented_operators.jl`,
while `residual_basis.jl` remains the geometry owner.

`HP-RG-PROTECT-ONEBODY-FN-01` approves that narrow source lane after the
successful audit. `augmented_operators.jl` owns private exact dense in-memory
transforms of `K`, per-center uncharged `U_A`, and assembled `H1` into the
protected fixed sector `F = [Z, M Qperp]`. `residual_basis.jl` may change only
if transform-ready geometry accessors or fields are missing. This lane does
not approve IDA/MWG interaction transforms, artifact support, public wiring,
matrix-action frameworks, Cr2 HF, residual default changes, or production
Hamiltonian claims.

`HP-RG-PROTECT-VEE-AUDIT-01` approves only an ignored measurement audit for an
in-memory Vee candidate in the same protected fixed sector. It may consume the
source-backed protected geometry, one-body helpers, and existing in-memory
interaction data to test whether transforming the `M = [G, R_compact]`
interaction through `F = [Z, M Qperp]` is finite, symmetric, and not
anomalously cheap for broad-`Z` directions. A bounded in-memory Cr2 HF replay
is allowed only after that Vee diagnostic gate passes. This ID does not approve
source edits, source-backed IDA/MWG interaction implementation, artifact
support, public wiring, production Hamiltonian workflow, screened-reference/
rho0 work, Vee scaling as the primary fix, rejected broad directions as MWG
residual channels, or Cr2 production claims.

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
