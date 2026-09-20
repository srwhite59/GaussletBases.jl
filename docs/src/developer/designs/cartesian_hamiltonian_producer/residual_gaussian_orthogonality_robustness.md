# Residual Gaussian Orthogonality And Cutoff Policy

This page is the canonical numerical-policy contract for the ordinary
owner-local Residual Gaussian basis. It distinguishes physical residual
retention, negative-metric failure, final-merge conditioning, `G-R`
orthogonality, and final residual identity. These are different checks and
must not share an implicit tolerance.

The basis algorithm itself is canonical in
[Residual Gaussian domain module](residual_gaussian_domain_module.md).

## Lifecycle And Precedence

| ID | Lifecycle | Durable result |
| --- | --- | --- |
| `HP-RG-ORTHO-FN-01` | Approved; bounded implementation | Stable Terminal Residual Construction; matrix-only premerge contract retained |
| `HP-RG-ORTHO-TEST-01` | Approved; bounded validation | Physical construction/operator checks; existing R3A snapshots unchanged |
| `HP-RG-IDTOL-FN-01` | Implemented historical; superseded | Former `identity_atol = 1e-8` default |
| `HP-RG-IDTOL-TEST-01` | Completed historical evidence | Be high-zeta identity-tolerance audit |
| `HP-RG-CUTOFF-FN-01` | Implemented historical; cutoff superseded | Former `5e-8` cutoff; current `identity_atol = 5e-8` originated here |
| `HP-RG-CUTOFF-TEST-01` | Completed historical evidence | Cr marginal-direction and H2 assertion validation |
| `HP-RG-CUTOFF-FN-02` | Implemented historical default; superseded by Pass 642 | Former `residual_occupation_cutoff = 1e-6`; robustness unchanged |
| `HP-RG-CUTOFF-TEST-02` | Completed evidence; active maintenance | H2 cutoff/provenance assertions and residual-only Cr2 audit |

Implementation commits are `76453cc39` for robust final identity,
`47e56593c` for the historical `1e-8` identity default, `f0b662dca` for the
historical `5e-8` cutoff/current identity value, and `1f7f04e56` for the
former `1e-6` cutoff. Manager-log Passes 131-132, 155-156, 171-172, and
183-185 preserve the evidence and policy transitions.

The IDTOL and CUTOFF-01 IDs no longer grant source or test work. They remain
addressable historical records so old artifacts, reports, and Git history can
be interpreted without treating their defaults as current.

## Current Defaults

Pass 643 accepts these shared defaults at b64499918. The
ordinary builder, terminal augmentation, collinear facade and protected-ladder
missing-key fallback must agree:

```text
residual_occupation_cutoff = 1e-8
tau_neg_abs                = 1e-12
tau_neg_rel                = 1e-12
tau_merge_abs              = 1e-12
tau_merge_rel              = 1e-12
orthogonality_atol         = 1e-10
identity_atol              = 5e-8
```

Steven's explicit policy supersedes the former 1e-6 default, not its recorded
evidence. The bounded implementation belongs to HP-COLLINEAR-PQS-RG-FN-01/TEST-01
under [cutoff forwarding](pqs_residual_gto_working_basis.md#Collinear-Residual-Cutoff-Forwarding).
The completed H10 selection comparison supports 1e-8 without another campaign.
Explicit 1e-6 comparisons and tighter 1e-10 requests remain valid inputs; all
negative-metric, merge, orthogonality, identity, representation and screening
checks still apply. This does not qualify Cr2, complete H10 operators or HF.

The explicit `1e-10` cutoff in
[Numerical-complete residual basis](numerical_complete_residual_basis.md) is a
separate internal opt-in. It is passed explicitly with injection and
compactness filtering disabled. It is not a production default, identity
tolerance, integral weight, or general conditioning knob.

## Metric And Merge Failure Rules

For any owner-local or final merge metric with eigenvalues `lambda`, define

```text
tau = max(tau_abs, tau_rel * max(maximum(lambda), 1)).
```

The following rules are binding:

1. `minimum(lambda) < -tau` is a materially negative metric and throws.
2. Owner-local retention is the strict physical test
   `lambda > residual_occupation_cutoff`.
3. The final merge additionally requires `minimum(lambda) > tau_merge`.
4. A zero or near-singular merge is a hard failure.
5. No eigenvalue may be floored, clamped upward, or retained through an
   alternate rank rule merely to complete construction.

Negative-metric tolerances decide whether the represented overlap geometry is
physically valid. The occupation cutoff decides which positive owner-local
residual directions belong to the ordinary production basis. Neither is the
final identity tolerance.

## Final Orthogonality And Identity

The matrix-only formula below remains binding outside the function-aware
terminal domain of [Stable Terminal Residual Construction](@ref). That section
replaces its numerical evaluation, not its physical tolerances, on that domain.

After owner-local selection, all owner blocks are concatenated and normalized
by one symmetric inverse square root of the explicitly symmetrized final merge
metric. The implementation then requires

```text
norm(T_G + X*T_A, Inf) <= orthogonality_atol.
```

It recomputes the symmetrized residual overlap `S_RR` and defines

```text
identity_error = maximum(abs, S_RR - I)
identity_scale = maximum(abs, S_RR).
```

The current scale-aware acceptance check is

```text
identity_error <= identity_atol *
                  (1 + max(1, identity_scale)).
```

This check may absorb only final floating-point cleanup error after positive
owner metrics, a healthy final merge, and strict `G-R` orthogonality. It must
not retain a direction below the occupation cutoff, hide a negative metric,
or replace the near-singular merge failure.

## Residual Premerge Gram Arithmetic

Pass 650 authorized only the premerge repair for GB-H10-Q8-RESIDUAL-20260918,
under HP-RG-ORTHO-FN-01/TEST-01. Baseline is
`aa2ac41c513635bb5a7208848551eeef0edd754b`. This is a construction-arithmetic
repair, not a cutoff, rank, conditioning policy or physical-basis admission.

Pass 651 accepts b974887980ef90bca6e753c049376d4548883dc0 and consumes this
implementation grant, including the narrow snapshot amendment below. The
following bounds preserve the completed contract, not an active execution task.
Fresh q5/q8 retained 270 directions naturally; complete rounded-input maximum
errors were 3.09e-9/4.45e-9 and selected physical errors 5.24e-9/1.86e-8.
Those physical checks are not a full metric certificate or consumer admission.
The amended R3A owner passed 464/464 plus 64/64; the original 463/464 failure
is retained without a claim about its historical cause. Full implementation
CI35425393601 and Docs35425393550 passed. No physical thresholds changed.

Change only the premerge `S_merge` evaluation in
`src/cartesian/cartesian_residual_gaussians/residual_basis.jl` function
`finalize_residual_gaussian_transform`, at most eight added source lines.
For the existing input transforms use

```text
C = X * T_A0
D = T_G0 + C
S_merge = _rg_sym(D' * D + T_A0' * S_AA * T_A0 - C' * C)
```

This identity is general: never drop `D'D`, including on injected paths.
Delete the replaced premerge four-term call. Preserve `_rg_sym`, inverse
square root, signs, return values/shapes, owner selection/order/occupations,
all thresholds and the final independent four-term `residual_gaussian_overlap`
assertion byte-for-byte. Preserve ordinary/injected callers and the accepted
basic integral repair. No helper, new type/API/metadata/cache, compensated
complement framework, clamp, extra cleanup or fallback is authorized.

At most 35 added lines in existing `test/misc/runtests.jl`: compact deterministic
high-cancellation and nonzero-D cases with independent high-precision Gram
evaluation, unchanged final identity gate, and negative/near-singular rejection.
Require accurate normalized Gram max error <=5e-8 in these fixtures. Explicitly
check nonzero-D contribution so an ordinary-only simplification cannot pass.
These local checks cover cancellation and the shared injected boundary that
an ordinary endpoint alone would miss. The sole additional test edit is the
snapshot-tolerance amendment below; no tracked large fixture. Preserve existing
misc, core, public Cartesian/residual-GTO and nested
R3A/supplemented owners, including its injected case. Run normal full three-job
source CI/Docs and ordinary package/docs/authority/self-test/Documenter/log/diff
checks. Do not repeat unrelated angular or paper calculations locally.

### Bounded q5 And q8 Acceptance

Reuse the four hash-frozen staged inputs listed in the task and diagnosis;
never rebuild a terminal basis or use the old 210-direction transform. Fresh
production construction must pass for q5 and q8 with the same ordered 270
ordinary s/p/d candidates and nuclei. Require natural 270 rank/27 per owner,
unchanged cutoff1e-8 and all existing merge/orthogonality/scale-aware identity
checks; never force that rank or loosen a failure. Report max-entry, induced
infinity row-sum and symmetry separately, actual thresholds, spectra and norms.

Independently evaluate the complete small rounded-input residual Gram and
selected physical subspaces using the actual base metric and d-capable analytic
inputs. Include original implicated directions q8[141,114], q5[114,60] and
new worst diagonal/off-diagonal directions. Require rounded-input and selected
physical max error <=5e-8 and row-sum <=1e-5, with all signed contributions
visible. Selected row sums are not full physical bounds: report coverage and
limitations, not complete physical certification or permission for consumers.
No alternate matrix may substitute for a failing production check.

Keep total task-specific diagnostic/acceptance work within 20 numerical minutes,
16GiB RSS and 2GiB additional scratch beyond staged inputs. Approximately88s
is already consumed by diagnosis plus independent compact design checks;
existing-owner tests and normal CI/setup wall time are separate. Use the
existing guarded runner, selected/tiled contractions and machine-local scratch;
no full dense final physical-metric rebuild. Record premerge and residual
construction time/allocations separately from physical validation/I/O. A material
cost regression requires review, not an unqualified speed claim. All frozen
inputs and hchain evidence remain read-only. No Hamiltonians, fields, HF,
q9/H20, correlation, release/version changes or consumer restart.

### Evidence And Failure Rule

Diagnosis `tmp/reviews/h10-q8-residual-diagnosis-2026-09-18.md`, SHA-256
`1e7e7495519d082222c89ecc676fa3cae915e703438644d50f46fcb2792e5233`,
reproduced q8 error1.26078e-7 against scaled threshold1.000000063e-7.
Accurate rounded inputs still fail at1.13263e-7; propagated premerge error
explains it. Formula-only two-quadratic error propagated through the ORIGINAL
transform is5.212e-9; it is not a repaired basis. Independent selected physical
error9.63739e-8 also fails. q5 passes. No cutoff/rank defect was established.
Independent design review verified the general identity, including nonzero D,
with ten compact 256-bit checks; no replacement q8 result has been claimed.

If this exact candidate fails original or independent acceptance, budget,
input identity or performance checks, make no implementation commit; preserve
evidence and report. Do not select another remedy, expand physical coverage
beyond resources, relax tolerances, change rank semantics or edit callers.
Broader architecture/API/numerical policy needs a separate decision. The bounded
cycle permits at most two in-scope correction rounds, not automatic expansion.
Hchain-doer remains paused even on success; no successor task is authorized.

### Occupation Snapshot Portability Amendment

Steven explicitly approved AMEND-GB-H10-Q8-RESIDUAL-20260918-2 after the
required R3A owner failed 463/464 at its minimum occupation snapshot. The
observed difference 1.6895772975145107e-14 exceeded its 1e-14 bound; this
original failure remains evidence, not a retroactive pass or established cause.

Under HP-RG-ORTHO-TEST-01, additionally permit exactly two `atol` substitutions
from 1e-14 to 1e-12 on the adjacent minimum/maximum residual_occupations
golden-value assertions in
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`, originally
lines 501-502, plus at most two comment lines explaining numerical portability.
Both golden values remain unchanged. At their scales, approximately 5.10e-4
and 1.22e-2, the absolute bound corresponds to relative errors about 2e-9
and 8.2e-11: still strict regression checks, not physical validity criteria.
No other assertion, production threshold, cutoff, rank or policy may change.

Preserve the source +3/-1 and misc-test +32 draft and completed q5/q8 evidence.
Rerun the affected R3A owner including its previously unreached 64-check section;
reuse unchanged successful checks, then require normal full implementation
CI/Docs and independent review. Do not repeat q5/q8 acceptance or perform a
baseline attribution campaign merely to explain these final digits. This narrow
user-approved exception supersedes the original test-edit freeze only here;
all other failure rules, budgets, exclusions and the two-round limit remain.

## Stable Terminal Residual Construction

Pass 652 authorizes H10-SUPPLEMENT-IMPLEMENTATION-CYCLE-20260920, baseline
67316863580b808d0ebe911f9fdf2e31906e46ed. Steven approved the reviewed design
and the task-local 240-added-source-line exception. This section supersedes
the historical four-term-assertion freeze only on the domain below. It grants
no HF, consumer restart, release or automatic successor. At most two in-scope
correction rounds; other scope, budget or scientific failures stop for review.

### Domain And Construction

Use one function-aware finalization for both non-injected terminal entrances:
`pqs_terminal_residual_gto_augmentation` and the supplement-taking
`pqs_terminal_residual_gto_augmented_hamiltonian`. Support their existing
PGDG primitive/stencil parent axes and Cartesian contracted supplements,
including s/p/d and higher nonnegative polynomial powers. Explicitly validate
that the evaluated axis representation is the actual parent representation.
Keep matrix-only, injected, protected and parent-backed composition paths
unchanged; no silently narrowed public signature or automatic rescue fallback.
If an ordinary supported terminal backend cannot use this representation,
stop and report the exact backend before committing; do not silently route it
to rounded arithmetic or add another evaluator.

Retain existing owner selection, ordering, cutoffs, signs and merge thresholds.
Separate selection from finalization without cloning it. A private call-local
finalizer is permitted, with the old matrix-only implementation as its default;
no public keyword, persistent callback or validation-bypass flag. Begin from
selected unnormalized T_G0/T_A0, never the old normalized/failing result.
Evaluate G*T_G0+A*T_A0 directly, preserving Gaussian tails and components
outside the finite parent. Use local terminal lifts and sparse accumulation
of contracted tensor coefficients without coefficient screening. Reproject
twice against the intended orthonormal G, recording exactly the same correction
P in T_G0. This is arithmetic cleanup, not a generalized metric solve.

For thin QR V=Q*R and small SVD R=L*Sigma*W', use the symmetric transform
U=W*diag(1/Sigma)*W'. Return (T_G0-P)*U and T_A0*U, then canonicalize signs.
Do not use inv(R) as the final gauge. Check Sigma^2 against the existing
absolute/relative merge gate; no floor, pruning, forced rank or changed
selection metadata. Preserve nonzero-D matrix-only behavior. Injected residuals
are orthogonal to their own fixed sector, so never reproject them against G
under this grant. Keep existing returned type/fields/orientation vocabulary.

### Finite Grid And Tail Freeze

Use positive Gauss-Legendre panels in each axis. For maximum axis angular
power l define T=max(12,sqrt(2*l+1)+10). Panel anchors for each parent primitive
or Gaussian key are center+t*sigma, t in (-T,-6,-2,0,2,6,T), with sigma equal
to its width or 1/sqrt(2*exponent). Include all centers, exponents and parent
primitive anchors, sorted uniquely; finite inputs and checked panel counts
are mandatory. No clipping to the parent box or screening small coefficients.
Existing evaluation/normalization conventions remain authoritative.

The finite schedule is exactly (order,tail)=(8,T),(12,T),(16,T),(16,T+4).
Construct once at (8,T); reevaluate the SAME returned coefficient functions
at all later grids, without renormalizing them or selecting the most favorable
result. Validate full RR and G-R on each grid. All must meet unchanged
scale-aware identity_atol and orthogonality_atol. Between consecutive grids
require RR max change <=identity_atol/10 and G-R max change
<=orthogonality_atol/10. Require RR-I and both cross orientations' induced
infinity row sums <=1e-5, and consecutive row-sum changes <=1e-6.
These are internal qualification checks, not replacements for the tighter
existing elementwise gates. No retry beyond this schedule or relaxed failure.

The final same-order tail enlargement isolates the tail check. Testing the
actual contracted, normalized returned functions accounts for amplification
and cancellation; raw-axis stabilization alone is insufficient. Angular-aware
anchors cover polynomial extent but are not universal certification. Independent
small analytic tests and saved-fixture selected analytic checks remain required.
If the bounded schedule cannot qualify an accepted input, stop with that input;
do not hard-code an angular cutoff, change its rank or start a grid framework.
Release construction tensors before validation, processing one grid at a time.

### Checked Memory Admission

Batch width b=min(8,nR), at least one. Before allocation compute all dimensions
and byte products with checked integer arithmetic, including panel/node counts
at the largest scheduled grid. Let N be the common axis-dimension product,
nP the parent-axis product, nK the Gaussian-key product, c_a axis column counts,
v_a largest grid node counts, and m the largest per-column intermediate across
all three tensor contractions/permutations (include N,nP,nK). Freeze the
conservative incremental estimate, in bytes:

```text
E = 512 MiB + 8*(4*N*nR + 12*b*m + 4*sum(v_a*c_a) +
                 6*nG*nR + 12*nR*nR + 2*nP*b + 2*nK*b)
```

Require E<=12 GiB and process high-water RSS at entry plus E<=16 GiB;
overflow/nonfinite dimensions or a failed bound throw before the large work.
This deliberately reserves prior temporaries and is not a guarantee about
unrelated allocations. Record measured RSS in qualification and stop above
16 GiB. No process scanning, public memory API or cache. Only one live
N-by-nR array, in-place QR and bounded batch buffers; no parent-by-final map,
dense parent metric, full tensor serialization or persistent representation.
Do not silently fall back or lower rank to fit memory. If this conservative
estimator excludes the saved fixture, report before allocating, not retune it.

### Operator And Test Acceptance

Operators use returned functions and full projection/normalization congruence.
Recompute moments and MWG through existing callers. Span preservation does
not justify rotating an old two-index IDA array. Preserve Galerkin nuclear
attraction, finite expansion and all operator source code. Test raw symmetry,
GR and RR overlap/kinetic/position/second moments and a small nuclear case
against independent analytic finite-expansion references. Keep existing
ordinary-fixture tolerances. For the selected cancellation fixture require
direct kinetic agreement <=1e-10 Ha and selected compact-assembly Ritz/occupied
perturbation <=1e-6 Ha including deltaH AND deltaS; neither is a full HF bound.
Invalid MWG widths or a failed operator test requiring operator-source changes
stop this task. No full-operator 1e-10 Ha or molecular-energy claim.

Test contracted translated tight/diffuse s/p/d and at least one existing
higher polynomial case, outside-parent span, rank failure, raw selected entry,
symmetric orientation, unchanged matrix-only/nonzero-D/injected behavior and
both terminal entrances. Physical validation replaces rounded cancellation
only on the repaired branch; preserve the latter as regression evidence, not
an obligatory rejecting gate or a suppressed exception. Report elementwise,
row-sum and symmetry measures separately.

### Exact Budget And Validation

Source additions including docstrings: residual_basis.jl<=30 and
pqs_terminal_residual_gto.jl<=210, total<=240; paths are the existing files
under src/cartesian/cartesian_residual_gaussians/ and
src/cartesian/cartesian_final_basis_realization/. No other source changes.
Existing tests: test/misc/runtests.jl<=35 and
test/driver_public/cartesian_residual_gto_mwg_system_runtests.jl<=100, total135.
Reader additions: docs/src/reference/export.md<=15 and
docs/src/manual/projected_q_shells.md<=20. No new files, APIs, result fields,
exports, caches, helpers outside those owners, dependency or workflow edits.
Delete/simplify replaced terminal logic rather than cloning selection.

One saved-H10 qualification across implementation/review/corrections:
60 cumulative numerical minutes,16 GiB peak RSS,8 GiB new local scratch.
Reuse the hash-frozen saved working basis and comparison artifacts in
/Users/srw/dmrgtmp/h10_supplement_design_20260920/. No basis rebuild or prototype
replay. Fresh raw selection must naturally retain270; check full span/refined
RR/cross and the four selected independent analytic/kinetic comparisons.
Charge failed work; separate construction, validation and I/O. No H10 full
operators/field/HF. Small nuclear/moment checks use small fixtures only.

Run existing misc, public Cartesian/collinear/residual, R3A/facade and relevant
injected/occupied-first owners, package load, docs_fast/docs, authority/self-test,
two deterministic renders, Documenter/log/diff and full three-job CI/Docs.
Use existing easy atom/diatomic fixtures for before/after time/allocation.
Stop for review if added construction time exceeds max(0.25*baseline,0.1s)
or peak RSS exceeds the above cap; no unqualified performance claim. No full
angular suite or repeated unrelated paper probes. Grant/closeout uses docs-only
checks and remote route; implementation uses normal full source CI.

Evidence: reviewed all270 report hash dd661c150e55f8e97079c9c70dc1715e441dc9eaa549b752c4dce1fd5bfdfcd1;
production design hash 61387b4d6756d5e5663ed65663e263892b60890614b15f82eeecebbde00064ef.
Construction11.409s, peak11.47GiB and selected operator limits are scratch
evidence, not timings of this symmetric production variant. Commit only after
the complete bounded acceptance succeeds. On a substantive failure preserve
the draft/evidence, do not broaden scope or manufacture a passing result.
Hchain-doer remains paused even after repository closeout.

## Cutoff History

The default sequence is historical evidence, not a menu of coequal policies:

```text
early owner-local path     occupation 1e-8; identity policy then current
HP-RG-IDTOL-FN-01          identity_atol 1e-8
HP-RG-CUTOFF-FN-01         occupation 5e-8; identity_atol 5e-8
HP-RG-CUTOFF-FN-02         occupation 1e-6; identity_atol remains 5e-8
Passes 642-643             standard occupation 1e-8 implemented; identity unchanged
```

The ORTHO pass addressed small final-identity overshoots with healthy spectra;
it did not change retained rank. The IDTOL pass admitted a Be high-zeta case
without dropping its real `6.151e-6` direction. CUTOFF-01 then explicitly
dropped a marginal Cr direction near `3.637e-8`. CUTOFF-02 generalized the
production cutoff to `1e-6` after residual-only Cr2 evidence found problematic
low-`H1_RR` sectors built from occupations around `1.27e-7` to `8.98e-7`.

The post-CUTOFF-02 audit dropped retained Cr2 counts from `68 + 68` to
`62 + 62` but did not eliminate every low residual-sector mode. That result is
historical measurement evidence, not authority for spectral pruning, another
cutoff change, or a Cr2 production claim.

## Active Residual-Sector Spectral Measurement

`HP-RG-SPECTRAL-AUDIT-01` remains an approved measurement-only follow-up to
`HP-RG-CUTOFF-FN-02`. It owns no tracked source, test, artifact, driver, or
workflow surface.

The exact optional local probe is:

```text
tmp/work/rg_spectral_cutoff1e6_audit.jl
```

It may report durable text/TSV evidence under the existing external CR2 run
directory, but those outputs are evidence rather than repository authority.
The audit must report:

- retained residual counts by owner;
- low eigenvalues of `K_RR`;
- low eigenvalues of `H1_RR = K_RR + sum_A Z_A U_A_RR`;
- owner weights and residual-occupation composition for low or flagged modes;
- comparison with an available one-center Cr residual baseline;
- whether low modes are dominated by the smallest retained occupations or by
  otherwise healthy retained directions.

Validation is package load, `git diff --check`, a residual-only Cr baseline
when available, and the current residual-only Cr2 fixture. No full HF or new
Hamiltonian artifact belongs to this lane.

The audit must not change production source, committed tests or fixtures,
artifacts or provenance, drivers, MWG/IDA, residual selection or merging,
cutoffs or tolerances, or add automatic pruning, kinetic guards, `H1_RR`
guards, dense `Vee`, HF, or solver workflow. If existing construction seams
cannot reconstruct `K_RR` and `H1_RR` cheaply enough, stop and report the exact
missing reusable seam. Do not add source instrumentation under this ID.

## Validation

Tracked maintenance coverage is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

It checks the current `1e-6` value both in memory and in
`supplement_provenance/occupation_cutoff`, along with owner-local metadata,
orthogonality, identity, exact augmented operators, final interaction, and the
bounded H2 endpoint.

`test/misc/runtests.jl` checks a synthetic positive near-null mode and rejects
a materially negative residual metric under the explicit numerical-complete
`1e-10` policy. That is opt-in policy coverage, not evidence that production
uses `1e-10`.

Historical N2, Be, Cr, and residual-only Cr2 measurements remain in the
manager log. Except for the exact active spectral probe above, machine-local
paths and old endpoint scalars are evidence rather than normative contracts.

## Non-Goals

This policy does not authorize:

- global candidate selection or global pivoted-Cholesky residual selection;
- eigenvalue flooring or alternate near-singular merge recovery;
- automatic width/zeta filtering or width-based conditioning repair;
- kinetic or `H1_RR` spectral guards;
- residual localization or Gaussian-array enrichment;
- injection, protected replacement, additive references, EGOI, or screened
  Hartree changes;
- MWG/IDA, raw-block, artifact, driver, public API, solver, or Cr2 workflow
  changes.

An explicit caller may override an already supported numerical keyword for a
separately authorized experiment. Such an override does not establish a new
default or change this precedence.
