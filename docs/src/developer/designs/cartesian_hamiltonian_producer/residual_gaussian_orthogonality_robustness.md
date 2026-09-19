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
| `HP-RG-ORTHO-FN-01` | Approved bounded premerge repair | Preserve implemented merge/rank/identity rules; change premerge evaluation only |
| `HP-RG-ORTHO-TEST-01` | Approved focused regression | Existing misc owner additions; unchanged H2/injected endpoint gates |
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

Pass 650 authorizes only the premerge repair for GB-H10-Q8-RESIDUAL-20260918,
under HP-RG-ORTHO-FN-01/TEST-01. Baseline is
`aa2ac41c513635bb5a7208848551eeef0edd754b`. This is a construction-arithmetic
repair, not a cutoff, rank, conditioning policy or physical-basis admission.

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
