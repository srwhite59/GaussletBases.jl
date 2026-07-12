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
| `HP-RG-ORTHO-FN-01` | Implemented current robustness | Symmetric final merge, hard near-singular failure, scale-aware identity check |
| `HP-RG-ORTHO-TEST-01` | Completed evidence; active maintenance | H2 endpoint plus historical N2 robustness evidence |
| `HP-RG-IDTOL-FN-01` | Implemented historical; superseded | Former `identity_atol = 1e-8` default |
| `HP-RG-IDTOL-TEST-01` | Completed historical evidence | Be high-zeta identity-tolerance audit |
| `HP-RG-CUTOFF-FN-01` | Implemented historical; cutoff superseded | Former `5e-8` cutoff; current `identity_atol = 5e-8` originated here |
| `HP-RG-CUTOFF-TEST-01` | Completed historical evidence | Cr marginal-direction and H2 assertion validation |
| `HP-RG-CUTOFF-FN-02` | Implemented current production policy | `residual_occupation_cutoff = 1e-6` |
| `HP-RG-CUTOFF-TEST-02` | Completed evidence; active maintenance | H2 cutoff/provenance assertions and residual-only Cr2 audit |

Implementation commits are `76453cc39` for robust final identity,
`47e56593c` for the historical `1e-8` identity default, `f0b662dca` for the
historical `5e-8` cutoff/current identity value, and `1f7f04e56` for the
current `1e-6` cutoff. Manager-log Passes 131-132, 155-156, 171-172, and
183-185 preserve the evidence and policy transitions.

The IDTOL and CUTOFF-01 IDs no longer grant source or test work. They remain
addressable historical records so old artifacts, reports, and Git history can
be interpreted without treating their defaults as current.

## Current Defaults

The ordinary production builder and terminal compatibility entry point use:

```text
residual_occupation_cutoff = 1e-6
tau_neg_abs                = 1e-12
tau_neg_rel                = 1e-12
tau_merge_abs              = 1e-12
tau_merge_rel              = 1e-12
orthogonality_atol         = 1e-10
identity_atol              = 5e-8
```

`HP-RG-CUTOFF-FN-02` owns the current production occupation cutoff. It changed
neither `identity_atol` nor the negative-metric, merge, or orthogonality
policies.

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

## Cutoff History

The default sequence is historical evidence, not a menu of coequal policies:

```text
early owner-local path     occupation 1e-8; identity policy then current
HP-RG-IDTOL-FN-01          identity_atol 1e-8
HP-RG-CUTOFF-FN-01         occupation 5e-8; identity_atol 5e-8
HP-RG-CUTOFF-FN-02         occupation 1e-6; identity_atol remains 5e-8
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
