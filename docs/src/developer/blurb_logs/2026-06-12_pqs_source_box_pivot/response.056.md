Pass 056 complete: ran the q=11 complete core/shell PQS He RHF scaling probe.

Source/test status:

- This pass stayed probe-only.
- No production source changed.
- No permanent test was added.
- Ignored probe artifacts:
  - `tmp/work/pqs_complete_core_shell_he_rhf_q11_probe.jl`
  - `tmp/work/pqs_complete_core_shell_he_rhf_q11_probe_summary.txt`

Actual dimensions:

```text
current_box                         = (1:13, 1:13, 1:13)
inner_box                           = (2:12, 2:12, 2:12)
raw_source_dims                     = (11, 11, 11)

expected support/core/shell         = 2197 / 1331 / 866
actual support/core/shell           = 2197 / 1331 / 866
expected shell retained/final dim   = 602 / 1933
actual shell retained/final dim     = 602 / 1933
final overlap identity error        = 8.427841757807641e-13
```

H1/J diagnostics:

```text
Z=1 H1 lowest energy                = -0.49927641075301576
Z=1 H1 error vs exact H             = +0.0007235892469842353

Z=2 H1 lowest energy                = -1.9949555260353655
Z=2 H1 error vs hydrogenic He+      = +0.005044473964634477
Z=2 H1 improvement vs q=5/pass053   = -0.14729360352923715
Z=2 H1 improvement vs q=7/pass054   = -0.034981508852791166
Z=2 H1 improvement vs q=9/pass055   = -0.008054778806425134

Z=2 H1 self-Coulomb J               = 1.2455645980701733
J delta vs q=5/pass053              = +0.22016430312620261
J delta vs q=7/pass054              = +0.0024491090684233274
J delta vs q=9/pass055              = +0.019489324702731237
```

RHF result:

```text
converged                           = true
iterations                          = 8
RHF one-electron energy             = -3.87801660119383
RHF electron-electron energy         = +1.0220690807649275
RHF total energy                    = -2.8559475204289022
He HF reference                     = -2.861679995612239
error vs He HF reference            = +0.005732475183336749
improvement vs q=5/pass053 total    = -0.13461023757573543
improvement vs q=7/pass054 total    = -0.045879470294499125
improvement vs q=9/pass055 total    = -0.006038358626971885

density trace, final                = 0.9999999999999989
electron count                      = 1.9999999999999978
density trace, pre-final            = 0.9999999999999936
pre-final electron-count proxy      = 1.9999999999999871
Fock symmetry error                 = 0.0
```

Timing:

```text
final-basis/support build            = 5.036738 s
one-electron build                   = 2.005209 s
density-interaction build            = 0.477374 s
RHF solve                            = 10.372039 s
total probe phases                   = 17.891361 s
```

Comparison ladder:

```text
q=5/pass053 final dimension          = 223
q=5/pass053 Z=2 H1                   = -1.8476619225061284
q=5/pass053 J                        = 1.0254002949439707
q=5/pass053 RHF total                = -2.7213372828531668

q=7/pass054 final dimension          = 561
q=7/pass054 Z=2 H1                   = -1.9599740171825744
q=7/pass054 J                        = 1.24311548900175
q=7/pass054 RHF total                = -2.810068050134403

q=9/pass055 final dimension          = 1115
q=9/pass055 Z=2 H1                   = -1.9869007472289404
q=9/pass055 J                        = 1.226075273367442
q=9/pass055 RHF total                = -2.8499091618019303

q=11/pass056 final dimension         = 1933
q=11/pass056 Z=2 H1                  = -1.9949555260353655
q=11/pass056 J                       = 1.2455645980701733
q=11/pass056 RHF total               = -2.8559475204289022
```

Decision:

- The final/pre-final RHF convention still appears accepted.
- q=11 gives the best physics result in the ladder: RHF is within about `0.00573 Ha` of the He HF reference, and Z=2 H1 is within about `0.00504 Ha` of `-2.0`.
- q=11 is still practical as a probe, but it is materially slower than q=9: measured phases are `17.89 s` vs q=9's `7.82 s`, mostly from RHF solve.
- Recommendation: q=9 is the better future compact gate candidate unless the gate is allowed to be around 20 seconds. q=11 should remain the stronger exploratory/reference probe. A manager fixture review should decide whether to promote q=9 or q=11.

Nonclaims:

- `signed_final_weight_density_used = false`
- `raw_no_division_density_used = false`
- `fixed_block_pair_data_authority_used = false`
- `gto_materialized = false`
- `driver_route_materialized = false`
- `exports_materialized = false`
- `artifacts_materialized = false`

Validation:

```text
julia --project=. tmp/work/pqs_complete_core_shell_he_rhf_q11_probe.jl
  PASS: probe completed and wrote summary artifact.

julia --project=. -e 'using GaussletBases; println("load ok")'
  PASS: load ok

git diff --check
  PASS
```

I did not run `test/nested/pqs_direct_retained_final_h1_runtests.jl` in this pass because no source changed.

Deletion/shrinkage report:

- No old production path became safe to delete from this probe-only scaling pass.
- Nothing tracked was deleted or simplified.
- No permanent test was added; the q=11 work remains an ignored `tmp/work` physics/performance probe.
- Old/fallback/oracle surfaces made less necessary: q=11 reinforces that the accepted final/pre-final density convention scales without signed final weights, raw no-division density, or fixed-block pair authority.
- Remaining stale/duplicate surface to retire next: the ladder probes should not all become permanent tests. If promoting a gate, choose one reviewed fixture and leave the rest as historical `tmp/work` artifacts.

-- repo-doer@macmini
