Pass 055 complete: ran the q=9 complete core/shell PQS He RHF scaling probe.

Source/test status:

- This pass stayed probe-only.
- No production source changed.
- No permanent test was added.
- Ignored probe artifacts:
  - `tmp/work/pqs_complete_core_shell_he_rhf_q9_probe.jl`
  - `tmp/work/pqs_complete_core_shell_he_rhf_q9_probe_summary.txt`

Actual dimensions:

```text
current_box                         = (1:11, 1:11, 1:11)
inner_box                           = (2:10, 2:10, 2:10)
raw_source_dims                     = (9, 9, 9)

expected support/core/shell         = 1331 / 729 / 602
actual support/core/shell           = 1331 / 729 / 602
expected shell retained/final dim   = 386 / 1115
actual shell retained/final dim     = 386 / 1115
final overlap identity error        = 2.9398462830787508e-12
```

H1/J diagnostics:

```text
Z=1 H1 lowest energy                = -0.4967467553707099
Z=1 H1 error vs exact H             = +0.0032532446292901085

Z=2 H1 lowest energy                = -1.9869007472289404
Z=2 H1 error vs hydrogenic He+      = +0.013099252771059611
Z=2 H1 improvement vs q=5/pass053   = -0.13923882472281202
Z=2 H1 improvement vs q=7/pass054   = -0.026926730046366032

Z=2 H1 self-Coulomb J               = 1.226075273367442
J delta vs q=5/pass053              = +0.20067497842347137
J delta vs q=7/pass054              = -0.017040215634307909
```

RHF result:

```text
converged                           = true
iterations                          = 8
RHF one-electron energy             = -3.8680809329777786
RHF electron-electron energy         = +1.018171771175848
RHF total energy                    = -2.8499091618019303
He HF reference                     = -2.861679995612239
error vs He HF reference            = +0.011770833810308634
improvement vs q=5/pass053 total    = -0.12857187894876354
improvement vs q=7/pass054 total    = -0.03984111166752724

density trace, final                = 1.0
electron count                      = 2.0
density trace, pre-final            = 1.0000000000000089
pre-final electron-count proxy      = 2.0000000000000178
Fock symmetry error                 = 0.0
```

Timing:

```text
final-basis/support build            = 3.207434 s
one-electron build                   = 0.820304 s
density-interaction build            = 0.271407 s
RHF solve                            = 3.524653 s
total probe phases                   = 7.823799 s
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
```

Decision:

- The final/pre-final RHF convention still appears accepted.
- The q=9 fixture continues improving H1 and RHF total, with RHF now within about `0.01177 Ha` of the He HF reference.
- The q=9 timing is still practical as an exploratory probe at about `7.82 s` for the measured phases.
- I would keep this as probe-only for now. It is close enough to be a candidate for a future smoke/physics gate, but the next step should be a fixture/design review before making it permanent, especially to decide whether q=9 is the right bounded gate or whether a slightly different fixture is more scientifically stable.

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
julia --project=. tmp/work/pqs_complete_core_shell_he_rhf_q9_probe.jl
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
- No permanent test was added; the q=9 work remains an ignored `tmp/work` physics/performance probe.
- Old/fallback/oracle surfaces made less necessary: q=9 further confirms the final/pre-final density convention without signed final weights, raw no-division density, or fixed-block pair authority.
- Remaining stale/duplicate surface to retire next: compact q=5/q=7 probe artifacts should remain exploratory; if a permanent gate is added later, avoid preserving all ladder probes as tests.

-- repo-doer@macmini
