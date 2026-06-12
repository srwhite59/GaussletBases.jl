Pass 054 complete: ran the larger q=7 complete core/shell PQS He RHF probe.

Source/test status:

- This pass stayed probe-only.
- No production source changed.
- No permanent test was added.
- Ignored probe artifacts:
  - `tmp/work/pqs_complete_core_shell_he_rhf_q7_probe.jl`
  - `tmp/work/pqs_complete_core_shell_he_rhf_q7_probe_summary.txt`

Actual dimensions:

```text
current_box                         = (1:9, 1:9, 1:9)
inner_box                           = (2:8, 2:8, 2:8)
raw_source_dims                     = (7, 7, 7)

expected support/core/shell         = 729 / 343 / 386
actual support/core/shell           = 729 / 343 / 386
expected shell retained/final dim   = 218 / 561
actual shell retained/final dim     = 218 / 561
final overlap identity error        = 2.8221522341276284e-13
```

H1/J diagnostics:

```text
Z=1 H1 lowest energy                = -0.4935104443467795
Z=1 H1 error vs exact H             = +0.006489555653220513

Z=2 H1 lowest energy                = -1.9599740171825744
Z=2 H1 error vs hydrogenic He+      = +0.04002598281742564
Z=2 H1 improvement vs pass 053      = -0.11231209467644598

Z=2 H1 self-Coulomb J               = 1.24311548900175
J change vs pass 053                = +0.21771519405777928
```

RHF result:

```text
converged                           = true
iterations                          = 8
RHF one-electron energy             = -3.7943561135826886
RHF electron-electron energy         = +0.9842880634482857
RHF total energy                    = -2.810068050134403
He HF reference                     = -2.861679995612239
error vs He HF reference            = +0.051611945477835874
improvement vs pass 053 total       = -0.0887307672812363

density trace, final                = 1.0
electron count                      = 2.0
density trace, pre-final            = 0.9999999999999991
pre-final electron-count proxy      = 1.9999999999999982
Fock symmetry error                 = 0.0
```

Timing:

```text
final-basis/support build            = 2.365947 s
one-electron build                   = 0.303753 s
density-interaction build            = 0.332236 s
RHF solve                            = 1.100219 s
total probe phases                   = 4.102156 s
```

Comparison to pass 053 compact fixture:

```text
pass 053 final dimension             = 223
pass 053 Z=2 H1                      = -1.8476619225061284
pass 053 H1 self-Coulomb J           = 1.0254002949439707
pass 053 RHF total                   = -2.7213372828531668

q=7 final dimension                  = 561
q=7 Z=2 H1                           = -1.9599740171825744
q=7 H1 self-Coulomb J                = 1.24311548900175
q=7 RHF total                        = -2.810068050134403
```

Decision:

- The final/pre-final RHF convention still appears accepted.
- The q=7 fixture improves the Z=2 one-electron energy and RHF total substantially, and the electron count/density traces remain coherent.
- This is still probe-only. It is a better physics point than pass 053, but not yet a permanent smoke gate because the fixture remains exploratory and the He HF error is still about `0.0516 Ha`.
- Recommended next direction is another larger/parameterized probe or a fixture-quality review, not a convention audit.

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
julia --project=. tmp/work/pqs_complete_core_shell_he_rhf_q7_probe.jl
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
- No permanent test was added; the q=7 work remains an ignored `tmp/work` physics/performance probe.
- Old/fallback/oracle surfaces made less necessary: the successful q=7 run further reduces the need to look at signed final-weight density, raw no-division density, or fixed-block pair data as anything other than rejected/oracle references.
- Remaining stale/duplicate surface to retire next: compact-box-only probe assumptions should not become acceptance language. If this route is promoted, use q=7 or a better reviewed fixture rather than the pass-053 compact smoke result as the scientific baseline.

-- repo-doer@macmini
