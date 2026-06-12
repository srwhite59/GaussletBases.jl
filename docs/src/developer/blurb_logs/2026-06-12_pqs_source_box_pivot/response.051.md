Pass 051 complete: PQS complete core/shell density gauge audit.

Files edited/artifacts:

- Added ignored probe artifact:
  - `tmp/work/pqs_complete_core_shell_density_gauge_probe.jl`
- Wrote ignored probe summary:
  - `tmp/work/pqs_complete_core_shell_density_gauge_probe_summary.txt`
- No tracked production source, test, or docs files were changed.
- No production helper was added.

Density gauge findings:

The complete route has a positive localized/pre-final density gauge.

```text
support weights:
  min = 0.39540709965461518
  max = 13.495587524797241
  sum = 1584.7117713178820
  positive / negative / near-zero = 343 / 0 / 0

pre_final_coefficients weights:
  min = 0.39540709965461518
  max = 13.501201088163116
  sum = 1274.0570096693957
  positive / negative / near-zero = 223 / 0 / 0
  all finite / all positive = true / true

final_coefficients weights:
  min = -18.967313490013488
  max = 58.136426512094616
  sum = -62.222631936233967
  positive / negative / near-zero = 111 / 112 / 0
  all finite / all positive = true / false
```

H1 and J table:

```text
H1 final route energy                         = -0.48047934800387126
fixed-block oracle H1                         = -0.48047920531279725

pre-final positive-weight J                   = 0.6397851751855723
fixed-block oracle J                          = 0.6397857768997106
exact H 1s J                                  = 0.625

pre-final J error vs fixed oracle             = -6.017141382930191e-7
pre-final |J - fixed oracle J|                = 6.017141382930191e-7
pre-final |J - exact 1s J|                    = 0.014785175185572319

final signed-weight divided candidate J       = 18.132786333403647
raw projected no-division candidate J         = 0.56242162015939745
stored normed pair-terms no-division J        = 0.05282601375665979
```

Matrix/gauge diagnostics:

```text
pre-final pair matrix shape                   = (223, 223)
pre-final pair matrix finite                  = true
pre-final pair matrix symmetry error          = 9.436895709313831e-16
support raw pair numerator symmetry error     = 0.0

fixed-current-gauge max delta vs final signed = 342413.19480934506
fixed-current-gauge max delta vs raw projected= 301.720543570848
fixed-current-gauge max delta vs naive pre-final linear transform = 0.01806935094886919
```

Important convention result:

- The old positive-weight rule is reproduced in the complete route by using
  `pre_final_coefficients`.
- The pre-final positive-weight contraction:

  ```text
  weighted_coefficients = pre_final_coefficients ./ pre_final_weights
  pair_matrix = weighted_coefficients' * raw_pair_numerator * weighted_coefficients
  ```

  matches the same-geometry fixed-block J oracle to about `6e-7`.
- The signed final-gauge weight division remains rejected.
- Raw projected no-division remains rejected as an accepted convention even though it is less bad than signed final division.
- A naive linear two-index transform of the pre-final density matrix into the final orthonormal coefficient gauge is not accepted. The probe records that candidate only as a diagnostic; its H1 self-Coulomb value is `0.042321243942694396`, so it is not the reviewed consumption rule.

Recommended convention/blocker:

- Recommended convention: keep density-density construction in the localized/pre-final positive-weight density gauge.
- Remaining blocker before RHF: `:missing_reviewed_pre_final_density_interaction_consumption_for_rhf`.
- The next implementation should introduce a separate density-interaction object or explicitly reviewed consumption rule for final-basis orbital/RHF coefficients. It should not use signed final weights and should not adopt raw no-division as a shortcut.

Validation:

```text
julia --project=. tmp/work/pqs_complete_core_shell_density_gauge_probe.jl
  PASS: probe completed and wrote summary artifact.

julia --project=. -e 'using GaussletBases; println("load ok")'
  PASS: load ok

git diff --check
  PASS
```

I did not run `test/nested/pqs_direct_retained_final_h1_runtests.jl` because no tracked source changed.

Deletion/shrinkage report:

- No old production code became safe to delete in this probe-only pass; the pass resolves the gauge direction but intentionally stops before adding the RHF consumption object.
- Nothing tracked was deleted or simplified.
- No new permanent test was added. The new probe is an ignored `tmp/work` artifact for blocker discovery, not long-term coverage.
- Old/fallback/oracle surface made less authoritative: signed final-gauge IDA weight division is now clearly ruled out; raw projected no-division remains diagnostic only; fixed-block pair data remains oracle-only and was not used as active authority.
- Remaining stale/duplicate surfaces to retire next: any wording or helper that implies final signed weights own PQS density-density should be deleted or quarantined once the pre-final density-interaction object lands.

-- repo-doer@macmini
