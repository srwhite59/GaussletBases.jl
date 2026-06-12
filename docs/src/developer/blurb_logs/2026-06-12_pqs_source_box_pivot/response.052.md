Pass 052 complete: added the complete PQS pre-final density-interaction seam and accepted the final-to-pre-final orbital consumption rule for the H1 self-Coulomb diagnostic.

Files edited:

- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  - added `pqs_complete_core_shell_pre_final_density_interaction`
  - added `pqs_complete_core_shell_pre_final_orbital_self_coulomb`
  - added private helper `_pqs_complete_core_shell_restricted_one_orbital_interaction_energy`
- `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - updated module boundary comments/docstring to include the narrow pre-final density-interaction seam.

Ignored probe artifacts:

- `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl`
- `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe_summary.txt`

Helper behavior:

- `pqs_complete_core_shell_pre_final_density_interaction(...)` builds the density interaction in the localized pre-final positive-weight gauge from:

  ```text
  pre_final_weights = pre_final_coefficients' * support_weights
  weighted_coefficients = pre_final_coefficients ./ pre_final_weights
  pair_matrix = weighted_coefficients' * raw_pair_numerator * weighted_coefficients
  ```

- It rejects/non-materializes the density matrix if pre-final weights are not finite and positive.
- It records `final_to_pre_final_coefficients = combined_lowdin_cleanup`.
- It does not use signed final IDA weights, raw no-division density, or fixed-block pair data as authority.
- `pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)` consumes a final-basis orbital by mapping:

  ```text
  c_prefinal = combined_lowdin_cleanup * c_final
  ```

  and then evaluates the restricted one-orbital direct-minus-exchange diagnostic in the pre-final density gauge.

Probe result:

```text
status                                  = accepted_pre_final_density_consumption_rule
final_dimension                         = 223
H1 energy                               = -0.48047934800387126

pre-final weight min/max/sum            = 0.39540709965461518 / 13.501201088163116 / 1274.0570096693957
pre-final positive/negative/near-zero   = 223 / 0 / 0
pre-final weights all positive          = true

pre-final pair matrix shape             = (223, 223)
pre-final pair matrix finite            = true
pre-final pair matrix symmetry error    = 9.436895709313831e-16

final-to-pre-final reconstruction error = 0.0
recorded reconstruction error           = 0.0
coefficient map                         = combined_lowdin_cleanup_times_final_coefficients
```

J diagnostic:

```text
H1 self-Coulomb J                       = 0.6397851751855723
fixed-block oracle J                    = 0.6397857768997106
delta vs fixed oracle                   = -6.017141382930191e-7
abs delta vs fixed oracle               = 6.017141382930191e-7
exact H 1s J reference                  = 0.625
abs error vs exact 1s                   = 0.014785175185572319

rejected signed final-weight J          = 18.132786333403647
rejected raw projected no-division J    = 0.56242162015939745
```

Final-to-pre-final consumption rule:

- Accepted for the H1 self-Coulomb diagnostic.
- The accepted rule is `c_prefinal = combined_lowdin_cleanup * c_final`.
- RHF is still not materialized in this pass. The new seam gives RHF a density-interaction object to consume next, but does not run or wire RHF.

Validation:

```text
julia --project=. tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl
  PASS: probe completed and wrote summary artifact.

julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
  PASS: 28/28, 2.9s

julia --project=. -e 'using GaussletBases; println("load ok")'
  PASS: load ok

git diff --check
  PASS
```

Deletion/shrinkage report:

- The signed final-gauge density-weight path is now less necessary and should remain rejected/quarantined; this pass did not delete it because it was only present as probe/diagnostic vocabulary, not as a tracked production helper.
- Raw projected no-division density remains diagnostic-only and was not promoted.
- Fixed-block pair data remains oracle-only; the production helper builds from route-owned complete support coefficients and raw pair numerator input.
- No tests were added. The focused `tmp/work` probe protects the new convention during this pass without adding permanent carrying cost.
- Nothing tracked was deleted. The pass added a new live seam needed by the next RHF step, so deletion would be premature.
- Remaining stale/duplicate surface to retire next: any probe/docs wording that implies final signed weights or raw no-division are viable density-density contracts should be removed once the manager accepts this seam.

-- repo-doer@macmini
