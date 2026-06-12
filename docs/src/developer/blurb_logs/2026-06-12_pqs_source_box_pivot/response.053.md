Pass 053 complete: ran the first complete core/shell PQS He RHF probe through the reviewed pre-final density-interaction seam.

Source/test status:

- This pass stayed probe-only.
- No production source changed in pass 053.
- No permanent test was added.
- Ignored probe artifacts:
  - `tmp/work/pqs_complete_core_shell_he_rhf_probe.jl`
  - `tmp/work/pqs_complete_core_shell_he_rhf_probe_summary.txt`

Final/pre-final Fock convention:

- Accepted for the exploratory probe.
- The occupied spatial orbital is solved in the orthonormal final basis.
- The electron-electron density response is evaluated in the pre-final localized positive-weight gauge using:

  ```text
  c_prefinal = combined_lowdin_cleanup * c_final
  rho_final  = c_final * c_final'
  rho_pre    = c_prefinal * c_prefinal'
  n_pre      = diag(rho_pre)
  G_pre      = 2 * diag(V_pre * n_pre) - rho_pre .* V_pre
  F_final    = H_final + L' * G_pre * L
  ```

- No dimension or convention conflict appeared in the probe.
- Signed final-weight density and raw no-division density remain rejected.

Probe result:

```text
status                              = materialized_exploratory_he_rhf
current_box                         = (1:7, 1:7, 1:7)
inner_box                           = (2:6, 2:6, 2:6)
final dimension                     = 223
support/core/shell support counts   = 343 / 125 / 218
shell final retained count          = 98
final overlap identity error        = 7.097873350991017e-14
nuclear charge                      = 2.0
```

One-electron and density diagnostics:

```text
Z=2 H1 lowest orbital energy        = -1.8476619225061284
error vs hydrogenic He+ -2.0        = 0.15233807749387163
H1 self-Coulomb J                   = 1.0254002949439707
density interaction status          = materialized_pqs_complete_core_shell_pre_final_density_interaction
density gauge                       = pre_final_localized_positive_weight
pre-final V shape                   = (223, 223)
pre-final V symmetry error          = 9.436895709313831e-16
```

RHF result:

```text
converged                           = true
iterations                          = 9
energy tolerance                    = 1e-11
RHF one-electron energy             = -3.6487559205762738
RHF electron-electron energy         = 0.9274186377231072
RHF total energy                    = -2.7213372828531668
He HF reference                     = -2.861679995612239
error vs He HF reference            = +0.14034271275907217
density trace, final                = 1.0000000000000009
electron count                      = 2.0000000000000018
density trace, pre-final            = 0.9999999999999978
pre-final electron-count proxy      = 1.9999999999999956
Fock symmetry error                 = 0.0
physically sane enough to consider promotion later = true
```

Interpretation:

- The RHF loop coherently consumes the final one-electron Hamiltonian and pre-final density interaction.
- The result is physically sane for this compact box: it converges, has positive electron-electron energy, the electron count is correct, and the total energy is negative.
- The energy is still high by about `0.14034 Ha` versus the He HF reference, so this should remain an exploratory physics probe rather than an acceptance test. The compact `1:7` / `2:6` PQS box is likely too coarse/small to judge final He quality.

Timing:

```text
final-basis/support build            = 1.949722 s
one-electron build                   = 0.209864 s
density-interaction build            = 0.240077 s
RHF solve                            = 0.434198 s
total probe phases                   = 2.833862 s
```

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
julia --project=. tmp/work/pqs_complete_core_shell_he_rhf_probe.jl
  PASS: probe completed and wrote summary artifact.

julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
  PASS: 28/28, 2.9s

julia --project=. -e 'using GaussletBases; println("load ok")'
  PASS: load ok

git diff --check
  PASS
```

Deletion/shrinkage report:

- No old production path became safe to delete from this probe-only pass.
- Nothing tracked was deleted or simplified.
- No new permanent test was added; the RHF work remains an ignored `tmp/work` physics probe.
- Old/fallback/oracle surfaces made less necessary: signed final-weight density and raw no-division density are further ruled out for RHF consumption; fixed-block pair data remains an oracle only.
- Remaining stale/duplicate surface to retire next: any probe or docs wording that frames the compact-box He RHF as accepted physics should be avoided or quarantined. If this path is promoted later, it should use a more meaningful fixture or explicitly mark the compact fixture as smoke/route validation only.

-- repo-doer@macmini
