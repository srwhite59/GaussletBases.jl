Pass 050 response: PQS complete core-shell final density J diagnostic

Summary:
- Created and ran the ignored density/J probe.
- Did not retain a production density-density helper because the intended final signed-weight division gives an unphysical H1 self-Coulomb value for the current complete PQS final basis.
- No tracked source, docs, or tests were changed in this pass.

Files edited:
- Ignored artifacts only:
  - `tmp/work/pqs_complete_core_shell_final_density_j_probe.jl`
  - `tmp/work/pqs_complete_core_shell_final_density_j_probe_summary.txt`

Helper names:
- Used existing helper:
  - `pqs_complete_core_shell_final_ida_weights`
- No new source helper was retained.
- The would-be helper `pqs_complete_core_shell_final_density_density_matrix` remains blocked by:
  - `:final_signed_ida_weight_division_gives_unphysical_j`

Probe route:
- Complete final basis:
  - `current_box = (1:7, 1:7, 1:7)`
  - `inner_box = (2:6, 2:6, 2:6)`
  - `raw_source_dims = (5, 5, 5)`
  - final dimension: `223`
- H1:
  - `-0.48047934800387126`
- Final IDA weights:
  - min: `-1.8967313490013488e+01`
  - max: `5.8136426512094616e+01`
  - sum: `-6.2222631936233967e+01`
  - positive count: `111`
  - negative count: `112`
  - near-zero count: `0`

Density/J diagnostic:
- Raw support numerator source:
  - `pgdg_intermediate.pair_factor_terms_raw`
- Candidate convention tested:
  - project raw numerator through complete final coefficients;
  - divide by `final_ida_weights[i] * final_ida_weights[j]` after projection.
- Candidate matrix checks:
  - shape: `(223, 223)`
  - finite: `true`
  - symmetry error: `1.2732925824820995e-11`
  - support raw numerator symmetry error: `0.0`
- H1 self-Coulomb contraction:
  - `:restricted_one_orbital_direct_minus_exchange`
  - same scalar convention used by the decomposed WL He H1/J diagnostic.
- Result:
  - J with raw numerator plus final signed-weight division: `18.132786333403647`
  - exact H 1s reference: `0.625`
  - absolute error: `17.507786333403647`
- This is not physically sane, so the final density-density route remains blocked.

Comparison variants:
- Raw numerator projected into the final basis with no final-weight division:
  - J = `0.56242162015939745`
  - absolute error vs `0.625` = `0.06257837984060255`
- Stored density-normalized `pair_factor_terms` projected with no final-weight division:
  - J = `0.05282601375665979`
- Same-geometry fixed-block oracle:
  - fixed H1 = `-0.48047920531279725`
  - fixed J = `0.6397857768997106`
  - fixed/current final-subspace identity error = `5.8875461168328425e-04`
  - fixed current-gauge density matrix max delta vs signed-weight-divided candidate = `3.4241319480934506e+05`
- Fixed-block pair data was used only as oracle/comparison, not route authority.

Interpretation/blocker:
- The current complete PQS final basis has signed final IDA weights. Blindly dividing the projected raw numerator by the signed final-weight outer product amplifies the interaction matrix and gives an unphysical J.
- The raw projected numerator without final-weight division is much closer to both the exact H reference and the fixed-block oracle, but adopting that would contradict the stated density convention without a reviewed rule.
- Exact blocker:
  - `:final_signed_ida_weight_division_gives_unphysical_j`
- Recommended next audit:
  - determine whether complete core/shell final density weights should be positive density-basis weights distinct from signed final function integrals, or whether the complete final-basis density interaction should use the raw projected numerator without the final-function integral division.

Validation:
- `julia --project=. tmp/work/pqs_complete_core_shell_final_density_j_probe.jl`
  - passed as a blocked diagnostic
  - note: this process precompiled `GaussletBases` after an intermediate source edit that was not retained.
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 28 checks in 2.9s
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
- `git diff --check`
  - passed
- `git status --short --branch`
  - clean before writing this response.

Deletion/shrinkage report:
- No old production code, tests, metadata, or compatibility path became obsolete because the density convention is still blocked.
- No tracked code was deleted or simplified.
- No permanent test was added. The new coverage is an ignored diagnostic probe only, because the candidate density-density route is not an accepted live contract.
- Remaining stale/duplicate surfaces to retire next:
  - any notes implying final signed function-integral weights are already validated as density weights for complete PQS density-density;
  - any future helper that uses fixed-block `pair_sum` as active authority instead of an oracle;
  - any probe that silently accepts raw projected/no-division J without first resolving the final density-weight convention.

-- repo-doer@macmini
