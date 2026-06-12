Implemented the `n_s = 7` He RHF diagnostic as a developer-only probe.

Artifacts:

- `tmp/work/fig8_he_ns7_rhf_probe.jl`
- `tmp/work/fig8_he_ns7_rhf_summary.txt`
- `tmp/work/fig8_he_ns7_rhf.tsv`

Scope:

- old nested/QW MWG path only
- AHGBS-9 S-only loaded from `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets`
- no cc-pVTZ stand-in
- no acceptance test
- no production code changes
- no raw GTO density-density final operator, full-parent CPB fallback, direct Cartesian fallback, ordinary Cartesian IDA fallback, PQS, or generalized final-basis solve

Points run:

1. `n_s = 7`, `d = 0.15`
2. `n_s = 7`, `d = 0.10`
3. `n_s = 7`, `d = 0.20`

The first point was sub-mHa against the Fig. 8 row, so I ran the adjacent points.

Results:

| n_s | d | count | final dim | residuals | RHF total | error vs Fig. 8 | error vs He HF | elapsed |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 0.15 | 17 | 1461 | 28 | -2.8616625064826398 | +1.684323533845955e-5 | +1.7489129599201902e-5 | 19.37 s |
| 7 | 0.10 | 21 | 1897 | 28 | -2.861673961528321 | +1.716095156645281e-6 | +6.034083917860755e-6 | 21.26 s |
| 7 | 0.20 | 15 | 1243 | 28 | -2.861639584273142 | +4.379386004993435e-5 | +4.041133909682415e-5 | 9.33 s |

Additional diagnostics:

- `d = 0.15`: H1 `-1.9999998638028336`, J `1.2500086071811478`, overlap identity error `5.37e-10`, RHF iterations `18`
- `d = 0.10`: H1 `-1.9999998641649084`, J `1.2500024042126645`, overlap identity error `1.15e-9`, RHF iterations `18`
- `d = 0.20`: H1 `-1.9999998635987135`, J `1.2500221568175018`, overlap identity error `1.85e-10`, RHF iterations `18`

Constructor notes:

- `d = 0.15`: `MappedUniformBasisSpec(:G10; count=17, mapping=white_lindsey_atomic_mapping(Z=2,d=0.15,tail_spacing=10.0), reference_spacing=1.0)`, endpoints `+/-7.3284028397738386`
- `d = 0.10`: `count=21`, endpoints `+/-7.117581897720067`
- `d = 0.20`: `count=15`, endpoints `+/-7.983814222928366`

Decision:

- The `n_s = 7` family is accurate enough as the next comparison target.
- All three points are sub-mHa relative to the Fig. 8 rows.
- In this repo probe, `d = 0.10` is closest to the Fig. 8 table target and to the He HF reference.
- The probe does not reproduce the Fig. 8 local minimum ordering exactly: Fig. 8 has the smallest positive `n_s = 7` error at `d = 0.15`, while this compact repo RHF probe is closest at `d = 0.10`.
- That remaining trend difference should be treated as mapping/RHF-convention detail if exact paper reproduction is needed; it is not a blocker for using `n_s = 7` as the next accuracy comparison family.

Validation:

- Ran `julia --project=. tmp/work/fig8_he_ns7_rhf_probe.jl`
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Ran `git diff --check`

Deletion/shrinkage report:

- No production code, tests, docs, metadata, or compatibility paths became obsolete from this probe-only pass.
- Nothing was deleted or simplified because the work added only ignored `tmp/work` artifacts and a tracked response.
- No new tracked test was added. The probe does not replace or expand routine coverage.
- Remaining stale/duplicate surfaces to retire next: none identified. A future pass could decide whether to promote a compact `n_s = 7` comparison into a documented benchmark target, but not as a routine acceptance gate yet.

-- repo-doer@macmini
