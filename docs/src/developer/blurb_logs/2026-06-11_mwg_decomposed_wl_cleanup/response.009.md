Implemented the H2 R = 4 old nested/QW restricted HF diagnostic as a developer-only probe.

Artifacts:

- `tmp/work/h2_r4_qw_hf_probe.jl`
- `tmp/work/h2_r4_qw_hf_summary.txt`
- `tmp/work/h2_r4_qw_hf.tsv`

Scope:

- old nested/QW route only
- `interaction_treatment = :ggt_nearest`
- H/cc-pVTZ molecular Gaussian supplement with `lmax = 1`
- restricted closed-shell HF only
- no ED
- no acceptance test
- no production code changes
- no decomposed-route changes, raw GTO density-density final operator, full-parent CPB fallback, direct Cartesian fallback, ordinary Cartesian IDA fallback, PQS, or generalized final-basis solve

Routes run:

| route | fixed block size | final dim | residuals | HF total | documented HF total | delta |
|---|---:|---:|---:|---:|---:|---:|
| default complete rectangular | `(1215, 463)` | 481 | 18 | -0.9109382643524664 | -0.910938264352 | -4.664046926450283e-13 |
| endcap/panel q=4,L=4 | `(1215, 443)` | 461 | 18 | -0.9109773150033322 | -0.910977315003 | -3.3217872896784684e-13 |

Fixture:

- H2, `R = 4.0` bohr
- nuclear coordinates: `((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))`
- `core_spacing = 0.5`
- `xmax_parallel = 6.0`
- `xmax_transverse = 4.0`
- `nside = 5`

Default complete-rectangular route:

- final overlap identity error: `2.435940338330056e-12`
- fixed overlap identity error: `1.3816121305989328e-14`
- residual owner set: `(1, 2)`
- residual width status: `all_nan`
- shared layer types: `(:_CartesianNestedCompleteShell3D, :_CartesianNestedCompleteShell3D)`
- shared layer columns: `(98, 114)`
- RHF one-electron energy: `-1.5634754903175414`
- RHF electron-electron energy: `0.40253722596507513`
- nuclear repulsion: `0.25`
- RHF iterations: `9`

Endcap/panel q=4,L=4 route:

- final overlap identity error: `1.0308420783644578e-12`
- fixed overlap identity error: `1.3816121305989328e-14`
- residual owner set: `(1, 2)`
- residual width status: `all_nan`
- shared layer types: `(:_CartesianNestedEndcapPanelShellLayer3D, :_CartesianNestedEndcapPanelShellLayer3D)`
- shared layer columns: `(96, 96)`
- RHF one-electron energy: `-1.5632903437804713`
- RHF electron-electron energy: `0.40231302877713915`
- nuclear repulsion: `0.25`
- RHF iterations: `9`

Timing:

| route | fixed/source | operators | HF solve | route total |
|---|---:|---:|---:|---:|
| default complete rectangular | 10.239579 s | 14.192019167 s | 0.370291542 s | 24.801889709 s |
| endcap/panel q=4,L=4 | 0.707715667 s | 8.027014334 s | 0.3231185 s | 9.057848501 s |

Decision:

- H2 R = 4 HF is a clean next physics target.
- The default complete-rectangular route reproduces the documented HF total to `4.7e-13 Ha`.
- The endcap/panel q=4,L=4 route also reproduces the documented HF total to `3.3e-13 Ha`.
- The probe confirms the expected dimensions, residual counts, residual ownership, and finite symmetric operator path.

Validation:

- Ran `julia --project=. tmp/work/h2_r4_qw_hf_probe.jl`
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Ran `git diff --check`

Deletion/shrinkage report:

- No production code, tests, docs, metadata, or compatibility paths became obsolete from this probe-only pass.
- Nothing was deleted or simplified because the task only established a developer diagnostic against existing documented rows.
- No new tracked test was added. The probe remains under `tmp/work` and does not replace or expand routine coverage.
- Remaining stale/duplicate surfaces to retire next: none identified in this pass. A future pass can decide whether to make a compact H2 HF comparison a documented benchmark target, but it should not become a routine long-running acceptance gate by default.

-- repo-doer@macmini
