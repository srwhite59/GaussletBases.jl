Implemented the Fig. 8 He RHF reproduction audit as a developer-only probe.

Artifacts:

- `tmp/work/fig8_he_rhf_target_reproduction_probe.jl`
- `tmp/work/fig8_he_rhf_target_reproduction_summary.txt`

Archive target:

- data table: `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets_fig8_heuhf_data.tsv`
- provenance note: `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets_fig8_heuhf_provenance.md`
- target row: `n_s = 5`, `d = 0.3`
- target energy: `-2.861543784624258`
- target signed error: `+1.362109879808e-04 Ha`

AHGBS-9 availability:

- The repo-vendored `data/legacy/BasisSets` does not contain `He AHGBS-9`.
- `/Users/srw/BasisSets` is absent on this machine.
- Your pointer was correct: `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets` contains `He AHGBS-9`.
- `/Users/srw/Library/CloudStorage/Dropbox/oldgaussletmodules/BasisSets` also contains `He AHGBS-9`.
- Loading through `legacy_atomic_gaussian_supplement("He", "AHGBS-9"; lmax = 0, basisfile = "/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets")` works.
- Loaded AHGBS-9 S-only shell/primitive count: `28 / 28`.

Constructor and mapping used:

- `MappedUniformBasisSpec(:G10; count = 11, mapping = white_lindsey_atomic_mapping(Z = 2, d = 0.3, tail_spacing = 10.0), reference_spacing = 1.0)`
- mapping: `AsinhMapping(a=0.38729833462074165, s=0.7745966692414834, tail_spacing=10.0)`
- `asinh_c = a * s = 0.3`
- physical endpoints from this repo constructor: `(-5.892850307983052, 5.892850307983052)`
- old nested fixed-block dimension before supplement: `419`

Legacy log convention mismatch recorded:

- Paper caption says `R_b = 7 bohr`.
- Legacy logs show `xmin = -8.0`, `xmax = 8.0`, `basradius = 8.0`.
- The detailed `He.5.3` log reports `doInvsqrt = true`, `rangeg = -8:8`, `nlet = 17`, and backbone coordinates about `+/-5.470267`.
- The current repo constructor above gives the matching basis count structure but not identical physical endpoints.

Old nested/QW MWG reproduction result:

- final dimension: `447`
- gausslet count: `419`
- residual count: `28`
- H1 lowest energy: `-1.9999998632985623`
- IDA self-Coulomb J: `1.2496940228276845`
- RHF one-electron energy: `-3.887511967391441`
- RHF electron-electron energy: `1.0254098228577182`
- RHF total: `-2.862102144533723`
- error relative to Fig. 8 `n_s = 5`, `d = 0.3` target: `-0.000558359909464734 Ha`
- RHF converged in `20` iterations.

Conclusion:

- The repo can now load AHGBS-9 S-only from the external GaussletModules `BasisSets` file and can build a 447-function old nested/QW MWG fixture.
- It does not currently reproduce the Fig. 8 `n_s = 5`, `d = 0.3` energy with this probe.
- The final dimension matches the Fig. 8/log count, so the earlier cc-pVTZ residual-count issue is resolved.
- The remaining mismatch is likely in one or both of:
  - exact legacy mapping/box/grid convention, since current endpoints differ from the detailed log;
  - RHF/update/energy convention between the compact repo probe and the original legacy run.
- This is not evidence for a decomposed-route-local bug; this pass used old nested/QW MWG as the comparison path.

Validation:

- Ran `julia --project=. tmp/work/fig8_he_rhf_target_reproduction_probe.jl`
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Ran `git diff --check`

Deletion/shrinkage report:

- No production code, tests, docs, metadata, or compatibility paths became obsolete from this audit-only pass.
- Nothing was deleted or simplified because the work was a source/mapping/reproduction audit and did not replace an active route.
- No new tracked test was added. The probe is developer-only under `tmp/work`.
- Remaining stale/duplicate surface to retire next: none identified. The next useful pass is a targeted legacy-mapping/RHF-convention audit if exact Fig. 8 reproduction is still required.

-- repo-doer@macmini
