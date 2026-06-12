Review result:

Accepted as a developer-only Fig. 8 reproduction audit. The pass found and
loaded the correct AHGBS-9 S-only supplement from the external
GaussletModules `BasisSets` file and constructed a 447-function old nested/QW
MWG fixture for the `n_s = 5`, `d = 0.3` target.

Main result:

The fixture now matches the Fig. 8 basis-count structure:

- gausslet count: `419`
- residual S-only AHGBS-9 directions: `28`
- final dimension: `447`

It does not yet reproduce the Fig. 8 energy:

- Fig. 8 target: `-2.861543784624258`
- repo reproduction: `-2.862102144533723`
- difference: `-5.58359909464734e-4 Ha`

This is too large to call the reproduction successful, and it is in the same
qualitative direction as the earlier side13 cc-pVTZ diagnostic: below the
external HF/Fig. 8 target.

Important narrowing:

The mismatch is no longer due to using cc-pVTZ or the wrong residual count.
The response records two more likely sources:

- exact legacy mapping/box/grid convention, since the repo constructor gives
  endpoints about `+/-5.89285` while the detailed legacy log reports backbone
  coordinates about `+/-5.470267`;
- RHF/update/energy convention, since the legacy `He.5.3` log converges the
  reported `uhfen` in four energy iterations while the compact repo probe RHF
  helper takes twenty iterations.

Additional manager check:

The legacy `He.5.3` log also reports:

- `Ntot = 447`
- `norm(O1 - eye(Ntot)) = 4.847810690997331e-11`
- `uhfen = -2.861543784624258`
- `(doside, dwidth, jflat, wi, gscalefac, polylim, lmaxadd, restrictedHF) =
  (5, 10.0, 39//2, 6.0, sqrt(2), 9, 0, true)`
- `doInvsqrt = true`, `rangeg = -8:8`, `nlet = 17`

Validation reviewed:

- `julia --project=. tmp/work/fig8_he_rhf_target_reproduction_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/fig8_he_rhf_target_reproduction_summary.txt`

Deletion/shrinkage review:

No old surface became obsolete. No test was added. This was the right choice:
the reproduction is still not exact enough to become an acceptance gate.

Commit/push:

Pending manager commit/push of the tracked response/review log.

Next target:

Run a targeted legacy mapping/RHF convention audit for `He.5.3`, centered on
the 447-function AHGBS-9 fixture. The next pass should determine whether the
remaining `0.558 mHa` mismatch comes from mapping/box construction, the compact
probe RHF helper, or residual MWG/IDA energy assembly.
