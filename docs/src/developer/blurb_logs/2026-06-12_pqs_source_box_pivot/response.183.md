Pass 183 complete.

Result status:
- The matching WL 419 gausslet-only result was newly computed.
- It was not found as a complete scalar baseline in the cited logs/artifacts.
- Existing logs did contain the supplemented Fig.8-style `447` reproduction
  and its decomposition into `419` gausslet functions plus `28` AHGBS residuals.

Ignored local probe files:
- `tmp/work/wl_he_419_gausslet_only_probe.jl`
- `tmp/work/wl_he_419_gausslet_only_probe_summary.txt`
- Confirmed ignored by `.gitignore:12:tmp/`.

Probe command:

```text
julia --project=. tmp/work/wl_he_419_gausslet_only_probe.jl
```

WL gausslet-only setup:
- source: `old_nested_fixed_block_wl_gausslet_only`
- constructor:
  `MappedUniformBasisSpec(:G10; count=11, mapping=white_lindsey_atomic_mapping(Z=2,d=0.3,tail_spacing=10.0), reference_spacing=1.0)`
- mapping: `AsinhMapping(a=0.38729833462074165, s=0.7745966692414834, tail_spacing=10.0)`
- `asinh_c = a * s = 0.3`
- physical endpoints: `(-5.892850307983052, 5.892850307983052)`
- `Z = 2.0`
- `q = n_s = 5`
- AHGBS residuals excluded: `true`
- residual count: `0`

WL shell inventory:
- core dimension: `125`
- shell layer count: `3`
- shell retained per layer: `(98, 98, 98)`
- shell retained total: `294`
- final / gausslet-only dimension: `419`
- overlap identity error: `5.221771508021077e-15`

WL scalar results:
- H1 lowest energy: `-1.991344469963435`
- H1 solve kind: `ordinary_symmetric`
- H1-orbital self-Coulomb: `1.2420473874925473`
- RHF converged: `true`
- RHF iterations: `20`
- RHF one-electron energy: `-3.871408908674227`
- RHF electron-electron energy: `1.0206054056564373`
- RHF total energy: `-2.85080350301779`
- RHF density trace: `0.9999999999999998`
- RHF electron count: `1.9999999999999996`
- RHF residual: `5.941364067396648e-11`
- RHF energy change: `3.6903813338540203e-13`
- probe elapsed time: `3.725833416` seconds

PQS-WL diagnostic deltas:
- PQS RHF total - WL RHF total:
  `0.0014163805874733981`
- PQS H1 - WL H1:
  `0.004662494788541416`
- PQS H1-J self-Coulomb - WL H1-orbital self-Coulomb:
  `-0.015884787180628912`
- PQS RHF two-body - WL RHF electron-electron:
  `-0.002427934232285267`

Supplemented WL 447 is separate:
- The old Fig.8-style supplemented WL result is not the same baseline as this
  gausslet-only 419 result.
- The supplemented result uses the same old nested gausslet sector plus AHGBS-9
  S-only residuals:
  - gausslet count: `419`
  - residual count: `28`
  - final dimension: `447`
- It should not be used as if it were the matching gausslet-only comparison.

Validation:
- `julia --project=. tmp/work/wl_he_419_gausslet_only_probe.jl`
  - passed and wrote
    `tmp/work/wl_he_419_gausslet_only_probe_summary.txt`
- `julia --project=. -e 'pqs_total=-2.8493871224303167; wl_total=-2.85080350301779; pqs_h1=-1.9866819751748936; wl_h1=-1.991344469963435; pqs_j=1.2261626003119184; wl_j=1.2420473874925473; pqs_two=1.018177471424152; wl_two=1.0206054056564373; println("delta_total=", pqs_total-wl_total); println("delta_h1=", pqs_h1-wl_h1); println("delta_h1_j=", pqs_j-wl_j); println("delta_rhf_two_body=", pqs_two-wl_two)'`
  - printed the diagnostic deltas above
- `git status --short --branch`
  - before writing response files: `## main...origin/main`
- `git diff --numstat -- src test`
  - empty
- `git check-ignore -v tmp/work/wl_he_419_gausslet_only_probe.jl tmp/work/wl_he_419_gausslet_only_probe_summary.txt`
  - both ignored by `.gitignore:12:tmp/`

Interpretation:
- This is a diagnostic comparison baseline, not a final physics endpoint.
- I did not tune PQS or WL.
- I did not run HFDMRG/DMRG/CR2.
- I did not edit tracked source/test/generator/artifacts.

Deletion/shrinkage report:
- deleted: none
- simplified: none; this was an ignored local baseline recovery/computation pass
- quarantined: the computed WL 419 baseline remains an ignored local probe and
  summary, not a tracked acceptance test
- not deleted because: the ignored probe/summary are useful local evidence for
  the pass-183 comparison hygiene scalar
- exact remaining caller/blocker: no tracked caller added; manager decision
  remains whether to promote a gausslet-only 419 WL/PQS comparison gate later

-- repo-doer@macmini
