# Atomic Core Compression vs Distortion Tradeoff

This pass compares two different follow-ons to the current trusted reduced
atomic anchor:

- direct tensor-style compression of the retained direct `5^3` core
- modest distortion / core-spacing adjustment with the trusted complete-shell
  source language left unchanged

## Scope

Atomic-only, on the active shared hybrid footing:

- `legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)`

The trusted reference point stays:

- corrected complete-shell hybrid at fixed dimension `589`

The outer complete-shell hierarchy is kept fixed for the direct core
compression variants.

## Variants Compared

Reference:

- current trusted corrected complete-shell anchor
- fixed dimension `589`

Direct core-compression variants:

- `5^3 -> 4^3`
- `5^3 -> 3^3`

using direct tensor-style core compression only:

- one `doside` contraction on the full `5`-point core interval in each axis
- then a tensor-product core basis of dimension `4^3 = 64` or `3^3 = 27`
- no inner shell construction

Modest distortion variants:

- keep the same complete-shell source language and same `lmax = 0` supplement
- compare `a = 0.20`, `0.25`, and `0.30`
- same `xmax = 10`
- same count-based He fixed-`a` family

So the distortion variants keep the trusted reduced fixed dimension `589`
while changing the near-origin/core-spacing choice modestly.

## Count-17 Primary Case

| Case | Fixed dim | Overlap error | `E1` | `⟨Vee⟩` | Projected fixed-only `⟨Vee⟩` shift | Ground capture | Avg first 4 capture | Weight min / max | Warmed time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Core `3^3` | `491` | `7.38e-12` | `-1.9709219606` | `1.2327605807` | `+2.72e-3` | `0.999431` | `0.996496` | `0.20548 / 3.25917` | `1.13 s` |
| Core `4^3` | `528` | `4.61e-12` | `-1.9709219606` | `1.2288892402` | `-1.79e-3` | `0.999431` | `0.996516` | `0.13043 / 3.25917` | `1.21 s` |
| Anchor `a = 0.25` | `589` | `4.04e-12` | `-1.9981842264` | `1.2489346081` | `+1.33e-4` | `0.999987` | `0.996655` | `0.06147 / 3.25917` | `1.89 s` |
| Anchor `a = 0.20` | `589` | `2.66e-12` | `-1.9984789266` | `1.2490214563` | `+1.22e-4` | `0.999990` | `0.996286` | `0.04763 / 3.19823` | `1.30 s` |
| Anchor `a = 0.30` | `589` | `2.22e-12` | `-1.9976857002` | `1.2486676292` | `+1.40e-4` | `0.999984` | `0.996906` | `0.07546 / 3.30558` | `1.56 s` |

Current trusted anchor geometry:

- direct core interval: `7:11`
- reference-center range: `(-2.0, 2.0)`
- physical-center range at `a = 0.25`: `(-0.3900208287, 0.3900208287)`

Core physical range under the modest distortion variants:

- `a = 0.20`: `(-0.3371644100, 0.3371644100)`
- `a = 0.25`: `(-0.3900208287, 0.3900208287)`
- `a = 0.30`: `(-0.4386856768, 0.4386856768)`

## Count-15 Robustness Case

| Case | Fixed dim | Overlap error | `E1` | `⟨Vee⟩` | Projected fixed-only `⟨Vee⟩` shift | Ground capture | Avg first 4 capture | Weight min / max | Warmed time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Core `3^3` | `491` | `1.29e-11` | `-1.9598828898` | `1.2297623630` | `+6.82e-3` | `0.998605` | `0.999092` | `0.31477 / 6.25532` | `0.68 s` |
| Core `4^3` | `528` | `6.75e-12` | `-1.9598828898` | `1.2196110549` | `-4.92e-3` | `0.998605` | `0.999169` | `0.17683 / 6.25532` | `0.73 s` |
| Anchor `a = 0.25` | `589` | `3.35e-12` | `-1.9970883635` | `1.2478983246` | `+6.69e-5` | `0.999996` | `0.999517` | `0.07750 / 6.25532` | `0.82 s` |
| Anchor `a = 0.20` | `589` | `9.89e-12` | `-1.9976624000` | `1.2475766734` | `+6.17e-5` | `0.999996` | `0.999586` | `0.06011 / 6.27886` | `0.81 s` |
| Anchor `a = 0.30` | `589` | `2.04e-12` | `-1.9965569658` | `1.2475239820` | `+6.96e-5` | `0.999995` | `0.999453` | `0.09505 / 6.22813` | `0.81 s` |

Current trusted anchor geometry:

- direct core interval: `6:10`
- reference-center range: `(-2.0, 2.0)`
- physical-center range at `a = 0.25`: `(-0.4896496380, 0.4896496380)`

Core physical range under the modest distortion variants:

- `a = 0.20`: `(-0.4269484801, 0.4269484801)`
- `a = 0.25`: `(-0.4896496380, 0.4896496380)`
- `a = 0.30`: `(-0.5469824910, 0.5469824910)`

## Interpretation

### Direct core compression does not help

Neither direct core-compression point is useful.

At both counts:

- `5^3 -> 4^3` badly degrades `E1`
- `5^3 -> 4^3` badly degrades `⟨Vee⟩`
- `5^3 -> 3^3` is also bad
- the projected fixed-only interaction transfer gets much worse
- the reduction is only
  - `589 -> 528`
  - `589 -> 491`

So the direct `5^3` core is not merely redundant fine resolution that can be
compressed away cheaply by a tensor-style local reduction.

### Modest distortion changes stay in the trusted regime

The modest `a = 0.20` and `a = 0.30` complete-shell anchors remain in the same
good regime as the current trusted `a = 0.25` anchor:

- overlap stays clean
- transformed weights stay finite and positive
- projected fixed-only interaction transfer stays small
- low-energy parent one-body capture stays excellent
- final `E1` and `⟨Vee⟩` remain close

Among these modest distortion variants:

- `a = 0.20` is the strongest “tighter core” candidate
- it improves `E1` slightly at both counts
- while keeping `⟨Vee⟩` and projected interaction transfer in the same trusted
  regime

So if the atomic line ever reopens one more tuning lever, the cleaner lever is
distortion/core-spacing, not `5^3` direct core compression.

## Recommendation

Plainly:

- `5^3 -> 4^3` does **not** help
- `5^3 -> 3^3` does **not** help
- a modest distortion adjustment does dominate those core-compression variants
  in the only practically relevant sense here:
  it preserves the trusted physical regime cleanly, while direct core
  compression does not

So the atomic roadmap should **stop atomic core-compression work** on this
line.

Practical conclusion:

- keep the current corrected complete-shell hybrid as the trusted reduced
  atomic anchor
- if a later atomic tuning pass is still desired, test a modest tighter map
  such as `a ≈ 0.20`
- otherwise do **not** spend more effort on `5^3` core compression
- move on to diatomics rather than continuing atomic hierarchy work on the
  present core-compression line
