# Atomic Hybrid Anchor Comparison

This pass settles the first trusted nonrecursive atomic anchor on the current
shared hybrid footing.

## Scope

Atomic-only, with the active current supplement path:

- `legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)`

used unchanged across:

- the unnested hybrid QW reference
- the shell-plus-core nested fixed block
- the corrected complete-shell nested fixed block

No `lmax = 1` physical claim is made here. That still requires a richer
nonseparable supplement model.

## Cases Compared

Primary benchmark:

- He
- fixed-`a` family with `a = 1/4`
- `xmax = 10`
- `count = 17`

Nearby robustness case:

- same setup, `count = 15`

In both cases the comparison uses the same shared named-basis supplement and
the same nearest/GGT QW consumer.

## Promotion Question

Can the corrected complete-shell line now be promoted from diagnostic-only to
the first trusted reduced nonrecursive atomic anchor beyond shell-plus-core?

The promotion standard for this pass was:

- no known construction bug on the path
- overlap and transformed fixed-block weights remain well behaved
- projected fixed-only interaction transfer stays close to the unnested hybrid
  reference
- low-energy parent one-body content is retained well
- final nearest/GGT `E1` and `⟨Vee⟩` stay close to the unnested hybrid
  reference
- the reduced line still gives a material dimension reduction

## Count-17 Results

Shared supplement:

- `LegacyAtomicGaussianSupplement(atom="He", basis="cc-pVTZ", lmax=0, shell_ls=[0], nshells=3, nactive_primitive=6, nactive_contracted=3, uncontracted=false)`

Results:

| Case | Fixed dim | Overlap error | `E1` | `⟨Vee⟩` | Projected fixed-only `⟨Vee⟩` shift | Ground capture | Avg first 4 capture | Weight min / max | Representative warmed time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Unnested hybrid QW | `4913` | `2.58e-12` | `-1.9985884400290912` | `1.248965557065537` | `0` | `1.0` | `1.0` | `0.06147 / 7.08141` | `~9.4 s` |
| Shell-plus-core | `1403` | `1.29e-12` | `-1.9985629071304736` | `1.2489021100178594` | `+1.51e-6` | `0.9999992615` | `0.9919606525` | `0.06147 / 2.97748` | `~3.5 s` |
| Corrected complete-shell | `589` | `4.04e-12` | `-1.9981842264017804` | `1.24893460807307` | `+1.33e-4` | `0.9999867458` | `0.9966550917` | `0.06147 / 3.25917` | `~1.6 s` |

Key count-17 comparisons:

- complete-shell vs unnested:
  - `ΔE1 = +4.04e-4`
  - `Δ⟨Vee⟩ = -3.09e-5`
- complete-shell vs shell-plus-core:
  - `Δ⟨Vee⟩ = +3.25e-5`
- dimension reduction:
  - `4913 -> 589`

## Count-15 Robustness Results

Results:

| Case | Fixed dim | Overlap error | `E1` | `⟨Vee⟩` | Projected fixed-only `⟨Vee⟩` shift | Ground capture | Avg first 4 capture | Weight min / max | Representative warmed time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Unnested hybrid QW | `3375` | `1.95e-12` | `-1.997177156084648` | `1.2478865431062496` | `0` | `1.0` | `1.0` | `0.07750 / 8.11068` | `~3.6 s` |
| Shell-plus-core | `1403` | `1.78e-12` | `-1.9971774480872162` | `1.2478972937836572` | `+8.84e-7` | `0.9999996055` | `0.9998654543` | `0.07750 / 5.58534` | `~3.0 s` |
| Corrected complete-shell | `589` | `3.35e-12` | `-1.9970883634888783` | `1.247898324597102` | `+6.69e-5` | `0.9999958918` | `0.9995172562` | `0.07750 / 6.25532` | `~0.9 s` |

Key count-15 comparisons:

- complete-shell vs unnested:
  - `ΔE1 = +8.88e-5`
  - `Δ⟨Vee⟩ = +1.18e-5`
- complete-shell vs shell-plus-core:
  - `Δ⟨Vee⟩ = +1.03e-6`
- dimension reduction:
  - `3375 -> 589`

## Promotion Decision

Yes. The corrected complete-shell hybrid line is now good enough to promote to
the first trusted reduced nonrecursive atomic anchor beyond shell-plus-core.

Why:

- the known construction bugs on that path are already fixed
- projected fixed-only interaction transfer is now close to the parent/unnested
  reference on both cases
- low-energy one-body capture is excellent on both cases
- final nearest/GGT `E1` and `⟨Vee⟩` stay in the same good physical regime as
  the unnested hybrid reference
- the reduction is material:
  - `4913 -> 589` at `count = 17`
  - `3375 -> 589` at `count = 15`
- runtime also drops materially in the same comparison

What remains true:

- shell-plus-core is still a trusted conservative comparison state
- but it is no longer the only trusted nonrecursive anchor

## Practical Trust Status

Trusted atomic nonrecursive anchors now are:

- unnested hybrid QW reference
- shell-plus-core hybrid
- corrected complete-shell hybrid

The corrected complete-shell line is the first trusted reduced anchor in that
set.

## Next Phase

With the atomic anchor now settled on the proper shared hybrid footing, the
next nesting phase should move from anchor establishment to hierarchy design:

- decide what the first post-nonrecursive nesting step should be
- keep the same atomic `lmax = 0` supplement footing
- leave diatomic placement/policy and true `lmax = 1` physical anchor work for
  later passes
