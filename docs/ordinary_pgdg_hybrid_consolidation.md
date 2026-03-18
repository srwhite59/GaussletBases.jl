# Ordinary PGDG Hybrid Consolidation

This note records the current practical reading of the ordinary PGDG branch.

## 1. What now looks good enough

The recent hydrogen and SHO checks support a clear practical conclusion:

- the ordinary PGDG backend is good enough in the friendly mild/core-supported
  regime
- the hard pure mapped small-`c` cases should still be read as stress tests
- the remaining raw kinetic-matrix mismatch is not, by itself, a reason to
  keep pushing backend refinements right now

So the ordinary branch can now be read this way:

- numerical mapped route remains the reference / validation path
- experimental PGDG-style analytic route is acceptable in the friendly hybrid
  regime
- the ordinary He solver should still wait until the practical basis-design
  questions are clearer

## 2. What should *not* be treated as settled

The current Asinh-based mapping family should not be treated as canonical.

At the moment it is only:

- the current working full-line map
- a narrow practical helper family
- a convenient way to reproduce the present hydrogen and SHO studies

It is **not** yet a claim that:

- `AsinhMapping` is the final best map
- the current `c,s` heuristics are settled
- `c` and `s` can be tuned independently without tradeoffs

The older gausslet work suggests a subtler picture:

- decreasing `c` often also required decreasing `s`
- for larger `Z`, both smaller `c` and smaller `s` may be needed
- completeness very near the nucleus and completeness at intermediate radius
  can pull in different directions

So the current mapping heuristics should be read as provisional.

## 3. Why the radial branch stays numerical

This does **not** reopen the radial PGDG idea.

The radial branch still has the harder analytic problem:

- half-line / `\Theta(r)` structure
- more awkward mapped primitive integrals
- less forgiving near-origin behavior

So the radial branch should stay numerical for now.

The practical PGDG-style analytic work remains on the ordinary Cartesian side.

## 4. The practical ordinary route

The practical ordinary-gausslet route is now understood to be:

- milder distortion
- larger `c` than the hard pure stress-test regime
- explicit core Gaussian augmentation
- overlap cleanup / localization on the combined 1D space

In other words, the intended route is hybrid:

- mapped ordinary gausslet backbone
- plus explicit core Gaussians

That is much closer to the practical White–Lindsey operating regime than the
hard pure mapped small-`c` tests.

## 5. Legacy basis-data and loader references

The relevant legacy sources are:

- basis-set data file:
  `~/BasisSets`
  which for the current setup is
  `/Users/srw/BasisSets`
  which currently resolves to
  `/Users/srw/Dropbox/GaussletModules/BasisSets`
- parser / file-format reference:
  `/Users/srw/Dropbox/GaussletModules/ReadBasis.jl`
- main hybrid Gaussian consumer reference:
  `/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl`
- older basis-construction / trimming reference:
  `/Users/srw/Dropbox/GaussletModules/Old/BasisGen.jl`

The natural later loader to adapt is a small modern equivalent of

```julia
getbasis(atom, basisname; maxl = ...)
```

from `ReadBasis.jl`, returning the same basic tuple structure

```julia
(l, zetas, coefficients)
```

rather than a wholesale port of the old modules.

## 6. Candidate H/He core-Gaussian sets to revisit later

From the current legacy `BasisSets` file, the first short candidate list to
revisit for later hybrid ordinary work is:

- H:
  - `cc-pVDZ`
  - `aug-cc-pVDZ`
  - `cc-pVTZ`
- He:
  - `cc-pVDZ`
  - `cc-pVTZ`
  - `cc-pVTZnc`

Two smaller special-purpose entries are also present and worth remembering as
reference points, but not as the first practical hybrid targets:

- `STO-3G`
- `He gauss8`

The old consumer path that actually injects these Gaussian sets into a hybrid
gausslet workflow is the `ReadBasis.getbasis(...)` +
`PureGaussianGausslet.contractedGTOs(...)` / `getGaussianbasis(...)` line,
with `Old/BasisGen.jl` as a secondary reference for trimming by width.

## 7. What the next expensive pass should be

The next heavy pass, when it is worth doing, should be about:

- practical hybrid basis design
- legacy-informed H/He core-Gaussian choices
- mapping-family and coupled `c,s` regime exploration

It should **not** be another backend-comparison campaign.
