# 1D Distorted-Gausslet PGDG Refinement Hierarchy

## Pseudocode

1. Start from one fixed one-dimensional distorted gausslet line.
   The central intermediate is not yet a three-dimensional Hamiltonian. It is a
   one-dimensional distorted-gausslet Coulomb data object.
   Code: current related infrastructure lives mainly in `src/ordinary_mapped_backends.jl`

2. For that one-dimensional distorted gausslet line, define the term-resolved
   Coulomb data against the Gaussian expansion of `1/r`.
   This includes the one-dimensional data needed for both:
   - nuclear-potential assembly
   - electron-electron / `Vee` assembly
   Code: current consumers include `src/ordinary_cartesian_ida.jl` and `src/ordinary_qiu_white_rg.jl`

3. In the pure distorted-gausslet route, build that one-dimensional Coulomb
   data directly in the distorted-gausslet representation.
   This is the conceptually clean reference limit.
   Code: not yet cleanly isolated as a separate repo path

4. In the PGDG surrogate route, represent the same one-dimensional
   distorted-gausslet line through a uniform Gaussian proxy line and build the
   Coulomb data there, then contract back to the working distorted-gausslet
   basis.
   Code: current proxy-based pieces live in `src/ordinary_mapped_backends.jl`

5. Introduce a refinement hierarchy for that PGDG surrogate route.
   The proxy line need not stop at the standard PGDG spacing. It can be
   refined systematically:
   - base PGDG level: proxy spacing `1/3`
   - one refinement level: proxy spacing `1/9`
   - two refinement levels: proxy spacing `1/27`
   - and so on
   Code: framing only for now; not yet implemented as a full hierarchy

6. At each refinement level, build the same one-dimensional Coulomb data
   object:
   - overlap / kinetic / position / `x^2`
   - nuclear Gaussian-factor terms
   - `Vee` pair-factor terms
   The refinement changes how faithfully the proxy line represents the original
   distorted gausslet line, not the role of the final one-dimensional data.
   Code: framing only for now; current bundle is only a first related step

7. Use that one-dimensional Coulomb data object as the input to the final
   Cartesian product assembly.
   Once the one-dimensional data are in hand, the final three-dimensional
   assembly is conceptually straightforward:
   - product assembly over `x`, `y`, `z`
   - short inner reduction over the Gaussian expansion terms
   Code: current term-first three-dimensional assembly lives in `src/ordinary_cartesian_ida.jl` and `src/ordinary_qiu_white_rg.jl`

8. Interpret the hierarchy operationally.
   In the intended first practical picture:
   - the `1/3` level is the base PGDG Gaussian line
   - the `1/9` level is obtained by one application of the stored `1/3`
     refinement mask
   - the `1/27` level is obtained by applying the same mask a second time
   This makes the hierarchy a repeated local refinement of one uniform Gaussian
   line rather than a collection of unrelated proxy constructions.

9. Distinguish two uses of the hierarchy.
   At the coarse `1/3` level, the PGDG line may act as a basis-realization
   layer and may still require basis cleanup or COMX-like postprocessing.
   At the finer `1/9`, `1/27`, ... levels, the refined PGDG line may instead
   be used only as an auxiliary integral-evaluation representation while the
   true basis remains the distorted gausslets.

10. Interpret the hierarchy as a bridge:
   - pure distorted-gausslet treatment at one end
   - standard PGDG at the coarse practical end
   - refined PGDG proxy levels in between
   This makes PGDG not only a practical surrogate, but a systematically
   refinable path toward the pure distorted-gausslet limit.

## References

- S. R. White and collaborators, legacy distorted-gausslet and hybrid code in `GaussletModules`
- Y. Qiu and S. R. White, "Hybrid gausslet/Gaussian basis sets"
- Repo algorithms page: `algorithms/qiu_white_residual_gaussian_route.md`

## What This Frames

This page is a provisional framing document for the repo's ordinary-branch
design.

It separates four ideas that had been partially mixed together in discussion:

- pure Qiu-White
- QW-PGDG
- the one-dimensional distorted-gausslet PGDG refinement hierarchy
- later COMX / localization / nesting lines

The most important conceptual point is that the one-dimensional
distorted-gausslet Coulomb data object is the real intermediate.

Once that one-dimensional object is available, the later three-dimensional
assembly is relatively simple and follows the familiar short-inner-loop
structure over the Gaussian expansion terms.

The next important distinction is that the hierarchy can be used in two
different ways:

- as a coarse basis-realization layer at the base `1/3` PGDG level
- as a progressively finer auxiliary integral-evaluation layer at `1/9`,
  `1/27`, ... while the true basis remains the distorted gausslets

Those two uses should be judged by different error criteria.

## Current Repo Status

This hierarchy is not yet fully implemented in the repo. The present page is
intended to stabilize the conceptual framework before more implementation work
grows around it.

Current state:

- the repo already has ordinary and Qiu-White-related proxy machinery in
  `src/ordinary_mapped_backends.jl`
- the repo already has efficient term-first three-dimensional assembly routes
  in `src/ordinary_cartesian_ida.jl` and `src/ordinary_qiu_white_rg.jl`
- the repo now has the first stored analytic ternary `1 -> 1/3` refinement
  mask and narrow local application helpers in
  `src/ordinary_pgdg_refinement_masks.jl`
- the repo does not yet expose the one-dimensional distorted-gausslet Coulomb
  data object as the main explicit internal abstraction
- the repo does not yet implement the full PGDG refinement hierarchy as
  full levels such as `1/3`, `1/9`, `1/27`, ... with all downstream ordinary
  consumers wired to it
- the repo does not yet separate the hierarchy explicitly into
  - coarse basis-realization use at `1/3`
  - finer auxiliary integral-evaluation use at `1/9`, `1/27`, ...

So this page should be read as a shared design vocabulary and provisional
algorithmic framing, not as a claim that the full hierarchy is already
available.

## Why This Matters

This framing helps keep several future lines separate:

- Pure Qiu-White:
  true distorted-gausslet line, no PGDG, no COMX
- QW-PGDG:
  Qiu-White-style hybrid/RG work using PGDG or proxy-Gaussian machinery for the
  hard one-dimensional distorted-gausslet Coulomb data
- 1D distorted-gausslet PGDG refinement hierarchy:
  a systematic family of proxy realizations of the one-dimensional
  distorted-gausslet Coulomb data, with both basis-realization and
  integral-evaluation roles
- COMX / later localization-contraction work:
  a separate later line, not part of the basic pure-QW or PGDG-refinement
  definitions

That separation should make later algorithm pages and implementation reviews
cleaner.

## Implementation Notes

If this hierarchy becomes an implemented internal object family, the repo
should eventually name that intermediate explicitly in code and tests.

The likely target is a shared internal one-dimensional bundle carrying the
term-resolved Coulomb data for a distorted-gausslet line, with multiple
construction backends:

- direct distorted-gausslet construction
- standard PGDG proxy construction at base spacing `1/3`
- refined PGDG proxy construction at finer proxy spacings such as `1/9`,
  `1/27`, ...

The likely first practical hierarchy uses one stored local refinement mask,
applied repeatedly to the base `1/3` Gaussian line.

For the current supporting math and numerics behind that first practical mask,
see:

- `docs/gaussian_refinement_analytic_mask_note.md`
- `docs/gaussian_refinement_analytic_mask_comparison.md`
- `docs/pgdg_refinement_hierarchy_first_mask.md`

That first mask is now settled and implemented narrowly, but the larger
hierarchy and downstream integration are still open, so this page still stops
mostly at the framing level.
