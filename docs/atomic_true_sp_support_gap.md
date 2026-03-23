# Atomic True-SP Support Gap

This note records the gap that existed before the explicit 3D atomic shell
route landed. It is now superseded for ordinary-QW and nested-QW by
`docs/atomic_true_sp_support_note.md`.

This pass checks whether the current atomic hybrid/QW path can honestly make
`lmax = 1` a true active physical supplement route.

## Historical result

At the time of this note, `lmax = 1` was still **not** a true active physical path in the then-present
ordinary-QW or nested-QW implementation.

The code now rejects that usage explicitly instead of silently collapsing it to
the centered `s` channel.

## Why the current route is not enough

The present active QW/nested supplement machinery is still built around one
centered separable 1D Gaussian channel:

- [`_qwrg_supplement_primitives_and_contraction(...)`](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/ordinary_qiu_white_rg.jl)
  returns one 1D primitive Gaussian list plus one contraction matrix
- [`_qwrg_raw_overlap_blocks(...)`](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/ordinary_qiu_white_rg.jl)
  and the related raw-block builders assemble 3D supplement blocks by cubing
  the same 1D channel across `x/y/z`
- [`hybrid_mapped_ordinary_basis(...)`](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/ordinary_hybrid.jl)
  likewise only admits a centered separable 1D supplement route

That is honest for active `s` support, but not for `p_x/p_y/p_z`.

True atomic `p` support needs explicit mixed-axis shell content, not “the same
1D centered Gaussian set cubed.”

## Legacy evidence

The legacy atomic line already did this in the right structural way:

- [`getbasis(...; maxl = lmaxadd)`](/Users/srw/Dropbox/GaussletModules/ReadBasis.jl)
  filtered named-basis shells by angular cutoff
- [`makeallcontractions(...)`](/Users/srw/Dropbox/GaussletModules/ReadBasis.jl)
  expanded each shell into explicit Cartesian `(l_x,l_y,l_z)` contractions
- [`getGaussianbasis(...)`](/Users/srw/Dropbox/GaussletModules/PureGaussianGaussletOld.jl)
  then kept those as explicit 3D Gaussian orbitals with separate `x/y/z`
  factors

So the missing piece in this repo is not the loader. It is the active
supplement representation and consumer algebra.

## What changed in that pass

The current active consumers now reject `LegacyAtomicGaussianSupplement` objects
that contain non-`s` shells:

- [`hybrid_mapped_ordinary_basis(...)`](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/ordinary_hybrid.jl)
- [`ordinary_cartesian_qiu_white_operators(...)`](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/ordinary_qiu_white_rg.jl)
  for both ordinary and nested fixed-block routes

This keeps the atomic line honest:

- `lmax = 0` is real and active
- `lmax = 1` is loader metadata only until the missing representation step is
  implemented

## Narrowest honest implementation path

To make true atomic `SP` support real without pretending, the next required
architectural step is:

- add an explicit atomic-centered 3D shell supplement object up to `lmax = 1`
- represent active supplement orbitals as contracted Cartesian shells
  `s, p_x, p_y, p_z`
- carry their separate `x/y/z` 1D factors explicitly in the active QW/nested
  raw-block assembly
- replace the current “same 1D supplement set cubed” assumption on the
  supplement side

That architectural step has now been implemented for ordinary-QW and nested-QW.
The one-dimensional hybrid builder remains honestly `s`-only.
