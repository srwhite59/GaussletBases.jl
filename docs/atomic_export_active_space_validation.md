# Atomic Export Active-Space Validation

This note records the next cross-repo validation step after the first sliced
export smoke check.

## 1. Why this came next

Schema compatibility is now established:

- the sliced atomic export loads through `HamIO.read_ham(...; validate=true)`
- the grouped keys and block shapes match the existing consumer expectation

So the next question was one layer deeper:

- do the existing active-space and block utilities accept the exported file
  cleanly?

This is still cross-repo validation work, not new solver development inside
`GaussletBases`.

## 2. Consumer path exercised

The deeper consumer path tested was:

- `work/slicedmrgutils/src/HamIO.jl`
- `work/slicedmrgutils/src/Ordering.jl`
- `work/slicedmrgutils/src/StageMap.jl`
- `work/slicedmrgutils/src/ActiveSpaceOps.jl`

The specific checks were:

1. write a fresh atomic sliced export
2. `read_ham(...; validate=true)`
3. `make_ordering(ham; within_slice=:from_file)`
4. `make_active_map(ham, ord, 1; lmax=1)`
5. `ActiveSpaceOps.dense_h1(ham, active_map)`
6. `ActiveSpaceOps.build_active_maps(...)`
7. `ActiveSpaceOps.build_dense_vblocks(...)`

## 3. Result

The active-space utilities accepted the export cleanly.

Observed on the current small atomic test:

- `read_ham(...; validate=true)` succeeded directly
- `active sites = 136`
- `active slices = 34`
- `dims_active[1] = 4`
- `dense_h1` reconstructed the producer-side sliced `H1` with infinity-norm
  error `0.0`
- `build_dense_vblocks` reproduced the producer-side density-density `Vee`
  exactly:
  - max diagonal difference `0.0`
  - max off-diagonal entry `0.0`

The metadata also remained explicit:

- `interaction_model = "density_density_ida"`

## 4. Convention check

No deeper convention mismatch showed up in this pass.

In particular:

- the slice ordering `l0_desc_mzigzag` was accepted as written
- the active-space map matched the exported within-slice labeling
- the pair-diagonal `Vblocks` convention was accepted cleanly by
  `ActiveSpaceOps.build_dense_vblocks(...)`

So the present sliced atomic export is semantically compatible with the
consumer-side active-space/block utilities tested here.

## 5. What this does and does not establish

This establishes:

- producer-side grouped export is readable
- slice/basis ordering is coherent downstream
- the current density-density / IDA model survives one more consumer layer

It does **not** add new solver work inside this package.

## 6. Next step

The next producer-side step should still defer ordinary export.

The better next external validation would be one slightly denser atomic
consumer check, for example through:

- `HFInitHFDMRGBridge.jl`
- or one small existing `slicedmrgutils` run path

Only after that should ordinary export be reconsidered.
