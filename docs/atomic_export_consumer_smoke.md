# Atomic Export Consumer Smoke Check

This note records the first external consumer-side smoke check for the atomic
Hamiltonian export path.

## 1. Why this came next

The package now has two atomic export layers:

- dense `fullida_dense_v1`
- grouped sliced/block atomic export

After those producer layers landed, the next clean step was not more in-package
solver work.

It was a narrow external validation through the existing consumer-side reader.

## 2. Consumer path tested

The consumer-side reference used for the smoke check was:

- `work/slicedmrgutils/src/HamIO.jl`

The specific external check was:

- write a small atomic sliced/block export file from `AtomicIDAOperators`
- load it with `read_ham(...; validate=true)`
- inspect the recovered slice counts, dimensions, labels, and metadata
- reconstruct a dense `H1` from the consumer-side block structure

## 3. Result

The external smoke check succeeded directly.

Observed:

- `read_ham(...; validate=true)` succeeded with no schema mismatch
- `nslices = 34`
- `dims[1] = 4`
- `within_slice = "l0_desc_mzigzag"`
- first label recovered as `r=1,l=0,m=0`
- `meta["interaction_model"] = "density_density_ida"`
- dense `H1` reconstructed from `ham.H1blocks` matched the producer-side
  sliced `H1` with infinity-norm error `0.0`

So the current sliced export is already compatible as written with the present
`HamIO.jl` reader shape.

The producer surface now also provides honest physical shell coordinates
through `layout/slice_coord`, while keeping the radial index available
separately for compatibility.

## 4. Interpretation

This is still producer/export work.

The point of the smoke check was:

- confirm that the grouped file shape is readable by the existing consumer
- confirm that the exported slice/basis metadata is coherent
- confirm that no tiny schema repair layer is needed immediately

It was not to add solver logic inside this package.

## 5. Next producer step

The next producer-side step should still not be ordinary export yet.

The better next move is a slightly denser atomic consumer validation outside
this package, for example:

- loading the sliced atomic file in a small existing `slicedmrgutils` driver
- checking a few downstream assumptions beyond `read_ham`

Only after that should ordinary export be reconsidered.
