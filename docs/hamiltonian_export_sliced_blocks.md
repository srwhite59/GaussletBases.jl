# Sliced Atomic Hamiltonian Export

This note records the second Hamiltonian export step for downstream solver
codes.

## 1. Why dense export came first

The dense bridge came first because it was the shortest honest export path for
the current atomic IDA model:

- the format already existed as `fullida_dense_v1`
- the surrounding ecosystem already had a consumer for it
- `AtomicIDAOperators` already had a coherent dense `H1` plus dense
  density-density `Vee` story

That made the dense bridge the fastest way to turn the package into a clean
Hamiltonian producer.

## 2. Why the next target is the sliced/block format

The next export target is the grouped sliced/block Hamiltonian shape already
consumed by `work/slicedmrgutils/src/HamIO.jl`.

That consumer expects grouped JLD2 data under:

- `layout/*`
- `basis/*`
- `ordering/*`
- `onebody/*`
- `twobody/*`

This differs conceptually from the dense bridge:

- the dense bridge is a single global `H1` plus global dense `Vee`
- the sliced export keeps the orbital grouping explicit
- slices, offsets, `(l,m)` labels, and block conventions become inspectable

## 3. What the present sliced export means physically

This export must still reflect the current model honestly.

For the present atomic IDA line, the interaction remains:

- density-density
- two-index IDA
- not a full four-index Coulomb tensor

So the sliced export does **not** pretend to be a full `(pq|rs)` Hamiltonian.
Its `twobody/Vblocks` encode the current density-density interaction in a
pair-diagonal block form, and the metadata says so explicitly.

## 4. First implementation scope

The first sliced export scope is:

- `AtomicIDAOperators` only
- radial-shell slices
- explicit `(l,m)` labels within each slice
- one-body blocks `H1blocks`, `Tblocks`, `Vnucblocks`
- density-density `Vblocks`

The ordinary Cartesian line should still come later. This pass stays atomic
only because the atomic branch already has the natural shell/channel structure
needed for a clean sliced export.

## 5. Producer only, not solver work

This is still producer/export work.

The goal is:

- write the grouped Hamiltonian file correctly
- preserve honest metadata about the present model
- make the file easy for external consumers to load

It is **not** to add more solver machinery inside this package.
