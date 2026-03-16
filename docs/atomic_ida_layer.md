# Static atomic IDA ingredients

This note describes the next atomic layer after the first one-electron `(l,m)` construction.

The goal here is **not** to solve the many-electron problem yet.

In other words, this pass should **not solve the many-electron problem**. It should only assemble the static ingredients cleanly.

The goal is to assemble the interacting atomic Hamiltonian ingredients cleanly and explicitly on top of the radial and angular layers that are now already in place.

## 1. The intended layering

The present atomic story should now be read in this order:

1. a radial basis and radial quadrature
2. radial one-body operators
3. explicit angular channels labeled by `(l,m)`
4. radial multipole operators
5. angular Gaunt and M-summed coupling factors
6. a combined interacting atomic data object

That is the right layer to build before any many-electron solve.

## 2. Why this is the next step

The package now already has:

- a clean radial substrate
- a shared primitive/contraction story under the radial operators
- an explicit one-electron atomic `(l,m)` layer for hydrogen

So the next question is no longer how to organize the one-electron atom.

The next question is:

**how should the static interacting atomic ingredients be represented so that they are explicit, inspectable, and ready for later He / IDA work?**

## 3. What is being assembled here

The interacting atomic layer should bundle, at minimum:

- the one-body atomic blocks
- the radial multipole tables
- the angular coupling data
- enough orbital and channel indexing information to identify what each block means

This is still a data-assembly step, not a solver step.

## 4. Radial part

The radial part already exists in `RadialAtomicOperators`.

That object gives:

- `ops.overlap`
- `ops.kinetic`
- `ops.nuclear`
- `centrifugal(ops, l)`
- `multipole(ops, L)`

So the interacting layer should not rebuild those ideas. It should sit directly on them.

## 5. Angular part

For the first interacting layer, the natural angular data are:

- explicit channels `(l,m)`
- Gaunt coefficients for coupling one-electron channels to multipole channels
- the corresponding M-summed angular kernels `Q_L`

Those angular kernels are the natural partners of the radial multipole tables in the current IDA-style factorization.

## 6. Why this is better than jumping straight to a solve

There are two reasons.

First, it keeps the physical ingredients visible:

- one-body radial terms
- radial multipoles
- angular selection rules
- channel and orbital indexing

Second, it makes later solver work much easier to trust. If the static interacting object is explicit and inspectable, then later HF, DMRG, or other many-electron code has a cleaner foundation.

## 7. What the first version should and should not do

The first version should:

- work for one radial basis and one `lmax`
- make the one-body and two-body ingredients explicit
- keep the indexing conventions simple
- make the selection rules easy to inspect

The first version should **not** yet:

- solve the many-electron Hamiltonian
- optimize performance aggressively
- introduce a large atom workflow layer

## 8. Recommended object structure

The clean object for this pass is a bundle in the spirit of:

- `AtomicIDAOperators`

containing:

- the one-body atomic operator bundle
- the radial multipole tables
- the angular Gaunt / kernel data
- the orbital/channel indexing metadata

The important point is that the object should expose the interacting structure, not hide it.

## 9. Orbital indexing

The simplest orbital convention for the first version is a **channel-major** ordering:

- order channels by the existing `YlmChannelSet`
- within each channel, order by radial basis index

That is consistent with the current one-body atomic block structure and easy to inspect in examples and tests.

## 10. Bottom line

This pass should establish the static interacting atomic structure on top of the now-clean radial-plus-angular architecture.

Only after that should the repository move on to an actual many-electron solve.

That next narrow step is now described in:

- [`atomic_ida_two_electron.md`](atomic_ida_two_electron.md)
