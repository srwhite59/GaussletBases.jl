# Atomic Mean-Field Supporting Notes

This page is the short synthesis page for the atomic mean-field supporting
notes.

For the current atomic branch status, read first:

- [`docs/current_atomic_branch.md`](current_atomic_branch.md)

## What this supporting chain is about

The current atomic mean-field line was built in a narrow sequence:

1. direct/Hartree term
2. exchange term
3. Fock-style algebraic combination
4. spin-aware UHF convention

Those notes are all still useful, but they are supporting notes rather than
the first pages a new reader should open.

## Recommended supporting-note order

Read them in this order:

1. [`docs/atomic_ida_direct.md`](atomic_ida_direct.md)
2. [`docs/atomic_ida_exchange.md`](atomic_ida_exchange.md)
3. [`docs/atomic_ida_fock.md`](atomic_ida_fock.md)
4. [`docs/atomic_ida_spin_fock.md`](atomic_ida_spin_fock.md)

Then, if you want the actual self-consistent consumer, read:

5. [`docs/atomic_ida_uhf.md`](atomic_ida_uhf.md)

## Current relevance

These notes are still current in the sense that they describe the present
implementation.

But they are not the best summary of the branch state. Their job is to explain
how the present direct / exchange / Fock / spin-aware line was assembled.

## What may be mergeable later

If the atomic line stays scientifically stable for a while, the most likely
future merge target is:

- `atomic_ida_direct.md`
- `atomic_ida_exchange.md`
- `atomic_ida_fock.md`
- `atomic_ida_spin_fock.md`

into one compact mean-field construction note.
