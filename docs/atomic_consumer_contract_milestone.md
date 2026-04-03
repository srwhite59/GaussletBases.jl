# Atomic Consumer-Contract Milestone

This note records the branch point reached after the recent atomic
producer/consumer integration pass.

## State at this point

- The atomic consumer contract is the main active producer-side milestone.
- A narrow bond-aligned heteronuclear diatomic `HeH+`-style line is already
  implemented and test-backed, but it remains checkpointed rather than
  promoted as a broad public workflow.
- The active producer-side milestone is commit `6415407`:
  `Add atomic HamV6 compatibility export contract`.

## What is now real on the producer side

`GaussletBases` now owns a clear atomic producer-side contract for downstream
HF / ED / DMRG-facing use.

That contract now includes:

- a public atomic IDA current-model interaction accessor
- an explicit HamV6 / HamIO compatibility export path
- consumer-facing top-level metadata including `meta["Z"]`
- the preserved native sliced export, kept intact as the package's honest
  native atomic sliced shape

## What was validated downstream

The contract was checked with a canonical atomic He ladder.

Validated pieces:

- `HamIO.read_ham`
- `hfdmrg_from_ham.jl`
- the atomic He HF bridge path
- internal consistency of the public current-model interaction accessor with
  dense export
- two-electron ED consistency on the same current model

The practical result was that no structural producer/consumer contract mismatch
showed up on that ladder.

## Practical conclusion

`GaussletBases` is now in a credible state as the preferred producer-side
foundation for the current atomic HF / ED / DMRG-facing downstream workflows.

This does **not** mean the package should absorb those consumer workflows.
It means the producer boundary is now explicit and tested enough that
downstream repos can depend on it more cleanly.

## Roadmap consequence

This is the point to stop pushing deeper atomic migration for the moment.

The intended next priority is:

- keep the heteronuclear diatomic line at narrow checkpoint scope rather than
  broadening it into a general molecule workflow
- keep the atomic consumer contract as the current proven branch point
- shift effort toward angular / Hooke-facing downstream support
