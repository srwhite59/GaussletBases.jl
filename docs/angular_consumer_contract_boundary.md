# Angular Consumer-Contract Boundary

This note records the current consumer-contract boundary reached after the
first angular exact-reference HamIO / HFDMRG-facing HF bridge pass.

## What already fits cleanly

The exact common low-`l` injected reference already fits the current
HamIO / HamV6 consumer language cleanly.

That currently valid producer-side bridge is:

- `angular_benchmark_exact_hamv6_payload(...)`
- `write_angular_benchmark_exact_hamv6_jld2(...)`

What it exports is intentionally narrow:

- the exact common low-`l` reference carried inside the angular benchmark line
- in the same grouped HamV6 / `HamIO` language already proven on the atomic
  radial line
- for the present HF consumer bridge path

This is an HF bridge boundary, not yet a true many-body DMRG bridge.

## What does not yet fit cleanly

The full mixed shell-local angular basis does **not** yet fit the current
HamIO / HamV6 consumer language cleanly.

That is not a numerical failure. The current shell-local angular line is
numerically sane through:

- one-electron benchmark
- HF-style benchmark
- small-ED benchmark

The mismatch is representational.

## Why the mismatch is representational

The current HamIO / HamV6 language assumes that each exported orbital can be
described by simple per-orbital labels such as:

- slice identity
- definite `l`
- definite `m`

The full mixed shell-local angular basis no longer has that simple structure.

Its orbitals are built from:

- injected exact low-`l` subspace pieces
- shell-local Gaussian / localized angular prototype pieces
- orthogonalized mixed combinations inside a shell

So a full orbital is no longer described honestly by one plain `l,m` tag.

## Likely missing representation layer

A fuller angular consumer contract will likely need a richer basis-description
layer than plain per-orbital `l,m`.

The minimum plausible ingredients look like:

- shell provenance
  - shell index / shell radius / shell order
- orbital kind
  - exact injected sector vs mixed remainder
- injected-subspace labeling
  - which exact low-`l` channels are represented exactly
- mixed-basis coefficient metadata
  - enough shell-local coefficient information to describe how each exported
    orbital sits inside the shell-local angular basis
- possibly a consumer capability flag
  - so downstream code can distinguish plain `l,m` orbitals from mixed angular
    orbitals without guessing

This note does **not** claim that all of those fields are final.
It only records the shape of the contract question.

## Practical branch-point conclusion

The current angular line has reached a useful but narrow downstream boundary:

- exact common low-`l` reference:
  compatible now with HamIO / HamV6 and the current HF bridge path
- full mixed shell-local angular basis:
  compatible now with the current in-memory dense HFDMRG-facing HF adapter
  path, but not yet honestly representable in the HamIO / HamV6 file language

So the next implementation question is a contract-design question, not a
broader downstream workflow push.

The next step should therefore be:

- a narrow angular export-design / consumer-contract formalization pass

and not yet:

- a larger downstream migration effort
- a Hooke-facing workflow push
- a claim of full mixed-basis HamIO compatibility
