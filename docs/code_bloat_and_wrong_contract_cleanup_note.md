# Code Bloat And Wrong-Contract Cleanup Note

This note records one general repo-side engineering rule.

It was sharpened during the Cartesian hybrid bundle/transfer work, but it is
meant to apply more broadly.

## Main Rule

There are two different costs of extra code:

- **volume cost**
- **conceptual cost**

Volume cost means:

- more lines
- more branches
- more names
- harder review
- harder later modification

Conceptual cost means:

- old code still appears to support a contract the repo no longer believes
- later code starts conforming to the wrong idea
- diagnostics and benchmarks reinforce the wrong mental model
- developers waste time debugging or optimizing a path that should not exist

The second cost is worse.

Small code that encodes the wrong contract is more dangerous than large code
that is conceptually clean.

## Cleanup Rule

When a contract changes:

- do not only patch the main implementation
- also remove or simplify the stale helpers, diagnostics, benchmarks, and
  docstrings that still encode the old idea

Do not keep semi-useless code around merely because it is small or because it
"might be useful later" if it preserves a conceptually wrong path.

If a path is no longer part of the intended repo contract, the default should
be:

- remove it
- or quarantine it clearly as debug-only / reference-only

## Current Trigger Example

The immediate example was the final-basis overlap handling on the Cartesian
hybrid bundle transfer line.

The standing repo policy for orthonormal final blocks is already recorded in:

- `docs/src/developer/numerical_contracts.md`
- `DESIGN.md`

That policy says:

- construct the final block to be orthonormal
- check that its overlap is `I +` small Float64 noise
- then treat the block as orthonormal

Accordingly:

- final-basis self-overlaps are diagnostic only
- they are not part of the normal downstream working model
- generalized-overlap downstream logic is not the intended final-basis path

So if code still:

- builds `S_AA` and `S_BB` on the normal final-basis transfer path
- solves generalized transfer equations on that path
- or benchmarks generalized-eigenproblem logic on that path

then that code is not just extra code. It is code that preserves the wrong
contract and should be removed or rewritten.

## Basis Bundle Refinement

One slight refinement is allowed on the bundle side.

The repo should distinguish between:

- **basis-defining core data**
- **basis-attached helpful sidecars**

Helpful sidecars such as factorized parent structure, shell/block grouping, or
cached axis tables may be worth carrying if they are genuinely useful
downstream.

But they must be treated as:

- auxiliary
- reproducible from the basis construction
- and not confused with the basis definition itself

This is different from keeping conceptually wrong mathematics alive in the
working path.

## Practical Rule Of Thumb

Use this priority order:

1. Remove code that preserves a wrong contract.
2. Simplify code that is only debug/reference baggage.
3. Only after that worry about raw line count.

If a feature adds code but supports a real active contract, it may be worth
keeping.

If a feature adds even a small amount of code while preserving a contract the
repo no longer intends to support, it should usually be deleted.
