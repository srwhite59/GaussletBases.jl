# Numerical Contracts

This page records a small internal engineering policy that is easy for machine
generated code to miss.

It is intentionally developer-facing rather than part of the main user manual.

## Orthonormal blocks

When a block is constructed to be orthonormal, the intended contract is:

1. build it to be orthonormal in the relevant metric
2. check that the resulting overlap is `I +` small unavoidable Float64 noise
3. then treat the block as orthonormal

Examples include:

- finalized PGDG / COMX-cleaned blocks
- nested fixed blocks
- any other internal block whose construction is explicitly meant to produce an
  orthonormal basis

## What not to do

Do not keep propagating a near-identity overlap matrix as if it were meaningful
mathematical data.

In particular:

- do not store `S = I + ε` by default when `ε` is just Float64 residue
- do not build downstream logic that keeps consulting such an `S`
- do not interpret tiny nonorthogonality at the `1e-12` to `1e-14` level as a
  real structural feature

## What to do instead

Use overlap matrices in this regime only for:

- construction
- validation
- assertions
- or one final cleanup step if the deviation from identity is not small enough

After that, the working representation should simply regard the block as having
identity overlap.

This is a coding and design policy, not a user-facing scientific statement.
