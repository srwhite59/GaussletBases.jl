# Round 008 Slice B Freeze

Status: Slice B frozen and bound

Date: 2026-06-20

Repo-manager froze Slice B after the successful Slice A terminal-basis
implementation and user approval to proceed to one-body operators.

## Newly Approved Slice B ID

```text
HP-FN-03
```

`HP-FN-03` is now implementation authority through
`docs/src/developer/cartesian_hamiltonian_producer_design.md` and the binding
rule in `AGENTS.md`.

## Already Approved Slice A IDs

```text
HP-OBJ-01
HP-OBJ-02
HP-FILE-01
HP-FN-00
HP-FN-01
HP-FN-02
HP-WIRE-01
```

## Still Candidate-Only

```text
HP-FN-04
HP-FN-05
```

These do not authorize IDA assembly, Hamiltonian artifact production, or driver
simplification.

## Freeze Boundary

Slice B authorizes final-basis one-body operator assembly only:

- kinetic matrix `K`;
- separated unit nuclear attraction matrices `U_A = -1/r_A`;
- one separable product-term helper accumulating into a caller-owned dense final
  matrix;
- no global support-space operator matrix;
- no atomic/diatomic branches;
- no IDA or artifact production.

Unlisted production surfaces require a prior docs-only amendment.
