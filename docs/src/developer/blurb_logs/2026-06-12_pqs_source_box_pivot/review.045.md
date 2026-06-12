Review 045: mechanically accepted, scientifically blocked.

The route-owned complete `223`-function H1 path now runs:

```text
complete core/surrounding-shell final basis
-> product-factor support operators
-> final one-body transfer
-> H = K + Z * V_center
-> ordinary symmetric eigensolve
```

The route boundaries are still clean: no fixed-block matrix authority, no
`_pqs_current_route_safe_term_matrices`, no generalized overlap solve, and no
IDA/RHF/GTO work.

However, the H / Z=1 result is not physically acceptable:

```text
H1 energy: -2.0638461028784776
exact H reference: -0.5
```

This is below the variational lower bound for a correctly represented
one-electron hydrogen Hamiltonian, so it is not a basis-quality issue. It is a
one-body operator/convention issue or a fixture/operator mismatch. The probe
records center charge `1.0`, so the next pass must audit the one-body
construction rather than relabeling the result as He+.

Next pass should not add features. It should isolate where the error enters:

- support-space generalized H1 using the same support overlap/Hamiltonian;
- final-basis H1 with kinetic/nuclear expectation decomposition;
- charge/sign convention for the uncharged by-center nuclear matrix;
- consistency of the 1D kinetic and Gaussian nuclear factor matrices with
  existing one-body kernels/oracles for the same geometry if possible.

Only a small convention fix should be implemented in that pass. If the issue is
not a small fix, report the exact blocker and keep this as non-acceptance
mechanical plumbing.

-- repo-manager@macmini
