Accepted pass 182.

The ignored local probe ran over the corrected 419-dimensional fixed-q PQS He
H1-J payload and produced the first private RHF diagnostic scalar:

```text
input contract = available_pqs_multilayer_complete_core_shell_rhf_input_contract
SCF status = materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload
converged = true
iterations = 7
total energy = -2.8493871224303167
one-body energy = -3.8675645938544685
two-body energy = 1.018177471424152
density trace = 1.9999999999999993
commutator residual = 3.2524884209517158e-9
idempotency error = 4.163336342344337e-17
```

This is a private Hamiltonian-validator result only. It is not route-driver
RHF, not a public solver surface, and not a physics endpoint.

The manager reran:

```text
julia --project=. tmp/work/pqs_he_419_rhf_diagnostic_probe.jl
```

and reproduced the same scalar results. `git diff --numstat -- src test` is
empty, as required for the no-tracked-source/test probe pass.

The next useful step is not more RHF machinery. It is comparison hygiene:
recover or compute the matching old WL gausslet-only q=5/n_s=5 fixed-shell He
RHF result, distinct from the supplemented AHGBS-9 final 447 result, so the PQS
419 scalar is compared to the correct baseline.

-- repo-manager@macmini
