# Same-Gauge IDA Nuclear External Potential

Status: approved source authority under `HP-PQS-IDA-NUCEXT-FN-01` and
validation authority under `HP-PQS-IDA-NUCEXT-TEST-01`.

This lane approves one missing primitive needed by screened-field audits:

```text
uN_IDA[A, i]
```

`uN_IDA[A, i]` is the point-nucleus external potential of the same normalized
final-row IDA density proxy used by `electron_electron_ida`. It is not a
Galerkin one-body matrix and not a substitute constructed from already
assembled nuclear-attraction operators.

## Required Definition

For final IDA row `i` and nucleus `A`,

```text
uN_IDA[A, i] =
    raw point-nucleus external-potential numerator for final IDA row i
    / same final IDA weight used by terminal electron_electron_ida
```

The normalization must be exactly the same final-row IDA weight convention
used in terminal `electron_electron_ida` assembly. For a nucleus with charge
`Z_A`, the sign convention is physical electron-nuclear attraction.

The helper may return a vector for a single nucleus or a by-center table. It
must report the final IDA weights used and must validate that those weights are
finite and positive.

## Owner

Preferred owner:

- `src/cartesian_final_basis_realization/pqs_terminal_ida.jl`

That file owns the terminal IDA proxy/weight convention. Do not make
`src/cartesian_ida_hamiltonian.jl` the owner; it stores finished Hamiltonians
and should not define the row-proxy convention.

Allowed narrow caller plumbing:

- PQS low-order/base materialization path only if needed for validation.

## Explicit Non-Goals

This helper is explicitly not:

- Galerkin `Vnuc`;
- `diag(Vnuc_G)`;
- row action from `Vnuc_G`;
- center evaluation `-Z/r_i`;
- screened-field `Delta_W`;
- `W_IDA`;
- `H1_eff`;
- constants or corrected Hamiltonian assembly;
- rho0/P0 revival;
- EGOI;
- artifact/public workflow/solver integration;
- Cr/Cr2 work.

## Validation

First validation target is H q5 with `core_spacing = 0.3`.

Approved validation for `HP-PQS-IDA-NUCEXT-TEST-01`:

- package load;
- focused ignored probe;
- `git diff --check`;
- final IDA weights used, with positivity and finiteness checks;
- `uN_IDA` finite/range by row class;
- direct-core rows compared to analytic Gaussian-proxy nuclear attraction
  where possible;
- diagnostic-only comparison to `diag(Vnuc_G)`, clearly marked not an
  acceptance criterion;
- diagnostic-only comparison to center `-Z/r`, clearly marked not an
  acceptance criterion;
- far-field behavior against `-Z/r` for localized rows.

Acceptance requires showing that the helper uses the same row proxy and
normalization convention as `electron_electron_ida` and validates on the tiny
H case. This lane must not form screened fields, corrected one-body matrices,
constants, artifacts, or solver inputs.
