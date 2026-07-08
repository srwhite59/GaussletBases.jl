# Screened-Vnuc Measurement Lane

Status: approved measurement-only authority under
`HP-PQS-SCREENED-VNUC-AUDIT-03`.

This lane continues the screened electron-nuclear measurement work now that the
same-gauge IDA nuclear external-potential primitive exists:

```text
uN_IDA[A, i]
```

The purpose is to measure whether a compact atom-local screening cloud can make
the electron-nuclear field benign when cancellation is formed in the same IDA
Coulomb gauge. This is not production Hamiltonian authority.

## Screened Sum

For an atom-local compact cloud represented by `q0`, the measurement object is:

```text
W_IDA = uN_IDA + V * q0
W_G   = Vnuc_G + J0_G
Delta_W = W_G - W_IDA
```

Here:

- `uN_IDA` is the point-nucleus external potential of the same normalized
  final-row IDA density proxy used by `electron_electron_ida`;
- `V * q0` is the IDA Coulomb field of the represented screening cloud;
- `Vnuc_G` is the Galerkin point-nucleus external potential for the comparison
  basis;
- `J0_G` is the Galerkin Coulomb field of the same compact cloud.

The important distinction is that the cancellation is first formed in the same
IDA Coulomb gauge. This lane must not correct bare `Vnuc` or bare `J0`
separately and then infer a screened result.

## Current Evidence

The `SCREENED-VNUC-IDA-PROMISE-AUDIT-02` H/Be/Be2 probe gave a positive static
signal:

- compact contracted-GTO clouds were represented cleanly:
  - H: `q0 = 0.9997902501 / 1`;
  - Be: `q0 = 3.9999512998 / 4`;
  - Be2: `q0 = 7.9997599506 / 8`;
- effective low one-body shifts were small:
  - H: `-2.49e-5 Ha`;
  - Be: `-6.19e-4 Ha`;
  - Be2: `-1.21e-3 Ha`;
- negative `Delta_W` modes were compact/core-local, not broad/protected or
  residual occupation incentives;
- Be2 radius `>= 2` `Delta_W` range was about `[-7.77e-3, 4.95e-3]`.

This is promising enough for continued measurement. It is not a production
claim.

## Allowed Measurement Work

Approved surfaces:

- ignored `tmp/work/*.jl` probes only;
- durable output tables under `/Users/srw/dmrgtmp`;
- H, Be, Be2, and Cr atom only;
- existing `uN_IDA` helper;
- compact contracted atom-local GTO clouds.

Allowed cloud choices:

- retained original compact `s1`/`s2` where available;
- minimal or STO-like contracted core shells;
- optional fake-RDM-selected compact shell directions when already available in
  probe code.

Allowed screening charge variants:

- neutral all-electron charge;
- closed-shell core charge, such as Cr `[Ar]`-like 18e;
- smaller shell-wise charges when physically meaningful.

Required reported objects include `W_IDA`, `W_G`, `Delta_W`, constants,
spectra, locality/radius bins, row-class decomposition, and orbital
expectations.

Optional Cr atom measurement is allowed only after H/Be/Be2 remain sane.

## Forbidden Work

This lane does not approve:

- tracked source edits;
- artifacts, schema changes, public workflow, or solver workflow;
- Cr2;
- production corrected Hamiltonians;
- EGOI changes or expansion;
- rho0/P0 affine-anchor revival;
- exact exchange;
- replacing `uN_IDA` with `diag(Vnuc_G)`, row action, exact `Vnuc_G`, or center
  metadata;
- changes to protected-localized `Vee`, residual selection, injection policy,
  or mapping defaults.

## Required Diagnostics

For each system and cloud, report:

- cloud definition, including source labels, exponents/widths if available,
  occupancy, and total charge;
- `q0` projected charge, projection loss, and range;
- `uN_IDA` final-weight sanity;
- `E_self_IDA = 0.5 * q0' * V * q0`;
- `E_self_G = 0.5 * (rho0 | rho0)_G`;
- `C_screen = E_self_IDA - E_self_G`;
- `W_IDA`, `W_G`, and `Delta_W` finite/symmetry checks;
- `Delta_W` diagonal, Frobenius, and spectral range;
- `Delta_W` locality/radius bins and row-class decomposition;
- low one-body spectra and low-mode sector makeup;
- occupied/reference orbital expectations of `Vnuc_G`, `J0_G`, `uN_IDA`,
  `V * q0`, and `Delta_W`;
- comparison against the previous `SCREENED-VNUC-IDA-PROMISE-AUDIT-02`
  numbers.

## Decision Rule

Continue only if:

- compact-GTO clouds remain well represented;
- `Delta_W` stays local/core-like;
- low one-body shifts remain modest;
- far bins decay;
- Cr atom, if run, does not introduce broad or valence-destabilizing modes.

Stop and report if:

- `q0` representation becomes poor;
- `Delta_W` becomes nonlocal or large in valence/far bins;
- low modes become broad/protected/residual occupation incentives;
- conclusions depend on correcting bare `Vnuc` or `J0` separately rather than
  the screened sum.

## Validation

Approved validation for this measurement lane:

- package load;
- ignored probe run;
- `git diff --check`;
- final git status.

A later source lane may be requested only after this measurement lane identifies
a stable cloud/charge convention.
