# Rho0 Row-Gauge Algebra Audit

Status: review packet for external/static analysis. This is not production
authority and not a source implementation proposal.

Authority context:

- `HP-RG-RHO0-GAL-AUDIT-01` allows ignored measurement probes for a
  rho0/Galerkin IDA correction over the protected-localized inherited-site
  baseline.
- The direct `C' V C` protected interaction transform is invalid and must not
  be reused.
- The current viable protected interaction baseline is localized injection
  with exact one-body operators and inherited pre-injection site-order IDA.

## Plain Result

The cheap Be2 row-gauge probe shows that the simple direct-`u0`
rho0/Galerkin correction is not ready for physical interpretation. The failure
appears before Cr2 and before HF.

The key same-basis check is:

```text
||(J_M - Diag(u0_M)) * w_M|| / ||J_M * w_M|| = 0.5626603070620595
```

That means the directly constructed `u0_M` is not in the same final IDA row
gauge as the Galerkin potential matrix `J_M` and the probe's reconstructed row
weights `w_M`. The issue is not primarily protected localization: the
localized `L` check using the same inherited `u0_M` has essentially the same
relative error, `0.561317734511763`.

Defining

```text
u_from_Jw = (J_M * w_M) ./ w_M
```

makes the row action zero by construction. That is useful as a diagnostic, but
it does not identify the canonical row-gauge object that should be used by the
rho0/Galerkin correction.

## Static-Audit Question

Find the correct final IDA row-gauge convention for the rho0/Galerkin
correction.

Specifically:

- Is `u0_M` being computed with the wrong density proxy, weights, or basis
  rows?
- Is `w_M` in the probe not the authoritative final IDA weight vector?
- Should the correction use a row-action object such as `(J*w)./w` instead of
  the direct proxy value?
- Is the current `J_M` built from the same raw pair numerator / final-weight
  convention as `electron_electron_ida`?
- Is there a missing constant, sign, or normalization convention in the
  intended rho0/Galerkin correction?

Do not infer Cr2 physics from the earlier rho0/Galerkin HF replay until this
row-gauge mismatch is resolved.

## Files

- `be2_rho0_row_gauge_algebra_probe.jl` - ignored probe copied verbatim from
  `tmp/work`.
- `summary.tsv` - top-level Be2 geometry and row-gauge summary.
- `row_gauge_checks.tsv` - main algebra checks for `M` and localized `L`.
- `weights.tsv` - weight summaries.
- `weight_compare.tsv` - `w_L_projected` versus `w_M`.
- `analytic_gaussian_checks.tsv` - normalized Gaussian charge/self-energy and
  center-potential checks.
- `linearity.tsv` - `N_screen` linearity checks.
- `occupied_orbital_shift_proxies.tsv` - non-physics shift proxies retained
  only for context.
- `stages.tsv` - timing/stage record.

## Key Evidence

From `summary.tsv`:

```text
case_label     Be2_row_gauge
alpha          8.0
base_dimension 431
compact_R      4
M_dimension    435
L_dimension    435
B_min          0.9525394022398224
B_lt_0p99      6
M_row_rel      0.5626603070620595
L_u0M_row_rel  0.561317734511763
L_projected_row_rel 1.1563262858691983e-16
weight_rel     0.001413654516807044
```

From `row_gauge_checks.tsv`:

```text
M_same_basis_direct_u0_M:
  row_rel = 0.5626603070620595
  row_abs = 395.3633410045879
  jw_norm = 702.6679082250965

M_same_basis_u_from_Jw:
  row_rel = 1.2281473375913627e-17

L_using_u0_M_and_w_M:
  row_rel = 0.561317734511763

L_using_u0_L_projected_and_w_L_projected:
  row_rel = 1.1563262858691983e-16
```

From `analytic_gaussian_checks.tsv`:

```text
charge error            -5.0959236830294685e-14
self-energy error       -1.9939605522267811e-13
center-potential error  -1.0505837000351903e-7
```

From `linearity.tsv`, `N_screen=0` is exactly zero and the correction is
linear in `N_screen` to roundoff. So the current suspect is row-gauge
construction, not gross normalization or nonlinearity.

## Interpretation

The earlier Cr2 rho0/Galerkin result should be considered invalid as physics.
The cheap Be2 probe already shows that the direct `u0_M` object does not match
the row action of the Galerkin potential in the same basis.

This packet is intended to let a reviewer inspect the probe code and tables
without depending on ignored `tmp/work` or machine-local `/Users/srw/dmrgtmp`
paths.

## Guardrails

This report does not approve:

- source edits;
- public API, driver, artifact, provenance, writer, or reader changes;
- Cr2 production claims;
- more HF with the current direct-`u0` rho0/Galerkin correction;
- `C' V C` interaction revival.

The next repo step should be a static/code audit of the IDA row-gauge
construction, followed by another cheap row-gauge probe only after a concrete
convention correction is identified.
