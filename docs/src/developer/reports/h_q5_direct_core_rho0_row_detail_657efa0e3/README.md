# H q5 Direct-Core Rho0 Row Detail

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
- Earlier full-`M` row-action checks mixed terminal weights with residual/MWG
  proxy weights. This packet avoids that category problem by using a tiny
  base-terminal H direct-core case only.

## Plain Result

The rho0/Galerkin blocker has been narrowed to a terminal direct-core IDA
proxy convention mismatch.

For H with `q=5`, `core_spacing=0.3`, and a one-Gaussian rho0 with
`alpha=8.0`, the direct-core terminal support weights already match the
authoritative parent product weights exactly. So the row-action mismatch is
not caused by using the wrong terminal weight vector.

The important split is:

```text
(J*w)/w    behaves like the point potential at the row center
u_direct   behaves like an equal-width Gaussian-smoothed potential
```

So the direct `u0` construction may be internally consistent with a
Gaussian-smoothed IDA density proxy, while the row-action diagnostic is asking
for point/constant-function behavior.

## Static-Audit Question

Identify the intended terminal direct-core IDA convention for rho0/Galerkin:

- Should direct-core `u0` represent the point/constant-function row-action
  object implied by `(J*w)/w`?
- Or should direct-core `u0` represent the Gaussian-smoothed density proxy
  currently produced by the MWG/direct-`u0` construction?
- If the smoothed proxy is correct, what is the right row-gauge diagnostic
  for rho0/Galerkin?
- If the point/constant-function proxy is correct, which source path is using
  the wrong density/proxy width for direct-core rows?

Do not infer Cr2 rho0/Galerkin physics until this convention is explicit.

## Files

- `h_q5_direct_core_rho0_row_detail.jl` - ignored probe copied from
  `tmp/work`.
- `summary.tsv` - top-level scalar diagnostics and conclusion.
- `direct_core_rows.tsv` - one row per direct-core terminal row, including
  point-potential and equal-width-smoothed comparisons.

The runtime PID file from `/Users/srw/dmrgtmp` is intentionally omitted.

## Key Evidence

From `summary.tsv`:

```text
alpha                         8.0
dimension                     517
direct_core_count             125
terminal_row_rel              0.035806146780649604
parent_product_row_rel        0.035806146780649604
parent_terminal_weight_maxabs 0.0
direct_core_residual_sq_fraction 0.9945083316193085
row_action_point_rel          0.001593816435896428
u_equal_width_average_rel     0.0003757432865293946
worst_row                     63
worst_state                   (7, 7, 7)
```

For worst row `63`, the row is at the direct-core center:

```text
terminal_support_weight       0.1591533792227156
parent_product_weight         0.1591533792227156
u_direct                      2.2566297073515225
row_action                    3.2019488345557066
analytic_point_potential      3.1915382432114616
analytic_equal_width_average  2.256758334191025
row_action_minus_point        0.010410591344244935
u_direct_minus_equal_width    -0.00012862683950265463
```

This row alone accounts for about `18%` of the full residual-square error,
and direct-core rows account for about `99.45%`.

## Interpretation

This packet supersedes the overly broad interpretation that direct `u0` is
simply in the wrong gauge. The better current interpretation is narrower:
there are two plausible row/proxy objects in play, and the diagnostic has not
yet established which one is the intended IDA convention for direct-core
rho0/Galerkin.

The earlier Cr2 rho0/Galerkin HF result should remain invalid as physics, but
this H direct-core result does not reject rho0/Galerkin as an idea.

## Guardrails

This report does not approve:

- source edits;
- public API, driver, artifact, provenance, writer, or reader changes;
- Cr2 production claims;
- more HF with the current rho0/Galerkin correction before the convention is
  explicit;
- `C' V C` interaction revival.

The next repo step should be a static/code audit of direct-core IDA/MWG row
proxy conventions, followed by another cheap row-gauge probe only after the
intended convention is stated.
