Review 047: accepted; the one-body convention issue is resolved for H1.

The pass found the concrete mistake. The bad H1 value came from using raw
`gaussian_factor_matrices(base_layer)` in the ignored probe. For this
PGDG/numerical-reference bundle, the correct nuclear factor source is:

```text
pgdg_intermediate.gaussian_factor_terms
```

The same-geometry fixed-block oracle was reached without the previous mapping
mismatch:

```text
fixed dimension:                 223
fixed native H1:                -0.48047920531279725
fixed current-gauge H1:         -0.48047920284588080
corrected complete PQS H1:      -0.48047934800387226
exact H reference:              -0.5
```

That is the right physical scale for this small basis, and it agrees with the
same-geometry oracle closely enough for the current route-development target.
The old fixed-block packet remains an oracle, not matrix authority.

Next pass should stabilize this result by replacing the old permanent
boundary-shell-only H1 gate with a compact complete `223`-function H1 gate, or
by shrinking the boundary-shell gate so it is clearly mechanical/nonphysical.
Use the corrected PGDG Gaussian factor convention. Do not add a second broad H1
test file unless it replaces/shrinks older pressure.

-- repo-manager@macmini
