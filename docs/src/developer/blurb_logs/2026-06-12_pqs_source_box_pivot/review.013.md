Accepted as a wiring/readiness probe.

Pass 013 built a one-center retained-source H1 diagnostic from existing CPBM
pieces:

```text
S = retained source overlap
T = retained source kinetic
V_unit = retained centered uncharged nuclear-by-center block
H_probe = T + Z * V_unit
```

The probe intentionally used a synthetic 3 x 3 x 3 self-pair with identity
source-axis transforms. It proves the retained-source blocks compose coherently
and that the charge application remains outside the source nuclear block. It is
not a physical PQS acceptance result.

Manager validation:

- `julia --project=. tmp/work/pqs_retained_source_h1_readiness_probe.jl`
- `git diff --check`

Important result:

The retained-source overlap was identity to roundoff for this synthetic
identity-transform fixture, so the diagnostic solve was ordinary rather than
generalized. The positive H1 value is only a wiring diagnostic and should not
be interpreted physically.

Next target:

Move from the synthetic 3 x 3 x 3 identity fixture to the real raw-product
source-box vocabulary: `CartesianRawProductSources.raw_product_box_plan(...)`
with `source_mode_dims = (5, 5, 5)` and
`pqs_boundary_product_mode_retained_rule(...)` retaining 98 boundary modes.
If only identity source-axis transforms are available, label that as a
diagnostic identity source-mode path and report the missing real PQS
source-axis transform builder as the next blocker.

-- repo-manager@macmini
