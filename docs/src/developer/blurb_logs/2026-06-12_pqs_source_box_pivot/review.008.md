Accepted.

Pass 008 correctly blocked on:

```julia
:missing_pqs_source_axis_gaussian_factor_source
```

The audit found that the active PQS source-pair record currently carries
source dimensions, ordering, retained rules, and summaries only. The raw
product source layer has `AxisSourceTransformFact`, but the default facts are
explicitly `:not_materialized` and carry no coefficient matrices or source-axis
representation data. Implementing Gaussian factor arrays from dimensions alone
would have been the wrong move.

The next contract gap is not the Gaussian contraction itself. It is the
plumbing for externally materialized source-axis transform/factor data:

```text
RawProductBoxPlan axis facts
-> retained-unit transform contract metadata
-> pair-block preflight record
-> later Gaussian factor source helper
```

No production code changed in this pass.

Manager validation:

- response reviewed against the live files;
- no source/test changes were present.

Next target:

Make the existing `AxisSourceTransformFact` placeholder usable for externally
materialized source-axis transforms, and plumb it to the PQS pair-block
preflight record as an explicit object or precise blocker. Do not implement
Gaussian factors yet.

-- repo-manager@macmini
