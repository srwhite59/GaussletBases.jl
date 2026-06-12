Review 037: accepted.

The direct retained-boundary algorithm shape now covers the one-body PQS
source terms needed for the current H1 seam:

- overlap;
- kinetic;
- by-center electron-nuclear.

The active retained nuclear wrappers now fill retained blocks directly from
retained source-mode tuples and Gaussian factor terms. The raw source-result
selector path remains as oracle/fallback over already materialized source
blocks.

Independent manager validation:

```text
julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

The CPBM contract passed, and the direct-vs-oracle comparisons are exact for
the covered fixture.

Next step: do a route-readiness probe for the explicit PQS final-basis H1 seam
using the direct retained overlap/kinetic/nuclear blocks. Do not immediately add
another permanent test. First confirm which old oracle/helper pressure the H1
path can replace or shrink.

-- repo-manager@macmini
