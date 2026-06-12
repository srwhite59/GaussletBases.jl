Purpose:

Identify or implement the PQS source-axis Gaussian factor source needed by the
new electron-nuclear block.

Why now:

Pass 007 implemented the 3D source-space electron-nuclear contraction when the
caller supplies term-first Gaussian factor arrays. The remaining missing piece
before an H/He+ one-electron diagnostic is to produce those arrays from the
actual PQS source-axis data carried by the source-box route.

Exact task:

Audit the current PQS source-box/retained-transform objects and determine
whether a CPBM helper can compute, for one source-pair record, term-first
Gaussian factor arrays:

```text
(nterms, left_axis_count, right_axis_count)
```

for x/y/z and a center location.

If the needed source-axis data exists, implement a narrow helper that produces
those arrays and feeds the pass-007 nuclear block.

If the route only carries dimensions/order and not enough axis representation
data, do not fake it. Report the exact blocker, preferably:

```text
:missing_pqs_source_axis_gaussian_factor_source
```

Code surfaces to inspect first:

```text
src/CartesianRawProductSources.jl
src/cartesian_retained_unit_transform_contracts.jl
src/cartesian_pair_block_materialization/preflight.jl
src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl
src/CartesianParentAxisFactors.jl
src/CartesianContractedParentMetrics.jl
src/CartesianCPBBlockProviders.jl
docs/src/developer/pqs_source_box_operator_framework.md
```

Implementation boundary if source-axis data exists:

- add only a narrow source-axis Gaussian factor helper;
- reuse existing analytic Gaussian factor conventions where possible;
- do not route through CCPM `_pqs_pqs_source_box_*nuclear_attraction*`
  wrappers;
- do not apply nuclear charge;
- do not assemble H1;
- do not add shell realization, Lowdin, IDA, Hamiltonian, driver adoption,
  exports, or artifacts.

Test policy:

If implemented, add one compact check in the existing PQS source-pair contract
test or a `tmp/work` probe if the shape is still exploratory. The check should
compare generated arrays to a direct existing Gaussian-factor helper for one
small source axis, not to CCPM physical nuclear wrappers.

If blocked, do not add tests. Just report the blocker and the missing data.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- if production code/tests changed, also run
  `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`

Report back:

- write `.agent_handoffs/response.008.md.tmp`, then atomically rename to
  `.agent_handoffs/response.008.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.008.md`;
- include inspected surfaces;
- include whether implementation happened or the exact blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
