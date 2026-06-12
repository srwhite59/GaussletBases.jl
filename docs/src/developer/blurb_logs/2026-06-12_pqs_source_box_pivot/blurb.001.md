Purpose:

Start the PQS pivot with a read-only restart audit. The goal is to identify the
first source-box-first implementation target without restarting the old broad
framework habit.

Why now:

The Be/GTO compile loop is paused with known nonblocking debt. The next physics
direction is PQS, but the route must start from the PQS source-box contract:

```text
raw product source box
-> retained source-mode rule/transform
-> source-space operator block
-> retained block by T_left' * O_source * T_right
```

Shell realization is oracle/adapter only. It must not become the route
algorithm.

Exact task:

Perform a read-only audit and write the first PQS restart status report. Do not
change production source, tests, docs, or examples in this pass except for the
required response file.

Read at least:

```text
AGENTS.md
docs/src/developer/root_live_baton_loop.md
docs/src/developer/pqs_source_box_operator_framework.md
docs/src/developer/raw_product_source_retained_transform_policy.md
docs/src/developer/projected_q_shell_policy.md
docs/src/developer/cartesian_route_retirement_ledger.md
```

Inspect code/tests enough to name the live surfaces for:

```text
CartesianRawProductSources
CartesianPairBlockMaterialization PQS/PQS raw source-space safe one-body blocks
current shell-realization bridge/readiness layers
raw product source box plan tests
PQS source-pair retained/final readiness tests
```

Audit questions:

1. What exact source files/functions already own raw product source boxes?
2. What exact source files/functions already materialize PQS/PQS raw
   source-space overlap/kinetic/position/x2 blocks?
3. What shell-realization or support-row surfaces exist today, and how should
   they be classified as oracle/adapter/debug rather than route authority?
4. What is the smallest first retained-rule object needed for a one-center
   source-box retained transform?
5. Can the first retained rule be the explicit boundary product selector:

   ```text
   source modes 5 x 5 x 5
   retained boundary modes 5^3 - 3^3 = 98
   ```

6. What exact implementation seam should the next blurb use for
   `PQSRetainedRule` / `PQSRetainedTransform`, if that remains the right next
   step?
7. What old code/test/docs are likely to become oracle-only or removable once
   retained source-mode block contraction works?

Trust boundary:

Read-only audit. No source edits, test edits, production docs edits, public API
changes, exports, shell realization implementation, Lowdin, IDA,
electron-nuclear, Hamiltonian, RHF, artifacts, PQS driver adoption, full-parent
fallback, or direct Cartesian fallback.

Test policy:

Do not add or run broad tests. A load check is optional. Use `rg`, file reads,
and existing docs/source inspection. If a tiny local probe is useful, put it
under `tmp/work` and explain why.

Validation:

- `git status --short --branch`
- `git diff --check`
- optional: `julia --project=. -e 'using GaussletBases; println("load ok")'`

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path is likely to become
  unnecessary after the first retained source-mode transform;
- what should stay as oracle/adapter/debug only;
- if nothing can be deleted yet, why this audit does not obsolete an existing
  surface.

Report back:

- write `.agent_handoffs/response.001.md.tmp`, then atomically rename to
  `.agent_handoffs/response.001.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.001.md`;
- include files/functions inspected;
- answer the audit questions;
- recommend the exact next implementation blurb;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
