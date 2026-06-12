Purpose:

Run a probe-first one-center PQS retained-source one-electron readiness check.

Why now:

The source-box pieces for one-body PQS are now connected:

```text
retained overlap block
retained kinetic block
retained centered uncharged nuclear-by-center block
```

Before adding a global route or final-basis machinery, check that those pieces
can assemble a coherent one-center retained-source H1 diagnostic.

Exact task:

Create a `tmp/work` probe, not a new long-term test, unless a truly missing
public seam is exposed.

Use existing CPBM helpers where possible:

- `pqs_source_pair_retained_overlap_block(...)`
- `pqs_source_pair_retained_kinetic_block(...)`
- `pqs_source_pair_retained_centered_electron_nuclear_by_center_block(...)`
- existing retained source-mode matrix helpers if they fit the one-self-pair
  path

The probe should:

- build one PQS retained self-pair with materialized source-axis transforms;
- use explicit analytic axis layers and overlap/kinetic factors;
- build `S`, `T`, and one uncharged by-center `V_unit`;
- form `H = T + Z * V_unit` in the probe only, as the Hamiltonian-stage charge
  application check;
- check finite/symmetric matrices and overlap positive definiteness if a solve
  is attempted;
- if solving, label it `:raw_retained_source_generalized_diagnostic` unless the
  retained source overlap is already identity;
- report the exact final-basis blocker, expected to be shell realization /
  Lowdin final-basis construction, if the overlap is nonidentity.

Do not make the raw generalized diagnostic an acceptance result. Do not add
generalized final-basis transfer logic.

Trust boundary:

No production route unless the probe exposes a very small missing seam and the
fix is clearly less than the probe scaffolding it replaces. No shell
realization, Lowdin cleanup, IDA, density-density, RHF, global driver, exports,
artifacts, full-parent fallback, direct Cartesian fallback, or old CCPM wrapper
adoption.

Test policy:

Do not add a test file for this pass. Use `tmp/work` artifacts for the probe.
If production code changes are unavoidable, add only compact coverage in the
existing CPBM contract file and explain why the artifact/test earned its cost.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test/artifact was added and why it earned its cost;
- any remaining stale or duplicate surfaces to retire next.

Validation:

- if production code/tests changed:
  `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- always run:
  `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- write `.agent_handoffs/response.013.md.tmp`, then atomically rename to
  `.agent_handoffs/response.013.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.013.md`;
- include probe artifact path and key matrix/solve diagnostics;
- include implementation or exact blocker;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
