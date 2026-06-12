# Successor Handoff: PQS Source-Box Pivot

Date: 2026-06-12
Role: repo-manager handoff for the GaussletBases PQS source-box pivot
Baseline implementation HEAD: `4090eae5 Connect PQS driver final H1 seam`

## Repository State

At handoff creation, the PQS implementation is at `4090eae5`, and `main` is
even with `origin/main` before this document is added. The last
committed/pushed implementation pass is:

```text
4090eae5 Connect PQS driver final H1 seam
```

Recent commits leading into the current state:

```text
39abc093 Build PQS driver source plan payload
9c013766 Select PQS shellification in driver shell stage
4a467ccb Log PQS H1 J driver input boundary
54f4f180 Add PQS H1 J driver diagnostic slot
0e4b2e43 Document PQS H1 J driver seam
cce93ab2 Add PQS multilayer H1 J payload
ed3fce43 Add PQS multilayer support density inputs
```

If this file has not yet been committed, the worktree should be dirty only
because of this handoff unless another agent or user change appears. If this
file has been committed, the worktree should be clean. Do not assume unrelated
dirty files are yours; inspect before editing.

## Baton Loop Status

The live baton loop is paused by explicit user request. Do not resume polling,
publish a new blurb, or start doer work until the user asks.

Relevant live/local files:

```text
.agent_handoffs/state.md
.agent_handoffs/blurb.095.md
.agent_handoffs/response.095.md
.agent_handoffs/review.095.md
```

Relevant tracked log directory:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/
```

The tracked pass-095 review is:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.095.md
```

If a future unattended baton loop is restarted, remember the current rules:

- the user must first paste/start the doer with clear startup instructions;
- avoid ambiguous wording such as "for pass N only" that might make doer stop
  after one pass;
- no UI escalation in unattended mode;
- if doer needs permission, it should write `.agent_handoffs/ATTENTION.md` and
  stop;
- manager polling should wait about one hour with roughly one-minute checks
  unless the user changes the rule;
- the old temporary ten-minute wait was not a standing rule and should not be
  cited as one.

## Current Technical State

The PQS source-box work has moved from raw metadata scaffolding to a real
driver-facing diagnostic seam. The current shape is:

```text
shellification/lowering-backed region plan
-> PQS multilayer source plan
-> complete core/shell final basis
-> complete core/shell H1 payload
-> compact H1/J diagnostic slot in cartesian_assembly/cartesian_report
```

Pass 095 connected the one-center `:pqs_source_box` driver path to the
complete core/shell final basis and H1 payload. The H1/J diagnostic now carries
the final dimension and H1 value while remaining blocked only on the density/J
inputs:

```text
:axis_weights
:raw_pair_factor_terms
```

The focused dry-run reported:

```text
final_dimension = 223
h1_energy = -5.6629907690725245
h1j_status = blocked on axis weights and raw pair terms
driver_route_materialized = false
```

This is not a production PQS route yet. It is an internal diagnostic route
seam. It makes no RHF, SCF, Fock, density iteration, GTO, export, artifact,
fixture promotion, explicit-box authority, or fixed-block authority claim.

Validation accepted for pass 095:

```text
doer focused one-center PQS source-box dry-run
doer test/nested/pqs_direct_retained_final_h1_runtests.jl
manager git diff --check
manager /Users/srw/Dropbox/codexhome/cjulia -e 'using GaussletBases; println("load ok")'
```

Plain `julia --project=.` hit the known Juliaup lockfile permission problem,
not a code failure. Prefer `cjulia` for routine manager load checks in this
environment.

## Recent Design Decisions

The driver spine is now the controlling architectural pressure. PQS helpers
should be pulled into:

```text
cartesian_system / cartesian_recipe
cartesian_parent
cartesian_shells
cartesian_units
cartesian_transforms
cartesian_pairs
cartesian_assembly
cartesian_report / cartesian_materialization
```

Do not grow a richer private PQS route while the driver remains a thin report
wrapper. The promotion path should be:

```text
private probe/oracle
-> module-owned machinery
-> driver-stage consumption
-> public route behavior
-> old code/test shrinkage
```

Shellification/lowering is the active authority for complete core/shell PQS
geometry. The explicit-box `pqs_multilayer_shell_source_plan(bundles, core_box,
outer_box; ...)` entry point is compatibility/probe bridge material only. New
route work should consume the shellification/lowering-backed region-plan entry
point.

`CartesianFinalBasisRealization` was introduced to keep final-basis realization
and final operator transfer out of `CartesianPairBlockMaterialization`. Do not
turn it into a broad PQS everything-module. Its appropriate scope is final
basis, final one-body transfer, narrow H1 seams, and narrow density diagnostics
that consume already-defined route inputs.

The support-space dense one-body helpers are acceptable for the current H1
diagnostic seam, but should not become the general PQS operator algorithm
without scale review. Direct retained-boundary kernels are the preferred
source-space operator shape where they apply.

The current H1/J density convention is the localized pre-final positive-weight
gauge. Do not substitute signed-final-weight division, raw no-division density,
or old density-normalized pair-factor authority without a reviewed convention
change.

## Current Auditor Inputs

The latest ChatGPT/auditor guidance accepted the recent direction but raised
three recurring warnings:

- `pqs_source_box_route_driver_helpers.jl` is still a large compatibility and
  report-field surface; avoid letting it become the owner of complete
  core/shell route logic.
- The ordinary one-center PQS driver dry-run now performs real final-basis and
  H1 construction, which explains longer test times. That should be explicitly
  requested or represented by a compact route object; it should not become an
  always-on cheap dry-run by accident.
- Recent testing is not badly bloated, but the focused H1 gate is near its
  limit and the large CPBM contract test remains inherited bloat pressure.

User guidance since pass 095:

- do not resume the baton loop while reviewing/writing this handoff;
- when asked questions during a future loop, do not stop automatically unless
  the user explicitly says to pause;
- use smaller tests while writing code, then a modest representative check once
  the seam works;
- avoid broad slow tests as routine validation;
- escalation requests stop unattended progress, so unattended blurbs must
  include the no-UI-escalation rule.

## Test And Runtime Policy

For the next PQS passes:

- use tiny smoke probes while editing driver seams;
- after the seam works, run one modest representative check;
- avoid broad route-driver/report/materialization harnesses unless the change
  directly touches those surfaces and the cost is explained;
- do not add a new permanent PQS test unless it replaces or shrinks older
  oracle pressure, protects a compact module contract, or guards a real
  scientific/workflow endpoint;
- do not grow `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  as a PQS development notebook.

## Intended Next Tasks

### 1. Design Cleanup Before More Physics

Make the complete core/shell final-basis and H1 construction explicitly
diagnostic or route-object-owned. The likely target is a compact internal route
object or payload that `pqs_source_box_route_driver_helpers.jl` calls and then
summarizes, rather than continuing to add local route logic and scalar report
fields in that file.

Success should be measured by reduced driver-helper sprawl and clearer cost
control. Do not add RHF or new physics in this cleanup.

### 2. Feed H1/J Density Inputs In A Narrow Pass

Once the diagnostic route object boundary is clear, provide the remaining
route-owned inputs for the existing H1/J diagnostic:

```text
axis/support weights
raw support pair numerator terms
```

The result should advance the driver H1/J payload from blocked-on-density-inputs
to materialized diagnostic status. It should still not add RHF, SCF, fixture
promotion, exports, or artifacts.

### 3. Review Fixture And Acceptance Policy Before RHF

Only after H1/J is driver-owned should RHF be considered. Before that, review
the physical fixture rule: box radius, central spacing `d`, distortion `s`,
`q`, and shell depth must move together. The side-13 and q-ladder results are
useful route/scaling evidence, not accepted convergence gates.

## Non-Goals For The Immediate Successor

Do not do these next unless the user explicitly changes direction:

```text
resume baton polling
publish blurb.096
add RHF/SCF/Fock iteration
add GTO work
promote q=9/q=11/side13 as an acceptance fixture
turn shell-support projection into the production operator path
add exports/artifacts
expand broad tests
rewrite all NamedTuples into structs
create a broad PQS everything-module
```

## Signoff

Prepared as a successor handoff by `repo-manager@macmini`.
