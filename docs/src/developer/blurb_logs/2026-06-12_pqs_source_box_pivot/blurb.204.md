Pass 204 - split `source_box_route_shadow.jl` into metric fallback, PQS pair-plan, and legacy fixture files

Role:
You are `repo-doer@macmini` implementing one bounded cleanup pass for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, and `BlurbStyle.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `ae1c4d99 Split contracted parent metrics core`
- `src/CartesianContractedParentMetrics.jl` is now a thin wrapper.
- The remaining hotspot is:
  `src/cartesian_contracted_parent_metrics/source_box_route_shadow.jl`
  currently about 8601 lines.

Goal:
Make `source_box_route_shadow.jl` stop being a mixed private junk drawer.
Split it into clearer behavior-preserving include files:

```text
src/cartesian_contracted_parent_metrics/product_staged_metric_fallbacks.jl
src/cartesian_contracted_parent_metrics/source_box_pair_shadow.jl
src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl
```

Keep `src/cartesian_contracted_parent_metrics/source_box_route_shadow.jl` only
as a very small compatibility include wrapper if that is the cleanest route, or
delete it if the parent wrapper includes the three new files directly. Do not
leave substantial implementation there.

Classification intent:

1. `product_staged_metric_fallbacks.jl`
   Keep code still needed by contracted-parent metric packet construction and
   resolved/staged metric fallback paths. Expected examples include, subject to
   actual dependency order:
   - staged axis utilities;
   - product-staged metric projection/fill helpers;
   - support-local retained entries;
   - staged unit entries;
   - resolved-payload low-order metric blocks;
   - fallback separable metric block helpers;
   - helpers directly used by `cartesian_contracted_parent_metric_packet(...)`
     or `cartesian_contracted_parent_metric_packet_dense_reference(...)`.

2. `source_box_pair_shadow.jl`
   Keep the potentially future-aligned PQS source-box pair-plan prototype
   machinery. Expected examples include:
   - `_pqs_raw_product_box_*` structural/operator plans;
   - PQS retained/boundary source-box selectors;
   - `_pqs_product_source_box_pair_plan`;
   - `_pqs_pqs_source_box_pair_plan`;
   - source-box block-from-1D-factor helpers that are part of the pair-plan
     prototype path.

3. `legacy_source_box_fixtures.jl`
   Put old private fixture/oracle route shadows here. Expected examples include:
   - product/doside source-box shadow/reference block family;
   - product/doside density-density, local Gaussian, and nuclear fixture blocks;
   - homonuclear raw-box geometry fixtures;
   - product slab/source-box fixture producers;
   - contact-cap, outer-mismatch, atom-box, and support-local oracle comparisons;
   - route-retained-unit fact audits that are not the metric fallback core and
     not the PQS pair-plan prototype.

Important:
- This is a classification split, not an algorithm pass.
- Preserve behavior and names.
- Do not add tests.
- Do not add public API.
- Do not touch H1/H1-J/RHF/PQS multilayer physics paths.
- Do not move code into public route modules.
- Do not delete PQS raw-product-box or source-pair plan prototype code in this
  pass.
- Do not delete metric-packet fallback code in this pass.

Deletion pressure:
- The active line-count rule still applies.
- For edits under `src`, `test`, `bin`, and the CR2 generator script, final
  tracked diff must be net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Acceptance condition:
`sum(deleted) > sum(added)` for those scoped files.

For this pass, a pure split may be line-neutral unless you also shrink trivial
wrapper/comment/blank-line overhead. Keep the split behavior-preserving, but do
not leave the pass net-positive. If you cannot satisfy the line budget without
deleting live endpoint/reference tests or changing behavior, write
`.agent_handoffs/ATTENTION.md` and stop.

Deletion-candidate audit:
Do not delete behavior yet, but report deletion candidates with caller status.
In particular audit callers of:

```text
_product_doside_source_box_pair_plan
_product_doside_source_box_reference_block
_product_doside_source_box_shadow_blocks
_product_doside_source_box_local_gaussian_sum_block
_product_doside_source_box_density_density_interaction_block
_product_doside_source_box_nuclear_attraction_by_center
```

Report whether these are called only by legacy fixture tests/probes, or by live
source/driver paths.

Validation:
Run the smallest checks that cover the split boundary:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. -e 'using Test, LinearAlgebra, SparseArrays, GaussletBases; include("test/nested/cartesian_contracted_parent_metric_packet_runtests.jl")'
```

Also run one focused legacy/private test if your caller audit shows the moved
legacy fixture or PQS pair-plan code has a short existing test. If the only
available test is expected to run longer than 60 seconds, do not run it by
default; explain why and report the exact test.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.204.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.204.md
```

Report:
- files created/removed;
- line count of each resulting include;
- which function families went into each include;
- source/test/bin scoped line budget added/deleted/net;
- deletion candidates and caller status;
- validation commands and results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Stop after writing the response. Manager will review, commit, push, and pause.

-- repo-manager@macmini
