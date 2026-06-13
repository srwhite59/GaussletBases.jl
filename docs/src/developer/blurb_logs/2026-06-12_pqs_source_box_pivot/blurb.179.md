Pass 179 - implement fixed-q multi-shell PQS source modes for the one-center He target.

Read this whole blurb before editing. The important correction from pass 178 is:

```text
The compact PQS 223 path is a seam/readiness route, not the physics target.
The first physics-facing PQS target is one-center He, q=5/n_s=5, matching the
old WL Fig.8-style fixed-shell inventory as closely as the current route can.
```

Target for this pass:

```text
He one-center source/final-basis inventory only
parent side: 11
core: fixed 5^3 = 125 support rows
shell layers: three complete shells around that core
PQS source mode shape per shell: fixed (5, 5, 5)
retained shell modes per layer: 5^3 - 3^3 = 98
gausslet-only final retained dimension: 125 + 3*98 = 419
AHGBS residuals: not in this pass
H1/H1-J/RHF/DMRG/HFDMRG/CR2: not in this pass
```

This is not Be2 or Cr2 work. Do not touch the CR2 artifact generator, CR2
handoff reports, Be2 routes, diatomic route files, HFDMRG, RHF, DMRG, exports,
artifacts, or public API.

Current suspected seam:

- `src/cartesian_terminal_lowering/region_contracts.jl`
  - `_pqs_complete_shell_contract(region, policy::PQSLowering)` currently records:

    ```julia
    q = policy.q
    source_mode_shape = CartesianCPB.shape(source)
    ```

  - Be careful: `CartesianCPB.shape(source)` is the support/source CPB box
    shape, which grows with shell extent. For fixed-q PQS retained modes, the
    source-mode shape should be the product-mode shape controlled by
    `policy.q`, i.e. `(q, q, q)` for this one-center complete-shell path.

- `src/pqs_multilayer_shell_region_plan.jl`
  - `PQSMultilayerShellLayerRegion` already carries `source_mode_shape`.
  - `_pqs_multilayer_contract_source_mode_shape(contract, raw_region)` already
    reads `contract.metadata.source_mode_shape`.

- `src/pqs_multilayer_shell_source_plan.jl`
  - `_pqs_multilayer_region_plan_layer_specs(region_plan)` currently does not
    pass `layer.source_mode_shape` through to the source-realization spec.
  - `_pqs_multilayer_realize_shell_source_plan(...)` currently derives:

    ```julia
    q_values = length.(inner_box)
    q = q_values[1]
    raw_source_dims = (q, q, q)
    selected_q = q
    ```

    That is the growing-q/q-ladder behavior. It must remain valid for the
    explicit-box compatibility bridge, but the shellification/lowering-backed
    entry point should consume the route-owned source-mode shape from the
    lowering contract.

Implementation intent:

1. For shellification/lowering-backed PQS source plans, carry
   `source_mode_shape` from each `PQSMultilayerShellLayerRegion` into the layer
   spec and use that as `raw_source_dims`.

2. For this first pass, support the fixed cubic source mode shape `(q, q, q)`
   only. If the shape is rectangular/non-cubic, block/throw with a clear
   label/message rather than inventing the rectangular policy here. The user
   explicitly noted that rectangular shells have a subtle long-side rule, but
   the atom-first target is cubic.

3. Preserve the explicit-box bridge behavior.
   `pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)` can
   still derive `q` from `inner_box` as before. Do not make explicit boxes
   route authority.

4. Record compact provenance in summaries/metadata if useful:

   ```text
   source_mode_shape_source = :terminal_lowering_contract
   fixed_source_mode_shape_used = true
   explicit_box_bridge still false for shellification-backed route
   ```

   Keep this compact. Do not create a new payload object or a field cloud.

5. Add/update the smallest focused test in
   `test/nested/pqs_direct_retained_final_h1_runtests.jl` or a nearby existing
   PQS source/final-basis test. Prefer modifying/shrinking the existing file
   over adding a new default test file.

   The test must prove the new one-center fixed-q inventory:

   ```text
   parent side/counts = (11, 11, 11)
   core_support_count = 125
   shell_layer_count = 3
   shell_support_count = (7^3-5^3) + (9^3-7^3) + (11^3-9^3) = 1206
   each shell record raw_source_dims == (5, 5, 5)
   each shell record retained_count == 98
   shell_final_retained_count = 294
   final_retained_count = 419
   final overlap identity error remains small
   ```

   Do not assert H1/J/RHF/physics energy in this pass unless it is already
   available without broadening the work. The goal here is the correct retained
   basis inventory.

Line-budget rule for this pass:

```text
git diff --numstat -- src test
```

must be net-negative:

```text
sum(deleted) > sum(added)
```

Docs/handoff logs do not count. Do not game this by deleting accepted
scientific endpoint checks. Good deletion candidates inside
`test/nested/pqs_direct_retained_final_h1_runtests.jl`:

- remove or shrink the explicit-box duplicate bridge construction/comparison
  (`explicit_box_plan`, `explicit_box_final_basis`, and equality assertions
  against the shellification-backed plan) if the new fixed-q
  shellification-backed source plan is the accepted route path;
- remove repeated nonclaim flag assertions that do not protect this seam;
- consider removing the support-nuclear axis-layer cross-check from this
  focused test if it is only preserving a now-stable helper convention.

Do not go hunting through unrelated old tests just to buy line budget. Use
directly related stale scaffold first. If this pass cannot be net-negative
while preserving the live contract, write `.agent_handoffs/ATTENTION.md`
explaining the exact blocker and stop.

Validation:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

If the focused test becomes much slower than before, report the elapsed time
and whether the new cost is source-plan/final-basis construction or an existing
H1 path. Do not add a broad suite.

Reporting requirements:

- exact files changed;
- observed fixed-q source/final-basis inventory;
- whether explicit-box bridge behavior was preserved;
- source/test line additions, deletions, and net;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Write the result to `.agent_handoffs/response.179.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.179.md
```

-- repo-manager@macmini
