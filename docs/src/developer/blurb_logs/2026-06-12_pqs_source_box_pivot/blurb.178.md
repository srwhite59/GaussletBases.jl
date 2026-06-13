Pass 178 - audit WL-matched multi-shell He q=5 target for PQS

Role: repo-doer@macmini

Task type: no-edit physics/route audit. Do not edit source, tests, docs,
tracked generators, or artifacts in this pass.

Purpose:

Reset the PQS comparison target before more implementation. The recent Be2
WL/PQS handoff artifact exposed a conceptual mismatch, but the deeper issue is
that we should not be trying Be2 first. The next serious PQS target should be
an atom:

```text
He atom
q = 5 / n_s = 5 first
WL-matched geometry and shellification policy
multi-shell construction, not a growing-core one-shell ladder
```

The important user correction is:

- The old PQS q-ladder was not crazy, but it was not the desired construction.
  It grew the core cube and used one boundary shell:

  ```text
  q=5   5^3 core  + boundary(5^3)  = 125 + 98  = 223
  q=7   7^3 core  + boundary(7^3)  = 343 + 218 = 561
  q=9   9^3 core  + boundary(9^3)  = 729 + 386 = 1115
  q=11  11^3 core + boundary(11^3) = 1331 + 602 = 1933
  ```

- The desired PQS construction is multiple shells around a fixed small local
  core, in the same spirit as the old WL He work. Start with q=5 / n_s=5 and
  match WL before moving to Be2.

- The old WL/QW He work did a very good job before PQS started. The relevant
  benchmark is not the tiny active 223-line gate by itself, but the old
  Fig. 8-style He validation:

  ```text
  n_s = 5, d = 0.3, AHGBS-9 S-only:
    final dimension 447
    RHF total about -2.862102144533723
    about -0.558 mHa from the Fig. 8 row

  n_s = 7, d = 0.10, AHGBS-9 S-only:
    final dimension 1897
    RHF total -2.861673961528321
    error about +1.716e-6 Ha vs Fig. 8 row
    error about +6.034e-6 Ha vs He HF reference
  ```

The goal of this pass is to identify the exact WL q=5/n_s=5 shell inventory
and the smallest PQS branch that would match it. Do not implement that branch
yet.

Why now:

The current CR2 Be2 handoff is not the right next physics step:

- PQS `221` is a compact source-box spine:

  ```text
  left PQS boundary 98 + midpoint slab 25 + right PQS boundary 98
  ```

  It is not an atom-local Be2 shellification with atom cubes.

- WL `2287` is a separate atom-growth/materializer route and is not the
  intended matched comparison.

Before CR2 runs HF on PQS, we need a clean atomic PQS Hamiltonian target. The
first serious comparison should be He q=5, WL-matched.

Current state / known surfaces:

Start by reading these specific surfaces:

```text
docs/src/developer/numerical_contracts.md
docs/src/developer/pqs_source_box_fixture_policy.md
docs/src/developer/pqs_near_term_final_basis_realization_plan.md
docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/summary.md
docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.007.md
docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.008.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.053.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.056.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.058.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.059.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.060.md
```

Then inspect the code paths, but do not edit them:

```text
test/nested/cartesian_wl_gausslet_he_atom_acceptance_runtests.jl
src/CartesianContractedParentMetrics.jl
src/cartesian_nested_diatomic.jl
src/cartesian_nested_faces.jl
src/cartesian_shellization_route.jl
src/pqs_multilayer_complete_core_shell_h1.jl
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Useful search anchors:

```text
_white_lindsey_low_order_materialized_seed_report
white_lindsey_shellification_decomposed_unit_pair_inventory
build_one_center_atomic_full_parent_shell_sequence
one_center_atomic_full_parent_fixed_block
pqs_multilayer_complete_core_shell_final_basis
pqs_multilayer_complete_core_shell_h1_payload
pqs_multilayer_complete_core_shell_h1_j_payload
_nested_diatomic_projected_q_shell_retained_count
_nested_projected_q_shell_boundary_comx_product_modes
_nested_projected_q_shell_layer
```

Exact audit questions:

1. What exactly was the good old WL He construction?

   Report, for the old WL/QW or decomposed WL q=5/n_s=5 target:

   - parent axis count / basis count;
   - mapping and spacing policy (`d`, `tail_spacing`, `gscalefac` or equivalent);
   - whether AHGBS-9 S-only supplement is part of the target or a separate
     higher-quality comparison layer;
   - core/direct retained count;
   - shell retained counts and shell count;
   - residual/supplement count if applicable;
   - final dimension;
   - H1, J/Vee, and RHF values already recorded;
   - whether the route is current decomposed WL, old nested/QW oracle, or both.

2. What is the correct first PQS q=5 atom target?

   Define the smallest target in plain terms. It should not be the old q-ladder.
   It should answer whether the first PQS target is:

   ```text
   5^3 core + one 98 shell only
   ```

   or:

   ```text
   fixed 5^3 core + multiple q=5 shells matching the WL n_s=5 inventory
   ```

   or some other WL-matched retained-region inventory. If the old WL q=5
   target includes GTO residual/supplement directions, separate the gausslet-only
   PQS target from the supplemented comparison target.

3. Can current PQS code express that target?

   Decide whether current code already has a route-owned way to build a
   multi-shell q=5 PQS He source plan/final basis/H1/Vee path, or whether it
   only has:

   - the old complete core/shell one-shell q-ladder;
   - the compact diatomic source-box spine;
   - source-box local pieces without a full atom multi-shell route;
   - probe-only or tmp/work code.

   Name the exact functions/files that are closest. If no proper path exists,
   say that directly.

4. What is the smallest next implementation seam?

   Do not code it. State the first safe implementation pass in one paragraph:

   - what object/helper would be added or adapted;
   - what old scaffold it should replace or make unnecessary;
   - what test/probe would validate it;
   - what line-deletion candidate can pay for the new surface.

Trust boundary:

- No code changes.
- No tracked docs changes.
- No generator changes.
- No artifact generation.
- No CR2, HFDMRG, HF, RHF, DMRG, Be2, Cr2, exports, HamV6, V6, Vblocks, or
  public API.
- Do not run broad tests.
- Do not run any command expected to exceed 60 seconds.
- If you need a local Julia probe, keep it under `tmp/work` only and explain
  why inspection alone was insufficient. Prefer read-only inspection.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule for future implementation:

This pass is no-edit, so the src/test/generator line budget should remain zero.
For the next implementation pass, assume the hard line-budget rule still
applies:

```text
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

must have more deleted than added lines unless the manager explicitly changes
the rule. While auditing, identify at least one deletion/shrinkage candidate
that can be bundled with the next implementation.

Decision rules:

- If the old WL q=5/n_s=5 construction is not reconstructable from current docs
  and code without running long probes, stop with the exact missing evidence.
- If current PQS can only do the growing-core one-shell family, do not stretch
  the wording. Say that the matching multi-shell PQS route does not exist yet.
- If a current PQS path looks close but uses shell/support-row oracle semantics
  instead of source-box-first PQS contraction, label that explicitly.
- If the first correct implementation step would require a scientific policy
  choice rather than engineering work, write that as the blocker.

Report back:

- Current git status.
- Whether any files changed. Expected: none.
- Concise inventory of the old WL He q=5/n_s=5 and n_s=7 validation status.
- Correct first PQS q=5 atom target, with expected retained counts/dimensions
  if known.
- Whether current PQS can express that target.
- Smallest next implementation seam.
- Deletion/shrinkage candidates for the next pass.
- Any stale assumptions corrected by this audit.

Write:

```text
.agent_handoffs/response.178.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.178.md
```

Sign:

```text
-- repo-doer@macmini
```

-- repo-manager@macmini
