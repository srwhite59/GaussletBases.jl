Pass 198 - audit H2 source-plan/final-basis inventory before H1-J.

Purpose:

Pass 197 materialized H2 H1, but the lowest H1 eigenvalue is positive:

```text
physics/h1_lowest = 0.14582426982296057
```

That is a physics warning. Do not continue to H1-J, density interaction, or RHF
until we understand whether the current 221-dimensional H2 basis is only a
source-box/boundary diagnostic rather than the intended atom-core-plus-shell
physics basis.

Task type:

Audit plus deletion. No H2 implementation.

Audit questions:

1. What exactly is in the current H2 source plan and final basis?

   Report:

   ```text
   parent axis counts
   source boxes
   retained unit keys/order
   retained unit kinds
   support/source dimensions
   retained counts per unit
   final dimension
   whether full 5^3 atom-core interiors are present
   whether only boundary COM/product modes are retained for atom boxes
   midpoint/product slab shape/count
   shell/multishell count, if any
   ```

2. Why is the dimension 221?

   Confirm the arithmetic. For example, current hints suggest:

   ```text
   pqs_left  = 98 boundary modes?
   pqs_right = 98 boundary modes?
   product/midpoint slab = 25?
   total = 221
   ```

   Verify from live route/source-plan objects, not just from memory.

3. Does this match the intended H2 physics target?

   Compare against:

   - old WL/QW default complete-rectangular H2 fixed-block row size `1215`;
   - old WL/QW default final dimension `481` with H/cc-pVTZ S/P supplement;
   - the user expectation that a physical atom-centered minimal target should
     include full `5^3` atom cores and then shells, not only boundary modes.

4. Interpret the positive H1.

   Do not fix it in this pass. Explain whether the positive lowest H1 is
   consistent with the current retained basis missing atom-core interior modes
   or other expected bound-state content.

5. Recommend the next target.

   Choose one:

   ```text
   A. Keep current 221-dimensional H2 as a route-smoke/source-box diagnostic
      only, and add explicit artifact labels later.

   B. Replace/extend the H2 driver target to a physical gausslet-only basis
      with full atom-core interiors plus shell layers before continuing to
      H1-J/RHF.

   C. Something else, with exact reason.
   ```

   Be specific about the next implementation seam if B is recommended.

6. Delete one safe stale scaffold.

   Requirements:

   - no include in default/integration/top-level runners;
   - no source caller;
   - not an explicit He/H2 driver endpoint;
   - not a WL H/H2 reference or scientific endpoint;
   - not needed by this audit.

   If no safe deletion exists, write `.agent_handoffs/ATTENTION.md` and stop.

Trust boundary:

- No H2 H1-J/density interaction.
- No private RHF for H2.
- No H2 source-plan/final-basis implementation changes.
- No supplemented WL/QW comparison.
- No supplement support.
- No Be2/Cr2 artifact work.
- No HFDMRG, DMRG, ECP, exports, public solver behavior, or Qiu-White
  correction work.
- Do not add this test to default runners.
- Preserve the visible staged driver style.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass must be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Expected implementation diff is deletion only. The audit response/log does not
count toward this source/test/bin budget.

Validation:

Use read-only inspection and focused commands only. Do not rerun broad suites.
At minimum:

```text
rg -n "<deleted_file_basename>|<deleted_testset_name>" test src bin
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Report back:

- inventory findings and arithmetic for 221;
- explanation of likely positive-H1 cause;
- recommended next target A/B/C;
- stale scaffold deleted and why safe;
- validation results;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini
