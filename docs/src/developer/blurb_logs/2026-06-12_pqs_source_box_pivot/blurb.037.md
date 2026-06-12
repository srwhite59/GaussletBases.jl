Purpose:
  Extend the direct retained-boundary PQS algorithm shape to by-center nuclear
  blocks, and tighten the retained-rule assertion introduced in pass 036.

Context:
  Pass 036 made retained overlap/kinetic avoid full raw source-space block
  materialization. The remaining one-body PQS retained helper that still uses
  the raw-source-block -> selector path is by-center electron-nuclear. The
  direct retained nuclear block can use the same retained source-mode tuples and
  centered Gaussian factor terms without filling the full raw source matrix.

Exact task:
  1. Tighten `_assert_pqs_source_record_retained_rule(...)` so it explicitly
     requires:

     `rule.retained_rule_kind === :boundary_comx_product_mode_selection`

     in addition to the current type, dimension, ordering, count, and selector
     checks.

  2. Add a direct retained-boundary by-center electron-nuclear helper that
     iterates retained mode tuples from `CRPS.retained_mode_indices(rule)` and
     fills only the retained block from supplied Gaussian factor terms and
     Coulomb coefficients.

  3. Route the active retained by-center nuclear wrappers through the direct
     retained construction when the record carries boundary source-mode
     retained rules:

     - `pqs_source_pair_retained_electron_nuclear_by_center_block(...)`
     - `pqs_source_pair_retained_centered_electron_nuclear_by_center_block(...)`

  4. Keep the existing raw-source-block -> retained selector path as oracle and
     fallback for already materialized source results. Do not delete
     `pqs_source_pair_retained_one_body_block(source_result, ...)`.

  5. Preserve the uncharged by-center convention. The direct retained nuclear
     result must record center key/index/location and nuclear charge, but must
     not apply charge or sum centers.

  6. Update comments/metadata so the active retained nuclear path does not
     claim full raw source-space materialization. Preserve clear nonclaims: no
     shell realization, Lowdin, final basis, Hamiltonian, IDA, RHF, driver,
     export, or artifact.

Test policy:
  Update the existing compact CPBM contract coverage only. Compare:

  - direct retained centered by-center nuclear block vs existing centered raw
    source block -> selector oracle;
  - retained shape/count for the existing boundary-rule fixture;
  - center metadata and uncharged by-center convention.

  Avoid broad metadata assertions.

Do not:
  - add shell realization, final basis, H1, IDA, density-density, RHF, driver
    wiring, exports/artifacts, or GTO changes;
  - add a new broad test file;
  - remove the raw-source selector oracle path.

Validation:
  - `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

Required response:
  - files edited;
  - names of direct retained nuclear helper(s) added or changed;
  - confirmation that active retained by-center nuclear avoids full raw source
    block materialization;
  - oracle comparison result and max difference;
  - validation run and result;
  - deletion/shrinkage report:
      - whether raw-source selector nuclear is now oracle/fallback only;
      - comments/metadata simplified;
      - test assertions added/removed;
      - remaining next target.

Continue the baton loop after writing `response.037.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini
