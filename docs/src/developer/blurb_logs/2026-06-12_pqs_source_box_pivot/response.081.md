Pass 081 complete.

Driver-stage mapping:

- `cartesian_system / recipe`
  - currently still test-local in `_pqs_h1_complete_fixture`;
  - owns center record, charge, Coulomb expansion, and temporary mapped-parent fixture choice.
- `cartesian_parent`
  - `_pqs_h1_test_bundle(7)`;
  - `_CartesianNestedAxisBundles3D(bundle7, bundle7, bundle7)`;
  - carries PGDG/intermediate axis data used by source planning and nuclear factors.
- `cartesian_shells`
  - `CartesianShellification.shellify(...)`;
  - one-center `OneCenterShellification(core_side = 5, q = 5)`.
- `cartesian_units`
  - `CartesianTerminalLowering.lower_terminal_regions(..., PQSLowering(q = 5))`;
  - `pqs_multilayer_shell_region_plan(...)`;
  - source CPBs and owned-support provenance come from shellification/lowering.
- `cartesian_transforms`
  - region-backed `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`;
  - `pqs_multilayer_complete_core_shell_final_basis(plan; ...)`;
  - shell descriptors, shell projection/Lowdin cleanup, collapsed shell sector, and final overlap identity.
- `cartesian_pairs`
  - H1-only support kinetic helper;
  - H1-only separated by-center electron-nuclear support helper;
  - axis-layer/origin-factor comparison remains an oracle/convention check.
- `cartesian_assembly`
  - still test-local:
    - transfer support kinetic and nuclear matrices through `CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix`;
    - build one-electron Hamiltonian;
    - solve H1.
- `cartesian_report / materialization`
  - currently test assertions only;
  - no route-owned compact H1 payload/report exists yet.

Recommended next pass:

- Add a narrow internal H1 assembly payload helper, not another source-plan helper.
- Candidate shape:
  - consume the region-backed source plan;
  - consume or build the complete core/shell final basis;
  - consume Coulomb expansion and center records;
  - use the existing H1-only support kinetic and separated by-center electron-nuclear helpers;
  - transfer one-body matrices through `CartesianFinalBasisRealization`;
  - return final kinetic, separated final by-center nuclear records, the one-electron Hamiltonian, the H1 solve, and a compact nonclaim summary.
- This can be implemented mechanically next pass. I do not see a design blocker for an H1-only assembly payload as long as it explicitly makes no IDA, density-density, RHF, driver, export, artifact, GTO, or fixture-rule claim.

Docs edited:

- Updated `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` with the driver-stage mapping and the H1 assembly payload recommendation.

Validation run:

- `git diff --check`
  - passed.

Deletion/shrinkage report:

- No source/test code was deleted because this was an audit/design pass.
- If the recommended H1 payload lands, old test/probe glue that can shrink:
  - manual support kinetic construction in the H1 test;
  - manual separated nuclear support construction in the H1 test;
  - manual final one-body transfer calls in the H1 test;
  - manual Hamiltonian/H1 solve assembly in the H1 test;
  - repeated nonclaim/status checks that can move into a compact payload summary.
- Remaining helper-vocabulary pressure:
  - support-space kinetic and nuclear helpers remain live H1 seam machinery for now;
  - the axis-layer/origin-factor nuclear comparison remains a small convention/oracle check;
  - explicit-box source planning remains a compact bridge comparison only.

-- repo-doer@macmini
