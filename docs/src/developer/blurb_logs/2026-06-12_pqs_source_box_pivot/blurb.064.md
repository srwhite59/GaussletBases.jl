Purpose:
  Fill the remaining direct retained PQS safe-term gap for position and x2
  moment blocks. This is an operator-shape cleanup, not a new physics fixture.

Context:
  The direct retained PQS path already has compact blocks for:

  - overlap: `pqs_source_pair_retained_overlap_block(...)`
  - kinetic: `pqs_source_pair_retained_kinetic_block(...)`
  - by-center nuclear:
    `pqs_source_pair_retained_electron_nuclear_by_center_block(...)`
    and the centered wrapper

  Local inspection shows no corresponding direct retained wrappers for
  position or x2 moment terms. Those terms matter for later moment/weight/MWG
  work and should not require materializing full raw source matrices when the
  retained mode tuples already determine the block entries.

Task:
  Add direct retained-boundary PQS moment blocks for position and x2 terms.
  Use the existing direct retained product-term machinery in
  `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`.

  Add narrow wrappers along these lines, following the existing naming and
  metadata style:

  - `pqs_source_pair_retained_position_x_block`
  - `pqs_source_pair_retained_position_y_block`
  - `pqs_source_pair_retained_position_z_block`
  - `pqs_source_pair_retained_x2_x_block`
  - `pqs_source_pair_retained_x2_y_block`
  - `pqs_source_pair_retained_x2_z_block`

  If the local naming convention suggests a more compact axis-dispatch helper,
  use it internally, but keep the public/internal surface easy to read and
  consistent with overlap/kinetic/nuclear.

Validation:
  Compare each new direct retained moment block against the existing
  raw-source-block then retained-selector path on the existing small synthetic
  PQS retained fixture.

  Keep test growth minimal:

  - one compact check group is enough;
  - dimensions and numerical equality are the important assertions;
  - do not assert every metadata/nonclaim field;
  - do not add a new physics test;
  - do not expand slow integration tests.

  If the only practical location is the existing CPBM contract test, keep the
  added assertions tightly scoped and consider shrinking nearby redundant
  metadata assertions if practical. Prefer not to add another standalone test
  file unless it replaces/shrinks older coverage.

Do not:
  - change the side13 PQS fixture rule;
  - add H1/RHF/IDA/density-density;
  - add driver wiring, exports, artifacts, or GTO paths;
  - promote any same-parent ladder point as acceptance;
  - study or codify a new `Z,d,s,ns` rule in this pass.

Docs:
  Update `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  only if needed to mark the direct retained position/x2 gap closed. Keep it
  concise.

Validation commands:
  - focused CPBM/retained PQS test that covers the new moment wrappers;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini
