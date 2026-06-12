Purpose:
  Reconcile the newly added direct retained PQS position/x2 wrappers with the
  existing generic retained one-body selector path, if this is a small local
  cleanup. Do not turn it into a broader API redesign.

Context:
  Pass 064 added direct retained wrappers for:

  - `pqs_source_pair_retained_position_x_block`
  - `pqs_source_pair_retained_position_y_block`
  - `pqs_source_pair_retained_position_z_block`
  - `pqs_source_pair_retained_x2_x_block`
  - `pqs_source_pair_retained_x2_y_block`
  - `pqs_source_pair_retained_x2_z_block`

  The manager removed the new export entries before committing because the
  blurb explicitly said not to add exports. Keep these wrappers qualified/internal
  unless a later reviewed API pass promotes them.

  Local inspection suggests the generic retained source one-body selector still
  only supports overlap/kinetic, for example:

  - `pqs_source_pair_retained_one_body_block(record, term; ...)`
  - `pqs_source_pair_retained_one_body_blocks(plan, term; ...)`
  - `pqs_retained_source_one_body_matrix(batch_result)`
  - `_supported_pqs_source_retained_safe_term_descriptor`
  - `_pqs_retained_source_matrix_supported_term`

Task:
  Audit this generic retained one-body selector seam.

  If the change is local and straightforward, add position/x2 support so callers
  can request retained moment blocks through the same generic retained selector
  path used for overlap/kinetic. Use the pass-064 direct wrappers/helpers as the
  implementation route; do not route through full raw source matrices except as
  oracle/reference.

  If the change requires a broader API decision, stop and report exactly why
  rather than forcing it.

Required behavior if implemented:

  - retained `:position_x`, `:position_y`, `:position_z`, `:x2_x`, `:x2_y`,
    `:x2_z` requests materialize direct retained blocks;
  - existing overlap/kinetic behavior is unchanged;
  - `pqs_retained_source_one_body_matrix(...)` accepts square retained
    self-pair moment matrices when requested;
  - no new exports;
  - no H1/RHF/IDA/density-density/driver/artifact changes;
  - no fixture-rule study.

Validation:
  Add or update only compact checks near the existing retained selector tests:

  - direct generic retained position/x2 block equals the named direct wrapper;
  - retained matrix wrapper accepts a square self-pair position/x2 result if
    this is in scope;
  - avoid metadata vocabulary assertions beyond dimensions, term, and equality.

  Run:

  - the focused CPBM retained/PQS contract test or the smallest test file that
    covers the edited seam;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Docs:
  Update `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  only if needed to say that the generic retained selector now covers
  position/x2 too. Keep it concise.

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini
