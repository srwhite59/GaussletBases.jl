Implemented the centered PQS source electron-nuclear composition helper.

Changes:

- Added/exported
  `pqs_source_pair_centered_electron_nuclear_by_center_block(record;
  axis_layers, coulomb_expansion, center_record)`.
- The helper calls
  `pqs_source_pair_centered_gaussian_factor_terms_1d(...)` and passes the
  generated source-mode Gaussian factors into the existing
  `pqs_source_pair_electron_nuclear_by_center_block(...)`.
- Added/exported the retained composition wrapper
  `pqs_source_pair_retained_centered_electron_nuclear_by_center_block(...)`,
  which composes the centered source block with the existing retained
  source-mode contraction.
- Preserved the current convention: by-center block, nuclear charge recorded
  but not applied, centers not summed, uncharged by-center source block. The
  wrappers do not call CCPM wrappers and do not build shell realization,
  Lowdin cleanup, IDA data, Hamiltonians, global routes, drivers, exports, or
  artifacts.
- Extended the existing CPBM source-pair contract test with compact checks that
  the centered source helper matches the supplied-factor path using the
  pass-011 centered factors, and that the retained wrapper matches
  `pqs_source_pair_retained_one_body_block(source_result)`.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old production path became obsolete in this pass. This is a convenience
  composition over the newly added centered factor generator and existing
  supplied-factor source block.
- No CCPM wrapper, fallback route, compatibility shim, global H1 route, or new
  metadata layer was added.
- The new test assertions are live-contract coverage inside the existing CPBM
  contract file; they do not create a new test suite and only check equivalence
  to the existing helper plus the live by-center nonclaim flags.
- Remaining stale/duplicate surfaces to retire next: old CCPM source-box
  electron-nuclear helper vocabulary should stay quarantined until a downstream
  PQS route consumes this centered source block path and no caller needs the old
  wrappers as an oracle.

-- repo-doer@macmini
