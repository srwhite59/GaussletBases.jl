Pass 047 response: PQS complete core-shell one-body oracle comparison

Summary:
- Built and ran `tmp/work/pqs_complete_core_shell_one_body_oracle_probe.jl`.
- Updated the ignored local H1 probe `tmp/work/pqs_complete_core_shell_h1_probe.jl` to use `bundle7.pgdg_intermediate.gaussian_factor_terms` instead of raw `gaussian_factor_matrices(base_layer)` for the centered nuclear factor source.
- No tracked production source or tests were changed.

Files edited/artifacts:
- Added ignored probe artifact:
  - `tmp/work/pqs_complete_core_shell_one_body_oracle_probe.jl`
  - `tmp/work/pqs_complete_core_shell_one_body_oracle_probe_summary.txt`
- Updated ignored probe artifact:
  - `tmp/work/pqs_complete_core_shell_h1_probe.jl`
  - `tmp/work/pqs_complete_core_shell_h1_probe_summary.txt`
- Changed function/object in ignored probe only:
  - `_support_nuclear_matrix(...)` now consumes `pgdg_intermediate.gaussian_factor_terms`.

Oracle reached:
- Same-mapping old nested fixed-block oracle reached with:
  - parent count 7
  - same `MappedUniformBasisSpec(:G10; count = 7, AsinhMapping(...))`
  - same Coulomb expansion
  - `one_center_atomic_full_parent_fixed_block(bundle7; expansion, nside = 5)`
- The fixed block is a trusted same-geometry convention oracle, but not exact final-basis authority for the new complete core/shell basis:
  - fixed dimension: 223
  - fixed overlap identity error: `3.7636560534792807e-14`
  - current/fixed final-subspace identity error after cross-overlap alignment: `5.8875461168328425e-04`
  - native fixed-block H1: `-4.8047920531279725e-01`
  - fixed-block H1 in current gauge: `-4.8047920284588080e-01`

Other oracle blockers:
- CCPM source-box nuclear helper family was not used as active authority. It is blocked for this exact comparison by the missing adapter from combined `core_then_shell` support rows to the product/PQS source-box unit objects those helpers require.
- CPBM PQS centered nuclear helper was not used as active authority. It requires a `PairBlockMaterializationRecord` for raw product source pairs, not the combined complete core/shell support row order.

Convention comparison:

| variant | support generalized H1 | final H1 | kinetic expectation | nuclear expectation | max delta vs current raw | sign/charge convention |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| current raw `gaussian_factor_matrices(base_layer)` | `-2.0639248059188007` | `-2.0638461028784740` | `1.1012644410060735` | `-3.1651105438845457` | `0.0` | negative uncharged electron-nuclear; charge applied at Hamiltonian assembly |
| `pgdg_intermediate.gaussian_factor_terms` | `-0.4805121123666641` | `-0.4804793480038722` | `0.5280430374324055` | `-1.0085223854362777` | support nuclear delta `0.7564844936233156`; final nuclear delta `0.28444711168458364` | negative uncharged electron-nuclear; charge applied at Hamiltonian assembly |
| same-mapping fixed-block oracle | final-only oracle | native `-0.48047920531279725`; current-gauge `-0.48047920284588080` | `0.5280493776478693` | `-1.0085285804937483` | final nuclear delta vs raw `0.28445322389723227`; vs PGDG `0.00016236200112270055` | fixed block stores positive Gaussian sum; negative sign applied for electron-nuclear |

Decision:
- This is a small local probe/operator-construction bug, not a final-basis transfer issue.
- The correct same-bundle nuclear factor source is `pgdg_intermediate.gaussian_factor_terms`.
- Raw `gaussian_factor_matrices(base_layer)` is the wrong convention for this PGDG/numerical-reference bundle and explains the overbound `-2.0638` value.
- No production source patch was needed because the route-owned final-basis helper consumes caller-supplied support operators; the bad convention lived in the ignored probe-side support operator construction.

Corrected H1 probe:
- `tmp/work/pqs_complete_core_shell_h1_probe.jl` now reports:
  - status: `materialized_pqs_complete_core_shell_final_h1_solve`
  - blocker: `nothing`
  - final dimension: 223
  - H1 energy: `-4.8047934800387126e-01`
  - error vs `-0.5`: `1.9520651996128735e-02`
  - solve kind: `ordinary_symmetric`
  - nuclear factor source: `pgdg_intermediate_gaussian_factor_terms`
  - raw base-layer Gaussian factor matrices used: `false`

Validation:
- `julia --project=. tmp/work/pqs_complete_core_shell_one_body_oracle_probe.jl` passed.
- `julia --project=. tmp/work/pqs_complete_core_shell_h1_probe.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- `git status --short --branch` showed no tracked source/test changes before writing this response.

Deletion/shrinkage report:
- No old production code, tests, metadata, or compatibility path became obsolete in this pass because the only incorrect path was an ignored probe-side support nuclear construction.
- Nothing tracked was deleted or simplified.
- The stale raw-base-layer probe convention is now isolated as a comparison variant in `pqs_complete_core_shell_one_body_oracle_probe`, not the corrected H1 probe path.
- No new permanent test was added; the new artifact is an ignored developer probe, not new live-contract coverage.
- Remaining stale/duplicate surface to retire next: any future complete core/shell one-body probe or handoff text that treats raw `gaussian_factor_matrices(base_layer)` as the PGDG nuclear-factor authority should be updated to `pgdg_intermediate.gaussian_factor_terms`.

-- repo-doer@macmini
