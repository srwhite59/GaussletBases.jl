Purpose:

Remove the first replacement-oriented Be S+P blocker:

```julia
:missing_driver_owned_decomposed_be_sp_fixture_wiring
```

Build a driver-owned decomposed atom+supplement construction seam so Be S+P can
reach the current combined-GTO/final-basis route without relying on H/H2+
fixture-local wiring.

Why now:

The old nested/QW Be S+P oracle is clean:

- q/ns `5 / 5`
- fixed dimension `615`
- raw S+P supplement orbitals `21`
- residual count `21`
- final dimension `636`
- RHF total `-14.574514244574694`

The current route already reports the low-level combined layout, one-electron
matrices, final-basis projection, residual MWG representation, and final-basis
density-density surfaces as present. The missing piece is not another old
oracle. It is a driver-owned seam that turns the Be q=5 parent/mapping,
shellification-backed decomposed inventory, center records, and authorized
S+P supplement into those existing route-global/final-basis pieces.

Added steering rule:

After the old oracle is clean, name the first new-route blocker whose removal
would make an old oracle-only surface less necessary. For this pass, that
blocker is the missing driver-owned decomposed Be S+P construction seam. The
old nested/QW route should remain an oracle comparator, not the architecture
parent.

Known surfaces:

- Current fixture-local wiring lives in
  `test/nested/cartesian_wl_gausslet_h_atom_acceptance_runtests.jl`, especially
  `_wl_decomposed_h_gto_supplement_acceptance_report`,
  `_wl_decomposed_h2plus_gto_supplement_acceptance_report`, and
  `_wl_final_basis_gto_one_electron_solve`.
- Current reusable route pieces include:
  - `white_lindsey_shellification_decomposed_unit_pair_inventory(...)`
  - `route_global_combined_gto_basis_layout(...)`
  - `route_global_mixed_gto_blocks_from_decomposed_units(...)`
  - `route_global_combined_gto_one_electron_matrices(...)`
  - `route_global_combined_gto_final_basis_projection(...)`
  - `route_global_residual_gto_mwg_representation(...)`
  - `route_global_combined_gto_final_basis_density_density_matrix(...)`
- The old oracle fixture used:
  - `legacy_atomic_gaussian_supplement("Be", "cc-pV5Z"; lmax = 1,
    basisfile = "/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets")`
  - `interaction_treatment = :mwg`
  - `residual_keep_policy = :near_null_only`

Exact task:

1. Add the narrowest driver-owned helper/seam needed to construct a decomposed
   one-center atom plus GTO supplement final-basis route for the Be S+P q=5
   fixture.
   - It may be private/internal.
   - It should consume the same Be q=5 mapping/parent setup and authorized
     supplement used by the oracle probe.
   - It should use shellification-backed decomposed inventory, not the
     low-order seed-only path.
   - It should call the existing route-global combined-GTO and residual-MWG
     pieces rather than duplicating them.

2. Retry the Be S+P current-route probe through that seam.
   - First aim for final-basis one-electron materialization.
   - If final-basis density-density is available, attempt the closed-shell RHF
     only through the final orthonormal basis.
   - Do not use raw generalized-overlap final solves.
   - Do not accept raw GTO density-density as final operator data.

3. Compare against the old oracle only at appropriate levels:
   - dimensions: old fixed `615`, old residual `21`, old final `636`
   - overlap/final-basis identity
   - one-body or RHF energy if the new route honestly reaches those stages
   - exact blocker if it cannot

4. Report whether this new seam replaces, shrinks, or makes obsolete any
   H/H2+ fixture-local GTO wiring. If not, explain why it is not removable yet
   and name the next removal condition.

Trust boundary:

Do not add public APIs, exports, driver defaults, PQS paths, ECP paths,
full-parent CPB fallbacks, direct Cartesian fallbacks, ordinary Cartesian IDA
fallbacks, raw GTO density-density final operators, generalized final-basis
solves, high-l Be, Be2, Cr, or H2 work.

Do not add tests by default. Use `tmp/work` probes for exploratory validation.
If production code is added, add at most one compact module-contract test only
if it protects a live reusable seam and either replaces/shrinks fixture-local
coverage or catches a non-obvious bug. Do not add helper-vocabulary or broad
metadata tests.

Decision rules:

- If the driver-owned seam can be added cleanly and used by the Be S+P probe,
  do it.
- If the seam would become a broad framework, stop and report the smaller seam
  proposal first.
- If Be S+P reaches final-basis one-electron but not RHF, report the exact
  electron-electron/final-density blocker.
- If Be S+P reaches RHF, report energy components and compare to the old oracle
  without claiming it is an acceptance test.
- If an old fixture-local helper becomes dead or duplicative, delete or shrink
  it in the same pass.

Artifacts:

- update or replace `tmp/work/be_atom_sp_qw_oracle_probe.jl` only if useful;
- create a clearly named `tmp/work` probe for the new current-route seam, for
  example `tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`;
- write a concise summary artifact for the new-route attempt.

Validation:

- run the relevant `tmp/work` probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if production code changes are made, run the most local existing test or one
  compact new seam test if justified.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- new helper/seam name and ownership;
- whether it used shellification-backed decomposed inventory;
- Be S+P current-route status and exact blocker or energies;
- comparison to old oracle dimensions/energies where meaningful;
- whether old H/H2+ fixture-local wiring was deleted, simplified, or still
  needed;
- validation run;
- deletion/shrinkage report.
