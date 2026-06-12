Purpose:

Advance Be S+P from final-basis one-electron materialization toward density-
density/RHF, but do it with phase attribution. Do not run a blind long
full-density job.

Why now:

The new private seam
`_white_lindsey_decomposed_atom_gto_final_basis_route(...)` removes the first
replacement-oriented blocker,
`:missing_driver_owned_decomposed_be_sp_fixture_wiring`. It reaches Be S+P
final-basis one-electron materialization with the same final dimension as the
old oracle:

- old oracle fixed/residual/final dimensions: `615 / 21 / 636`
- new route retained/supplement/final dimensions: `615 / 21 / 636`
- final overlap identity error: about `1.01e-10`
- no full-parent CPB, direct Cartesian, ordinary IDA, raw-GTO-final-density, or
  generalized-final-solve fallback

The one-electron route took about `148` seconds. The first uninstrumented
full-density attempt reached the density stage and was stopped because it did
not provide useful phase attribution. The next pass should make the expensive
parts visible before attempting RHF.

Exact task:

1. Add or use narrow timing attribution for the private Be S+P seam.
   Attribute at least:
   - parent-axis setup;
   - shellification/decomposed inventory;
   - overlap;
   - kinetic;
   - electron-nuclear by center;
   - mixed GTO blocks;
   - combined one-electron assembly;
   - final-basis projection;
   - residual moment matrices;
   - residual MWG representation;
   - gausslet density-density;
   - final-basis density-density assembly;
   - RHF solve, if reached.

2. Run the Be S+P current-route probe with `build_density_density = true`.
   Use the same fixture as the old oracle:
   - Be, `Z = 4`, all-electron;
   - q/ns `5 / 5`;
   - `legacy_atomic_gaussian_supplement("Be", "cc-pV5Z"; lmax = 1,
     basisfile = "/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets")`;
   - `interaction_treatment = :mwg`;
   - `residual_keep_policy = :near_null_only`.

3. If final-basis density-density materializes, attempt closed-shell RHF only
   in the final orthonormal basis.
   - Do not use raw generalized-overlap final solves.
   - Do not accept raw GTO density-density as final operator data.
   - Do not substitute old nested/QW matrices into the new route.

4. Compare against the old oracle only after the new route honestly reaches the
   matching stage:
   - old oracle RHF one-electron `-19.06620047058102`;
   - old oracle electron-electron `4.491686226006327`;
   - old oracle total `-14.574514244574694`.

Trust boundary:

Keep this private/internal. Do not add public APIs, exports, driver defaults,
PQS, ECP, high-l Be, Be2, Cr, H2, full-parent CPB fallback, direct Cartesian
fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or
generalized final-basis solves.

Do not add tests by default. Use `tmp/work` probes for this pass. If a tiny
production timing hook is needed, keep it private and report its carrying cost;
do not add broad timing framework code.

Decision rules:

- If density-density materializes and RHF converges, report the energy
  components, delta from old oracle, and phase timing.
- If density-density blocks, report the exact blocker and the last completed
  phase.
- If runtime is dominated by one obvious phase, stop after recording that
  phase and recommend the next optimization/replacement target.
- If phase attribution requires a broad framework, stop and propose the
  smaller instrumentation seam first.
- If the new route makes any H/H2+ fixture-local GTO setup dead or duplicative,
  delete or shrink it in the same pass. If not, state the next removal
  condition.

Artifacts:

- update or replace `tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`;
- write a concise summary artifact with phase timings;
- optional TSV for phase rows.

Validation:

- run the Be S+P probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if production source changes are made, run the most local relevant existing
  test unless the run is clearly too expensive, and explain any skipped test.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- whether timing attribution was added to production seam or only the probe;
- phase timing table;
- density-density status and exact blocker or success;
- RHF status, energies, and delta from old oracle if reached;
- whether any H/H2+ fixture-local wiring was deleted, simplified, or still
  needed;
- validation run;
- deletion/shrinkage report.
