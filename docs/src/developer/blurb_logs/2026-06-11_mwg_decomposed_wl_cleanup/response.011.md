Reran the Be atom S+P old nested/QW oracle probe with the authorized GaussletModules basis source.

Artifacts:
- `tmp/work/be_atom_sp_qw_oracle_probe.jl`
- `tmp/work/be_atom_sp_qw_oracle_summary.txt`

Exact basisfile used:
- `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets`

Old nested/QW oracle result:
- status: `:materialized_be_sp_old_nested_qw_oracle`
- atom/Z: Be / 4.0, all-electron
- q/ns: 5 / 5
- parent side count: 15
- physical outer radius: 10.056710949734484
- fixed dimension: 615
- raw S+P supplement orbital count: 21
- residual count: 21
- discarded supplement count: 0
- final dimension: 636
- final overlap identity error: 1.01090247284219e-11
- residual centers finite: true
- residual widths finite and positive: true

RHF result:
- converged: true
- iterations: 27
- one-electron energy: -19.06620047058102
- electron-electron energy: 4.491686226006327
- total energy: -14.574514244574694
- density trace: 2.0
- electron count: 4.0

Timings:
- supplement load: 0.005382541 s
- supplement 3D representation: 5.4583e-5 s
- fixed-block construction: 2.206159209 s
- QW/MWG operator construction: 3.673783375 s
- restricted closed-shell RHF: 2.032742834 s

Current decomposed/final-basis Be S+P audit:
- status: `:blocked_current_decomposed_final_basis_be_sp_audit`
- blocker: `:missing_driver_owned_decomposed_be_sp_fixture_wiring`
- lower-level surfaces present:
  - combined layout: true
  - combined one-electron matrices: true
  - final-basis projection: true
  - residual MWG representation: true
  - final-basis density-density matrix: true
- raw GTO density-density accepted as final operator: false
- generalized-overlap final solve used: false
- full-parent CPB fallback used: false
- ordinary Cartesian IDA fallback used: false

Interpretation:
The old q=5 Be S+P QW/MWG oracle is now clean. The current route has the relevant low-level combined-GTO and residual-MWG surfaces, but there is not yet a driver-owned decomposed Be S+P construction seam that turns the old q=5 parent/mapping plus authorized Be supplement into the current shellification-backed final-basis route. I did not patch around that with a fixture-local fallback.

Validation:
- Ran `julia --project=. tmp/work/be_atom_sp_qw_oracle_probe.jl`
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Ran `git diff --check`

Deletion/shrinkage report:
- No production code, tests, docs, metadata, or compatibility paths became obsolete.
- Nothing was deleted or simplified; this was a `tmp/work` oracle/readiness probe.
- No tracked test was added.
- Remaining stale/duplicate surface to retire next: fixture-local H/H2+ GTO final-basis wiring should be generalized into a driver-usable decomposed atom+supplement construction seam before Be S+P current-route RHF can be attempted.

-- repo-doer@macmini
