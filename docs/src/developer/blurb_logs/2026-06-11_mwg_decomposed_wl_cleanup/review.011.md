Review result:

Accepted as a clean Be atom S+P old nested/QW oracle plus current-route
readiness audit. The doer used the explicitly authorized GaussletModules
`BasisSets` file, built the old oracle first, and did not patch around the
new-route blocker.

Old oracle result:

- atom/Z: Be / `4.0`, all-electron
- q/ns: `5 / 5`
- parent side count: `15`
- physical outer radius: `10.056710949734484`
- fixed dimension: `615`
- raw S+P supplement orbitals: `21`
- residual count: `21`
- discarded supplement count: `0`
- final dimension: `636`
- final overlap identity error: `1.01090247284219e-11`
- RHF total: `-14.574514244574694`
- RHF one-electron/electron-electron:
  `-19.06620047058102` / `4.491686226006327`
- RHF converged in `27` iterations

Current-route blocker:

The lower-level combined-GTO and residual-MWG surfaces are present, but the
newer route does not yet have a driver-owned decomposed Be S+P construction
seam. The exact reported blocker is:

```julia
:missing_driver_owned_decomposed_be_sp_fixture_wiring
```

This is the first replacement-oriented blocker: removing it should make the
old nested/QW route less necessary as route authority and leave it as oracle
comparison only.

Validation reviewed:

- `julia --project=. tmp/work/be_atom_sp_qw_oracle_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/be_atom_sp_qw_oracle_summary.txt`

Deletion/shrinkage review:

No old surface became obsolete yet. No test was added. The next pass should
target the fixture-local H/H2+ GTO final-basis wiring as the surface to
generalize or shrink, rather than adding another oracle-only benchmark.

Commit/push:

Pending manager commit/push of this tracked response/review log and the next
published blurb.

Next target:

Add the driver-owned decomposed atom+supplement construction seam needed for
Be S+P. The pass should explicitly report whether any H/H2+ fixture-local GTO
wiring was deleted, simplified, or made obsolete.
