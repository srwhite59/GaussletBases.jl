Review result:

Accepted. The pass shrinks the He acceptance test to the live scientific
contract and removes broad report-field/timing/helper-vocabulary assertion
pressure.

Corrections made:

Removed the remaining full-report `println`, which was leftover audit output
and worked against the cleanup goal.

Validation:

- `julia --project=. test/nested/cartesian_wl_gausslet_he_atom_acceptance_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Commit/push:

Committed as `88ae9dbd Shrink He WL acceptance assertions`; push handled by
repo-manager after review.

Next target:

Run a side13 He + GTO final-basis RHF diagnostic using the now-materialized
final-basis density-density matrix.
