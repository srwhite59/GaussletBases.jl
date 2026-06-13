Pass 108 manager review

Accepted.

The implementation adds a private H1-Aufbau initial-density payload and stays
inside scope. It validates the RHF input contract and H1 payload, symmetrizes
the H1 matrix before diagonalization, occupies the lowest closed-shell
orbitals, and returns a spin-summed final-density matrix with compact summary
metadata. It does not call the one-step Fock helper, run SCF, compute RHF
energy, wire the route driver, or add report fields.

Validation repeated by manager:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl`
  passed: 27/27.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed.
- `git diff --check` passed.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the future SCF loop no longer needs to own H1-Aufbau
  initialization logic.
- quarantined: initial density remains private diagnostic-only, with no SCF or
  convergence claim.
- not deleted because: no existing initial-density path was being replaced.
- exact remaining caller/blocker: no route-driver caller. Next seam is a
  private closed-shell RHF iteration payload that consumes the input contract,
  H1 payload, density interaction, and initial-density payload.

-- repo-manager@macmini
