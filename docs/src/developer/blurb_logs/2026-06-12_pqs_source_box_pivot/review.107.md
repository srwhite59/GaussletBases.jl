Pass 107 manager review

Accepted.

The change is test-only and adds the requested convention guard. The one-step
two-body energy for a closed-shell one-orbital final density is now checked
against the existing
`pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)` diagnostic on the
same synthetic density interaction. The observed comparison is 2.0 versus 2.0,
and the test preserves the density-convention and no-SCF/nonconvergence
assertions.

Validation repeated by manager:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`
  passed: 38/38.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed.
- `git diff --check` passed.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: none directly; this guards the private RHF one-step convention.
- quarantined: RHF remains private diagnostic-only.
- not deleted because: the one-step test is now the narrow active convention
  guard.
- exact remaining caller/blocker: no route-driver caller. Next useful seam is a
  private closed-shell initial-density payload before a full SCF loop.

-- repo-manager@macmini
