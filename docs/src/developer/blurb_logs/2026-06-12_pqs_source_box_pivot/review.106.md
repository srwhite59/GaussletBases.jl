Pass 106 manager review

Accepted.

The implementation adds a private one-step complete core/shell RHF diagnostic
payload over an externally supplied spin-summed final-basis density. It stays
inside the pass boundary: no SCF loop, no convergence claim, no route-driver
wiring, no report aliases, no exports/artifacts, and no fixture promotion.

The contraction is tied to the existing pre-final positive-weight density
interaction and the existing one-orbital direct-minus-exchange convention. The
payload keeps large matrices out of summary/metadata and carries the required
nonclaims. The test is small and synthetic.

Validation repeated by manager:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`
  passed: 31/31.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed.
- `git diff --check` passed.

Follow-up before SCF:

- Add one tiny convention cross-check that constructs a closed-shell
  one-orbital density and asserts the one-step payload's two-body energy agrees
  with `pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)` on the same
  density interaction/orbital. This should be test-only unless it exposes a
  typo.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the one-step density contraction now has a separate private seam
  instead of pressure to put Fock/energy behavior into H1/J.
- quarantined: RHF remains private diagnostic/prototype; the payload explicitly
  reports no SCF and no convergence.
- not deleted because: this was the first one-step contraction seam and no
  older RHF implementation surface exists.
- exact remaining caller/blocker: no route-driver caller. Before a private SCF
  loop, the one-step convention should be cross-checked against the existing
  one-orbital self-Coulomb helper.

-- repo-manager@macmini
