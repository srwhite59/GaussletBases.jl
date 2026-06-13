Pass 111 manager review

Accepted.

The audit identifies a real contract issue that should be corrected before a
real compact PQS SCF probe: the SCF loop returns a post-diagonalization
`final_density`, while `final_one_step_payload` and `final_total_energy` are
computed from the density at the start of the converged iteration. That is close
under the density tolerance, but the payload labels are not tight enough for a
diagnostic result.

Required next corrective pass:

- On converged return, recompute the one-step payload using the returned
  `final_density`.
- Return that recomputed payload as `final_one_step_payload`.
- Keep iteration records as convergence history.
- Add summary metadata making the ordinary final-basis metric and final
  one-step consistency explicit.
- Add/update the focused SCF test to assert
  `scf.final_one_step_payload.final_density ≈ scf.final_density`.

Validation/status:

- No code changes in pass 111.
- No Julia commands.
- Reported git status: clean and even with origin/main.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: next pass should simplify the SCF payload semantics before real
  probing.
- quarantined: SCF remains private diagnostic/prototype.
- not deleted because: no implementation occurred.
- exact remaining caller/blocker: real compact PQS probe remains blocked by the
  final one-step/final-density consistency fix.

-- repo-manager@macmini
