Pass 112 manager review

Accepted.

The corrective pass fixes the SCF payload consistency issue found in pass 111.
On converged return, the helper now recomputes the one-step payload using the
returned `final_density`, returns that payload as `final_one_step_payload`, and
uses its energy for `summary.final_total_energy`. Iteration records remain
convergence history. The summary now explicitly records the ordinary
orthonormal final-basis metric and the final one-step density match.

Validation repeated by manager:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  passed: 35/35.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed.
- `git diff --check` passed.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: final SCF payload semantics now point energy and final one-step
  diagnostics at the returned final density.
- quarantined: SCF remains private diagnostic/prototype behavior; no route,
  report, public API, export, or artifact behavior was added.
- not deleted because: this was a corrective pass on an active private helper.
- exact remaining caller/blocker: no route-driver caller. The next step can be
  a local ignored `tmp/work` compact real PQS SCF probe.

-- repo-manager@macmini
