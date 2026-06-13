Pass 105 manager review

Accepted.

The implementation adds the intended private complete core/shell RHF input
contract boundary and stays inside the pass scope. It validates route-owned
source plan, final basis, H1 payload, density inputs, Coulomb expansion,
explicit electron count, closed-shell occupation policy, and explicit fixture
role. It does not build a Fock matrix, run SCF, compute RHF energy, wire driver
or report fields, add public API, or promote a physics fixture.

The test is appropriately small and synthetic. It catches the available
closed-shell route-smoke contract plus missing electron count, odd/open-shell
input, missing fixture role, and missing H1 payload blockers.

Validation repeated by manager:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl`
  passed: 30/30.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed.
- `git diff --check` passed.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: RHF now has a separate private input-contract boundary, reducing
  pressure to grow H1/J helpers or report aliases into RHF/SCF behavior.
- quarantined: RHF remains private diagnostic/prototype; H1/J remains private
  diagnostic; fixture roles remain explicit.
- not deleted because: this was the first RHF boundary pass and there was no
  prior RHF implementation surface to remove.
- exact remaining caller/blocker: no production caller. Next seam should be a
  private one-step Fock/energy builder over an externally supplied final-basis
  density, still without SCF.

-- repo-manager@macmini
