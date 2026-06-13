Pass 127 manager review

Accepted and pushed.

Commit:

- `460d428e Default private Fock DIIS history to 8`

The diff is exactly the intended private-control default update:

- `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)` now
  resolves omitted Fock-DIIS `max_history` to 8 instead of 6.
- Explicit caller-supplied `max_history` remains unchanged.
- The focused SCF test expectation was updated from 6 to 8.

Validation reported by doer:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  passed: 91/91
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed
- `git diff --check` passed

Manager checks:

- `git show --check HEAD` passed.
- Commit touched only the intended source/test files.
- Commit was pushed to `origin/main`.

Decision:

- Next pass should be local ignored validation that the compact route-smoke
  probe converges when `max_history` is omitted and the new private default is
  used.
- Do not route-wire RHF, promote fixtures, loosen tolerances, or treat this as
  serious HF. Serious HF remains an `hfdmrg`/CR2 downstream boundary.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: callers can now use the confirmed private Fock-DIIS default
  instead of explicitly setting history 8 for the route-smoke diagnostic.
- quarantined: no new ignored artifacts from pass 127.
- not deleted because: this was a two-line private-control default update with
  no stale tracked path to retire.
- exact remaining caller/blocker: route-driver RHF remains intentionally
  unwired; the new default still needs one local compact route-smoke
  confirmation with omitted `max_history`.

-- repo-manager@macmini
