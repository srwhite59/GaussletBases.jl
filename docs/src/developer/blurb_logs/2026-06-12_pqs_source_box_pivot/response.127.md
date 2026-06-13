Pass 127 response - default private Fock-DIIS history to 8

Commit:
- `460d428e Default private Fock DIIS history to 8`

Diff summary:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
  - changed the private Fock-DIIS default history resolution from `6` to `8`
    only when `mixing_kind === :fock_diis` and `max_history === nothing`.
  - explicit caller-supplied `max_history` remains unchanged.
- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  - updated the focused expectation from `fock_diis_control.max_history == 6`
    to `== 8`.

Validation:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  - Passed: `91 / 91`
  - Note: run precompiled `GaussletBases`.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`
- `git diff --check`
  - Passed.

Git status after commit:
```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:
- deleted: none.
- simplified: private Fock-DIIS default now matches the confirmed compact
  route-smoke control setting.
- quarantined: no new ignored artifacts in this pass.
- not deleted because: no stale tracked path was introduced or retired by this
  two-line default change.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)` now
  defaults `:fock_diis` to `max_history = 8`; route-driver integration remains
  intentionally blocked and no serious-HF/public claim is made.

-- repo-doer@macmini
