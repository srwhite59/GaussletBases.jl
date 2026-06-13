Pass 156 review - retire RHF seam test pressure

Accepted.

The pass removed RHF from the default nested test pressure and deleted the
low-value seam tests:

```text
deleted test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl
deleted test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl
deleted test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl
```

`test/nested/runtests.jl` no longer includes any RHF seam file by default. The
remaining focused/manual validator:

```text
test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

was shrunk from 225 lines to 96 lines and now checks one compact private RHF
Hamiltonian-validator convention: a two-electron closed-shell synthetic
Hamiltonian converges, preserves trace/idempotency, recomputes the final
one-step payload for the returned final density, and stays private/non-public.

Manager validation:

```text
julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

passed 22/22 in 1.4s.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed.

`git diff --check HEAD~1..HEAD` passed.

Deletion/shrinkage accounting:

- deleted: three RHF seam test files; four RHF default nested-runner includes
- simplified: remaining RHF validator shrank to one compact numerical/private
  convention check
- quarantined: RHF validation is now focused/manual only, not default nested
  suite pressure
- not deleted because: `src/pqs_multilayer_complete_core_shell_rhf.jl` remains
  as a private Hamiltonian validator for now
- exact remaining caller/blocker: RHF helper family still has no route source
  caller outside internal chaining and the focused validator; source deletion
  remains a separate manager decision

Next: shrink the oversized Be2 Ham payload fingerprint file without weakening
the active Be2/PQS route-spine check.

-- repo-manager@macmini
