Pass 156 response - retire RHF seam-test pressure from default suite

Implemented and committed:

```text
9bdd7942 Retire RHF seam test pressure
```

Files deleted:

- `test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`

Runner includes removed:

- Removed all four RHF includes from `test/nested/runtests.jl`, including the
  remaining SCF validator. The compact RHF validator is now focused/manual
  only, not default nested-suite pressure.

Remaining RHF validator:

- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
- line count: 96
- check count: 22 `@test` assertions
- retained behavior: two-electron closed-shell synthetic density trace,
  convergence, residual/idempotency checks, final one-step recomputation for
  returned final density, and private/non-public/export/artifact nonclaim.
- removed behavior: broad control-payload inventory, DIIS option inventory,
  blocked-path seam checks, repeated status vocabulary, and repeated nonclaim
  flag checks.

Validation:

```text
julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

passed:

```text
PQS complete core-shell RHF SCF payload | 22 pass | 1.4s
```

The command also precompiled `GaussletBases` first:

```text
57002.4 ms  ✓ GaussletBases
1 dependency successfully precompiled in 57 seconds. 50 already precompiled.
```

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed:

```text
load ok
```

```text
rg -n "pqs_multilayer_complete_core_shell_rhf_(input_contract|initial_density|one_step|scf)_runtests" test/nested/runtests.jl test/runtests.jl
```

returned no matches, as expected.

```text
git diff --check
```

passed with no output.

Git status after commit:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: three RHF seam test files; four RHF default nested-runner includes
- simplified: the remaining RHF SCF validator was reduced from 225 lines to 96
  lines and from broad seam/status inventory to one compact private-Hamiltonian
  convention check
- quarantined: RHF validation is now focused/manual only; it is no longer part
  of the default nested runner
- not deleted because: `src/pqs_multilayer_complete_core_shell_rhf.jl` remains
  untouched by instruction; the compact SCF validator remains because RHF is
  still a private Hamiltonian validator
- exact remaining caller/blocker: the RHF helper family still has no route
  source caller outside its own internal chaining and the focused SCF validator;
  source deletion remains blocked on a separate manager decision about retiring
  or retaining RHF as a private Hamiltonian validator

-- repo-doer@macmini
