Pass 155 review - test-retirement audit

Accepted. The audit stayed read-only and identified the highest-value deletion
target clearly.

The key finding is that the four private RHF seam tests are still included in
the default nested runner:

```text
test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl
test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl
test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl
test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

Together they are about 640 lines of default test pressure for a private RHF
validator path. That is the wrong carrying-cost balance now that the repo has
HFDMRG for serious HF/DMRG work and GaussletBases RHF is only a Hamiltonian
sanity validator.

The Be2 focused Ham payload fingerprint is also oversized, but it is not in the
default nested runner and it protects the current active Be2/PQS route spine.
It should be shrunk, but not in the same pass as the RHF deletion work.

Decision:

- pass 156 should retire RHF seam-test pressure first;
- remove RHF tests from the default nested runner;
- delete the input-contract, initial-density, and one-step seam test files;
- shrink the SCF test file to one compact private-Hamiltonian validator, or
  replace it with one smaller validator file if cleaner;
- do not delete RHF source helpers yet;
- do not touch the Be2 fingerprint file in pass 156.

Validation should be narrow: run the remaining compact RHF validator, check the
runner no longer includes deleted files, load the package, and run
`git diff --check`. Do not run the full nested default runner as a deletion-pass
gate unless there is a concrete reason and the runtime is called out first.

Deletion/shrinkage accounting:

- deleted: none by pass 155 itself
- simplified: deletion plan narrowed to RHF default test pressure first
- quarantined: none
- not deleted because: pass 155 was audit-only
- exact remaining caller/blocker: RHF source helpers have no source callers
  outside their own file/tests, but source deletion is blocked by the decision
  to keep RHF as a private Hamiltonian validator for now; Be2 route helper
  tests are still blocked by active `:missing_diatomic_hamiltonian_consumer_contract`

-- repo-manager@macmini
