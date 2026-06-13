Pass 155 response - no-edit test-retirement audit for PQS/RHF/Be2 route work

Read-only audit only. I inspected the requested startup/policy docs, the current
fixture policy under `docs/src/developer/pqs_source_box_fixture_policy.md`
(`fixture_role_policy_2026-06-12.md` was not present under that exact name),
review 154, `test/runtests.jl`, `test/nested/runtests.jl`, the four RHF files,
and the Be2 Ham payload fingerprint file.

Default-runner inventory:

- `test/runtests.jl` includes `test/nested/runtests.jl` when the nested group is selected.
- `test/nested/runtests.jl:28-31` includes all four RHF seam files by default:
  `pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl`,
  `pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl`,
  `pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`, and
  `pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`.
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  is not in the default nested runner. It is a focused route-driver file.

Deletion candidates:

- `test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl`:
  delete from default pressure, and likely delete the file if pass 156 keeps one
  compact RHF validator elsewhere/in-place. It mostly protects private staged
  object names, missing-input blockers, and route-smoke nonclaim flags. The
  still-live bug is closed-shell input convention handling, but that should be
  covered once in the compact validator rather than as a standing seam file.
- `test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl`:
  delete. The useful bug is H1-Aufbau spin-summed density construction, but it
  is synthetic and already exercised through the SCF path. The compact
  replacement should assert density trace/idempotency once.
- `test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`:
  delete as a separate seam. The useful bug is the restricted direct-minus-
  exchange convention and one-step energy accounting; keep that convention in
  one compact RHF validator if RHF remains at all. Delete repeated blockers and
  nonpromotion flags.
- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`:
  keep only as the file to shrink in pass 156, or replace it with one smaller
  validator file. Delete most control/status/blocker inventory, DIIS default
  field assertions, repeated residual summary vocabulary, and repeated
  public/export/artifact false checks. The compact live check is: a tiny
  closed-shell synthetic Hamiltonian converges, preserves trace/idempotency,
  uses the intended ordinary final-basis commutator residual, and reports
  private/non-public once.
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  blocked/no-probe testset (`Be2 PQS Ham payload readiness fingerprint`, about
  lines 99-425): delete or reduce to at most one no-probe boundary smoke. It
  mainly preserves transitional blocked payload vocabulary and repeated
  nonclaim flags.
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  probe-enabled testset (`Be2 PQS probe-enabled Ham readiness fingerprint`,
  about lines 427-947): keep, but shrink hard. Keep one source-plan semantic
  smoke, one final-basis semantic smoke, and one H1/Ham handoff readiness check.
  Delete duplicated object-kind/status aliases, available/missing object
  inventories at every stage, repeated materialization false flags, and repeated
  public/export/artifact/RHF/WL/H1-J false checks.

Keep and why:

- Keep a compact RHF/private-Hamiltonian validator only if `src/pqs_multilayer_complete_core_shell_rhf.jl`
  remains as a live private validator. It catches a real convention bug that a
  downstream Be2 route shape fingerprint would not localize: closed-shell
  density trace/idempotency, Fock construction, and energy accounting.
- Keep a compact Be2 probe-enabled route fingerprint because it exercises the
  active source-plan -> final-basis -> H1 -> Ham-input/handoff route assembly
  in `cartesian_assembly(...)`. This is the current live Be2 PQS route contract.
- Keep the single readiness blocker
  `:missing_diatomic_hamiltonian_consumer_contract` temporarily. It is the exact
  remaining route blocker before adding a real Hamiltonian consumer contract.

Temporary keeps and exact blockers:

- RHF validator: temporary keep until manager either deletes
  `src/pqs_multilayer_complete_core_shell_rhf.jl` or wires a reviewed private
  Hamiltonian-validator consumer. Exact blocker:
  `RHF is private validator only; no route-owned Be2/PQS consumer currently
  calls the RHF helper family`.
- Be2 readiness blocker assertion: temporary keep until the Hamiltonian
  consumer contract exists. Exact blocker:
  `:missing_diatomic_hamiltonian_consumer_contract`.

Proposed pass 156 deletion plan:

1. Remove the four RHF includes from `test/nested/runtests.jl`.
2. Delete three RHF seam files:
   `pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl`,
   `pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl`, and
   `pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`.
3. Shrink `pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl` to one
   focused private validator, or delete it too and add one smaller replacement
   smoke only if manager wants a renamed file. Expected shrinkage if kept
   in-place: remove roughly 400-500 RHF test lines and remove four default
   nested includes.
4. Shrink `pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
   by deleting the blocked/no-probe field inventory and replacing the
   probe-enabled assertions with compact semantic checks. Expected shrinkage:
   several hundred lines, while preserving the route assembly check.

Source helper deletion pairing:

- RHF source deletion can be paired only if manager decides RHF is no longer a
  live private validator. Caller check found no source callers of the main RHF
  helper family outside `src/pqs_multilayer_complete_core_shell_rhf.jl`; current
  callers are tests and internal helper chaining.
- Do not delete Be2 diatomic Ham-input/handoff/readiness helpers in pass 156.
  `src/pqs_source_box_route_driver_helpers.jl` still calls them from
  `cartesian_assembly(...)`, and they protect the live Be2 route handoff seam.

Validation:

- Read-only inspection commands only.
- `git status --short --branch`:

```text
## main...origin/main
```

Deletion/shrinkage report:

- deleted: none in this audit
- simplified: none in this audit
- quarantined: none in this audit
- not deleted because: pass 155 is inventory/decision only
- exact remaining caller/blocker: RHF helpers have no source callers outside
  their own file and tests, but deletion is blocked by the decision whether RHF
  remains as a private Hamiltonian validator; Be2 route helpers are still called
  by `cartesian_assembly(...)`, and the route readiness blocker is
  `:missing_diatomic_hamiltonian_consumer_contract`

-- repo-doer@macmini
