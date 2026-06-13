Pass 154 review - split diatomic PQS route helpers

Accepted with one manager-side mechanical repair.

Doer implemented the requested behavior-preserving file split. The private
diatomic complete-core/shell route payload subsystem moved out of
`src/pqs_source_box_route_driver_helpers.jl` into:

```text
src/pqs_source_box_diatomic_complete_core_shell.jl
```

The split moved the diatomic readiness, support-window, raw-box,
source-realization, source-plan, final-basis, H1, Ham-input, and Hamiltonian
handoff payload definitions plus their direct helpers. The broad route-driver
glue and `cartesian_assembly(...)` remain in
`pqs_source_box_route_driver_helpers.jl`.

The commit is now:

```text
b3b09d9a Split diatomic PQS route helpers
```

I amended the original doer commit only to remove a trailing blank line at EOF
in the new file so the committed range passes `git diff --check`. No behavior,
names, statuses, blockers, summaries, or assertions were changed by that amend.

Validation observed:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed after precompile.

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

passed:

```text
Be2 PQS Ham payload readiness fingerprint: 225/225
Be2 PQS probe-enabled Ham readiness fingerprint: 350/350
```

`git diff --check HEAD~1..HEAD` now passes.

Deletion/shrinkage accounting:

- deleted: no behavior deleted; 3,465 lines were removed from the oversized
  route-driver helper file by moving the private diatomic subsystem into its
  own file
- simplified: `pqs_source_box_route_driver_helpers.jl` no longer carries the
  full diatomic complete-core/shell payload subsystem
- quarantined: diatomic payloads remain private/internal; dense `Vee`,
  export/HamV6, CR2/HFDMRG readiness, RHF/SCF, WL, and H1-J remain nonclaimed
- not deleted because: this was a behavior-preserving split and compatibility
  report surfaces remain live
- exact remaining caller/blocker: `cartesian_assembly(...)` still calls the
  moved private helpers; route readiness remains blocked at
  `:missing_diatomic_hamiltonian_consumer_contract`

Next: audit test carrying cost before adding more route features.

-- repo-manager@macmini
