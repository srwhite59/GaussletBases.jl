Pass 167 review - accepted

Reviewed implementation commit:

```text
9b6d4789 Add Be2 CR2 inspection artifact
```

Verdict: accepted.

The revised implementation fixes the two issues found in the pass-166 draft:

- the White-Lindsey route now has schema-shaped placeholder groups for
  `route`, `readiness`, `system`, `final_basis`, `one_body`, `two_body`, and
  `validation`;
- the focused Be2 test again guards the compact PQS route-spine semantics:
  source-plan object kind, support order, route retained order, final dimension,
  and final support row order.

The new private writer is intentionally not a public API or solver/export path.
It writes a read-only JLD2 inspection bundle and TSV fingerprint from plain
arrays and compact metadata:

- final one-body Hamiltonian;
- pre-final pair matrix;
- final-to-pre-final coefficients;
- pre-final and support weights;
- support raw pair numerator;
- nuclear/electron/spin metadata;
- density gauge and raw pair convention;
- CR2 readiness flags with solver/export still false.

Validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

passed 25/25 in 50.1s.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed.

```text
git diff --check
```

passed.

Line budget:

```text
src/test added:   250
src/test deleted: 276
net:              -26
```

The budget was paid by deleting the stale non-default private synthetic RHF SCF
validator:

```text
test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

`rg` found no runner references to that deleted file.

Remaining blocker:

CR2 solver/export handoff remains blocked by
`:missing_cr2_solver_handoff_format`; downstream solver/HamV6/HFDMRG readiness
remains blocked by `:missing_hfdmrg_density_density_contract` and related
downstream contract work.

-- repo-manager@macmini
