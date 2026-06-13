Pass 161 review - thin diatomic Hamiltonian consumer contract

Accepted.

This pass met the new corrective line-budget rule:

```text
git diff --numstat HEAD~1..HEAD -- src test
45  142  src/pqs_source_box_diatomic_complete_core_shell.jl
9   12   test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Totals:

```text
src/test added:   54
src/test deleted: 154
net:             -100
```

The consumer contract is now much thinner. It keeps the live contract fields:
status/blocker, route family, source handoff reference/status, readiness,
available/missing objects, summary, and metadata. The copied scalar fields that
belong to the handoff were removed from the consumer struct and constructor.

The new downstream nonclaim helper centralizes the downstream-missing objects
and false readiness flags. The route behavior stayed intact:

```text
readiness blocker = :missing_hfdmrg_density_density_contract
consumer contract available
private_inspector_ready = true
downstream readiness flags remain false
```

Manager validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

passed 55/55 in 47.7s.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed.

`git diff --check HEAD~1..HEAD` passed.

Deletion/shrinkage accounting:

- deleted: copied consumer scalar fields, scalar-copy builder code, repeated
  downstream false-flag lists, and copied-field test assertions
- simplified: consumer now points at handoff authority and carries compact
  readiness/nonclaim summary
- quarantined: consumer remains private inspect-only; HFDMRG density-density,
  sliced integrals, HamV6/export, CR2 format, dense `Vee`, H1/J, and RHF/SCF
  remain unavailable
- not deleted because: handoff still owns the route objects and compact scalar
  inspection summary; readiness still owns the active downstream blocker
- exact remaining caller/blocker: readiness remains blocked on
  `:missing_hfdmrg_density_density_contract`

Next: decide which downstream contract is worth implementing first. Do not add
another source/test payload without the line-negative rule or an explicit
exception.

-- repo-manager@macmini
