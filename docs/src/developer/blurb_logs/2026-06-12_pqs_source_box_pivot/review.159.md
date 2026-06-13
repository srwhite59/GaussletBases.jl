Pass 159 review - private Hamiltonian consumer contract

Accepted, with immediate carrying-cost caution.

The pass added a private inspect-only Hamiltonian consumer contract and, more
importantly, replaced the active readiness blocker:

```text
old: :missing_diatomic_hamiltonian_consumer_contract
new: :missing_hfdmrg_density_density_contract
```

That satisfies the hard gate in the blurb. The route no longer says there is no
consumer contract. It now says the private inspector contract exists, while
real downstream contracts for HFDMRG density-density `H,V`, sliced integrals,
HamV6, and CR2 handoff format remain missing.

Manager validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

passed 63/63 in 47.8s.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed.

`git diff --check HEAD~1..HEAD` passed.

The caution: this pass added 288 lines, mostly in the diatomic complete-core/
shell file. It bought a real state transition, but it also creates a new
handoff/consumer/readiness duplication surface. Do not add another downstream
payload before auditing whether this new contract can be made thinner or
whether repeated flags/fields can be consolidated.

Deletion/shrinkage accounting:

- deleted: old active readiness pressure on
  `:missing_diatomic_hamiltonian_consumer_contract`
- simplified: readiness now separates private inspectability from downstream
  HFDMRG/CR2/HamV6 readiness
- quarantined: consumer contract is private inspect-only and temporary; dense
  `Vee`, sliced integrals, H1/J, RHF/SCF, WL, public API, export/HamV6, CR2,
  HFDMRG, and artifacts remain unavailable
- not deleted because: the existing handoff remains the route-owned source of
  H1, pre-final density interaction, convention labels, and Be2 metadata
- exact remaining caller/blocker: readiness now blocks on
  `:missing_hfdmrg_density_density_contract`; source file carrying cost and
  handoff/consumer/readiness duplication need audit before downstream format
  implementation

-- repo-manager@macmini
