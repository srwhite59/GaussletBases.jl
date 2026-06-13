Pass 157 review - shrink Be2 PQS Ham fingerprint

Accepted.

The focused Be2/PQS Ham payload fingerprint was reduced from 947 lines to 195
lines, and from 575 assertions to 49 assertions. The no-probe blocked testset
was deleted entirely. The remaining probe-enabled test now checks the live
route spine:

```text
source plan -> final basis -> H1 -> Ham input -> Hamiltonian handoff readiness
```

The test keeps the important semantic and numerical facts:

- source plan object kind is `:pqs_diatomic_complete_core_shell_source_plan`
  and not the old one-center source-plan kind;
- support counts and retained dimensions are the expected Be2 compact route
  values;
- final basis materializes with dimension 221 and `:core_then_shell` support
  row order;
- H1 materializes, is finite/symmetric, and has the reviewed lowest energy
  `-0.27746109235228694`;
- Ham input carries the positive-weight density gauge and raw-numerator pair
  convention;
- private handoff carries H1, density interaction, nuclear metadata, electron
  count, and spin sector;
- readiness remains blocked only on
  `:missing_diatomic_hamiltonian_consumer_contract`.

Manager validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

passed 49/49 in 47.6s.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed.

`git diff --check HEAD~1..HEAD` passed.

Deletion/shrinkage accounting:

- deleted: no-probe blocked testset; repeated field inventories; repeated
  nonclaim flags; old one-center blocked Ham payload assertions
- simplified: remaining Be2 focused test is now a compact route-spine
  semantic/numerical fingerprint
- quarantined: Be2 check remains focused/manual, not default nested-suite
  pressure; H1/J, dense `Vee`, RHF, WL, CR2/HFDMRG readiness, public API,
  export/HamV6, artifacts, and physics endpoint claims remain out of scope
- not deleted because: the compact probe-enabled test is still the active
  Be2/PQS route-spine gate
- exact remaining caller/blocker: route readiness remains blocked at
  `:missing_diatomic_hamiltonian_consumer_contract`

Next: audit the missing Hamiltonian consumer contract against downstream CR2
and HFDMRG expectations before implementing it.

-- repo-manager@macmini
