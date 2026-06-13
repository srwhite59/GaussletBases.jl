Pass 157 response - shrink Be2 PQS Ham payload fingerprint test

Implemented and committed:

```text
bba678e2 Shrink Be2 PQS Ham fingerprint test
```

File changed:

- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Shrinkage:

- old line count: 947
- new line count: 195
- diff stat: 27 insertions, 779 deletions
- old focused pass count from the pre-shrink file: 575 total
  (`225` no-probe + `350` probe-enabled)
- new focused pass count: 49

No-probe testset:

- Deleted the no-probe blocked testset entirely.
- I did not keep a no-probe boundary smoke because the retained live contract
  is the probe-enabled source-plan -> final-basis -> H1 -> Ham-input/handoff
  route spine. The no-probe path mostly preserved missing-object and nonclaim
  vocabulary for stages that now exist in the probe-enabled route.

Semantic groups remaining:

- assembly/route identity;
- source-plan semantic smoke: object kind, old source-plan nonclaim,
  core/shell support counts, shell retained count, precleanup dimension,
  shell coefficient shape/finite values, support order, route retained order,
  and disjoint core/shell support indices;
- final-basis semantic smoke: materialized payload, final dimension, support
  row order, precleanup dimension, and old source-plan object-kind nonclaim;
- H1 numerical/semantic check: materialized payload, final dimension, finite
  symmetric one-body Hamiltonian, finite lowest energy, and reviewed H1 energy
  `-0.27746109235228694` with `atol = 1.0e-12`;
- Ham input / handoff readiness: density gauge, raw pair-factor convention,
  support weight count, pre-final pair shape, private inspect-only handoff,
  carried one-body Hamiltonian and density interaction references, nuclear
  charges/coordinates/repulsion, electron count, spin sector, and readiness
  blocked only on `:missing_diatomic_hamiltonian_consumer_contract`.

Validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

passed:

```text
Be2 PQS probe-enabled Ham readiness fingerprint | 49 pass | 47.9s
```

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed:

```text
load ok
```

```text
git diff --check
```

passed with no output.

Git status after commit:

```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:

- deleted: no-probe blocked testset; repeated available/missing inventories;
  repeated materialization false flags; repeated public/export/artifact/RHF/WL/
  H1-J nonclaim checks; old one-center complete-core/shell Ham payload blocked
  assertions; route-shadow blocker assertions superseded by materialized
  source-plan/final-basis/H1/Ham-input/handoff payloads
- simplified: probe-enabled testset now checks compact route-spine semantics
  instead of every intermediate private payload field
- quarantined: remaining Be2/PQS check is still focused/manual, not part of the
  default nested runner; H1/J, dense `Vee`, RHF, WL, CR2/HFDMRG, public API,
  export/HamV6, artifacts, and physics endpoint claims remain out of scope
- not deleted because: the probe-enabled Be2 route-spine check remains the
  active focused gate for source-plan -> final-basis -> H1 -> Ham handoff
  readiness
- exact remaining caller/blocker: route readiness remains blocked at
  `:missing_diatomic_hamiltonian_consumer_contract`

-- repo-doer@macmini
