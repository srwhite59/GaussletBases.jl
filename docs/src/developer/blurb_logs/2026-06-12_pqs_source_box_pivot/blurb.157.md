Pass 157 - shrink Be2 PQS Ham payload fingerprint test

Role: repo-doer@macmini

Task type: test deletion/shrinkage implementation.

Purpose:

Shrink the oversized Be2/PQS Ham payload fingerprint test so it protects the
live route spine instead of preserving every transitional private payload field.
This is the second test-carrying-cost cleanup after retiring RHF seam pressure.

Governing policy:

- `AGENTS.md` test scope/deletion policy: tests are code with runtime,
  maintenance, and conceptual cost; development scaffolding tests should be
  deleted once superseded.
- `docs/src/developer/pqs_source_box_operator_framework.md`: source-box-first
  PQS is the algorithmic frame; shell/support-row contraction remains
  oracle/debug.
- `docs/src/developer/pqs_source_box_fixture_policy.md`: compact route facts are
  route-smoke/convention diagnostics, not physics endpoints.

Current state:

`test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
is about 947 lines. It is not in the default nested runner, but it is the
focused active Be2/PQS route-spine check. It currently has:

- a large no-probe blocked testset that mostly preserves old missing-object and
  nonclaim vocabulary;
- a large probe-enabled testset that checks every intermediate payload field.

Exact task:

Shrink only:

```text
test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Do not touch source files.

Delete the blocked/no-probe testset:

```text
@testset "Be2 PQS Ham payload readiness fingerprint"
```

unless you find one no-probe assertion that is still a live route contract and
not covered by the probe-enabled route. If you keep any no-probe check, keep it
to a tiny boundary smoke and explain why.

Shrink the probe-enabled testset to compact semantic checks for the active
route spine:

1. Assembly/route identity:
   - `assembly.object_kind == :cartesian_assembly`
   - route family is `:pqs_source_box`

2. Source-plan semantic smoke:
   - source plan exists and has object kind
     `:pqs_diatomic_complete_core_shell_source_plan`
   - it does not claim `:pqs_multilayer_shell_source_plan`
   - core/support counts are `25`, `250`, shell retained count `196`
   - precleanup/final retained dimension is `221`
   - shell coefficient matrix shape is `(250, 196)`
   - support order is `(:product, :pqs_left, :pqs_right)`
   - route retained order is `(:pqs_left, :pqs_right, :product)`
   - core and shell support indices are disjoint

3. Final-basis semantic smoke:
   - final basis materializes
   - final dimension is `221`
   - support row order is `:core_then_shell`
   - old one-center source-plan object kind is false

4. H1 numerical/semantic check:
   - H1 payload materializes
   - final dimension is `221`
   - final one-body Hamiltonian matrix is finite and symmetric/hermitian within
     a tight tolerance
   - lowest H1 energy is finite; if the existing value is stable in the test,
     assert it against the reviewed current fingerprint
     `-0.27746109235228694` with a reasonable tight tolerance

5. Ham input / handoff readiness check:
   - Ham input payload materializes
   - density gauge is `:pre_final_localized_positive_weight`
   - raw pair-factor convention is `:raw_numerator`
   - support weight count is `275`
   - pre-final pair matrix shape is `(221, 221)`
   - Hamiltonian handoff payload materializes as private inspect-only
   - one-body Hamiltonian and density interaction references are carried
   - nuclear charges `(4.0, 4.0)`, coordinates
     `((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))`, nuclear repulsion `4.0`,
     electron count `8`, and spin sector `:closed_shell_singlet` are recorded
   - readiness remains blocked only on
     `:missing_diatomic_hamiltonian_consumer_contract`

Delete:

- repeated available/missing object inventories at every stage;
- repeated materialization false flags;
- repeated public/export/artifact/RHF/WL/H1-J false checks;
- old one-center `complete_core_shell_ham_payload` blocked-input assertions
  unless one remains essential to a live contract;
- route-shadow blocker assertions that were only useful before source-plan,
  final-basis, H1, Ham-input, and handoff payloads existed.

Trust boundary:

- Test-only deletion/shrinkage.
- No production/source changes.
- No new tests unless a much smaller helper inside the same file replaces
  deleted scaffolding.
- No route-driver behavior change.
- No Hamiltonian consumer implementation, dense `Vee`, H1-J, RHF, WL, CR2,
  HFDMRG, public API, export/HamV6, artifact, fixture promotion, or physics
  endpoint work.

Validation:

Run the focused Be2 file even though it can take about one minute. This is the
representative gate for the active Be2/PQS route spine after deleting hundreds
of assertions:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

Do not run the full nested default suite for this shrink pass unless a concrete
ambiguity appears and you explain why first.

Decision rules:

- If the compact probe-enabled route no longer passes, stop and report the
  first failing semantic fact; do not add back the whole old assertion wall.
- If the reviewed H1 energy value is not stable, keep only finite/symmetric H1
  checks and report the observed value.
- If the file cannot be shrunk substantially while preserving the route-spine
  check, stop with the exact blocker instead of making cosmetic edits.

Report back:

- old and new line count;
- old and new test count if easy to report;
- whether the no-probe testset was deleted or reduced;
- which semantic groups remain;
- validation commands/results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Commit if validation passes, with a message like:

```text
Shrink Be2 PQS Ham fingerprint test
```

-- repo-manager@macmini
