Pass 137 - audit diatomic PQS source-plan producer contract

Purpose:

Do a no-edit contract audit for the missing Be2/diatomic PQS complete-core/shell
source-plan producer. The next implementation should be based on a clear object
contract, not on forcing the diatomic route skeleton through the one-center
shellification producer.

Why now:

Pass 136 showed that `probe_parent_axis_construction = :auto` makes the parent
axis-bundle object available for the focused Be2/PQS route. The remaining
readiness blocker is therefore:

```text
:missing_diatomic_complete_core_shell_source_plan_producer
```

This is the main seam between current Be2/PQS route skeletons and a usable
private Ham payload for the medium-term Be2 WL/PQS comparison that CR2 can later
exercise externally.

Audit surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`
  - diatomic readiness payload from pass 135
  - one-center complete-core/shell source-plan payload
  - final-basis/H1/H1-J/Ham payload handoff
- `src/pqs_source_box_route_driver_skeletons.jl`
- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_pqs_product_source_box_route_skeleton`
- `src/pqs_multilayer_shell_region_plan.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- Focused Be2 fingerprint test:
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- Governing docs:
  - `docs/src/developer/pqs_source_box_operator_framework.md`
  - `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  - `docs/src/developer/pqs_source_box_fixture_policy.md`

Questions to answer:

1. What exact fields of `pqs_multilayer_shell_source_plan` are consumed by:
   - `pqs_multilayer_complete_core_shell_final_basis`
   - `pqs_multilayer_complete_core_shell_h1_payload`
   - density-input/H1-J helpers?
2. Which of those fields correspond to universal route concepts that a
   diatomic source-box producer can provide, and which are one-center
   shellification-specific?
3. What structured Be2/PQS objects already exist in the route skeleton:
   source boxes, retained units/ranges, pair entries, retained dimension,
   route-axis counts, parent axis bundle, center metadata, Coulomb expansion
   inputs?
4. Should the next producer return:
   - a new diatomic source-plan object;
   - a common abstract/duck-typed source-plan shape consumed by existing final
     basis/H1 helpers;
   - or a blocked/readiness payload plus a later final-basis adapter?
5. What is the smallest implementation pass that would advance the private Ham
   payload without making shell/support-row contraction route authority?
6. What existing tests can be reused, and what new test pressure can be avoided?

Deliverable:

Write a precise recommendation for the next coding pass. Include:

- proposed object/helper names;
- required input objects;
- required output fields or summary fields;
- explicit blockers;
- whether the first implementation should materialize anything or only return a
  structured blocked/available source-plan payload;
- focused validation command(s);
- deletion/shrinkage forecast.

Trust boundary:

- No source edits.
- No test edits.
- No commits except the handoff response.
- No final-basis/H1/H1-J/Ham materialization.
- No RHF/SCF work.
- No WL payload, public API, exports, artifacts, hfdmrg, or CR2 execution.
- No scalar report-field clouds.
- Do not force one-center shellification/support-row semantics onto the
  diatomic source-box route.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not request interactive command approval during unattended baton work. If
  approval would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- Read-only inspection should be enough.
- Use `rg`, `sed`, `git show`, and similar read-only commands.
- Do not run tests unless a short focused probe is necessary to answer the
  contract question.
- If any command is expected to exceed 60 seconds, report why instead of
  running it by default.

Report back:

- Files/helpers inspected.
- Existing source-plan fields and consumers.
- Existing Be2 route-skeleton facts relevant to a producer.
- Recommended source-plan producer contract.
- Smallest next implementation pass and validation.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
