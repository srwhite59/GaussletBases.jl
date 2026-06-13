Pass 163 - add CR2 read-only inspection view line-negatively

Role: repo-doer@macmini

Task type: source/test implementation with strict line-negative budget.

Purpose:

Expose the current Be2/PQS private Hamiltonian handoff as a compact CR2
read-only inspection view, without claiming CR2 solver/export readiness. This
should help the future CR2 agent compare Be2 WL and PQS route metadata/H1 facts
without needing a solver-ready `V` contract yet.

Hard constraints:

- Do not add a new payload type.
- Do not add a new test file.
- Do not make CR2, HFDMRG, HamV6, export, or artifact readiness true.
- Do not change the active readiness blocker away from:

```text
:missing_hfdmrg_density_density_contract
```

Line-budget rule:

This pass edits `src/` and `test/`, so it must be net-line-negative across
tracked `src/` + `test` files.

Measure:

```text
git diff --numstat -- src test
```

Acceptance condition:

```text
sum(deleted) > sum(added)
```

Target:

- at least 25 deleted lines;
- no more than 20 added lines;
- net reduction of at least 5 lines.

If this cannot be done safely, write `.agent_handoffs/ATTENTION.md` and stop.

Implementation seam:

Use existing objects in:

```text
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Add a compact read-only inspection view as a summary/readiness field on the
existing consumer contract and/or readiness summary. Do not create
`_PQS...CR2...Payload` or any new struct.

Suggested compact fields:

```text
cr2_read_only_inspector_ready = true
cr2_solver_ready = false
cr2_export_ready = false
cr2_handoff_blocker = :missing_cr2_solver_handoff_format
two_body_representation_kind = :pre_final_density_interaction
density_gauge = :pre_final_localized_positive_weight
raw_pair_factor_convention = :raw_numerator
```

Keep these as compact summary/readiness labels derived from the existing
handoff. Do not copy large matrices or object references into new places.

Pay for the additions by deleting/consolidating remaining duplicated nonclaim
flags in the handoff layer:

- handoff conventions/summary/metadata still repeat false flags for CR2/HFDMRG/
  HamV6/export readiness;
- use `_pqs_source_box_route_driver_diatomic_downstream_ham_nonclaims()` where
  helpful;
- if a handoff-level false flag is not tested and is already represented by
  consumer/readiness, delete it.

Test update:

Only update the existing compact test:

```text
test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Add or adjust minimal assertions:

- CR2 read-only inspector is ready;
- CR2 solver/export readiness remains false;
- CR2 handoff blocker is `:missing_cr2_solver_handoff_format`;
- two-body representation is explicitly `:pre_final_density_interaction`;
- density gauge and raw pair convention are explicitly named;
- overall readiness blocker remains `:missing_hfdmrg_density_density_contract`.

Remove any old assertion that only preserved now-deleted duplicate flags.

Forbidden:

- No new payload/struct.
- No new test file.
- No CR2 run.
- No HFDMRG run.
- No solver-ready `H,V` claim.
- No dense `Vee`, final-space `V`, `V6`, `Vblocks`.
- No H1/J materialization.
- No RHF/SCF/DIIS.
- No WL comparison run.
- No public API.
- No export/HamV6/JLD2/artifact writing.
- No docs/handoff line-budget gaming.

Validation:

Run:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

Decision rules:

- If a CR2 read-only view requires a new payload type, stop and write
  `ATTENTION.md`; that would violate this pass.
- If line budget is not negative across `src` + `test`, do not commit; write
  `ATTENTION.md`.
- If the readiness blocker changes away from
  `:missing_hfdmrg_density_density_contract`, stop and report.

Report back:

- where the CR2 read-only view lives;
- status/blocker labels added;
- confirmation that CR2 solver/export readiness remains false;
- `git diff --numstat -- src test` totals:
  - src/test added:
  - src/test deleted:
  - net:
- validation commands/results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Commit if validation passes and line budget is negative, with a message like:

```text
Add CR2 read-only Hamiltonian inspection view
```

-- repo-manager@macmini
