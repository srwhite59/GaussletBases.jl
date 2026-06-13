Pass 167 - revise CR2 inspection artifact schema before commit

Role: repo-doer@macmini

Continue from the current uncommitted pass-166 working tree. Repo-manager has
not accepted or committed pass 166 yet.

The pass-166 implementation is close and validates, but it needs one correction
before acceptance:

- the White-Lindsey route group is too thin for the "same schema" requirement;
- the focused Be2 test lost all compact source-plan/final-basis semantic
  guards;
- line budget should not be met by over-compressing the new artifact writer.

Keep the current private writer direction:

- private JLD2 + TSV inspection artifact;
- PQS route populated from the current Be2 handoff;
- White-Lindsey unavailable/not-applicable placeholders;
- no public API/export;
- no HamV6, dense final-space `V`, dense four-index `Vee`, `V6`, `Vblocks`,
  solver bundle, HF/RHF/DMRG, H1/J promotion, Qiu-White correction
  implementation, CR2 run, HFDMRG run, or downstream solver.

Required corrections:

1. Make the White-Lindsey route group schema-shaped.

   It does not need WL arrays yet, but it should have the same top-level group
   families as PQS where practical:

   - `routes/white_lindsey/route`
   - `routes/white_lindsey/readiness`
   - `routes/white_lindsey/system`
   - `routes/white_lindsey/final_basis`
   - `routes/white_lindsey/one_body`
   - `routes/white_lindsey/two_body`
   - `routes/white_lindsey/validation`

   Use explicit `:unavailable`, `:not_applicable`, `nothing`, empty arrays, or
   compact blocker labels. Do not synthesize PQS pre-final/source-box data for
   WL.

2. Restore a compact semantic guard for the PQS route spine.

   Do not bring back the old 50+ assertion wall. Add only a few checks or
   artifact keys that protect the live semantics:

   - source plan object kind is `:pqs_diatomic_complete_core_shell_source_plan`;
   - support order is `(:product, :pqs_left, :pqs_right)`;
   - route retained order is `(:pqs_left, :pqs_right, :product)`;
   - final basis dimension is `221`;
   - final support row order is `:core_then_shell`.

   These may be direct test checks or compact artifact fields. Keep the H1
   finite/symmetric/lowest-energy checks.

3. Do not make schema construction hard to read just to satisfy line budget.

   Compact argument listings are fine. Removing structural newlines from
   schema/control-flow blocks is not. Prefer clear grouped NamedTuples/Dict
   entries even if they cost a few lines.

Line-budget permission:

- The line-budget rule still applies: final tracked `src` + `test` diff must be
  net-negative by `git diff --numstat -- src test`.
- You may pay for the budget by deleting stale source/test pressure outside the
  active Be2 file.
- Explicitly allowed deletion candidate:

  ```text
  test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
  ```

  This is a remaining private synthetic RHF validator. RHF is not the product,
  CR2/HFDMRG own serious HF, and this file is not included in the default test
  runner. If deleting it is needed to keep the artifact writer readable and
  line-negative, delete it.

Do not delete:

- H1 dimension/finite/symmetry/lowest-energy checks;
- the focused Be2 artifact smoke;
- scientific endpoint tests;
- public export tests unrelated to this lane.

Validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

If you delete the RHF SCF test, also report whether any runner included it. Do
not run a broad suite just to prove deletion of a non-default stale test.

Report:

- final status and commit SHA if you commit; otherwise say uncommitted;
- validation results;
- exact line-budget numbers;
- whether the WL group now has same top-level schema families;
- whether the compact PQS route-spine semantic guard is present;
- readability notes if you intentionally kept compact argument formatting;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
