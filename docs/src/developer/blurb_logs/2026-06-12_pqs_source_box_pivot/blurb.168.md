Pass 168 - private Be2 PQS CR2 artifact generator path

Role: repo-doer@macmini

Implement the next CR2-requested step: a repo-owned private artifact-generation
path under ignored `tmp/work`, so CR2 does not need to call underscored writer
helpers directly.

Context:

CR2 accepted the pass-167 JLD2+TSV schema as the right first private read-only
Be2/PQS inspection artifact seam. CR2 requested next:

1. add a private durable artifact-generation path;
2. fill producer/fixture/backend provenance;
3. later populate White-Lindsey under the same schema;
4. later add supplement/residual-GTO and correction metadata placeholders.

This pass should do only the first step and cheap provenance needed by that
step.

Allowed implementation surface:

- A tracked script under `tmp/work/`, if the repo's ignore rules allow tracking
  explicit files there with `git add -f`; or a small tracked script under a
  repo-owned script/dev path if tracking under ignored `tmp/` is too awkward.
- Existing private helpers in
  `src/pqs_source_box_diatomic_complete_core_shell.jl`.
- Existing focused test
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  only if needed for a small smoke.

Decision rule:

- Prefer a tracked generator script that writes ignored outputs to:

  ```text
  tmp/work/be2_wl_pqs_cr2_inspection_artifact/
  ```

- The generated JLD2/TSV outputs must remain ignored and must not be committed.
- If tracking a file under ignored `tmp/` is awkward, use a tracked location
  like `dev/` or `scripts/` only after checking existing repo conventions.

Generator behavior:

- Build the same probe-enabled Be2/PQS assembly used by the focused test.
- Call the private pass-167 payload/writer.
- Write:

  ```text
  tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
  tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
  ```

- Print both paths and a short readiness summary.
- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense final-space `V`,
  dense four-index `Vee`, `V6`, `Vblocks`, solver/export code, or Qiu-White
  corrections.

Cheap provenance to add if straightforward:

- producer repo commit;
- dirty-state marker;
- generator entrypoint/path;
- generated timestamp;
- fixture basics already known in the generator: atom symbols, q, n_s,
  reference/tail spacing, parent axis counts, map backend, route kind.

Keep this provenance compact. Do not add a large field cloud to route reports.
If adding provenance to the private payload would make the pass grow, put it in
the generator's payload call or artifact metadata, not in unrelated route
readiness summaries.

Line-budget rule:

- If you edit tracked `src/` or `test/`, final `src` + `test` diff must be
  net-negative by `git diff --numstat -- src test`.
- If this pass only adds/edits a tracked script outside `src` and `test`, report
  that the source/test line-budget rule is not triggered.
- Do not delete scientific endpoint tests or public export tests to pay for
  this.

Validation:

Run the generator once with:

```text
julia --project=. <generator path>
```

Then verify:

```text
ls -l tmp/work/be2_wl_pqs_cr2_inspection_artifact
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git status --short --branch
```

If generated outputs appear in `git status`, stop and fix ignore/tracking
behavior before reporting success.

Report:

- generator path;
- generated output paths;
- top-level JLD2/fingerprint confirmation;
- producer/fixture provenance fields included;
- whether source/test line-budget was triggered;
- validation results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
