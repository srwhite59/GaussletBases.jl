Pass 206 - shrink the PQS projected q-shell local-layer integration test

Role:
You are `repo-doer@macmini` implementing one bounded cleanup pass for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, and `BlurbStyle.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `86e9ef9c Extract low order materialization`
- The target file is:
  `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
- It is about 6476 lines and begins:
  `# Integration/slow test. Do not include in default nested runner.`
- It is valuable historically, but it now carries too much development-era
  route-shadow vocabulary and repeated nonclaim assertions.

Goal:
Shrink the giant integration test without weakening the live mathematical
contracts.

This is a test-retirement/shrink pass, not a source refactor.

Preserve live checks for:
- cubic q=5/L=5 PQS geometry and retained-mode counts, especially boundary COMX
  product-mode count 98;
- rectangular q=5/L=7 analog, especially retained count 130;
- raw product-box plan contracts, boundary selector, source mode ordering, and
  explicit raw product-box reference comparisons;
- PQS/PQS source-box self/cross block comparisons where they validate formulas;
- density-density normalization semantics, including density-normalized vs
  raw-weighted pair factors;
- center-resolved nuclear attraction sign/charge convention and finite/symmetry
  checks.

Remove or collapse development-era scaffolding where safe:
- repeated `private_shadow_only`, `fixture_only`, `production_supported == false`
  assertions;
- repeated `qwhamiltonian_consumes == false`,
  `public_default_consumes == false`, `packet_adoption == false`,
  `fixed_block_routing == false`, and similar no-go flags;
- assertions whose only role is preserving old route-shadow vocabulary rather
  than detecting a numerical/formula bug;
- the block explicitly labeled:
  `Legacy smoke only: CPB provider tests own detailed one-body product math.`
  This block calls `_product_doside_source_box_shadow_blocks`; delete it if no
  nearby live mathematical check depends on it.

Do not:
- split this into several new test files in this pass;
- add a new test unless it replaces substantially more old lines in the same
  pass;
- edit source code unless an exact caller search proves a helper is now dead and
  the deletion is smaller/risk-free;
- delete live formula/reference comparisons for raw product-box, PQS/PQS,
  density-density, or nuclear-attraction behavior;
- touch H1/H1-J/RHF, driver, artifact, WL, CR2, or source module code.

Line-count rule:
The active line-count rule still applies. For edits under `src`, `test`, `bin`,
and the CR2 generator script, final tracked diff must be net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Acceptance condition:
`sum(deleted) > sum(added)` for those scoped files.

Expected shape:
- This pass should be substantially line-negative in the target test.
- If you cannot remove at least a meaningful chunk of scaffold assertions
  without compromising live contracts, write `.agent_handoffs/ATTENTION.md` and
  stop rather than making a tiny cosmetic edit.

Validation:
Because this pass edits the slow integration test itself, running that test is
the correct validation even if it may exceed 60 seconds. Explain/report elapsed
time. Prefer Julia-level timing, for example:

```sh
julia --project=. -e 't = @elapsed include("test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl"); println("elapsed_s=", t)'
```

If the test has harness assumptions, include the minimal imports needed, similar
to:

```sh
julia --project=. -e 'using Test, LinearAlgebra, SparseArrays, GaussletBases; t = @elapsed include("test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl"); println("elapsed_s=", t)'
```

Also run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.206.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.206.md
```

Report:
- lines removed/added in the target test;
- which assertion families were removed/collapsed;
- which live mathematical contracts were preserved;
- whether `_product_doside_source_box_shadow_blocks` still has any test/source
  pressure after this pass;
- validation commands and elapsed time;
- source/test/bin scoped line budget added/deleted/net;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Stop after writing the response. Manager will review, commit, push, and pause.

-- repo-manager@macmini
