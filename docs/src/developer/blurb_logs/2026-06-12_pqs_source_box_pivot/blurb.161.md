Pass 161 - shrink diatomic Hamiltonian consumer duplication

Role: repo-doer@macmini

Task type: source/test shrink implementation.

Purpose:

Reduce the carrying cost introduced by the private diatomic Hamiltonian
consumer contract while preserving the state transition from pass 159:

```text
readiness blocker = :missing_hfdmrg_density_density_contract
```

The consumer contract should remain private inspect-only, but it should stop
copying handoff-owned scalar fields that have no live downstream caller.

Hard line-budget rule:

This pass edits `src/` and `test/`, so it must be net-line-negative across
tracked `src/` + `test/` files.

Measure at the end with:

```text
git diff --numstat -- src test
```

Acceptance condition:

```text
sum(deleted) > sum(added)
```

Target:

- at least 70 deleted lines;
- no more than 30 added lines;
- preferred net reduction: 65+ lines.

Do not count docs/handoff files. Do not satisfy this by deleting accepted
scientific endpoint tests, deleting only whitespace/comments, or moving code to
untracked/tmp files. If this cannot be done safely, write
`.agent_handoffs/ATTENTION.md` with the exact blocker and stop.

Exact task:

1. In `src/pqs_source_box_diatomic_complete_core_shell.jl`, thin
   `_PQSDiatomicCompleteCoreShellHamiltonianConsumerContractPayload`.

   Keep only fields needed for the live contract:

   - `status`
   - `blocker`
   - `route_family`
   - `source_handoff`
   - `source_handoff_status`
   - `readiness`
   - `available_objects`
   - `missing_objects`
   - `summary`
   - `metadata`

   Delete copied scalar fields such as final dimension, one-body reference,
   one-body status, representation kind/status, density gauge, raw pair
   convention, support counts, pair shapes, nuclear metadata, electron count,
   and spin sector from the consumer struct and constructor call. These already
   live in the handoff and handoff summary.

2. Delete consumer-builder scalar-copy code that exists only to populate those
   copied fields. Use the handoff summary directly if a compact summary still
   needs one or two facts.

3. Consolidate repeated downstream nonclaim flags.

   Add at most one small private helper or local bundle, for example:

   ```julia
   _pqs_source_box_route_driver_diatomic_downstream_ham_nonclaims()
   ```

   It should centralize the false readiness flags and downstream missing
   objects:

   - `hfdmrg_density_density_ready = false`
   - `hfdmrg_sliced_ready = false`
   - `hamv6_export_ready = false`
   - `cr2_ready = false`
   - `public_api = false`
   - `exports_materialized = false`
   - `artifacts_materialized = false`
   - missing objects for HFDMRG density-density, sliced integrals, HamV6, and
     CR2 format

   Keep it tiny. If the helper costs more lines than it saves, use a local
   named tuple instead.

4. Keep readiness behavior unchanged:

   - consumer contract available on the probe-enabled Be2 path;
   - readiness blocker remains `:missing_hfdmrg_density_density_contract`;
   - `:diatomic_hamiltonian_consumer_contract` remains available, not missing;
   - no downstream readiness flag becomes true.

5. In
   `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`,
   remove assertions that only preserve copied consumer scalar fields. Keep only:

   - consumer payload available;
   - `private_inspector_ready`;
   - source handoff identity;
   - one compact check that downstream readiness flags are false;
   - readiness blocker replacement.

Forbidden:

- No new payloads.
- No new test file.
- No dense `Vee`.
- No final-space `V`.
- No `V6`/`Vblocks`.
- No H1/J materialization.
- No RHF/SCF/DIIS.
- No WL comparison.
- No CR2 or HFDMRG run.
- No public API.
- No export/HamV6/JLD2/artifact writing.
- No docs-only line-budget gaming.

Validation:

Run:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

The focused Be2 test may take about one minute; it is the representative gate
for this shrink because the pass changes the route readiness/consumer contract.

Decision rules:

- If the readiness blocker changes away from
  `:missing_hfdmrg_density_density_contract`, stop and report.
- If the line-budget condition is not met, do not commit; write
  `.agent_handoffs/ATTENTION.md`.
- If a copied scalar field appears to have a real source caller outside tests,
  preserve it and report the caller; otherwise delete it from the consumer.

Report back:

- fields removed from the consumer payload;
- downstream nonclaim consolidation approach;
- old and new Be2 test assertion count if easy;
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

Commit if validation passes and the line budget is negative, with a message like:

```text
Thin diatomic Hamiltonian consumer contract
```

-- repo-manager@macmini
