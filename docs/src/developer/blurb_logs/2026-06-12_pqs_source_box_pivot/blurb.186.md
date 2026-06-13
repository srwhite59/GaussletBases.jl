Pass 186 - add compact WL baseline deltas to the PQS He gate, paid for by deleting stale term-storage scaffolding.

Purpose:

Pass 185 fixed the tracked PQS He fixture so it uses the WL parent mapping and
now asserts:

```text
PQS final dimension = 419
PQS H1 lowest = -1.991334820314074
PQS H1-J self-Coulomb = 1.2420423900074902
```

Pass 183 recovered the matching old WL gausslet-only 419 baseline:

```text
WL H1 lowest = -1.991344469963435
WL H1-orbital self-Coulomb = 1.2420473874925473
WL RHF total = -2.85080350301779
```

The next tracked gate should make the comparison hygiene explicit without
turning this into a broad WL reconstruction test.

Task type:

Implementation cleanup. Add a tiny scalar comparison to the focused PQS He test
and delete stale integration-only test scaffolding.

Primary live physics target:

```text
He q=5/n_s=5 gausslet-only, same WL parent lattice and same 419 shell inventory.
```

Exact task:

1. In:

   ```text
   test/nested/pqs_direct_retained_final_h1_runtests.jl
   ```

   add compact WL-baseline delta assertions near the existing H1 and H1-J
   fingerprints. Use the pass-183 WL constants directly; do not reconstruct the
   WL route in this test.

   Suggested values:

   ```text
   WL H1 lowest = -1.991344469963435
   WL H1 self-Coulomb = 1.2420473874925473
   PQS H1 - WL H1 = 9.649649361120893e-6
   PQS H1-J self-Coulomb - WL self-Coulomb = -4.997485057112172e-6
   ```

   Keep it compact. The point is to preserve the matched-WL comparison signal,
   not to create a second WL fixture in this file.

2. Delete the stale integration-only test file:

   ```text
   test/nested/one_center_atomic_compact_fixed_block_term_storage_runtests.jl
   ```

   and remove its include from:

   ```text
   test/nested/integration_runtests.jl
   ```

   Why this deletion is allowed:

   - It is an integration/slow scaffolding test, not a default physics gate.
   - Its `!hasproperty(..., :term_storage/:gaussian_terms/:pair_terms)` checks
     are already covered in
     `test/nested/cartesian_nested_fixed_block_qw_pgdg_adapter_runtests.jl`
     near the top of that file.
   - Its remaining nested glass-box fixed-dimension checks are not the active
     He/PQS Hamiltonian contract.
   - It is preserving old fixed-block storage vocabulary more than it is
     protecting current physics.

3. Do not delete accepted endpoint checks or the current focused PQS He gate.

Trust boundary:

- Do not edit production source.
- Do not add RHF to tracked tests.
- Do not reconstruct WL inside the PQS test.
- Do not run Be2/Cr2/HFDMRG/DMRG/artifact generators.
- Do not touch supplemented WL 447 behavior.
- Do not preserve the deleted test through an adapter or renamed include.

Line-budget rule:

This pass must be net-negative in tracked source/test/generator lines:

```text
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Acceptance condition:

```text
total deleted > total added
```

This pass should be comfortably net-negative because it deletes the stale
49-line integration test.

Validation:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Do not run the slow integration suite as a routine gate for this deletion.

Report:

- exact WL constants/deltas added to the focused PQS test;
- deleted file and removed include;
- validation output;
- source/test/generator line additions/deletions/net;
- deletion/shrinkage report:
  - deleted;
  - simplified;
  - quarantined;
  - not deleted because;
  - exact remaining caller/blocker.

Commit and push if validation passes.

Write the result to `.agent_handoffs/response.186.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.186.md
```

-- repo-manager@macmini
