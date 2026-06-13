Pass 184 - isolate whether the PQS/WL 419 gap is a parent-lattice mismatch.

Purpose:

Pass 183 recovered the matching old WL gausslet-only He baseline:

```text
WL source = old_nested_fixed_block_wl_gausslet_only
Z = 2
d = 0.3
q = n_s = 5
final dimension = 419
RHF total = -2.85080350301779
H1 lowest = -1.991344469963435
```

The current corrected PQS 419 diagnostic is:

```text
PQS source = shellification/lowering-backed fixed-q complete-core/shell path
final dimension = 419
RHF total = -2.8493871224303167
H1 lowest = -1.9866819751748936
```

The raw diagnostic delta is therefore:

```text
PQS RHF total - WL RHF total = +0.0014163805874733981 Ha
PQS H1 - WL H1 = +0.004662494788541416 Ha
```

Do not interpret that as a true PQS-vs-WL contraction result yet. I checked the
current fixture source, and the parent lattices/mappings are probably not the
same:

```text
test/nested/pqs_direct_retained_final_h1_runtests.jl
  _pqs_h1_test_bundle(count)
    uses AsinhMapping(a = 0.25, s = asinh(xmax / 0.25) / (...), tail_spacing = 10.0)
    with xmax = 8.0

pass-183 WL probe
  uses white_lindsey_atomic_mapping(Z = 2, d = 0.3, tail_spacing = 10.0)
  mapping = AsinhMapping(a = 0.38729833462074165,
                         s = 0.7745966692414834,
                         tail_spacing = 10.0)
```

The intended comparison is:

```text
same parent lattice / same shellification inventory
WL and PQS diverge only in the contraction/source representation after that
```

So this pass should answer one narrow question:

```text
If PQS uses the same one-center parent lattice/mapping as the old WL q=5/n_s=5
gausslet-only baseline, what are the PQS 419 H1/H1-J/private RHF diagnostics?
```

Task type:

This is a no-tracked-source/test diagnostic probe pass. Use ignored
`tmp/work` scripts and summaries. Do not edit production source, tracked tests,
tracked generators, docs, or artifacts except the response file.

Exact task:

1. Build a local ignored PQS probe that mirrors the current corrected PQS 419
   path, but replaces `_pqs_h1_test_bundle(11)`'s parent axis construction with
   the WL mapping:

   ```julia
   build_basis(MappedUniformBasisSpec(
       :G10;
       count = 11,
       mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.3, tail_spacing = 10.0),
       reference_spacing = 1.0,
   ))
   ```

   Keep the corrected fixed-q shellification/lowering path:

   ```text
   OneCenterShellification(core_side = 5, q = 5)
   PQSLowering(q = 5)
   source mode shape = (5, 5, 5)
   three shell layers
   retained per shell = 98
   final dimension = 419
   ```

2. Compute the same diagnostics already used for the current PQS pass:

   ```text
   parent mapping summary
   parent physical endpoints
   source/final-basis inventory
   final overlap identity error
   H1 lowest energy
   H1 symmetry/finiteness
   H1-J self-Coulomb
   private RHF diagnostic total
   private RHF one-body and two-body components
   RHF convergence/iteration/residual/trace/idempotency
   elapsed time
   ```

3. Compare the WL-mapped PQS result to both:

   ```text
   old current-PQS custom-mapping result
   old WL 419 gausslet-only result from pass 183
   ```

   Report deltas:

   ```text
   WL-mapped PQS RHF total - WL RHF total
   WL-mapped PQS H1 - WL H1
   WL-mapped PQS H1-J self-Coulomb - WL H1 self-Coulomb
   WL-mapped PQS RHF two-body - WL RHF two-body
   WL-mapped PQS values - current custom-mapping PQS values
   ```

Decision rules:

- If the WL-mapped PQS result moves close to WL, then the pass-183 gap was
  mostly a parent lattice/mapping mismatch. Report that; do not change tracked
  tests yet.
- If the WL-mapped PQS result remains about 1 mHa or worse above WL, then the
  next investigation is likely the PQS contraction/source representation or
  one-body construction, not just the parent lattice. Report which scalar
  component dominates.
- If the probe cannot be built without tracked source changes, write
  `.agent_handoffs/ATTENTION.md` and stop. Do not patch production code in this
  pass.

Trust boundary:

- Do not tune WL or PQS.
- Do not add or update tracked tests.
- Do not edit `src/`.
- Do not run CR2, HFDMRG, DMRG, Be2, or artifact generators.
- Do not use the supplemented WL 447 result as the comparison baseline.
- Do not treat private RHF as a product solver; it is only a Hamiltonian
  diagnostic.

Line-budget/diff rule:

- `git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  should be empty.
- If tracked `src`, `test`, or tracked generator changes become necessary,
  stop with `ATTENTION.md`.

Suggested artifacts:

- Create ignored files similar to:

  ```text
  tmp/work/pqs_he_419_wl_parent_mapping_probe.jl
  tmp/work/pqs_he_419_wl_parent_mapping_probe_summary.txt
  ```

Validation:

```text
julia --project=. tmp/work/pqs_he_419_wl_parent_mapping_probe.jl
git status --short --branch
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Report:

- probe path and summary path;
- exact parent mapping parameters and endpoints;
- source/final-basis inventory;
- H1/H1-J/RHF diagnostics;
- deltas versus pass-183 WL 419;
- deltas versus current custom-mapping PQS 419;
- whether the likely next issue is parent mapping, PQS contraction/source
  representation, or one-body construction;
- no tracked source/test/generator changes.

Write the result to `.agent_handoffs/response.184.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.184.md
```

-- repo-manager@macmini
