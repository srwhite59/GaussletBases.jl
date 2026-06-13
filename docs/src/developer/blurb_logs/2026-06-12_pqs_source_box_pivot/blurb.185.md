Pass 185 - align the tracked PQS He 419 fixture to the WL parent mapping.

Purpose:

Pass 184 showed that the mHa-scale PQS/WL gap was mostly not a PQS contraction
failure. It was a parent-lattice/mapping mismatch:

```text
old tracked PQS fixture:
  AsinhMapping(a = 0.25, s = asinh(xmax / 0.25) / (...), tail_spacing = 10.0)

old WL 419 baseline:
  white_lindsey_atomic_mapping(Z = 2.0, d = 0.3, tail_spacing = 10.0)
```

When PQS used the same WL parent mapping, the 419-dimensional diagnostics were:

```text
PQS H1 lowest = -1.991334820314074
PQS H1-J self-Coulomb = 1.2420423900074902
PQS private RHF diagnostic total = -2.850817886618113
```

and the deltas versus WL 419 were small:

```text
RHF total delta = -1.4383600322798173e-5 Ha
H1 delta = +9.649649361120893e-6 Ha
H1-J self-Coulomb delta = -4.997485057112172e-6 Ha
```

The tracked focused PQS He gate should therefore use the WL parent mapping.
That makes the tracked fixture match the intended comparison rule:

```text
same parent lattice / same shellification inventory
WL and PQS diverge only after that
```

Task type:

Implementation cleanup in the focused tracked test only.

Primary file:

```text
test/nested/pqs_direct_retained_final_h1_runtests.jl
```

Exact task:

1. Change `_pqs_h1_test_bundle(count::Int)` so it builds its axis basis with:

   ```julia
   white_lindsey_atomic_mapping(Z = 2.0, d = 0.3, tail_spacing = 10.0)
   ```

   Keep:

   ```text
   MappedUniformBasisSpec(:G10; count = 11, reference_spacing = 1.0)
   OneCenterShellification(core_side = 5, q = 5)
   PQSLowering(q = 5)
   fixed source mode shape = (5, 5, 5)
   final retained dimension = 419
   ```

2. Do not preserve the old custom `xmax = 8.0`, `a = 0.25` mapping path in this
   focused fixture. Delete those local setup lines rather than keeping an
   option or compatibility branch.

3. Replace loose numerical smoke assertions with compact scalar fingerprints:

   ```text
   h1.lowest_energy ≈ -1.991334820314074
   h1_j_summary.self_coulomb ≈ 1.2420423900074902
   ```

   Use tolerances tight enough to catch drift but not brittle to harmless
   floating-order noise; `atol = 1.0e-10, rtol = 0.0` is fine unless you see a
   concrete reason to loosen it.

4. Keep the inventory checks that protect the corrected shellification:

   ```text
   shell layer count = 3
   source dims per shell = (5,5,5)
   retained per shell = 98
   final dimension = 419
   final overlap identity error small
   ```

5. Shrink low-value metadata/nonclaim assertions in this file to pay for the
   new scalar fingerprints. Good deletion candidates:

   ```text
   repeated exports/artifacts/nonclaim checks in both H1 and H1-J sections
   redundant finite/positive checks now superseded by exact scalar fingerprints
   old custom mapping setup lines
   ```

   Preserve checks that catch real numerical or convention bugs: matrix
   finiteness, symmetry, final dimension, source dims, retained counts, density
   gauge, and density-interaction materialization.

Trust boundary:

- Do not edit production source in this pass.
- Do not add RHF to the tracked test.
- Do not edit Be2/Cr2 artifacts or generators.
- Do not run HFDMRG/DMRG/CR2.
- Do not add a second fixture or compatibility option for the old custom
  mapping.
- Do not compare to supplemented WL 447.

Line-budget rule:

This implementation pass must be net-negative in tracked source/test/generator
lines:

```text
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Acceptance condition:

```text
total deleted > total added
```

If you cannot make the pass net-negative without deleting live scientific checks,
write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Report:

- exact mapping now used by the tracked fixture;
- final dimension and retained inventory;
- H1 and H1-J scalar fingerprints;
- tests run and pass counts;
- source/test/generator line additions/deletions/net;
- deletion/shrinkage report:
  - deleted;
  - simplified;
  - quarantined;
  - not deleted because;
  - exact remaining caller/blocker.

Commit and push if validation passes.

Write the result to `.agent_handoffs/response.185.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.185.md
```

-- repo-manager@macmini
