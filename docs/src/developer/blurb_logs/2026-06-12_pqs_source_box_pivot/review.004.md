Accepted.

The pass stayed properly narrow. It added a retained PQS source-mode one-body
selector and batch selector for the two requested terms:

```text
:overlap
:kinetic
```

The implementation delegates through the existing raw source-space one-body
selector and then applies the retained source-mode contraction from pass 003.
Position and x2 were intentionally left out, which is the right choice for
this pass because the next consumer only needs overlap and kinetic.

The test extension is compact enough: it checks selector equivalence to the
term-specific retained helpers, batch ready/skipped counts, retained block
space, and the no-shell/no-Lowdin boundary. It does not add a new file or a
new broad metadata suite.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Build the smallest retained PQS source-mode global one-body matrix layer:
one retained PQS unit, one self-pair, overlap and kinetic only. This should
consume the retained source-mode batch result directly and place it into a
dense retained matrix with dimension equal to the retained source-mode count.
Do not generalize to multi-unit placement, shell realization, electron-nuclear,
IDA, Hamiltonian assembly, or drivers in this pass.

-- repo-manager@macmini
