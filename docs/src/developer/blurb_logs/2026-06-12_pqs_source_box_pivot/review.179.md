Accepted pass 179.

This is the right correction for the atom-first PQS target. The implementation
now separates the shell support box from the fixed PQS source-mode shape:

- PQS lowering records `source_mode_shape = (q, q, q)` and preserves the
  physical support box shape separately as `source_box_shape`.
- Region-backed multi-layer PQS source planning carries each layer's
  `source_mode_shape` into realization.
- Source realization uses that route-owned shape when present, while the
  explicit-box bridge still falls back to `length.(inner_box)`.

The observed inventory is the one we needed before doing any physics:

```text
parent side = 11
core support = 125
shell layers = 3
shell support = 1206
shell raw_source_dims = (5,5,5) for all layers
shell retained_count = 98 for all layers
shell final retained = 294
final retained = 419
```

Validation was focused and passed:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

The source/test diff is net-negative by 94 lines, with directly related
scaffold removed from the old compact H1 test.

One caveat for the next pass: this file is now an inventory/final-basis gate,
not an H1 gate. That is acceptable for pass 179 because H1 was explicitly out
of scope, but it should not linger as a misleading "final_h1" test. The next
engineering step should either materialize H1 on this 419-dimensional fixed-q
basis or rename/reshape the test around its actual contract after H1 coverage
is restored elsewhere.

-- repo-manager@macmini
