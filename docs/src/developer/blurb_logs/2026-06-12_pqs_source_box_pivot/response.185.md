Pass 185 complete.

Files changed:

- `test/nested/pqs_direct_retained_final_h1_runtests.jl`

Exact mapping now used by the tracked fixture:

```julia
white_lindsey_atomic_mapping(Z = 2.0, d = 0.3, tail_spacing = 10.0)
```

Final dimension and retained inventory:

- current box: `(1:11, 1:11, 1:11)`
- core box: `(4:8, 4:8, 4:8)`
- raw source dims: `(5, 5, 5)`
- shell layer count: `3`
- retained per shell: `98`
- core support count: `125`
- shell support count: `1206`
- shell final retained count: `294`
- final retained dimension: `419`
- final overlap identity error remains checked as `< 1.0e-10`

H1 and H1-J scalar fingerprints now asserted:

- `h1.lowest_energy ≈ -1.991334820314074` with `atol = 1.0e-10, rtol = 0.0`
- `h1_j_summary.self_coulomb ≈ 1.2420423900074902` with `atol = 1.0e-10, rtol = 0.0`

Validation:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
```

passed:

```text
PQS fixed-q complete core-shell inventory gate | 32 pass, 32 total, 6.4s
```

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed with `load ok`.

```text
git diff --check
```

passed with no whitespace errors.

```text
git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

reported:

```text
6	28	test/nested/pqs_direct_retained_final_h1_runtests.jl
```

Line budget:

- added: `6`
- deleted: `28`
- net: `-22`

Git status at validation point:

```text
## main...origin/main
 M test/nested/pqs_direct_retained_final_h1_runtests.jl
```

Deletion/shrinkage report:

- deleted: old custom `xmax = 8.0`, `a = 0.25` mapping setup; redundant H1/H1-J metadata nonclaim checks; redundant loose finite/positive scalar checks superseded by exact scalar fingerprints.
- simplified: focused He fixture now uses the WL parent mapping directly, with compact physics-level H1 and H1-J fingerprints.
- quarantined: none.
- not deleted because: inventory, matrix finiteness/symmetry, final dimension, density gauge, and density-interaction materialization checks still protect live scientific and convention contracts.
- exact remaining caller/blocker: focused tracked gate still directly constructs H1 and H1-J payloads because RHF is intentionally not part of this tracked test and the route still needs these diagnostics as the current He seam validation.

Commit/push: pending immediately after this response write per the pass blurb.

-- repo-doer@macmini
