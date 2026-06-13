Pass 186 complete.

Files changed:

- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `test/nested/integration_runtests.jl`
- deleted `test/nested/one_center_atomic_compact_fixed_block_term_storage_runtests.jl`

WL constants and deltas added to the focused PQS He gate:

- WL H1 lowest: `-1.991344469963435`
- WL H1 self-Coulomb: `1.2420473874925473`
- asserted `PQS H1 - WL H1 ≈ 9.649649361120893e-6`
- asserted `PQS H1-J self-Coulomb - WL self-Coulomb ≈ -4.997485057112172e-6`
- both delta checks use `atol = 1.0e-10, rtol = 0.0`

Deleted file and removed include:

- deleted stale slow integration scaffold `test/nested/one_center_atomic_compact_fixed_block_term_storage_runtests.jl`
- removed its include from `test/nested/integration_runtests.jl`

Validation:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
```

passed:

```text
PQS fixed-q complete core-shell inventory gate | 34 pass, 34 total, 6.4s
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
0	1	test/nested/integration_runtests.jl
0	49	test/nested/one_center_atomic_compact_fixed_block_term_storage_runtests.jl
5	0	test/nested/pqs_direct_retained_final_h1_runtests.jl
```

Line budget:

- added: `5`
- deleted: `50`
- net: `-45`

Git status at validation point:

```text
## main...origin/main
 M test/nested/integration_runtests.jl
 D test/nested/one_center_atomic_compact_fixed_block_term_storage_runtests.jl
 M test/nested/pqs_direct_retained_final_h1_runtests.jl
```

Deletion/shrinkage report:

- deleted: stale slow integration-only term-storage scaffold and its integration runner include.
- simplified: focused PQS He gate now records the matched-WL H1 and H1-J comparison signal with two direct scalar delta assertions, without reconstructing WL.
- quarantined: none.
- not deleted because: the focused PQS He gate remains the active 419-dimensional seam validation; the QW-PGDG adapter test still carries the live `term_storage`/`gaussian_terms`/`pair_terms` absence checks.
- exact remaining caller/blocker: no remaining include or tracked caller for `one_center_atomic_compact_fixed_block_term_storage_runtests.jl`; no blocker for deleting that file.

Commit/push: pending immediately after this response write per the pass blurb.

-- repo-doer@macmini
