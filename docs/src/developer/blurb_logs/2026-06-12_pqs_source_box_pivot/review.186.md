Accepted pass 186.

Manager inspected the diff and reran:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
```

which passed:

```text
PQS fixed-q complete core-shell inventory gate | 34 pass, 34 total, 6.3s
```

The focused PQS He gate now carries the matched-WL comparison signal directly:

```text
WL H1 lowest = -1.991344469963435
PQS H1 - WL H1 = +9.649649361120893e-6

WL H1 self-Coulomb = 1.2420473874925473
PQS H1-J self-Coulomb - WL self-Coulomb = -4.997485057112172e-6
```

Doer also deleted:

```text
test/nested/one_center_atomic_compact_fixed_block_term_storage_runtests.jl
```

and removed its slow integration-runner include. This was the right deletion:
the deleted file mostly preserved old fixed-block storage vocabulary, while the
live absence checks are already covered in the QW/PGDG adapter test.

Tracked source/test/generator line budget:

```text
5 added, 50 deleted, net -45
```

The tracked q=5 He gate now says the important thing: with the same parent
lattice and 419-dimensional shell inventory, PQS and old WL agree at the
10-microhartree H1/H1-J scale. The remaining private RHF diagnostic from pass
184 is PQS lower than WL by about 14 microhartree, but RHF remains diagnostic,
not a product solver.

-- repo-manager@macmini
