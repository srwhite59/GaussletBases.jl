Accepted pass 180.

This restores the H1 gate on the corrected 419-dimensional fixed-q one-center
PQS basis and deletes the old standalone compact 223 Ham payload scaffold.

Observed fixed-q He H1 smoke:

```text
final dimension = 419
lowest H1 energy = -1.9866819751748936
matrix finite = true
symmetry max abs error = 2.3092638912203256e-14
solve kind = ordinary_symmetric
generalized-overlap solve = false
```

Validation passed:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

The source/test line budget is net-negative by 80 lines:

```text
35 added, 115 deleted
```

The deleted file was an old compact 223 route-driver Ham payload test that was
not in the default nested runner and had no live source/test caller. Removing
it is consistent with the current atom-first target and with the repo policy
against preserving old route-shadow scaffold.

Next step should be H1-J/density interaction on the same 419 basis, still
private and diagnostic. Do not move to RHF or Be2/Cr2 until the 419 H1-J
convention is checked.

-- repo-manager@macmini
