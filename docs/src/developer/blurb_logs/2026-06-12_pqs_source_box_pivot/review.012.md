Accepted.

Pass 012 added the intended composition only:

```text
centered Gaussian source-mode factors
-> existing uncharged by-center source electron-nuclear block
-> optional existing retained source-mode contraction
```

The wrappers preserve the current source-box convention: by-center matrix,
nuclear charge recorded but not applied, centers not summed, and no shell
realization, Lowdin cleanup, IDA, Hamiltonian, driver, export, or artifact
claim. The implementation did not call the old CCPM source-box wrappers.

The test coverage remains in the existing CPBM source-pair contract file. It
checks that the centered helper matches the already-tested supplied-factor path
and that the retained wrapper matches `pqs_source_pair_retained_one_body_block`.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Run a probe-first one-center PQS retained-source H1 readiness check using the
existing retained overlap, kinetic, and centered nuclear helpers. This should
be a `tmp/work` probe unless a missing public seam is exposed. Any generalized
solve must be labeled raw retained-source diagnostic only, not final acceptance.

-- repo-manager@macmini
