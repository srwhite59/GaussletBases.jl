Accepted with a manager fix.

Pass 023 added the intended CPBM helper:

```text
pqs_source_shell_final_one_body_from_boundary_matrix(...)
```

It consumes a `pqs_retained_source_one_body_matrix(...)`-style result in
`:retained_pqs_source_modes`, not an arbitrary raw source matrix, and computes:

```text
O_final = L' * O_boundary * L
```

for overlap and kinetic only. It does not claim electron-nuclear, H1,
Hamiltonian assembly, IDA, density-density, RHF, driver, export, or artifact
behavior.

Manager fix:

The first implementation normalized `:source_overlap` and `:source_kinetic`,
but the live `pqs_retained_source_one_body_matrix(...)` path uses:

```text
:retained_source_overlap
:retained_source_kinetic
```

I patched the term normalizer and compact test to use the live retained-source
labels while still accepting the older source labels as aliases.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next blocker remains:

```text
:missing_pqs_shell_boundary_electron_nuclear_operator_source
```

Before extending the final-boundary helper to nuclear or assembling H1, the next
pass should audit the existing PQS retained centered electron-nuclear by-center
path on the real projected-q-shell fixture. It must prove the retained-source
nuclear block equals the shell-projected boundary nuclear operator, by center,
with charge recorded but not applied and centers not summed.

-- repo-manager@macmini
