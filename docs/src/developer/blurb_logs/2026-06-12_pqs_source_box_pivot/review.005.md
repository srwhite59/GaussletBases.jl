Accepted with a manager cleanup.

The pass added the first one-unit retained PQS source-mode global one-body
matrix helper. It accepts retained source-mode overlap/kinetic batch results,
requires exactly one materialized self-pair block, rejects skipped ready
self-pairs, and returns the retained block as the dense matrix. That is the
right scope for this stage:

```text
retained source-mode self-pair block
-> one-unit retained source-mode matrix
```

It does not perform multi-unit placement, shell realization, Lowdin cleanup,
electron-nuclear, IDA, Hamiltonian assembly, or driver adoption.

During review I trimmed one result-shape issue before committing. The helper
initially returned a flat set of nonclaim booleans for electron-nuclear,
density-density, IDA, Hamiltonian, driver, exports, and artifacts. That was
more field-cloud-like than needed, so the committed form keeps explicit
shell/Lowdin drift guards and groups the broader nonclaims in
`nonclaim_capabilities`.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Audit the electron-nuclear source-factor seam before implementing it. PQS H/He+
one-electron diagnostics need nuclear attraction, but that is the first
non-safe one-body term and should reuse existing Coulomb Gaussian expansion
conventions rather than invent a new route.

-- repo-manager@macmini
