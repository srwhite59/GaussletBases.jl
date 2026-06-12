Accepted.

This pass implemented the intended first PQS electron-nuclear kernel: a raw
product-source by-center block from caller-supplied term-first Gaussian factor
arrays. The sign and charge convention is correct for the new route:

```text
sum_t (-coefficients[t]) * Gx[t] * Gy[t] * Gz[t]
```

The result records the center and nuclear charge, but does not apply charge and
does not sum centers. The retained wrapper is only the existing retained
source-mode selector contraction, which is exactly the right scope for this
pass.

The common raw product-source metadata extraction was simplified and reused by
the new nuclear helper. No old CCPM helper was adopted as route authority, and
no source-axis Gaussian generator, global matrix, H1 diagnostic, IDA,
Hamiltonian, or driver path was added.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Identify and, if possible, implement the module-owned source of term-first 1D
Gaussian factor arrays for PQS source modes. If the current PQS route only
carries source dimensions and not enough axis data to compute centered Gaussian
factors, block precisely instead of inventing a fallback.

-- repo-manager@macmini
