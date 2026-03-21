The current convergence check for the base QW-PGDG route is blocked by
residual-space robustness logic, not by runtime.

That is now the narrow target:

- the `count = 9` base QW-PGDG nearest/GGT path is already in a good physical
  regime
- runtime is already short enough that `count = 11` and `count = 13` are
  reasonable probe cases
- but the higher-count runs were being blocked by two robustness problems:
  - the residual keep threshold was too aggressive
  - the `:ggt_nearest` route was still coupled to unconditional MWG residual
    width extraction

The intended fix is small:

1. lower the residual keep threshold to a true null-mode scale
2. let `:ggt_nearest` depend only on residual centers
3. leave MWG width extraction only on the `:mwg` branch

The current implementation uses a simple scale-aware residual keep rule:

- absolute floor `1e-8`
- relative floor `0.1 * maximum(residual_overlap_eigenvalues)`

so the retained residual space is filtered by

- `max(1e-8, 0.1 * λmax)`

rather than by a single hard absolute cutoff.

The physical interpretation is unchanged. Residual Gaussians in the PGDG route
are expected to be subtle and lightly occupied, so the keep criterion should
remove genuine null modes, not physically meaningful residual directions at the
`1e-3` to `1e-2` level.

This is therefore a robustness pass, not a new architecture pass. The goal is
simply to unlock the intended `count = 9 / 11 / 13` convergence check on the
existing base QW-PGDG nearest/GGT path.
