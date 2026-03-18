# Ordinary PGDG Supporting Notes

This page is the short synthesis page for the ordinary PGDG supporting-note
chain.

For the current ordinary branch status, read first:

- [`docs/current_ordinary_branch.md`](current_ordinary_branch.md)

## What this supporting chain is about

The present ordinary PGDG interpretation was built through a sequence of
development notes:

1. deciding that the Coulomb expansion length was not the real bottleneck
2. testing COMX/localization
3. refining the primitive proxy
4. checking the distortion regime
5. pivoting to the current backend split

Those notes are still valuable because they explain why the current ordinary
branch is interpreted the way it is.

## Recommended supporting-note order

Read them in this order:

1. [`docs/ordinary_pgdg_decision.md`](ordinary_pgdg_decision.md)
2. [`docs/ordinary_pgdg_comx.md`](ordinary_pgdg_comx.md)
3. [`docs/ordinary_pgdg_proxy_refinement.md`](ordinary_pgdg_proxy_refinement.md)
4. [`docs/ordinary_pgdg_distortion_regime.md`](ordinary_pgdg_distortion_regime.md)
5. [`docs/ordinary_pgdg_backend_pivot.md`](ordinary_pgdg_backend_pivot.md)

Then, if you want the later branch-state follow-ons, read:

6. [`docs/ordinary_pgdg_localized_backend.md`](ordinary_pgdg_localized_backend.md)
7. [`docs/ordinary_pgdg_one_body_fidelity.md`](ordinary_pgdg_one_body_fidelity.md)
8. [`docs/ordinary_pgdg_hybrid_regime.md`](ordinary_pgdg_hybrid_regime.md)
9. [`docs/ordinary_sho_spectral_test.md`](ordinary_sho_spectral_test.md)
10. [`docs/ordinary_pgdg_hybrid_consolidation.md`](ordinary_pgdg_hybrid_consolidation.md)

## Current relevance

These notes are mostly supporting development notes, not current workflow
pages.

Their value is that they preserve the decisions that led to the current
ordinary-branch interpretation:

- friendly hybrid regime is the practical target
- strong/small-`c` pure mapped cases are stress tests
- `AsinhMapping` and coupled `c,s` heuristics are still provisional

## What may be mergeable later

If the ordinary branch stays stable for a while, the most likely future merge
targets are:

- the early development chain from `ordinary_pgdg_decision.md` through
  `ordinary_pgdg_backend_pivot.md`
- the later branch-state chain from `ordinary_pgdg_localized_backend.md`
  through `ordinary_pgdg_hybrid_consolidation.md`
