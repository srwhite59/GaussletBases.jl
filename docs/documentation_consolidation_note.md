# Documentation Consolidation Note

This note records the purpose of the current consolidation pass.

## 1. What duplication still remains

The recent structure pass fixed the main navigation problem, but some
duplication still remains.

The main overlaps are now:

- `README.md`, `docs/index.md`, and `docs/architecture.md` still all describe
  the branch structure at different lengths
- the atomic mean-field notes still read as a flat stack unless a reader
  already knows which ones are current and which ones are supporting
- the ordinary PGDG notes still read as a long development chain unless a
  reader already knows the current branch interpretation

So the remaining problem is less about missing pages and more about how the
supporting notes are grouped.

## 2. What now looks mergeable

The two clearest note chains to compress are:

- the atomic mean-field chain
  - `atomic_ida_direct.md`
  - `atomic_ida_exchange.md`
  - `atomic_ida_fock.md`
  - `atomic_ida_spin_fock.md`
- the ordinary PGDG development chain
  - `ordinary_pgdg_decision.md`
  - `ordinary_pgdg_comx.md`
  - `ordinary_pgdg_proxy_refinement.md`
  - `ordinary_pgdg_distortion_regime.md`
  - `ordinary_pgdg_backend_pivot.md`

Those notes are still useful, but they no longer need to compete directly with
the current branch-status pages.

## 3. What this pass should do

This pass should reduce repetition without deleting development history:

- keep the current top-level skeleton
- add one small synthesis page for the atomic mean-field chain
- add one small synthesis page for the ordinary PGDG chain
- point the current branch-status pages to those synthesis pages
- trim one remaining overlap point in `docs/architecture.md`
- add lightweight supporting-note status lines where they help

## 4. What this pass should not do

This pass should not:

- redesign the scientific story
- archive large numbers of notes
- broaden the public claims of the ordinary or atomic branches
- loosen the current wording discipline around the radial numerical path,
  ordinary PGDG experimental status, or provisional mapping heuristics

## 5. Intended result

After this pass, the docs should feel stable enough to leave mostly alone until
the science changes again:

- top-level pages for navigation and current status
- synthesis pages for the main supporting-note chains
- individual development notes retained, but clearly demoted
