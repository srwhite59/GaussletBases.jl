# One-Center Atomic Nested Shell Contract Note

This note supersedes the earlier one-center `ns7/ns9` localization diagnosis
and the earlier compressed-shell atomic default.

## Correct Atomic Contract

For the one-center atomic Cartesian nested line, the intended legacy/W&L shell
contract is:

- full parent coverage
- `working_box = (1:n, 1:n, 1:n)`
- complete shell layers peeled until the inner direct cube reaches `nside`
- shell increment
  - `ns^3 - (ns - 2)^3`
- face retains
  - `(ns - 2) × (ns - 2)` on all six faces
- edge retains
  - `ns - 2` on all twelve edges
- corners
  - direct carry-through on all eight corners

So for `ns = 7`, every complete shell layer must retain exactly:

- faces: `6 * 5^2 = 150`
- edges: `12 * 5 = 60`
- corners: `8`
- total increment: `150 + 60 + 8 = 218`

This is the canonical one-center atomic shell contract now used by the repo.

## What Was Wrong Before

Two earlier ideas were wrong for the supported one-center atomic path:

- a central inferred working-box diagnostic fixture
- compressed atomic shell defaults like `(4,3)` / `3`

The central-box fixture was the wrong atomic object because it did not preserve
full parent coverage. The compressed atomic shell defaults were a mistaken
reading of the legacy/W&L shell contract.

The generic internal rectangular-shell builder still supports explicit retained
counts for non-atomic uses, but that compressed contract is no longer the
canonical or default supported behavior for:

- `build_one_center_atomic_full_parent_shell_sequence(...)`
- `one_center_atomic_full_parent_fixed_block(...)`

## Canonical Repo-Native Path

The supported one-center atomic helpers are now:

- `build_one_center_atomic_full_parent_shell_sequence(...)`
- `one_center_atomic_full_parent_fixed_block(...)`

These helpers always build on the full parent cube and always use the
legacy/W&L complete-shell retention counts for the chosen `nside`.

The old local central-box diagnostic helper is now explicitly quarantined as:

- `wrong_central_box_atomic_fixture`

and should not be reused for one-center atomic contraction diagnosis.

## Repo-Native Structure Diagnostics

The repo now exposes the atomic shell/core structure directly through:

- `one_center_atomic_nested_structure_diagnostics(...)`
- `one_center_atomic_nested_structure_report(...)`

That diagnostics path reports, at minimum:

- parent side count
- working-box side count
- `nside`
- core side count
- number of shell layers
- expected shell increment
- expected and actual face / edge / corner retained counts
- total expected gausslet count
- total actual gausslet count

This makes basis-size discrepancies separable into:

- shell contract
- working-box choice
- supplement choice

## Structural Anchors

The structural contract is now pinned directly in repo tests:

- `ns = 5`
  - increment `= 5^3 - 3^3 = 98`
- `ns = 7`
  - increment `= 7^3 - 5^3 = 218`
- `27^3` working box with `ns = 7`
  - shell layers `= 10`
  - total gausslet count `= 343 + 10 * 218 = 2523`

The `27^3` count is intentionally a shell-structure statement, independent of
any later supplement choice.

## Diagnosis Boundary After Hardening

The earlier one-center `ns7/ns9` bug hunt was aimed at the wrong object because
it used a central-box fixture. That diagnosis is superseded.

Future one-center atomic contraction work should use only the full-parent
helpers above, together with the structure diagnostics, so the shell contract
and working-box choice stay explicit.
