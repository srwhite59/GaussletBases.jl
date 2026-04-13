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
- `build_one_center_atomic_legacy_profile_shell_sequence(...)`
- `one_center_atomic_legacy_profile_fixed_block(...)`

These helpers always build on the full parent cube and always use the
legacy/W&L complete-shell retention counts for the chosen `nside`.

The contract split is now explicit:

- modern canonical path
  - full parent coverage
  - `working_box = (1:n, 1:n, 1:n)`
- legacy-profile reproduction path
  - explicit inner working box such as `(2:28, 2:28, 2:28)` on a `29^3`
    parent lattice
  - same exact shell increment and face/edge/corner counts

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

That distinction is now repo-native for the common `ns = 7` comparison:

- parent `29^3`, full-parent working box
  - total gausslets `= 2741`
- parent `29^3`, legacy-profile inner `27^3` working box
  - total gausslets `= 2523`

For the literal one-center Ne legacy-profile comparison point

- `Z = 10`
- `d = 0.03`
- parent side count `= 29`
- legacy-profile working box `= (2:28, 2:28, 2:28)`
- `nside = 7`
- supplement `= repo-v6z-sp`, `lmax = 1`

the repo now reproduces the structural gausslet profile exactly:

- legacy-profile gausslet count `= 2523`
- raw supplement orbital count `= 25`

But the old default supplement merge path still did **not** reproduce the
legacy final basis count:

- expected legacy total basis count `= 2548`
- previous default merged basis count `= 2531`

So the shell contract and working-box contract were aligned, and the remaining
one-center Ne discrepancy had moved to the supplement merge / residual-space
layer rather than the nested shell-retention layer.

## Legacy-Profile Ne Residual-Space Boundary

The remaining `2548 - 2531 = 17` count gap is now localized to the residual
supplement keep rule inside:

- `_qwrg_residual_space(...)`

For the anchored one-center Ne legacy-profile case:

- parent side `= 29`
- working box `= (2:28, 2:28, 2:28)`
- `nside = 7`
- supplement `repo-v6z-sp`, `lmax = 1`

the supplement-side structure is:

- raw supplement orbital count `= 25`
- supplement overlap numerical rank `= 25`
- residualized supplement overlap numerical rank before keep `= 25`

So the reduction does **not** happen in the raw supplement loader or in the
raw overlap assembly. It happens at the explicit residual keep stage:

- `keep_tol = max(1e-8, 0.1 * maximum(residual_overlap_eigenvalues))`
- anchored Ne value: `keep_tol = 1.9458934805275673e-4`

On that case:

- kept residual count `= 8`
- discarded residual count `= 17`
- maximum discarded residual eigenvalue
  - `1.9270157729463577e-4`
- minimum kept residual eigenvalue
  - `1.9955105831311726e-4`

The discarded directions are therefore not numerically null in the usual
sense:

- residual numerical-null threshold used for diagnosis
  - `1e-12`
- discarded directions above that null threshold
  - `17 / 17`

So the repo is not losing directions because the residualized supplement space
collapses to rank `8`. The current repo keeps only `8` because the present
relative keep rule drops every residualized direction below `10%` of the
largest residual overlap eigenvalue.

That diagnosis is now implemented directly in the repo.

The downstream naming has now been cleaned up:

- `nested_profile = :legacy_profile`
  controls the inner working-box geometry contract
- `residual_keep_policy = :near_null_only`
  controls supplement residual-direction retention

Those are independent contracts and no longer share the same preferred public
name.

The canonical atomic supplement keep rule is now:

- `residual_keep_policy = :near_null_only`
- keep if orthogonalized residual-overlap eigenvalue `> 1e-8`

The old name:

- `residual_keep_policy = :legacy_profile`

remains accepted only as a compatibility alias for `:near_null_only`.

The retained residual block is now explicitly stabilized after selection:

- if `R` is the selected raw residual coefficient block, the repo forms
  `Skeep = R' * Sraw * R`
- then replaces `R` by
  `R * inv(sqrt(Symmetric(Skeep)))`

This keeps all selected near-null-surviving directions while making the
retained residual block orthonormal again to numerical precision before:

- `raw_to_final`
- residual center extraction
- residual width extraction
- downstream one-body / interaction transforms

That closes the later `ns = 9` nested atomic failure boundary: the problem was
not another keep-policy bug, but the fact that the kept full-rank residual
block was no longer orthonormal enough for downstream center extraction after
selection.

On the anchored Ne legacy-profile case:

- kept residual count `= 25`
- total basis count `= 2548`
- low one-body ladder
  - `[-49.99990524677932, -12.499971049252997, -12.49997104925268, -12.499971049252615, -12.49997099214756, -5.555472605898918, -5.555472605898094, -5.5554726058980854]`

So the remaining one-center Ne basis-count gap was indeed due to aggressive
residual truncation, not shell law, working-box choice, or raw supplement rank
loss. The later `ns = 9` atomic blocker was a residual-block orthonormality
issue after selection, and is now addressed by explicit stabilization rather
than reintroducing aggressive pruning.

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
