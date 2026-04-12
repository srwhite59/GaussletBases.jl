# Homonuclear Diatomic Policy Correction

This note records the first correction to the landed bond-aligned diatomic box
policy on the Cartesian QW line.

The relevant commits are:

- `8a414a3` `Preserve midpoint slab in homonuclear diatomic split`
- `8433f91` `Preserve odd shared-shell retains on symmetric homonuclear intervals`

Together, these commits correct the first homonuclear `H2` geometry line in two
ways:

- odd bond-axis working intervals now keep the midpoint as a shared
  `nx × ny × 1` slab rather than assigning that row to one child
- symmetric tangential shared-shell intervals now keep an odd retained local
  side count, so the localized side space retains a near-zero center

These are intentionally narrow corrections:

- they apply to the current bond-aligned homonuclear diatomic route
- they do not settle heteronuclear midpoint handling
- they do not settle the general adaptive local-side-count rule

## Why The Correction Was Needed

The first visualization/debug pass on `H2` at `R = 1.4` showed that the old
midpoint split on an odd bond-axis interval produced asymmetric child boxes
after one shared shell:

- left child `7 × 7 × 6`
- right child `7 × 7 × 5`

That let only one child continue through another shell-contraction stage and
produced an avoidable homonuclear asymmetry in the compressed fixed-center
picture.

The subsequent local `doside` / `COMX` trace then showed that the remaining
shared-shell missing-center behavior came from even retained local side spaces
on symmetric intervals:

- interval `2:8`
- retained count `4`
- localized centers
  `[-2.6740914029, -0.9276287253, 0.9276287253, 2.6740914029]`
- no near-zero center

Odd retained local count `3` on the same kind of symmetric interval restored
the near-zero center.

## What Is Now Landed

### 1. Shared midpoint slab

For a bond-aligned homonuclear diatomic whose bond-axis working interval has
odd length:

- reserve the midpoint as a shared `nx × ny × 1` slab
- split the remainder into equal left/right child boxes
- keep the slab as a direct shared region between the child subtrees

For the representative `H2`, `R = 1.4` case, the corrected post-shared-shell
geometry is now:

- left child `(2:8, 2:8, 2:6)` = `7 × 7 × 5`
- shared slab `(2:8, 2:8, 7:7)` = `7 × 7 × 1`
- right child `(2:8, 2:8, 8:12)` = `7 × 7 × 5`

The left/right child physical widths are now symmetric.

### 2. Odd shared-shell retains on symmetric tangential intervals

In the current bond-aligned homonuclear shared shell only:

- if a tangential local interval is symmetric about zero
- and the provisional face retain count is even
- reduce that retain count by one

This is confined to the shared shell:

- not edges
- not child/core local retains
- not heteronuclear cases
- not chains

It is a narrow homonuclear correction, not a general adaptive local-retain
formula.

## H2 Validation

The correction was confirmed at two bond lengths.

### `R = 1.4`

With the corrected homonuclear shared-shell behavior:

- fixed dimension: `655 -> 637`
- overlap `||S-I||∞`: `5.995e-15 -> 5.995e-15`
- projected fixed-only capture:
  `0.9999972527022969 -> 0.999997252702297`
- projected fixed-only energy:
  `-1.279762295532455 -> -1.279762295532455`
- projected fixed-only `⟨Vee⟩`:
  `0.7772582435002245 -> 0.777258329001004`
- hybrid nearest/GGT `E1`:
  `-1.2837079015541708 -> -1.2837079015541661`
- hybrid nearest/GGT `⟨Vee⟩`:
  `0.7802897903866277 -> 0.7802898627470577`

The representative nested `xz` projection became visibly more centered.

### `R = 2.0`

The same narrow correction again behaved well:

- fixed dimension: `597 -> 579`
- overlap `||S-I||∞`:
  `1.3100631690576847e-14 -> 1.3100631690576847e-14`
- projected fixed-only capture:
  `0.9981550127852054 -> 0.998155012785205`
- projected fixed-only energy:
  `-1.0688846589800904 -> -1.0688846589800904`
- projected fixed-only `⟨Vee⟩`:
  `0.6596694587811625 -> 0.659669666366456`
- hybrid nearest/GGT `E1`:
  `-1.1016387068355626 -> -1.1016387068355642`
- hybrid nearest/GGT `⟨Vee⟩`:
  `0.66188159470718 -> 0.6618817788926099`

Again, the traced shared-shell contractions regained the near-zero center and
the geometry picture improved.

## What Remains Out Of Scope

This correction does not yet settle:

- heteronuclear midpoint handling
- linear-chain midpoint policy
- a general adaptive local-side-count rule
- edge or child/core odd-retain rules
- arbitrary non-linear geometries

## Roadmap Consequence

The first bond-aligned homonuclear diatomic line is now corrected at the two
policy points that the `H2` debug artifacts exposed:

- symmetric midpoint handling
- symmetric shared-shell center retention

So the next geometry-family work should move outward from this corrected
homonuclear baseline rather than reopening those two settled corrections.
