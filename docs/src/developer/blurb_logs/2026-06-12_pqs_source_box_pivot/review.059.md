Review 059: accepted; repeated one-cell shell composition is the right first
multi-layer PQS seam.

The audit found that the one-cell guard in `cartesian_nested_faces.jl` is real,
but it does not force a broad rewrite. The viable route-owned shape is:

```text
core:    (4:10)^3
shell 1: (3:11)^3 / (4:10)^3
shell 2: (2:12)^3 / (3:11)^3
shell 3: (1:13)^3 / (2:12)^3
```

The prototype showed disjoint shell supports, full parent coverage, shell
isometry errors around `1e-14` to `3e-14`, and a collapsed shell sector that
works with the existing complete core/shell final-basis helper:

```text
core support:       343
shell support:      1854
shell retained:     1206
final dimension:    1549
final overlap err:  5.51e-13
```

This is a good place to implement one narrow route seam. Keep the driver spine
in mind: this belongs in the shells/transforms/final-basis route lifecycle, not
as a separate private physics lane. Also keep the spacing rule provisional: PQS
may not need the same q/ns intuition as WL, and any `Z`, `d`, `s`, radius, and
shell-depth rule needs later review rather than a new study right now.

-- repo-manager@macmini
