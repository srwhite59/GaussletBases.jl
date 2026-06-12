Review 063: accepted; stop same-parent fixture probing.

The side-13 parent-only ladder did exactly what it needed to do. With parent
count 13, `AsinhMapping(c = 0.1, s = 1.0, tail_spacing = 10.0)`, and outer box
`(1:13)^3`, varying the direct core side and shell layer count gave:

```text
core side 5   final dim 1429   RHF -2.8366493020529360
core side 7   final dim 1549   RHF -2.8372556463894707
core side 9   final dim 1717   RHF -2.8373268318655693
core side 11  final dim 1933   RHF -2.8373293921422116
```

The results are monotone but effectively plateaued after core side 7. This says
the side-13/core7 split was not a special artifact, but it also says the
remaining `~0.02435 Ha` gap is not fixed by repartitioning the same parent box
between direct core and shell layers.

Do not continue same-parent core/shell ladders. The physical fixture question is
now the parent/mapping rule: radius, central spacing `d`, distortion `s`, `Z`,
direct core size, and shell depth together. The user explicitly noted that a
new rule may need more study than we want right now, so keep that caveat in the
docs but do not launch a broad fixture-rule study in the baton loop.

The next useful implementation step is a contained operator-shape gap: direct
retained PQS safe-term blocks exist for overlap, kinetic, and by-center nuclear,
but position/x2 moment blocks are still missing from the direct retained route.
Add those rather than doing more physical extrapolation.

Deletion/shrinkage: no permanent tests or source changed in pass 063. The
result demotes same-parent repartitioning probes to reference evidence only;
they should not become acceptance coverage.

-- repo-manager@macmini
