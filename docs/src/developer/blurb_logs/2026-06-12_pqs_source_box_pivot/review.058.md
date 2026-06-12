Review 058: accepted; the next blocker is multi-layer PQS shell/source
construction, not more q-only scaling.

The audit aligned PQS with the WL fixture logic. The q ladder used fixed
`R ~= 8`, one surrounding shell, and changing central spacing through q. That is
route/scaling evidence, not a physical convergence ladder.

The attempted WL-aligned side13 fixture:

```text
parent count: 13
mapping:      AsinhMapping(c = 0.1, s = 1.0, tail = 10)
current box:  (1:13)^3
inner box:    (4:10)^3
raw dims:     (7,7,7)
shell layers: 3
```

blocked on:

```text
ArgumentError: projected q-shell requires a one-cell raw boundary around a strict inner box
```

That means the next work should not be another q ladder or gate promotion. It
should identify the minimal route-owned way to represent multi-layer PQS shells
and then rerun final-basis/H1 before any RHF interpretation.

-- repo-manager@macmini
