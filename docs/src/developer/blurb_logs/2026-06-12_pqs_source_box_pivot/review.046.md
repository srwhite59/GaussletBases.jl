Review 046: accepted; the bad value is a support-operator convention blocker.

The audit did the right separation. The final-basis transfer is not the source
of the overbinding:

```text
support generalized H1:  -2.0639248059188007
final H1:               -2.0638461028784740
support overlap cond:    1.000000000000083
```

The lowest final vector decomposition is also internally consistent:

```text
kinetic:                 1.1012644410060735
Z=1 nuclear:            -3.1651105438845457
total:                  -2.0638461028784723
```

This rules out the easy mistakes: not a Z=2 label, not double charge
application, not the wrong final sign, not a generalized final solve, and not
the complete-basis Lowdin transfer. The active blocker is:

```text
:blocked_nonacceptance_h1_operator_convention
```

The next useful pass is a same-geometry one-body oracle audit. Avoid comparing
against the default old seed if its mapping differs; instead try to construct
an old/fixed-block or existing one-body route using the exact basis and
geometry from the probe. If that cannot be done, report the exact missing
adapter. Also compare the raw `gaussian_factor_matrices(base_layer)` nuclear
construction against existing PGDG/source-box nuclear factor paths, because the
large negative nuclear expectation is already present on support rows.

Do not add more PQS features until this one-body convention is resolved or
precisely classified as a larger missing kernel.

-- repo-manager@macmini
