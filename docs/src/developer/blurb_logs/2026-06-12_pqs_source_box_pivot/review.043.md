Review 043: accepted; the real first blocker is final-basis independence.

Pass 043 usefully prevented a bad implementation. Combining:

```text
125 raw/core source modes
+ 98 Lowdin-cleaned boundary modes from the same 5 x 5 x 5 source space
```

produces a `(223, 223)` overlap matrix with rank `125`, not a `223`-function
final basis. That is the right reason to stop. The issue is not yet
core-shell operator placement; the issue is that the shell sector used in the
attempt was not an independent surrounding shell.

The current exact blocker is:

```text
:missing_route_owned_shell_projection_lowdin_for_combined_basis
```

For the intended physical fixture, the shell must be the surrounding layer
outside the direct core. The old White-Lindsey seed shows the target geometry
only as an oracle/reference:

```text
working box:          (1:7, 1:7, 1:7)
direct core:          (2:6, 2:6, 2:6), count 5^3 = 125
surrounding support:  7^3 - 5^3 = 218 support points
retained shell funcs: 98
total retained dim:   125 + 98 = 223
```

Next pass should target that route-owned complete-shell final-basis
realization. The projected q-shell machinery may already have most of the
needed pieces if it can build a surrounding-shell descriptor over
`current_box=(1:7)^3`, `inner_box=(2:6)^3`, and a `5 x 5 x 5` source-mode
transform. If it cannot, report the exact smaller blocker.

Do not use the old fixed-block packet as active matrix authority, and do not
advance to operator placement, IDA, density-density, RHF, or GTO until the
independent `223`-function final basis is materialized or the exact missing
geometry/Lowdin piece is identified.

-- repo-manager@macmini
