Review 041: accepted, and continue with the complete core-plus-shell H1 target.

Pass 041 did the right audit-only thing. It did not promote the boundary-shell
H1 fixture, and it identified the next physical basis target cleanly:

```text
5 x 5 x 5 inner core:     125 modes
one PQS boundary shell:    98 modes
combined complete basis:  223 modes
```

The existing production surfaces cannot yet express that combined direct-core
plus boundary-shell final basis. The exact blocker is:

```text
:missing_combined_direct_core_boundary_shell_final_basis_realization
```

The all-source `125`-mode fallback is correctly blocked and should not become
the preferred target. It would require a deliberate all-source retained-rule
contract, and it still would not answer the current user-level intent of a
`5 x 5 x 5` inner core plus one surrounding shell.

Next pass should be a narrow implementation pass for the missing complete H1
route. The acceptable outcomes are:

- materialize the combined `223`-function one-center final-basis H1 path; or
- stop with a smaller exact blocker inside that route, such as missing
  route-owned shell projection/Lowdin data or missing combined core/shell
  operator-block placement.

Do not move to final IDA weights, density-density, RHF, GTO, or driver work
until this complete one-center H1 path exists or a reviewed alternative
complete-basis contract is chosen.

-- repo-manager@macmini
