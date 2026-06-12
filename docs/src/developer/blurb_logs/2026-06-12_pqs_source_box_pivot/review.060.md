Review 060: accepted as the first route-owned multi-layer PQS shell source
seam.

The implementation stayed narrow: `pqs_multilayer_shell_source_plan(...)`
builds repeated legal one-cell shell descriptors, validates support
disjointness and coverage, and collapses the shell layers into one shell sector
for the existing complete core/shell final-basis helper. It did not generalize
the final-basis helper, add RHF, add driver wiring, or make fixed-block matrices
authority.

The side13 final-basis/H1 smoke is the right next milestone:

```text
shell layers:       3
final dimension:    1549
overlap identity:   5.51e-13
Z=1 H1:            -0.494223730383033
Z=2 H1:            -1.975561823201342
```

This is close to the earlier WL side13 H1 scale, so the next test should follow
the WL progression: H1 first, then a pre-RHF self-Coulomb/J diagnostic, then
RHF only if the density convention is sane. Keep the spacing rule provisional;
PQS may not need the same q/ns values as WL, and the `Z`, `d`, `s`, radius, and
shell-depth relation can wait for a more detailed study.

-- repo-manager@macmini
