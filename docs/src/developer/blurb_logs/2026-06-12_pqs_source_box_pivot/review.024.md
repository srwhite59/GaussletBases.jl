Accepted.

Pass 024 proved that the existing CPBM retained centered electron-nuclear
by-center path is a valid retained boundary operator source for the real
projected-q-shell fixture:

```text
retained_source_nuclear(center) == P' * V_shell_support(center) * P
```

to roundoff for both:

- origin center `(0.0, 0.0, 0.0)`;
- off-origin center `(0.15, -0.2, 0.25)`.

The probe preserved the important nuclear convention:

```text
nuclear charge recorded = true
nuclear charge applied = false
centers summed = false
uncharged by-center convention = true
```

and did not call `_pqs_current_route_safe_term_matrices(...)` or materialize
H1, Hamiltonian assembly, IDA, density-density, RHF, driver, export, or
artifact behavior.

Manager validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- reviewed `tmp/work/pqs_retained_nuclear_boundary_audit_probe_output.txt`

Next implementation target:

Add a by-center-specific final-basis nuclear helper that consumes the existing
retained-source nuclear boundary block, transforms it with:

```text
V_final(center) = L' * V_boundary(center) * L
```

and preserves center metadata with charge still unapplied and centers still
separated. Do not assemble H1 yet.

-- repo-manager@macmini
