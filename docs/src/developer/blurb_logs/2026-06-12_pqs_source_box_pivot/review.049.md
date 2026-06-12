Review 049: accepted; final IDA weights are now a route-owned seam.

The helper:

```text
pqs_complete_core_shell_final_ida_weights
```

projects support-row integral weights into the actual complete final basis:

```text
final_ida_weights = final_coefficients' * support_weights
```

The same-geometry fixed-block oracle comparison is strong after gauge
alignment:

```text
max weight delta: ~6.9e-13
```

The current final-gauge weights are signed:

```text
min: -18.967313490013488
max:  58.136426512094616
positive count: 111
negative count: 112
near-zero count: 0
```

Signed weights are not automatically a bug for an orthogonal final gauge, but
they make the next convention check important. Density-density must use the
correct raw-numerator projection first, then divide at the final retained
weight boundary. Do not derive another weight convention in the density pass.

Next pass should implement or probe only the first final-basis
density-density/self-Coulomb diagnostic, not RHF.

-- repo-manager@macmini
