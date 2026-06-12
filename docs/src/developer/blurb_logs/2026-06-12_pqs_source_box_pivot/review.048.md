Review 048: accepted; the durable H1 gate now protects the live route.

The old permanent `98`-function boundary-shell H1 gate has been replaced by a
compact complete core/surrounding-shell H1 gate:

```text
current_box:            (1:7, 1:7, 1:7)
inner_box:              (2:6, 2:6, 2:6)
core support count:     125
shell support count:    218
shell retained count:    98
final dimension:        223
nuclear factor source:  pgdg_intermediate.gaussian_factor_terms
H1:                    -0.48047934800387226
fixed-block oracle:    -0.48047920531279725
oracle delta:           1.43e-7
```

The test file shrank substantially and removed the synthetic retained-unit /
pair-materialization scaffolding that existed only to support the old
boundary-shell-only gate. This is the right test-diet pattern: replace the
transitional surface with the live physical route instead of adding a second
test file.

Next step should be final IDA weights for the complete final basis. Do not go
straight to density-density or RHF. First materialize and validate the actual
integral weights of the `223` final functions:

```text
support integral weights
-> final_coefficients' * support_weights
-> final IDA weights
```

Compare with the same-geometry fixed-block weights if possible, but keep the
old fixed block oracle-only.

-- repo-manager@macmini
