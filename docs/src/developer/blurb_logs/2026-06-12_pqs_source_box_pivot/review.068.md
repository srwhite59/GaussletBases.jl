Review 068: accepted.

This pass did the intended shrinkage. The H1 fixture now builds its final basis
through:

```julia
plan = pqs_multilayer_shell_source_plan(...)
final_basis = pqs_multilayer_complete_core_shell_final_basis(plan; ...)
```

and removed the manual staged layer/descriptor/shell-plan construction plus
manual final-basis overlap assembly from the test fixture. The test remains a
compact H1 gate and still checks the scientific/workflow facts that matter.

One small remaining duplicate is visible: the test-local
`_pqs_h1_support_product_matrix` now overlaps with the private production helper
`_pqs_multilayer_support_product_matrix` added in pass 067. A small follow-up
can remove that duplicate if it does not make the test less clear.

-- repo-manager@macmini
