Accepted as the source-transform readiness closeout.

Pass 017 reran the 5 x 5 x 5 retained-source H1 probe using the new
repo-owned source-axis transform fact builder:

```text
pqs_source_axis_transform_facts_from_pgdg_axes(...)
```

The probe no longer calls the old nested source-box helper. It reproduces the
pass-015 old-transform oracle H1 and overlap diagnostics to roundoff:

```text
lowest H1 diagnostic = 0.0320561000473788
h1_delta_from_pass015 = 0.0
overlap_identity_error_delta_from_pass015 = 0.0
```

The source-axis transform blocker is cleared for this retained-source probe.
The H1 solve remains diagnostic only, because it is still before PQS shell
realization / Lowdin final-basis construction.

Manager validation:

- `julia --project=. tmp/work/pqs_source_box_repo_transform_h1_readiness_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next blocker:

```text
:missing_pqs_shell_realization_lowdin_final_basis_construction
```

Next target:

Do a read-only audit/target-card pass for the shell-realization/Lowdin
boundary before writing code. The key risk is accidentally turning shell-row
support contraction into the PQS algorithm. The audit should identify the
minimal final-basis object, transform direction, old oracle surfaces, and the
first implementation step.

-- repo-manager@macmini
