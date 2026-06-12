Accepted.

The source change is small and on-target. The one-center `:pqs_source_box`
driver path now selects terminal shellification instead of treating
shellification as not applicable, and one-center terminal shellification uses
`CartesianShellification.OneCenterShellification`. This advances the first
missing owner from pass 092 without adding lowering/source-plan/final-basis/H1
construction.

Current status from the focused smoke:

- `low_order_shellization_source = :terminal_cartesian_shellification_geometry`;
- `low_order_terminal_shellification_selected = true`;
- terminal lowering plan is available;
- H1/J remains blocked on the later complete core/shell inputs.

Validation accepted:

- doer's focused dry-run smoke;
- manager rerun of a compact one-center dry-run smoke;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Validation caveat:

The manager rerun of the compact dry-run smoke passed, but took roughly
90 seconds before producing output. That is too slow for routine per-pass
validation. Future passes should keep using the smallest possible smoke and
should avoid broad route-driver/report/materialization files unless the edit
genuinely requires them.

Deletion/shrinkage review:

- No old code became obsolete yet because this pass only enables shell-stage
  selection.
- No tests were added, which is correct.
- The H1 gate was not expanded, which preserves the current test diet.
- The next stale surface to shrink remains fixture-local H1/J reconstruction
  after driver-owned region/source/final-basis inputs exist.

Next pass should target the next missing owner: lower/translate the selected
terminal shellification into `pqs_multilayer_shell_region_plan(...)` and
`pqs_multilayer_shell_source_plan(...)`, without entering final-basis or H1/J
unless that source-plan seam is already complete.

-- repo-manager@macmini
