Review 084: accepted as audit; defer implementation.

The audit correctly maps the existing PQS complete core/shell density helpers:

- `pqs_complete_core_shell_final_ida_weights(...)` is final-weight diagnostic
  data, not the current H1/J density-interaction boundary.
- `pqs_complete_core_shell_pre_final_density_interaction(...)` owns the current
  convention: raw support-pair numerator first, positive support-row weights in
  `core_then_shell` order, then localized pre-final positive-weight division.
- `pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)` owns final H1
  orbital consumption by mapping final orbital coefficients into the pre-final
  density gauge.

The docs update is useful because it replaces stale "future H1 payload" text
and records the H1/J convention precisely. However, do not proceed directly to
H1/J implementation. The current PQS multi-layer implementation has too many
seams concentrated in `src/pqs_multilayer_shell_source_plan.jl`; the next
implementation pass should be file/module boundary cleanup before more physics
features are added.

Manager validation:

- `git diff --check` passed.

Deletion/shrinkage:

- no source/test deletion expected in this audit;
- the next shrinkage target is organizational: split the broad PQS multi-layer
  file into clearer seams while preserving behavior.

-- repo-manager@macmini
