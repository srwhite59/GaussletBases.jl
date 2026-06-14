# ATTENTION - pass 234 line-budget blocker

Pass 234 is blocked by the explicit line-budget/deletion rule:

```text
Keep src/ + test/ + bin/ net-negative.
If honest net-negative implementation is blocked, write ATTENTION.md and stop.
```

I attempted the bounded implementation for the independent H2 PQS
support-region fingerprint, then reverted my source edits before writing this
file.

Observed implementation shape:

- Compact support-plan helper near
  `_pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload(...)`.
- Used `CartesianShellification.raw_terminal_geometry(...)` from parent axis
  bundle centers and route atom locations.
- Grouped:
  - two `:atom_local_core` regions plus one `:midpoint_slab` as
    `:atom_contact_core`;
  - `:shared_molecular_shell` regions outside-in as
    `:shared_shell_1`, `:shared_shell_2`.
- Generated the expected counts `(275, 578, 362)` in the support order.

Line-budget result from the compact version before revert:

```text
111     0       src/pqs_source_box_diatomic_complete_core_shell.jl
8       0       src/pqs_source_box_route_driver_reporting.jl
0       10      src/pqs_source_box_route_driver_skeletons.jl
```

Net result: `119` added / `10` deleted, net `+109`.

I did not find a safe same-surface deletion large enough to pay for the helper:

- The fake-PQS source-plan candidate/reference path is still wired into
  assembly/reporting and should not be deleted in this pass.
- Support-window helpers are still called by the complete-core-shell diagnostic
  path.
- No accepted endpoint tests or fake-PQS guard fields were touched.

Current git status before this file was written:

```text
## main...origin/main
```

Exact blocker:

```julia
:line_budget_net_negative_blocked_without_unsafe_deletion
```

Recommended manager decision:

- approve a small positive-line exception for this support-region authority
  pass; or
- provide an explicit deletion target large enough to keep the implementation
  net-negative.

-- repo-doer@macmini
