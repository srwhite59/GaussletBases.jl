Pass 208 review: accepted.

The pass did not fake a physical H2 source plan. It added a private physical
gausslet source-plan payload and deliberately reports the current sharper
blocker:

```text
status = :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
blocker = :missing_atom_contact_core_support_rows
```

The payload carries the reviewed H2 463 target inventory:

- support order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- retained order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- support counts `(275, 578, 362)`
- retained counts `(251, 98, 114)`
- expected final dimension `463`
- `retained_atom_core_interiors = true`
- `source_plan_role = :atom_contact_core_plus_pqs_shared_shells`
- `supplement_policy = :none`

The H2 221 diagnostic plan was not reused, and no placeholder support rows,
coefficient matrices, final basis, H1, H1-J, RHF, public API, export, or CR2
behavior was added.

The artifact/report behavior is appropriately narrow:

- `cartesian_assembly` now carries
  `diatomic_physical_gausslet_source_plan_payload`.
- `physical_gausslet_target_summary` exposes compact source-plan status,
  blocker, object kind, and missing objects.
- The H2 physical target artifact now reports
  `:missing_atom_contact_core_support_rows` instead of the broader
  `:missing_physical_gausslet_source_plan`.

The deletion/shrinkage side is acceptable. The pass removed 109 lines of stale
receipt/record diagnostic assertion pressure from the old high-order
source-construction integration test while preserving the basis size, operator
shape, finite-output, symmetry, overlap, and route-dimension checks. It also
removed the non-contract H2 target artifact path print.

Validation:

- Doer ran package load; passed.
- Doer ran the H2 463 artifact-readiness test; 36/36 passed, about 62.8 s.
- Manager reran package load; passed.
- `git diff --check` passed.

Source/test/bin line budget is net negative:

```text
added:   121
deleted: 124
net:      -3
```

Next useful step: audit/construct the missing atom-contact-core support rows and
shared-shell coefficient/contraction source for the physical 463 source plan.
Do not proceed to final basis or H1 until that source-plan construction is
real.

-- repo-manager@macmini
