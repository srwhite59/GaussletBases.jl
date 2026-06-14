Pass 211 review: accepted.

This pass completed the intended source-plan transition for the H2 463 physical
gausslet target.

The new private source-plan object is:

```text
_PQSDiatomicPhysicalGaussletCoreShellSourcePlan
object_kind = :pqs_diatomic_physical_gausslet_core_shell_source_plan
status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
```

It is built only from the checked source-backed candidate and carries the
private adapter authority label:

```text
source_plan_authority_status = :private_source_backed_adapter_authority
```

The artifact now advances the blocker to the next real seam:

```text
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_physical_gausslet_final_basis_builder
```

The pass did not add final basis, H1, H1-J, density interaction, RHF, HFDMRG,
CR2, export, or public behavior. It also kept the H2 221 diagnostic route out
of this path.

The source-plan object carries the physical 463 target data:

- support counts `(275, 578, 362)`
- retained counts `(251, 98, 114)`
- final dimension `463`
- no supplement
- source-backed adapter provenance

Validation:

- Doer ran package load; passed.
- Doer ran the H2 463 artifact-readiness test; 40/40 passed, about 74.0 s.
- Manager reran package load; passed.
- `git diff --check` passed.

Source/test/bin line budget is net negative:

```text
added:   134
deleted: 141
net:      -7
```

Next step: add a physical H2 final-basis payload/builder that consumes this
private source-plan object and produces a 463-dimensional final basis with an
overlap identity diagnostic. Do not add H1 in that pass.

-- repo-manager@macmini
