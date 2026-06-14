# Pass 240 response - independent H2 PQS numerical source-plan audit

No files were edited for the audit.

## Recommended Next Numerical Object

Implement the shared-shell shell-realization coefficient payload first, not a
full combined source-plan payload.

Recommended private boundary name:

```text
_PQSIndependentH2SharedShellRealizationPayload
```

or a helper-returned internal NamedTuple with the same role.

This payload should materialize, per shared shell:

```text
raw_source_plan::CartesianRawProductSources.RawProductBoxPlan
retained_rule::CartesianRawProductSources.PQSBoundaryProductModeRetainedRule
shell_support_indices
shell_support_states
shell_projection
lowdin_cleanup
shell_final_coefficients = shell_projection * lowdin_cleanup
realized_overlap / identity error summary
```

It should not yet assemble the complete core/shell final basis, H1, H1-J, RHF,
supplements, CR2, exports, or public API.

The atom-contact core does not need a numerical identity matrix before
final-basis assembly. It needs only support indices/states and a descriptor that
the retained source modes are direct identity-like modes. The existing complete
core/shell final-basis helper builds the core identity block internally once it
has the core support rows and overlap block.

## Exact Implementation Seam And Owner File

Primary owner file:

```text
src/pqs_source_box_diatomic_complete_core_shell.jl
```

Suggested seam:

```text
_pqs_source_box_route_driver_independent_h2_shared_shell_realization_payload(...)
```

called from:

```text
_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload(...)
```

after the descriptor is available and before promoting the independent source
plan to numerical availability.

The eventual next combined source-plan object should become available only
after both shared-shell realization payloads are available and after the
atom-contact support rows are recorded. Until then, keep:

```text
source_plan.status = :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
```

## Existing Functions To Reuse

Safe to reuse without fake/WL authority:

- `CartesianRawProductSources.raw_product_box_plan(...)`
  - Use with the route-owned shared-shell source CPB and
    `source_mode_dims = (5, 5, 5)`.
  - This is metadata-only unless axis transform facts are explicitly
    materialized.

- `CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(...)`
  - Already gives the independent q=5 boundary COMX retained selector with
    retained count `98`.

- `CartesianFinalBasisRealization.pqs_source_shell_realization_final_basis(...)`
  - Reusable after the route supplies shell support rows, shell overlap,
    shell projection, and Lowdin cleanup.
  - It is the right shell-realization contract for one shared shell.

- `CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(...)`
  - Reuse later, after the source plan has core support rows and combined
    shared-shell final coefficients.

- `CartesianTerminalLowering` and `CartesianRetainedUnits` metadata paths
  - `:direct_*_identity_cpb` retained units are useful for atom-contact identity
    descriptors.
  - `:pqs_filled_source_cpb` retained units correctly mark the shared shells as
    `:shell_projection_lowdin` with retained rule
    `:pqs_boundary_comx_product_modes`.

Reusable with caution:

- `_nested_projected_q_shell_layer(...)`
  - It currently performs the full support-row shell projection and Lowdin
    cleanup for a box shell.
  - It is source-box/PQS machinery, but it is old diagnostic/private machinery
    and carries timing labels plus broader projected-shell vocabulary. A next
    implementation can reuse the mathematical pieces if the route-owned
    support boxes are passed explicitly and no fake/WL source-backed data is
    imported.

- `CartesianContractedParentMetrics._pqs_shell_realization_plan(...)`
  - It is the current shell projection plus Lowdin adapter.
  - It should remain an internal adapter, not become route authority by itself.

Do not reuse for independent authority:

- `bond_aligned_diatomic_nested_fixed_source(...)`
- the fake-PQS source-backed candidate payload;
- WL/QW fixed-source retained transforms or coefficient matrices;
- the old left/product/right complete-core-shell source-plan producer as the
  independent physical route authority.

## Blockers / Missing Primitives

The current descriptor records counts but does not yet retain enough geometry
to build the realization coefficients. The next implementation needs to carry
or recover:

```text
atom-contact support indices/states
shared_shell_1 support indices/states
shared_shell_2 support indices/states
shared-shell raw outer boxes / source CPBs
shared-shell inner exclusion boxes
parent axis metrics for shell overlap
per-shell shell_projection and Lowdin cleanup
```

The support-region helper currently computes terminal regions locally and then
mostly discards them. The implementation pass should either return compact
support row descriptors from
`_pqs_source_box_route_driver_independent_h2_support_region_plan(...)` or build
a small sibling payload that owns those rows.

Recommended blocker until the shared-shell coefficients exist:

```text
:missing_independent_pqs_shared_shell_realization_coefficients
```

After that payload exists, the next blocker should narrow to:

```text
:missing_independent_pqs_complete_core_shell_source_plan_assembly
```

or, if the same pass also assembles the source plan:

```text
:missing_independent_pqs_complete_core_shell_final_basis
```

## Per-Unit Decision

`:atom_contact_core`:

- No dense identity/source coefficient matrix should be materialized yet.
- Keep direct support rows and identity-like source-mode descriptor.
- Final-basis assembly can use direct core support plus overlap blocks.

`:shared_shell_1` and `:shared_shell_2`:

- Raw product source plan: metadata/source CPB plus q=5 source-mode ordering.
- Boundary retained selector: `pqs_boundary_product_mode_retained_rule`, count
  `98`.
- Actual coefficient-producing seam: shell support projection plus Lowdin
  cleanup.
- Final retained unit transform: the shell realization result, not the raw
  retained selector alone.

## Deletion / Shrink Candidates For Next Implementation

deleted:
- none in this audit.

simplified:
- next implementation can replace the broad
  `:missing_independent_pqs_source_plan_numerical_materialization` blocker with
  `:missing_independent_pqs_shared_shell_realization_coefficients`;
- support/descriptor report aliases that only restate the same blocker can be
  collapsed once the shell-realization payload exposes a compact summary.

quarantined:
- fake-PQS H2 463 source-backed reproduction remains separate;
- old left/product/right complete-core-shell source realization remains old
  route machinery, not independent H2 physical authority.

not deleted because:
- no broad test deletion is appropriate after pass 239;
- fake-PQS endpoint coverage still protects the WL/QW reproduction;
- existing raw product source, final-basis realization, terminal lowering, and
  retained-unit contract tests protect live module contracts.

exact remaining caller/blocker:
- `:missing_independent_pqs_shared_shell_realization_coefficients`

## Whether Next Pass Should Implement

Yes, implement the shared-shell realization coefficient payload next, but keep
it narrow:

1. carry compact support-row/source-box descriptors for the two shared shells;
2. materialize q=5 raw product source plans and boundary retained rules;
3. build shell projection plus Lowdin cleanup for each shared shell;
4. report per-shell retained count `98`, coefficient shape, identity error, and
   fake/source-backed usage flags;
5. keep final basis, H1, H1-J, RHF, supplements, CR2, exports, and public API
   blocked.

## Forbidden Paths Confirmed Avoided

No Julia commands were run. Inspection was read-only only.

The audit did not use or recommend fake-PQS/WL coefficient matrices,
fixed-source retained transforms, final basis, H1, H1-J, RHF, supplements, CR2,
exports, or public API.

## Git Status

```text
## main...origin/main
```

-- repo-doer@macmini
