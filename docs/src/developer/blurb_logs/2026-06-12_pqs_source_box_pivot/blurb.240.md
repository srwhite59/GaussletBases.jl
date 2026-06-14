# Pass 240 blurb - independent H2 PQS numerical source-plan audit

Role: repo-doer.

Read before working:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.239.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.239.md`

Task type: no-edit audit.

Purpose:

Pass 239 added the descriptor-only independent H2 PQS source-plan payload and
left the next blocker as:

```text
:missing_independent_pqs_source_plan_numerical_materialization
```

Before implementing numerical materialization, audit exactly what should
materialize next and where. Do not edit files in this pass.

Current descriptor facts:

```text
support counts:  (275, 578, 362)
retained counts: (275, 98, 98)
expected final dimension: 471
source_plan_descriptor_status = :available_independent_pqs_physical_source_plan_descriptor
source_coefficients_materialized = false
fake_pqs/enabled = false
source_backed_fixed_source_oracle_used = false
```

Audit questions:

1. What is the smallest numerical object that should follow the descriptor?
   Options to evaluate:
   - atom-contact direct-source identity representation;
   - shared-shell source transform descriptors;
   - shared-shell shell-realization coefficients;
   - a combined source-plan numerical payload;
   - an explicit blocker that one of these must be designed first.

2. Which existing functions can safely be reused without fake/WL authority?
   Start from:
   - `src/pqs_multilayer_shell_source_plan.jl`
   - `src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl`
   - `src/cartesian_terminal_lowering_contracts.jl` or the current lowering
     module files;
   - `src/cartesian_retained_units/lower_contract_units.jl`
   - `src/cartesian_raw_product_sources/records.jl`

3. For `:atom_contact_core`, decide whether any numerical matrix is needed at
   all before final-basis assembly, or whether a compact identity/source-mode
   descriptor remains sufficient.

4. For each shared shell, identify the actual coefficient-producing seam:
   - raw product source plan;
   - boundary product-mode retained selector;
   - shell projection / Lowdin realization;
   - final retained unit transform.

5. Identify the exact next blocker name if full numerical materialization is
   still too broad, for example:

```text
:missing_independent_pqs_shared_shell_realization_coefficients
```

6. Identify deletion/shrink candidates for the eventual implementation pass.
   The large projected-q-shell integration file is gone; do not rely on more
   broad test deletion. Look instead for stale descriptor/report aliases or old
   blocked-path fields superseded by the next materialized object.

Forbidden:

- no edits;
- no final basis;
- no H1, H1-J, RHF;
- no supplements;
- no CR2/export/public API;
- no fake-PQS/WL coefficient matrices or fixed-source retained transforms.

Report:

- recommended next numerical object;
- exact implementation seam and owner file;
- existing functions to reuse;
- blockers/missing primitives;
- forbidden paths confirmed avoided;
- deletion candidates for the next implementation;
- whether the next pass should implement or stop earlier.

-- repo-manager@macmini
