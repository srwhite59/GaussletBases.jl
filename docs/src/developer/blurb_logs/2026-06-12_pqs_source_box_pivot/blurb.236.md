# Pass 236 blurb - audit independent H2 PQS retained-rule/source-plan seam

Role: repo-doer.

Read before starting:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.234.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.234.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.235.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.235.md`

Task:

No-edit audit. Identify the smallest honest next seam for independent H2 PQS
retained-rule/source-plan authority now that the support regions are generated
from shellification geometry.

Current generated support units:

```text
:atom_contact_core => 275
:shared_shell_1    => 578
:shared_shell_2    => 362
```

Current blockers:

```julia
:missing_independent_pqs_atom_contact_core_retained_rule
:missing_independent_pqs_shared_shell_2_retained_rule
:missing_independent_pqs_physical_source_plan_materializer
```

Audit questions:

1. What should the independent PQS retained rule be for `:atom_contact_core`?
   Is the old WL/fake retained count `251` physically/algorithmically
   justified by PQS, or should independent PQS retain a different count?
2. What should happen for `:shared_shell_1`? The q=5 boundary product-mode rule
   plausibly explains retained count `98`; name the exact existing function or
   object that should generate it for the H2 support region.
3. What should happen for `:shared_shell_2`? The old/fake count `114` is not
   explained by the standard q=5 boundary count. Is there an independent PQS
   rule that can generate `114`, or should the independent route expose a
   different retained count/blocker?
4. Should the next implementation pass build a partial source plan
   (`shared_shell_1` only), or define the atom-contact-core retained rule first?
5. Which existing source objects should own the source-plan contract:
   `CartesianRawProductSources`, `CartesianTerminalLowering`, PQS multilayer
   source plans, or a new private H2 route source-plan object?
6. What deletion/shrink candidates remain after pass 235, if the next
   implementation pass needs line budget?

Surfaces to inspect:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_multilayer_shell_region_plan.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `src/cartesian_raw_product_sources/`
- `src/cartesian_terminal_lowering/`
- `src/cartesian_shellification/`
- `src/cartesian_retained_units/`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`

Forbidden:

- no source/test/doc edits;
- no final basis, H1, H1-J, RHF, supplements, CR2, export, public API;
- no fake-PQS/WL coefficient matrices or fixed-source retained transforms;
- do not claim retained counts `(251, 98, 114)` are independent PQS results
  unless the audit finds exact PQS authority for each.

Report:

- direct answers to each audit question;
- exact functions/types/files that can be reused;
- exact missing function/type if absent;
- recommended next implementation pass boundary;
- deletion/shrink candidates for that pass;
- `git status --short --branch`;
- confirmation that no files were edited.

-- repo-manager@macmini
