Pass 264 - retire legacy source-box fixture builders unblocked by metadata export deletion

Context:
- Current HEAD should include
  `09105671 Retire current-route metadata export stack`.
- Pass 243 preserved three retained-unit fixture builders because
  `current_route_metadata_export.jl` still called them:
  - `_pqs_contact_cap_product_doside_unit`
  - `_pqs_outer_mismatch_product_doside_units`
  - `_pqs_atom_box_support_dense_units`
- Pass 263 deleted `current_route_metadata_export.jl` and proved there are no
  remaining `src/test/bin` callers for that old metadata stack.
- The broader atom-growth/contact-cap/outer-mismatch route concepts are still
  live elsewhere. This pass is only about the old
  `cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl` fixture
  builders that were retained solely for the deleted metadata-export path.

Task:
Do a caller-audited deletion of the now-unblocked legacy source-box fixture
builders in:

```text
src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl
```

Candidate functions:
- `_pqs_contact_cap_product_doside_unit`
- `_pqs_outer_mismatch_product_doside_units`
- `_pqs_atom_box_support_dense_units`
- directly private helpers used only by those functions, such as
  `_pqs_contact_cap_parent_coefficient_matrix`, if the caller audit proves they
  have no other live callers.

Decision rule:
1. Audit callers with `rg` before deletion.
2. If callers are limited to the candidate family itself and historical docs/logs,
   delete the functions and any private helpers used only by them.
3. If any active source or test caller remains, stop and report an exact deletion
   card. Do not add compatibility wrappers.

Strict exclusions:
- Do not touch current atom-growth/contact-cap/outer-mismatch route code outside
  `legacy_source_box_fixtures.jl`.
- Do not edit independent H2 PQS route/source-plan/final-basis/H1/H1-J/RHF/
  support-partition/provider code.
- Do not touch provider blocks, supplement values, CR2/export, HamV6, or public
  API.
- Do not delete unrelated source-box pair-plan or density/nuclear convention
  helpers in the same pass.
- Do not chase historical docs/log references unless a small retirement note is
  clearly helpful.

Validation:
- `git diff --check`.
- Package load:
  `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Caller proof after deletion:
  run `rg` checks for the deleted function/helper names in `src test bin`.

Line budget:
- Scoped `src + test + bin` must be net-negative.

Report:
- Caller audit results.
- Functions deleted and functions preserved.
- Scoped line count for `src + test + bin`.
- Validation commands.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
