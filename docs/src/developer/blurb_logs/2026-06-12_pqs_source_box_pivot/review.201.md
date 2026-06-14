Pass 201 review

Accepted.

This implements exactly the right first surface for the physical H2 target:
route-owned inventory and driver artifact labels, with no source-plan,
final-basis, H1, H1-J, density interaction, RHF, supplement, WL comparison, or
export behavior.

Accepted target inventory:

- Route kind is
  `:bond_aligned_diatomic_physical_gausslet_core_shell_pqs`.
- Parent axis counts are `(x = 9, y = 9, z = 15)`.
- Support and retained order is
  `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`.
- Support counts are `(275, 578, 362)`.
- Retained counts are `(251, 98, 114)`.
- Expected final dimension is `463`.
- `retained_atom_core_interiors = true`.
- `source_plan_role = :atom_contact_core_plus_pqs_shared_shells`.
- `supplement_policy = :none`.

Accepted implementation boundary:

- The inventory is explicitly from the reviewed pass-200 WL/QW gausslet-only
  contract.
- The old supplemented route is not constructed.
- The current 221-dimensional diagnostic route remains separate and unchanged.
- The new physical route writes target inventory only and blocks on
  `:missing_physical_gausslet_source_plan`.
- The focused artifact test confirms final basis, H1, H1-J, and private RHF are
  not materialized.

Validation reviewed and rerun:

- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
  passed, 37/37, 1m05.4s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Deleted-test live search had no `src/test/bin` hits.
- Scoped line budget including the two new untracked test/input files is
  `509 added / 1656 deleted`, net `-1147`.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_lw_plan_batch_runtests.jl`.
- Deleted `test/nested/cartesian_pair_block_one_body_accessors_contract_runtests.jl`.
- The default mixed one-body consumer smoke remains and was updated to remove
  the stale pointer to the deleted accessor scaffold.

Next manager decision:

Pause the H2 implementation lane for the requested cleanup discussion. The next
published blurb should audit/extract/delete the private-global-overlap subtree
from `src/pqs_source_box_route_driver_helpers.jl`, not continue H2 source-plan
implementation.

-- repo-manager@macmini
