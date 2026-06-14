Pass 207 review: accepted.

The H2 221 driver path is now named for what it is: a source-box diagnostic,
not a gausslet-only physics endpoint. The ambiguous input/test names were
retired in favor of:

- `test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`

The renamed test now asserts the role and blocker fields before scalar/H1
details:

- `route/artifact_role = :source_box_diagnostic`
- `physics/endpoint_ready = false`
- `physics/endpoint_blocker = :retained_atom_core_interiors_missing`
- `basis/retained_atom_core_interiors = false`
- `basis/source_plan_role = :boundary_source_box_diagnostic`
- `comparison/ready = false`
- `comparison/blocker = :supplemented_reference_not_comparable_to_gausslet_only`

The new endpoint manifest is compact and useful. It separates:

- He 419 physical endpoint checks.
- H2 221 source-box diagnostic, endpoint-blocked.
- H2 463 physical gausslet target inventory, source-plan-blocked.

This is the right taxonomy before any physical H2 source-plan work.

Manager validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Doer ran the renamed H2 221 diagnostic artifact test: 52/52 passed,
  elapsed about 86.5 seconds.
- Doer ran the H2 463 target inventory artifact test: 37/37 passed,
  elapsed about 61.6 seconds.
- I did not rerun the two driver tests because doer already ran them and both
  are above the 60 second threshold.

Line budget:

- Source/test/bin rename-aware budget is net negative by one line.
- This pass removed an old non-contract diagnostic print rather than adding new
  source/test pressure.

No driver/source-plan/H1-J/RHF behavior was added. The next physics blocker is
still the H2 463 physical target source-plan producer:
`:missing_physical_gausslet_source_plan`.

-- repo-manager@macmini
