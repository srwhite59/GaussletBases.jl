Pass 187 review: accepted as the driver-owned He 419 PQS seam audit.

The audit answered the right question and did not edit tracked implementation
surfaces. The important conclusion is that the current driver already has the
right visible laboratory-script spine, but cannot yet express the exact
White-Lindsey parent mapping used by the accepted He q=5/n_s=5 PQS gate.

Key accepted findings:

- Preserve `bin/cartesian_ham_builder.jl` as visible staged construction:
  system, recipe, parent, shells, units, transforms, pairs, assembly, report,
  materialization, save. Do not replace this with an opaque `run_driver(config)`.
- The next driver input should use current driver naming where possible:
  `parent_axis_counts = (x = 11, y = 11, z = 11)`,
  `parent_axis_probe_family = :G10` / clearer alias `parent_axis_family = :G10`,
  `q = 5`, `n_s = 5`, `tail_spacing = 10.0`.
- The missing source seam is the one-center parent mapping rule:
  `_cartesian_one_center_parent_basis_object` currently hard-codes
  `IdentityMapping()`. The driver needs explicit parent mapping inputs for
  `:white_lindsey_atomic_mapping` with `Z = 2.0`, `d = 0.3`, and
  `tail_spacing = 10.0`.
- H1 and H1-J should remain owned by `cartesian_assembly`, using the existing
  complete-core/shell diagnostic route payload. Materialization/save should
  publish compact artifact facts, not reconstruct route internals.
- The replacement test should be driver-shaped:
  input file -> `bin/cartesian_ham_builder.jl` -> output artifact -> read and
  check artifact values.
- Once the driver-owned artifact test covers the same physics gate, the direct
  hand-built test `test/nested/pqs_direct_retained_final_h1_runtests.jl` should
  be deleted rather than preserved as scaffold pressure.

One caution for the next blurb: the audit listed
`src/pqs_source_box_route_driver_reporting.jl`; this file exists in the live
repo, but the doer should still use live hooks rather than invent a parallel
reporting surface if the exact artifact hook belongs in the helper/save path.

Recommended next pass:

- Implement the He 419 PQS driver-owned artifact gate.
- Add the minimal parent mapping inputs.
- Add `test/driver_inputs/he_pqs_q5_wlmap.jl`.
- Add a thin driver-output-reading test.
- Delete the direct hand-built He PQS test once the new test passes.
- Keep source/test/bin/generator line budget net-negative.

Validation reviewed:

```text
git status --short --branch
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

No implementation changes were made in pass 187.

-- repo-manager@macmini
