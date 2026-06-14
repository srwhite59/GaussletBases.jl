Pass 221 - wire H2 no-supplement WL reference values into artifact

Role:
You are `repo-doer@macmini` implementing one bounded comparison-artifact wiring
pass for GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`,
and the unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `23e4f361 Publish PQS blurb 220`
- Pass 220 locally reproduced the matching no-supplement WL/old-QW H2 463 fixed
  block and scalar reference values.
- The active blocker is:
  `:wl_h2_gausslet_only_reference_values_not_yet_wired_to_driver_artifact`.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. H2 221 remains
source-box diagnostic only.

Reviewed WL reference values from pass 220:

```text
wl_h1_lowest = -0.7946609179724673
wl_h1_orbital_self_coulomb = 0.45696639804337047
wl_rhf_one_electron_energy = -1.5611571934181985
wl_rhf_electron_electron_energy = 0.40220533775308426
wl_rhf_electronic_energy = -1.1589518556651142
wl_rhf_nuclear_repulsion = 0.25
wl_rhf_total_with_nuclear_repulsion = -0.9089518556651142
```

Purpose:
Wire these reviewed no-supplement WL values into the existing H2 physical driver
artifact as compact comparison fields and assert small PQS-vs-WL deltas. This
turns the current H2 q5 gausslet-only endpoint into a reviewed WL/PQS comparison
endpoint, still explicitly no-supplement and private/diagnostic.

Task:
1. Add the reviewed WL no-supplement values to the H2 physical driver input or
   a compact local comparison helper, following the existing visible-driver
   style. Prefer the input file if that keeps the route script transparent.
2. Extend the physical H2 artifact/report comparison fields compactly:

```text
comparison/wl_h1_lowest
comparison/delta_h1
comparison/wl_h1_self_coulomb
comparison/delta_h1_j
comparison/wl_rhf_electronic_energy
comparison/delta_rhf_electronic_energy
comparison/wl_rhf_nuclear_repulsion
comparison/pqs_rhf_total_with_nuclear_repulsion
comparison/wl_rhf_total_with_nuclear_repulsion
comparison/delta_rhf_total_with_nuclear_repulsion
```

If you can reuse a smaller existing field set without ambiguity, do so. Do not
use an ambiguous molecular `wl_rhf_total` unless you also label whether it is
electronic-only or nuclear-repulsion-included.

Expected endpoint transition:

```text
comparison/ready = true
comparison/blocker = nothing
physics/endpoint_ready = true
physics/endpoint_blocker = nothing
comparison/reference_label = "WL/QW H2 R=4 gausslet-only 463"
```

Keep the no-matrix candidate fields from pass 219. Keep
`comparison/old_supplemented_wl_qw_scalar_references_blocked = true`.

Expected numerical assertions:

Use tight but not brittle tolerances. Suggested checks:

```text
abs(delta_h1) < 1e-10
abs(delta_h1_j) < 1e-10
abs(delta_rhf_electronic_energy) < 1e-9
abs(delta_rhf_total_with_nuclear_repulsion) < 1e-9
```

The PQS electronic RHF value from pass 217 was about
`-1.1589518556683855`, so the electronic delta should be around
`-3.3e-12`.

Do not:
- recompute WL values in the driver;
- add WL matrices to the artifact;
- add GTO/MWG supplement;
- compare to supplemented WL/QW H2 values;
- alter the H2 221 diagnostic route except possibly shrinking stale assertions
  if needed for the line budget;
- add public API, export, HamV6, CR2, HFDMRG, or DMRG behavior;
- pass final-basis self-overlap as downstream working data;
- hide driver behavior behind an opaque wrapper;
- add a new test file.

Test expectation:
Update the existing H2 physical endpoint test to assert comparison readiness,
the WL values, and the deltas. Do not add a broad/default-runner test.

Line-count rule:
The active source/test/bin line-count rule applies. Final tracked diff under
`src`, `test`, `bin`, and the CR2 generator script must be net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Candidate shrink surface:
- `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`
  still carries many parent-materialization and nonclaim assertions for the
  old 221 diagnostic route. If source/test additions need budget, shrink that
  test toward its live contract:
  - role is `:source_box_diagnostic`;
  - endpoint ready false with `:retained_atom_core_interiors_missing`;
  - final dimension 221;
  - H1 finite/symmetric;
  - no comparison/RHF fields.
- Preserve active He/H2 endpoint tests.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

If you shrink the 221 diagnostic test, also run it:

```sh
julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl
```

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.221.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.221.md
```

Report:
- comparison readiness and endpoint readiness before/after;
- exact comparison fields added;
- WL constants used;
- PQS-vs-WL deltas;
- confirmation supplemented WL/QW scalar references remain blocked;
- source/test/bin scoped added/deleted/net line count;
- validation commands and elapsed times;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini
