Pass 227 - materialize H2 supplement representation request

Role:
You are `repo-doer@macmini` implementing one bounded representation seam for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the
unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should include pass 226.
- H2 R=4 q=n_s=5 gausslet-only PQS/WL physical endpoint is accepted and
  comparison-ready at final dimension 463.
- Pass 226 added `_PQSDiatomicPhysicalGaussletSupplementRequestPayload`.
- For `supplement_policy = :mwg_residual_gto`, the current first blocker is
  `:missing_gto_supplement_representation`.

Architectural guardrail:
WL and PQS should share the same physical support/shell decomposition. For this
H2 target:

```text
support counts  = (275, 578, 362)
retained counts = (251, 98, 114)
retained order  = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
gausslet final dimension = 463
```

The GTO/MWG supplement policy sits above the WL/PQS retained-transform
distinction. Do not treat PQS as needing a different supplement theory.

Task:
Materialize a route-owned, matrix-free GTO supplement representation for the
existing H2 physical supplement request.

Use the existing route-neutral machinery rather than inventing a new
representation:

```julia
legacy_bond_aligned_diatomic_gaussian_supplement(
    "H", "cc-pVTZ", nuclei; lmax = 1
)

basis_representation(supplement)
```

where `nuclei` are the H2 centers `(0,0,-2)` and `(0,0,2)` from the route
request. The resulting object should be a
`CartesianGaussianShellSupplementRepresentation3D` or a clearly equivalent
existing representation object.

Expected facts:

```text
fixture_label = :h2_r4_physical_gausslet_q5
basis_name = "H/cc-pVTZ"
lmax = 1
uncontracted = false
atom symbols = ("H", "H")
nuclear charges = (1, 1)
bond axis/length = :z, 4.0
expected raw supplement orbital count = 18
```

Implementation boundary:
- It is acceptable to add a small private representation payload or to extend
  the request payload if that keeps the object boundary clearer.
- Carry only compact metadata in summaries/artifacts.
- Do not copy orbital exponent/coefficient arrays into the artifact.
- Do not build provider blocks or matrices.

Decision rule:
- For `supplement_policy = :none`, representation status remains
  `:not_requested`.
- For `supplement_policy = :mwg_residual_gto`, create the representation if the
  route request is available.
- If representation construction succeeds:
  - request status should become available, or otherwise clearly report that
    the representation fact is available;
  - representation status should become something like
    `:available_pqs_physical_gausslet_gto_supplement_representation`;
  - preflight blocker should advance to
    `:missing_provider_gto_supplement_blocks`.
- If construction fails because an existing basis/representation contract is
  missing, do not fake metadata. Leave the request blocked with the exact
  blocker and report the missing contract.

Artifact/reporting:
Add or update compact fields only. A small `supplement_representation` group is
acceptable:

```text
status
blocker
object_kind
basis_name
lmax
atom_symbols
center_count
orbital_count
matrices_materialized = false
provider_blocks_materialized = false
```

Keep `supplement_request` compact. Do not add a field cloud to the route report.

Do not:
- build GTO/GTO, mixed gausslet/GTO, MWG residual, or density-density matrices;
- add supplemented H2 scalar values;
- change accepted no-supplement H2 H1/H1-J/RHF values or WL deltas;
- add public API/export/HamV6/CR2/HFDMRG/DMRG behavior;
- pass final-basis self-overlap as downstream working data;
- add a broad test file;
- revive component-smoke/CR2 sidecar vocabulary.

Line-count rule:
This pass must be net-negative under source/test/bin. Measure:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Candidate shrink surfaces:
- Keep shrinking stale assertion pressure in
  `test/nested/pqs_source_box_route_driver_report_runtests.jl` if it only
  preserves transitional route-report field names.
- In the H2 endpoint test, replace repeated request/preflight field checks with
  one compact request/representation/preflight fingerprint.
- Do not delete accepted He/H2 endpoint numerical checks.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 physical endpoint test if it is modified. It is expected to exceed
60 seconds because it rebuilds the accepted H2 endpoint:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

If representation construction uses the existing GTO supplement code, also run
the smallest existing focused representation/source test you identify, unless it
is clearly a slow integration test.

Always run:

```sh
git diff --check
git diff --cached --check
```

Response file:
Write `.agent_handoffs/response.227.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.227.md
```

Report:
- representation payload/status behavior for `:none` and `:mwg_residual_gto`;
- actual representation object kind and orbital count;
- resulting request/preflight blocker transition;
- artifact fields written;
- confirmation no matrices or supplemented values were built;
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
