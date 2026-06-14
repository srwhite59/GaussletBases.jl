Pass 212 - add the H2 physical gausslet final-basis seam, or block precisely

Role:
You are `repo-doer@macmini` implementing one bounded final-basis seam for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the
unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `67c880e9 Add H2 physical source plan object`
- The H2 physical source plan is now private route-owned and available:
  - object kind `:pqs_diatomic_physical_gausslet_core_shell_source_plan`
  - status `:available_pqs_diatomic_physical_gausslet_core_shell_source_plan`
  - authority `:private_source_backed_adapter_authority`
  - support counts `(275, 578, 362)`
  - retained counts `(251, 98, 114)`
  - final dimension `463`
- The current blocker is:
  `:missing_physical_gausslet_final_basis_builder`

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. This is the current
PQS physics target being developed.

Architectural guardrail:
Treat WL and PQS as sharing the same physical shell/support decomposition
wherever possible. For this target the common decomposition is:

```text
atom_contact_core support: 275   retained: 251
shared_shell_1 support:    578   retained: 98
shared_shell_2 support:    362   retained: 114
total retained dimension:  463
```

The route vocabulary should describe this shared physical shell/support plan.
PQS should differ in retained transform and operator-construction path, not by
inventing a different physical basis vocabulary. Do not encode this as a
PQS-only physical object if a shared shell/support label is available.

Future supplement guardrail:
Do not implement GTO/MWG supplement in this pass. The intended later model is:
build full/intermediate gausslet and mixed GTO operator matrices, then apply the
chosen final-basis transform. PQS does not need a separate supplement theory,
but the physical 463 gausslet-only final basis must exist first.

Purpose:
Advance the H2 physical route from an available source-plan object to either:

```text
:available_pqs_physical_gausslet_final_basis
```

or a precise final-basis blocker that explains exactly what convention/data is
missing.

Important implementation caution:
Do not blindly reuse
`CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(...)` unless
you prove its assumptions match this physical source plan.

That helper appears to assume a direct identity core:

```text
core final columns == core support rows
```

but this H2 physical route has:

```text
atom_contact_core support rows = 275
atom_contact_core retained cols = 251
```

So this pass likely needs either:

1. a small private physical final-basis builder for contracted core plus
   contracted shared shells, or
2. a precise blocked payload if the source-plan coefficient matrices are not
   yet shaped/owned correctly for that builder.

Exact task:
Add a private H2 physical final-basis payload/builder path that consumes
`_PQSDiatomicPhysicalGaussletCoreShellSourcePlan`.

It should:
- verify the core and shared-shell coefficient matrix shapes against the support
  and retained counts;
- assemble the final retained transform in support order
  `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
- use the source plan's route-owned support rows/states/metrics rather than the
  H2 221 source-box diagnostic route;
- produce an orthonormal final-basis object or compact payload with:
  - status;
  - blocker;
  - final dimension `463`;
  - support counts `(275, 578, 362)`;
  - retained counts `(251, 98, 114)`;
  - final-overlap identity error as a diagnostic scalar;
  - transform/source-plan provenance;
  - H1/H1-J/density/RHF/export/public flags false.

If the final basis materializes, expected transition:

```text
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_physical_gausslet_h1_builder
route/final_basis_status = :available_pqs_physical_gausslet_final_basis
basis/final_dimension = 463
```

If it cannot materialize cleanly, keep the source plan available and return a
more specific blocker than the current generic one. Examples:

```text
:physical_gausslet_core_coefficient_shape_mismatch
:physical_gausslet_shared_shell_coefficient_shape_mismatch
:missing_physical_gausslet_support_overlap
:physical_gausslet_final_overlap_not_identity
```

Final-basis overlap policy:
For the final working basis, do not export or preserve a working `S` matrix for
downstream generalized-overlap use. Self-overlap is diagnostic only. The output
should record an identity-error diagnostic, not route consumers toward a
generalized-overlap solver.

Do not:
- implement H1;
- implement H1-J or density interaction;
- implement RHF/HF/DMRG/HFDMRG;
- add GTO/MWG supplement;
- compare to supplemented WL/QW H2 references;
- use or mutate the H2 221 source-box diagnostic plan;
- add CR2/export/public API work;
- add broad tests or default-runner includes.

Line-count rule:
The active source/test/bin line-count rule still applies. For edits under
`src`, `test`, `bin`, and the CR2 generator script, final tracked diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs and blurb logs do not count. If source lines are added, pay for them by
deleting stale source/test pressure in the same pass. Prefer deleting
transitional metadata/report-vocabulary checks over touching endpoint or
physics checks.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 463 artifact-readiness test with Julia-level timing:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

If the test runtime remains around 75 seconds, report that. If it grows
substantially, report the likely phase. Do not start a broader suite unless a
failure requires a narrowly related check.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.212.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.212.md
```

Report:
- whether the final-basis seam materialized or blocked;
- final-basis status/blocker;
- support/retained counts and final dimension;
- final-overlap identity-error diagnostic if materialized;
- artifact fields changed;
- source/test/bin scoped line budget added/deleted/net;
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
