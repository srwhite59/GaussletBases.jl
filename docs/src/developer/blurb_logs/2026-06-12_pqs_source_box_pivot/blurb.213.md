Pass 213 - materialize H2 physical gausslet H1, or block precisely

Role:
You are `repo-doer@macmini` implementing one bounded H1 seam for GaussletBases.
Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the unattended baton
rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `ad9e1836 Add H2 physical gausslet final basis`
- The H2 physical q5 gausslet-only source plan and final basis are available:
  - source-plan status:
    `:available_pqs_diatomic_physical_gausslet_core_shell_source_plan`
  - final-basis status:
    `:available_pqs_physical_gausslet_final_basis`
  - support counts `(275, 578, 362)`
  - retained counts `(251, 98, 114)`
  - final dimension `463`
  - final overlap identity error about `1.6e-13`
- The current blocker is:
  `:missing_physical_gausslet_h1_builder`

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. This is the active
PQS physics target. The H2 221 route remains source-box diagnostic only.

Architecture guardrail:
WL and PQS should share the physical shell/support decomposition. For this H2
target, keep the route vocabulary tied to:

```text
atom_contact_core + shared_shell_1 + shared_shell_2
```

PQS-specific work belongs in the retained transform/operator path, not in a
different physical basis vocabulary.

Task:
Add a private H2 physical one-electron/H1 payload path that consumes:

```text
_PQSDiatomicPhysicalGaussletCoreShellSourcePlan
_PQSDiatomicPhysicalGaussletFinalBasisPayload
```

The payload should materialize, or block precisely, the following:

```text
support kinetic matrix
support electron-nuclear by-center matrices
final kinetic matrix
final electron-nuclear by-center matrices
final H1 Hamiltonian
ordinary final-basis H1 solve
```

Use the physical support order:

```text
(:atom_contact_core, :shared_shell_1, :shared_shell_2)
```

and the physical final-basis coefficients from
`:available_pqs_physical_gausslet_final_basis`.

Reuse allowed:
- You may reuse the low-level support product matrix helpers and centered
  Gaussian factor-term helpers already used by the diatomic complete-core/shell
  H1 path.
- You may reuse the lower linear algebra pattern of
  `final_operator = C' * O_support * C`.
- You may reuse `pqs_complete_core_shell_final_one_electron_hamiltonian(...)`
  and `pqs_complete_core_shell_final_h1_solve(...)` only if their input object
  contracts genuinely match or if you provide a small private physical wrapper
  that preserves the correct physical object/status labels.

Do not force old labels:
- Do not require final-basis status
  `:available_pqs_complete_core_shell_final_basis`.
- Do not require source-plan object kind
  `:pqs_diatomic_complete_core_shell_source_plan`.
- Do not relabel the physical 463 path as the old diagnostic complete-core/shell
  path.

Expected success transition:

```text
route/h1_status = :materialized_pqs_physical_gausslet_h1_solve
route/h1_materialized = true
physics/h1_lowest = finite Float64
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_physical_gausslet_h1_j_builder
```

H1 sanity rule:
For this physical H2 target, the H1 lowest eigenvalue should be finite and
negative. If the H1 solve materializes but the lowest eigenvalue is nonnegative,
do not call the H1 route accepted; return a precise blocker such as:

```text
:physical_gausslet_h1_lowest_energy_nonnegative
```

and report the scalar.

If H1 cannot materialize, keep the final basis available and return a precise
blocker. Examples:

```text
:missing_physical_gausslet_h1_inputs
:physical_gausslet_support_kinetic_shape_mismatch
:physical_gausslet_support_electron_nuclear_shape_mismatch
:physical_gausslet_final_one_body_shape_mismatch
:physical_gausslet_h1_hamiltonian_not_symmetric
```

Artifact/test expectations:
- Turn on `run_h1 = true` in
  `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`.
- Extend the driver artifact writer only enough to expose compact H1 facts:
  - route H1 status/materialized flag;
  - final dimension;
  - `physics/h1_lowest`;
  - H1 symmetry/finiteness diagnostics if already compactly available.
- Update
  `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
  to assert the H1 route facts and the negative finite H1 sanity check.
- Do not add a tight H1 reference baseline yet. Report the observed H1 lowest
  value for review.

Do not:
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
deleting stale source/test pressure in the same pass. Do not delete the active
He or H2 driver endpoint tests.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 463 driver-artifact test with Julia-level timing:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

This test was about 74 seconds after pass 212. If it grows substantially,
report the likely phase. Do not start a broader suite unless a failure requires
a narrowly related check.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.213.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.213.md
```

Report:
- H1 materialized or blocked;
- H1 status/blocker;
- final dimension;
- H1 lowest energy and symmetry/finiteness diagnostics;
- support/final one-body statuses;
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
