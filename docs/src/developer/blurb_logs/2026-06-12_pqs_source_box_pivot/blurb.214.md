Pass 214 - materialize H2 physical density/H1-J input seam, or block precisely

Role:
You are `repo-doer@macmini` implementing one bounded density/H1-J input seam
for GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and
the unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `970190a1 Add H2 physical gausslet H1`
- The H2 physical q5 gausslet-only source plan, final basis, and H1 are
  available:
  - final dimension `463`
  - H1 lowest energy about `-0.7946609179724462`
  - H1 Hamiltonian finite and symmetric
- The current blocker is:
  `:missing_physical_gausslet_h1_j_builder`

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. This is the active
PQS physics target. H2 221 remains source-box diagnostic only.

Purpose:
Build the next private route-owned density-interaction/H1-J input seam for the
physical 463 H2 path. The goal is to materialize the support density inputs and
pre-final density interaction, or return a precise blocker. This is not RHF,
not solver-ready HF, not GTO/MWG supplement, and not public export.

Task:
Add a private physical H2 density/H1-J input payload that consumes:

```text
_PQSDiatomicPhysicalGaussletCoreShellSourcePlan
_PQSDiatomicPhysicalGaussletFinalBasisPayload
_PQSDiatomicPhysicalGaussletH1Payload
```

It should materialize, or block precisely:

```text
density provenance from the parent/source axis bundles
support weights in physical support order
support raw pair numerator matrix
pre-final density interaction matrix
compact H1-J/self-Coulomb diagnostic if the existing convention is unambiguous
```

Use the physical support order:

```text
(:atom_contact_core, :shared_shell_1, :shared_shell_2)
```

and the physical final-basis coefficients/Lowdin cleanup from
`:available_pqs_physical_gausslet_final_basis`.

Allowed math/conventions:
- The density gauge should remain:
  `:pre_final_localized_positive_weight`
- Raw pair-factor convention should remain:
  `:raw_numerator`
- Support weights and raw pair factors should come from
  `CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(...)`
  or the same structured provenance used by the existing diatomic complete
  core/shell density path.
- The pre-final density interaction should follow the existing convention:

```text
pre_final_weights = pre_final_coefficients' * support_weights
weighted_coefficients = pre_final_coefficients ./ pre_final_weights
pre_final_pair_matrix = weighted_coefficients' * raw_pair_numerator * weighted_coefficients
final_to_pre_final_coefficients = combined Lowdin cleanup
```

Important caution:
Do not blindly call
`CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction(...)`
unless you prove its validator accepts the physical final-basis contract. It is
documented around the old direct-core complete-core/shell final basis. It is OK
to reuse the math in a small private physical helper with physical labels.

Expected success transition:

```text
route/h1_j_materialized = true
route/h1_j_status = :materialized_pqs_physical_gausslet_h1_j_payload
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_physical_gausslet_rhf_or_solver_contract
density_interaction/status = materialized
density_interaction/density_gauge = :pre_final_localized_positive_weight
density_interaction/raw_pair_factor_convention = :raw_numerator
```

If you materialize a compact self-Coulomb/H1-J scalar, report it, but do not add
a tight baseline yet. If the density interaction materializes but the scalar
diagnostic convention is not clear, leave H1-J blocked on a precise scalar
diagnostic blocker while preserving the density input payload.

Expected blocker examples:

```text
:missing_physical_gausslet_density_provenance
:missing_physical_gausslet_axis_weights
:missing_physical_gausslet_raw_pair_factor_terms
:physical_gausslet_support_weight_nonpositive
:physical_gausslet_raw_pair_numerator_shape_mismatch
:physical_gausslet_pre_final_density_interaction_blocked
:physical_gausslet_h1_j_scalar_convention_unreviewed
```

Artifact/test expectations:
- Turn on `run_h1_j = true` in
  `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl` only if this pass
  materializes the H1-J/density seam. If you return a precise blocker before
  materialization, keep the request state honest.
- Extend the artifact with compact density/H1-J fields only:
  - density interaction status;
  - density gauge;
  - raw pair factor convention;
  - pre-final pair matrix shape;
  - support weight count;
  - support raw pair shape;
  - H1-J/self-Coulomb scalar if reviewed/materialized.
- Update the H2 physical driver endpoint test to assert those compact fields.
- Do not store or require final-basis self-overlap as downstream working data.

Do not:
- implement RHF/HF/DMRG/HFDMRG;
- add GTO/MWG supplement;
- compare to supplemented WL/QW H2 references;
- use or mutate the H2 221 source-box diagnostic plan;
- add CR2/export/public API work;
- add a final dense Ham/export contract;
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

This test was about 76 seconds after pass 213. If it grows substantially,
report the likely phase. Do not start a broader suite unless a failure requires
a narrowly related check.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.214.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.214.md
```

Report:
- density/H1-J seam materialized or blocked;
- statuses/blocker;
- density gauge and raw pair-factor convention;
- support weight count and positivity;
- support raw pair shape;
- pre-final pair matrix shape/finiteness/symmetry;
- any self-Coulomb/H1-J scalar if materialized;
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
