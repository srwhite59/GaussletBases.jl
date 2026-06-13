Pass 148 - audit diatomic H1 consumer seam

Purpose:

Do a no-edit audit for the next Be2/PQS blocker:

```text
:missing_diatomic_complete_core_shell_h1_consumer
```

The route now has a private diatomic source plan and private final basis. Before
adding H1 materialization, check exactly how the existing one-center H1 path
gets support kinetic and electron-nuclear matrices, and what has to change for
the private diatomic source-plan object.

Read first:

- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/pqs_multilayer_support_one_body.jl`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_source_box_fixture_policy.md`

Audit questions:

1. What exact plan fields/properties do these functions need?
   - `pqs_multilayer_support_kinetic_matrix(plan)`
   - `pqs_multilayer_support_electron_nuclear_by_center_matrices(plan; ...)`
   - `pqs_complete_core_shell_final_one_body_matrix(final_basis, support_matrix)`
   - `pqs_multilayer_complete_core_shell_h1_payload(plan; final_basis, ...)`

2. Which requirements are already satisfied by
   `_PQSDiatomicCompleteCoreShellSourcePlan` and
   `_PQSDiatomicCompleteCoreShellFinalBasisPayload`?

3. Which requirements are missing or ambiguous?
   Pay particular attention to:
   - source-plan object-kind guards
   - support row order labels
   - `bundles` and `metrics` convention
   - nuclear center record / charge / coordinate source
   - electron-nuclear sign convention
   - final-basis transfer shape and ordering
   - whether H1 can reuse support-one-body helpers safely

4. What is the smallest safe implementation seam?
   Options to evaluate:
   - broaden support-one-body guard to accept both the old one-center source
     plan and the new private diatomic source plan, if the required fields and
     support-order convention are genuinely compatible;
   - add private diatomic support-one-body wrappers that call lower helpers
     without changing the old public/internal one-center helpers;
   - add a blocked H1 payload with a sharper missing convention/fact if reuse is
     not yet justified.

5. What should the next H1 payload report if implemented?
   Include proposed compact fields such as:
   - status/blocker
   - source-plan status
   - final-basis status
   - support kinetic status
   - support electron-nuclear status
   - final H1 status
   - final dimension
   - one-body energy only if already naturally available
   - nonclaim flags for H1-J/Ham/RHF/public/export/artifact

Trust boundary:

- No code edits.
- No commits.
- No H1/H1-J/Ham materialization.
- No support-density changes.
- No RHF/SCF/DIIS work.
- No WL payload, public API, exports, artifacts, hfdmrg, or CR2 execution.
- Do not promote shell/support-row contraction, raw product-box probes, or old
  WL adapter paths to route authority.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not ask for interactive approval during unattended baton work. If approval
  would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- Read-only audit only.
- `git status --short --branch`
- No Julia commands are required unless you need a short local introspection
  probe. If you do run one, report why and keep it local/ignored.

Report back:

- Existing H1/support-one-body consumer contract, with file/line references.
- Which requirements the diatomic source plan/final basis already satisfy.
- Which requirements remain missing or ambiguous.
- Recommended smallest implementation pass.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
