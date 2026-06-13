Pass 132 - no-edit Be2 WL/PQS driver comparison readiness audit

Purpose:

Audit what is needed for the medium-term goal: the CR2 agent should be able to
try a Be2 comparison between White-Lindsey and PQS driver payloads from
GaussletBases.

Do not implement anything in this pass.

Why now:

Pass 131 added the first private PQS complete core/shell Ham payload seam, but
it is currently validated on a compact one-center route-smoke. The user wants
the driver to become usable enough that CR2 can compare WL and PQS, perhaps for
Be2. We need to identify the exact readiness gap before coding.

Working target:

```text
GaussletBases driver
-> comparable WL and PQS Be2 payloads
-> explicit Hamiltonian/final-basis conventions
-> CR2 agent can run an external comparison/probe
```

Boundary reminders:

- Source-box-first PQS remains the algorithmic framing.
- Shell/support-row data is diagnostic/oracle support, not route authority.
- Retained diagnostic weights are not IDA/quadrature weights.
- RHF is private validator only, not the product.
- Serious HF belongs to `codexhome/work/hfdmrg`.
- CR2 downstream validation is done by the CR2 agent after the driver line is
  ready.

Read these surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`
  - new `complete_core_shell_ham_payload`
  - route dry-run/stage functions
  - current WL low-order/materialization/preflight paths
- `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`
- `test/nested/pqs_source_box_route_driver_report_runtests.jl`
  - read only; it is broad/slow, do not choose it as a routine validation gate
- `test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `src/fullida_dense_export.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- CR2 context only, no execution:
  - `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/AGENTS.md`
  - `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/answers.md`

Audit questions:

1. Current driver status:
   - What can the driver produce today for one-center PQS?
   - What can it produce today for Be2/diatomic PQS?
   - What can it produce today for WL/low-order Be2?
   - Which paths produce object-carrying payloads versus only report aliases or
     materialization statuses?

2. Comparable payload target:
   - What would be the smallest comparable WL/PQS payload pair for a Be2 CR2
     probe?
   - Should the first comparison be H1-only, H1 plus density-interaction
     surrogate, or a fuller Ham payload?
   - Which conventions must be identical or explicitly labeled: final basis
     ordering, center metadata, one-body terms, Coulomb expansion, density
     gauge, pair-factor convention, IDA/MWG status, and export status?

3. Gaps:
   - Is the new PQS Ham payload one-center-only in practice?
   - Is there a route-owned diatomic/Be2 PQS final-basis/H1/Ham payload seam?
   - Is there an analogous WL Ham payload seam, or only materialization/export
     status?
   - What exact object is missing before CR2 can consume both sides?

4. Smallest next implementation:
   - If the next step is a PQS diatomic/Be2 Ham payload audit/fingerprint, say
     so.
   - If the next step is a WL-side object-carrying Ham payload, say so.
   - If the next step is a driver option/request object to emit private payloads
     for CR2, say so.
   - Do not recommend export/public API unless a private payload cannot answer
     the comparison question.

Trust boundary:

- No file edits.
- No Julia commands required unless a tiny load/read-only check is necessary.
- No implementation.
- No broad route-driver report test.
- No hfdmrg or CR2 execution.
- No artifact/export writing.
- No RHF/SCF work.
- No fixture promotion.

Decision rules:

- If the Be2 comparison is blocked mainly by PQS being one-center-only, identify
  the smallest PQS Be2 payload/fingerprint pass.
- If it is blocked mainly by lack of a WL object-carrying payload, identify the
  smallest WL payload pass.
- If both sides exist but are not exposed in a CR2-friendly way, recommend a
  private handoff/fingerprint pass, not public API.
- Stop and report if the only path requires broad export/artifact plumbing.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Report back:

- Plain-language driver status.
- One-center PQS status.
- Be2/diatomic PQS status.
- Be2 WL status.
- Exact blockers for CR2 Be2 WL/PQS comparison.
- Recommended next 2-3 passes.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
