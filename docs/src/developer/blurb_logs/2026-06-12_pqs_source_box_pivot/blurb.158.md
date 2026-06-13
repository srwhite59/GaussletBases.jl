Pass 158 - audit Be2/PQS Hamiltonian consumer contract for downstream use

Role: repo-doer@macmini

Task type: no-edit audit. Do not edit files, do not commit.

Purpose:

Define the smallest honest Hamiltonian consumer contract needed to move the
Be2/PQS route from private inspect-only handoff toward downstream CR2/HFDMRG
inspection, without implementing export, dense `Vee`, HamV6, or production
workflow behavior yet.

Why now:

After pass 157, the active Be2/PQS route spine is compactly tested and blocked
only on:

```text
:missing_diatomic_hamiltonian_consumer_contract
```

The repo now has a private handoff carrying H1, density interaction, pair
representation/provenance, nuclear metadata, electron count, and spin sector.
Before adding a consumer object, inspect what downstream consumers actually
need.

Read first:

- `AGENTS.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_source_box_fixture_policy.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.157.md`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`

Downstream read-only inspection targets:

- `~/Dropbox/codexhome/work/hfdmrg`
- `~/Dropbox/codexhome/work/cr2`

Inspect only enough to answer the contract question. Prefer `rg` and small file
reads. Do not edit anything outside this repo. If either downstream directory is
not readable in this environment, report that exact blocker; do not request
interactive approval.

Audit questions:

1. What does the current Be2/PQS private Hamiltonian handoff already provide?
   Summarize the concrete available pieces:
   - one-body Hamiltonian/final dimension;
   - density interaction / pre-final pair matrix;
   - final-to-pre-final coefficients;
   - support weights and raw pair numerator;
   - raw pair-factor convention and density gauge;
   - nuclear charges/coordinates/repulsion;
   - electron count and spin sector;
   - explicit nonclaims.

2. What does HFDMRG or CR2 expect for a Hamiltonian input today?
   Look for scripts, loaders, HamV6/JLD2 conventions, one- and two-body matrix
   names, nuclear metadata, electron-count/spin conventions, and comparison
   workflows.

3. Is the current Be2/PQS handoff enough for a read-only inspector/comparison
   probe, or does it need a new compact consumer object first?

4. If a consumer object is needed, what should be the smallest private
   contract?
   Candidate shape might be a private route-owned object that says:
   - `status`
   - `blocker`
   - `source_handoff_status`
   - final dimension
   - one-body Hamiltonian reference/status
   - two-body representation kind/status
   - density/pair convention labels
   - nuclear/electron/spin metadata
   - downstream readiness flags
   - missing export/format blocker

5. What should remain blocked?
   Be explicit about whether blockers should be:
   - `:missing_hamiltonian_consumer_contract`
   - `:missing_dense_two_body_integrals`
   - `:missing_hamv6_export_contract`
   - `:missing_cr2_handoff_format`
   - or something else.

Trust boundary:

- No source edits.
- No tests.
- No export/HamV6/JLD2 writing.
- No CR2/HFDMRG runs.
- No dense `Vee`.
- No H1/J materialization.
- No RHF/SCF/DIIS.
- No WL comparison run.
- No public API or artifact behavior.

Decision rules:

- If downstream expects dense four-index integrals and the current handoff only
  provides a factored/density-interaction representation, report that mismatch
  directly. Do not propose smuggling dense data into summaries.
- If the current handoff is already enough for a private read-only inspector,
  propose the smallest inspector/consumer pass and the exact nonclaims.
- If the next implementation would require choosing an export format, stop at
  the format-decision boundary rather than implementing.

Validation:

Read-only audit only. Run:

```text
git status --short --branch
```

Do not run Julia tests or downstream solver scripts in this audit.

Report back:

- current Be2/PQS handoff inventory;
- downstream HFDMRG/CR2 input expectations, with file citations;
- gap analysis;
- recommended next implementation pass, or blocker if a design decision is
  needed first;
- exact proposed status/blocker vocabulary;
- git status;
- deletion/shrinkage report:
  - deleted: none in this audit
  - simplified: none in this audit
  - quarantined: none in this audit
  - not deleted because: audit only
  - exact remaining caller/blocker:

-- repo-manager@macmini
