Pass 130 - no-edit Hamiltonian constructor payload audit

Purpose:

Pivot from RHF feature work to the Hamiltonian constructor boundary. Audit what
the one-center PQS source-box route should produce as a clean route-owned Ham
payload for downstream consumers, without implementing it yet.

Why now:

The private RHF diagnostic has done its useful job as a Hamiltonian stress test:
it forced density/pair-factor conventions to be explicit, exercised the H1/J
data path, and caught the convergence-gating bug. Further RHF capability is a
distraction. Serious HF/DMRG belongs to `codexhome/work/hfdmrg`, and downstream
Cr2 testing should be done by the CR2 agent once GaussletBases exposes a clean
Hamiltonian line.

The product direction should be:

```text
PQS route
-> final basis
-> H1 + electron-electron data
-> clean Ham payload / downstream handoff boundary
```

RHF remains a private validator, not the product.

Read these surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`
  - complete core/shell diagnostic route payload
  - H1/J diagnostic payload and report fields
  - density input payload and pair-factor provenance
- `src/pqs_multilayer_complete_core_shell_h1.jl`
  - final basis and H1/H1J numerical payload helpers
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
  - only as diagnostic/validator context, not as product direction
- `src/pqs_source_box_route_driver_reporting.jl`
  - current report/artifact vocabulary around `cartesian_ham_builder`
- `src/fullida_dense_export.jl`
  - existing dense export conventions and warnings
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- `docs/src/developer/pqs_source_box_fixture_policy.md`
- Context only, no execution:
  - `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hfdmrg/README.md`
  - `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/AGENTS.md`
  - `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/answers.md`

Audit questions:

1. What should the first private PQS Hamiltonian constructor payload contain?
   Consider:
   - final basis object and compact final-basis summary;
   - final H1 matrix/payload;
   - electron-electron representation available today, likely through H1/J
     density interaction inputs rather than a full four-index tensor;
   - Coulomb expansion provenance;
   - source plan / retained transform provenance;
   - center/nuclear metadata needed by downstream users;
   - density gauge and pair-factor convention labels;
   - dimensions and ordering labels;
   - nonclaim flags.

2. What should it explicitly not claim?
   Confirm:
   - no production HF solver;
   - no RHF product surface;
   - no public API yet;
   - no export/artifact writing in the first payload pass;
   - no IDA/MWG semantic promotion;
   - no CR2 science acceptance;
   - no serious-HF claim.

3. What existing objects can be carried rather than flattened?
   Avoid scalar report-field clouds. Prefer a compact payload object carrying
   route-owned source payload, final basis, H1 payload, density inputs, H1/J
   diagnostic payload or density interaction object, and summaries.

4. What is missing for a useful downstream handoff?
   Distinguish:
   - already route-owned objects;
   - objects available only in diagnostic helpers;
   - objects available only as report aliases;
   - truly missing Hamiltonian data.

5. What should the smallest implementation pass be?
   It should probably add a private Ham payload seam, not export/public API. If
   that is wrong, say why.

Trust boundary:

- No file edits.
- No Julia commands required unless a tiny load/read-only check is necessary.
- No RHF route wiring.
- No new SCF/DIIS work.
- No public API.
- No export/artifact implementation.
- No hfdmrg or CR2 execution.
- No fixture promotion.

Decision rules:

- If the payload seam is clear, recommend the smallest implementation pass with
  exact files and one focused validation target.
- If the route currently lacks electron-electron data suitable for a Ham
  payload, report the exact missing object and recommend the smallest data
  contract pass.
- If existing export/report code would tempt a broad public surface, recommend
  staying private and summary-only first.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Report back:

- Source inventory with relevant file/function names.
- Proposed private Ham payload name and fields.
- Required convention labels and nonclaims.
- Missing data, if any.
- Smallest next implementation pass.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
