Pass 169 - audit White-Lindsey population seam for CR2 inspection artifact

Role: repo-doer@macmini

Task type: read-only audit. Do not edit source, tests, docs, generator scripts,
or artifacts in this pass.

Purpose:

CR2 accepted the private Be2/PQS JLD2+TSV generator path. The remaining CR2
ask is narrower:

- populate the White-Lindsey route under the same schema;
- add richer fixture/backend provenance;
- add supplement/residual-GTO and correction/EGOI metadata placeholders;
- still do no solver/export/HFDMRG/DMRG/RHF work.

Before implementation, audit the existing White-Lindsey surfaces and decide
which one, if any, can honestly populate the Be2 comparison artifact. The main
risk is accidentally labeling a one-center, seed, oracle, or unsupported
supplement path as a Be2-comparable WL route.

Context:

- Current clean baseline is `a7f07000`.
- Generator accepted by CR2:
  `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- Current generated outputs are ignored and should remain ignored:
  `tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2`
  `tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv`
- Current CR2 review:
  `/Users/srw/Dropbox/codexhome/work/cr2/reports/be2_pqs_cr2_inspection_artifact_schema_review_2026-06-13.md`

Read first:

- `AGENTS.md`
- `BlurbStyle.md`
- the CR2 review above
- the current generator script
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/cartesian_bundle_export.jl`
- `src/cartesian_shellization_route.jl`
- `src/cartesian_atom_growth_route_driver_helpers.jl`
- relevant WL tests found by `rg`, especially:
  - `test/nested/pqs_source_box_route_driver_report_runtests.jl`
  - `test/nested/white_lindsey_materialized_seed_runtests.jl`
  - `test/nested/cartesian_ham_builder_diatomic_config_smoke_runtests.jl`

Audit questions:

1. Which current WL path can produce Be2-relevant arrays for the CR2 schema:
   route-configured low-order WL ham bundle, diatomic atom-growth route, old
   materialized seed, ordinary QW operators, or none?
2. Does that path produce a diatomic Be2 object comparable to the PQS route, or
   only one-center/seed/oracle data? Do not fake comparability.
3. Which exact existing fields could populate these schema groups:
   `routes/white_lindsey/system`, `final_basis`, `one_body`, `two_body`,
   `validation`, and `metadata`?
4. What two-body representation should WL report if populated? If it is a
   final-basis density-density interaction matrix, say that. Do not invent PQS
   pre-final density-interaction fields for WL.
5. What status/blocker labels should remain if WL cannot yet be populated
   honestly?
6. What richer fixture/backend provenance can be added without route report
   field clouds?
7. What supplement/residual-GTO, correction/EGOI, and Qiu-White atom-local HF
   placeholders can be added now as explicit unavailable metadata?
8. What is the smallest implementation seam for the next pass?
9. What stale source/test surfaces are deletion or shrinkage candidates that
   could pay for the next implementation under the line-budget rule?

Known caution:

`src/ordinary_qw_raw_blocks.jl` currently says ordinary QW operators do not yet
support true active `l > 0` atomic supplements. Do not claim residual-GTO,
supplement, or Qiu-White correction readiness unless the audited path really
has the needed atom-local data and convention.

Hard boundaries:

- No file edits.
- Do not generate new artifacts.
- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense final-space `V`,
  dense four-index `Vee`, `V6`, `Vblocks`, or solver/export code.
- Do not promote WL oracle/debug paths to production route authority.
- Do not call underscored helpers from CR2; this audit is for GaussletBases'
  repo-owned generator path.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule for the next implementation:

- This audit should not edit `src/` or `test/`.
- For any later pass that edits tracked `src/` or `test/`, final `src` + `test`
  diff must be net-negative by `git diff --numstat -- src test`, unless the
  manager explicitly grants an exception.
- Do not satisfy that by deleting scientific endpoint tests or moving code to
  untracked/tmp files.

Report back:

- recommended WL source path, or precise reason no honest path exists yet;
- availability/blocker table for the schema groups above;
- proposed WL schema labels and representation-kind names;
- provenance/placeholders that should be added next;
- smallest implementation seam;
- deletion/shrinkage candidates for the next implementation;
- commands run;
- `git status --short --branch`;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
