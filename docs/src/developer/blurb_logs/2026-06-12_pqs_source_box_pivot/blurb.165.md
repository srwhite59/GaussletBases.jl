Pass 165 - map CR2 inspection artifact schema, no edits

Role: repo-doer@macmini

This is a no-edit audit/schema pass. Do not change source, tests, docs, runtime
behavior, artifacts, or handoff files except your response file.

Goal:

CR2 answered that the first Be2 WL-vs-PQS handoff should be a compact
CR2-specific JLD2 inspection artifact plus a small JSON-or-TSV fingerprint.
Before adding any writer, map the requested fields to what GaussletBases can
produce today for PQS and for White-Lindsey, and define the smallest private
artifact schema.

Read:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`
- `~/Dropbox/codexhome/work/cr2/reports/be2_wl_pqs_downstream_handoff_request_2026-06-13.md`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- existing JLD2/export examples as needed, especially:
  - `src/fullida_dense_export.jl`
  - `src/sliced_ham_export.jl`
  - `src/cartesian_bundle_export.jl`

Also inspect only enough existing White-Lindsey Be2 surfaces to answer whether
the same schema can be filled now for `:white_lindsey`. Use `rg` first. Do not
run a broad test suite.

Questions to answer:

1. What exact top-level artifact schema should the private CR2 inspection JLD2
   use? Prefer plain arrays, strings/symbol labels, numbers, vectors/tuples, and
   NamedTuple-like metadata. Do not store private `_PQS...` structs.
2. For `:pqs_source_box`, which requested fields are already available from
   the current handoff?
3. For `:white_lindsey`, which requested fields are already available in an
   analogous route/build path, and which are missing or blocked?
4. Should the companion fingerprint be TSV or JSON for the first pass? Since
   JLD2 is already a dependency and JSON is not in the main project, TSV is
   probably the conservative default unless you find an existing JSON writer in
   main src.
5. What is the smallest first implementation seam?
6. What old readiness/fingerprint assertions or scaffolding could be deleted or
   shrunk in the implementation pass to satisfy the source/test line-budget
   rule?

Schema should cover, at least as availability/blocker fields:

- route label/family/kind and repo commit/dirty marker;
- system metadata: atoms, charges, Cartesian coordinates, nuclear repulsion,
  electron count, spin convention, bond length/axis, units;
- fixture metadata: q, n_s, spacing, parent axis counts/dimension, physical
  box/extents, backend/fallback/provenance flags;
- final-basis metadata: final dimension/order, retained/source labels, range or
  owner labels when authoritative, explicit unavailable/mixed labels otherwise;
- pre-final/support metadata: pre-final dimension, support row order,
  support/pre-final weights, and retained diagnostic weight nonclaim;
- one-body arrays/metadata: final one-body H, low H1 spectrum, and availability
  for overlap/kinetic/electron-nuclear-by-center;
- two-body inspection data: current pre-final density interaction,
  pre-final pair matrix, final-to-pre-final coefficients, pre-final weights,
  support weights, support raw pair numerator, raw pair-factor provenance;
- ordering contracts for each matrix;
- validation scalars: dimensions, finite/symmetry defects, density gauge, raw
  pair convention, readiness flags, blockers.

Hard boundaries:

- No source/test/docs edits.
- No artifact writing.
- No public API.
- No HamV6, dense final-space `V`, dense four-index `Vee`, `V6`, `Vblocks`,
  solver bundle, HF/RHF/DMRG, H1/J promotion, Qiu-White correction
  implementation, exports, or artifacts.
- Do not call CR2/HFDMRG or run downstream solvers.
- Do not ask the user for permission from inside the unattended loop. If you
  hit a permission/external-state blocker, write `.agent_handoffs/ATTENTION.md`
  and stop.

Line-budget rule for the later implementation pass:

- Any pass editing tracked `src/` or `test/` must be net-line-negative across
  `src` + `test` by `git diff --numstat -- src test`.
- This audit should identify deletion/shrinkage candidates so the next code
  pass can pay for itself.

Report back with:

- proposed JLD2 group/key schema;
- proposed TSV fingerprint columns;
- PQS field availability table;
- White-Lindsey field availability table;
- smallest first implementation seam;
- deletion/shrinkage candidates for the implementation pass;
- commands run and git status.

-- repo-manager@macmini
