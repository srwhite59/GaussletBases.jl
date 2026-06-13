Pass 170 - populate WL route in private CR2 inspection generator

Role: repo-doer@macmini

Task type: implementation, generator-only.

Purpose:

CR2 accepted the Be2/PQS inspection artifact generator and pass 169 identified
the smallest honest White-Lindsey population seam. Populate
`routes/white_lindsey` in the private generator using the route-configured
diatomic atom-growth WL path. This should give CR2 a read-only WL/PQS artifact
with both routes labeled honestly.

Allowed implementation surface:

- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- If absolutely necessary, tiny local helper code in the same generator
  directory.

Do not edit tracked `src/`, `test/`, or docs in this pass.

Existing source/test surfaces to mirror:

- `test/nested/cartesian_ham_builder_diatomic_config_smoke_runtests.jl`
  shows the accepted route-configured diatomic WL config shape:
  - `route_family = :white_lindsey_low_order`
  - `atom_symbols = ("Be", "Be")`
  - `nuclear_charges = (4, 4)`
  - `atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))`
  - `parent_axis_counts = (x = 9, y = 7, z = 9)`
  - `q = 5`
  - `n_s = 5`
  - `materializer_backend = :pgdg_localized_experimental`
  - `materializer_nside = 5`
  - `low_order_shellization_policy = :atom_growth_complete_rectangular`
  - `route_configured_diatomic_ham_interaction_treatment = :ggt_nearest`
  - `save_basis_artifact = true`
  - `save_ham_artifact = true`
  - `white_lindsey_expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false)`
- The same test checks the expected route ham bundle shape and labels:
  - `basis/final_dimension == 1271`
  - `ham/model_kind == "ordinary_cartesian_operators"`
  - `ham/interaction_treatment == "ggt_nearest"`
  - `ham/overlap`, `ham/one_body_hamiltonian`, and `ham/interaction_matrix`
    have shape `(1271, 1271)`.

Implementation guidance:

- Add a generator-local WL helper that runs or invokes the existing
  route-configured diatomic atom-growth materialization.
- Prefer an in-process/mktempdir implementation, similar to the existing
  diatomic config smoke, rather than committing intermediate WL basis/ham
  artifacts.
- Temporary basis/ham/report files created to obtain the WL arrays must stay in
  `mktempdir()` or ignored `tmp/work` output space and must not be committed.
- Do not depend on test helper functions. It is okay to duplicate a small config
  writer inside the generator if that is the narrowest route-owned path.

WL schema requirements:

- Change `routes/white_lindsey/status` from unavailable to an available
  route-populated status only when the route materialization and ham bundle are
  actually available.
- Label the route as final-basis ordinary Cartesian/Qiu-White WL data, not PQS
  source-box data.
- Populate compact groups from existing bundle values:
  - `system`: atoms, charges, locations, bond axis if available, route family,
    fixture basics.
  - `final_basis`: final dimension, final integral weights, labels, centers,
    basis kind/materialization kind.
  - `one_body`: overlap and one-body Hamiltonian, plus kinetic/nuclear split if
    already present in the ham bundle.
  - `two_body`: `interaction_matrix` only as a final-basis density-density
    interaction matrix with representation labels such as
    `representation_kind = :final_basis_density_density_matrix` and
    `interaction_treatment = :ggt_nearest`.
  - `validation`: shape, finite, symmetry, weight positivity, nuclear metadata,
    and route/ham adapter statuses.
  - `metadata`: compact route-configured materialization provenance.
- If a field is not available, mark it unavailable with a status/blocker. Do
  not synthesize PQS-style pre-final fields for WL.

Placeholders to add or preserve:

- supplement/residual-GTO: unavailable;
- Qiu-White atom-local HF inputs: unavailable;
- correction/EGOI/stationary/cusp metadata: unavailable;
- MWG/IDA route-configured diatomic ham support: pending/unavailable;
- solver/HFDMRG/RHF/DMRG/HamV6/export readiness: false.

Hard boundaries:

- Do not edit `src/` or `test/`.
- Do not run CR2, HFDMRG, HF, RHF, DMRG, H1/J, HamV6, dense four-index `Vee`,
  `V6`, `Vblocks`, or solver code.
- Do not promote old WL seed/one-center/debug paths into the artifact.
- Do not claim supplement/residual-GTO or Qiu-White correction readiness.
- Do not request UI escalation. In unattended baton mode, if a required command
  needs permission, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

- Because this pass should edit only the tracked generator outside `src/` and
  `test/`, the `src`/`test` line-budget rule should not be triggered.
- If you find you must edit tracked `src/` or `test/`, stop and write
  `.agent_handoffs/ATTENTION.md` unless the final `src` + `test` diff is
  net-negative by `git diff --numstat -- src test`.

Validation:

Run the generator once:

```text
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

This may take longer than 60 seconds because it builds the PQS inspection route
and the WL route-configured diatomic ham bundle. That is acceptable for this
pass because the generated artifact is the deliverable.

Then read back the JLD2/TSV enough to confirm:

- producer commit and dirty flag;
- PQS route still ready;
- WL route status now available, or precise blocker if not;
- WL final dimension;
- WL overlap/H1/two-body shape labels;
- WL two-body representation kind and interaction treatment;
- solver/export flags remain false for the CR2 inspection artifact.

Also run:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git status --short --branch
```

If generated JLD2/TSV or temporary WL basis/ham files appear in `git status`,
fix the generator/ignore behavior before reporting success.

Report:

- files changed;
- whether `src`/`test` line-budget was triggered;
- WL source path actually used;
- generated output paths;
- readback summary for PQS and WL;
- any unavailable placeholders and blockers;
- validation commands/results;
- `git status --short --branch`;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
