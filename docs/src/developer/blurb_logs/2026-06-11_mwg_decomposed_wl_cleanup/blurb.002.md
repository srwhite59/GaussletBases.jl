Purpose:

Implement the residual MWG density-density kernel for the final-basis
`V_gR` / `V_RR` blocks needed by He + GTO RHF.

Why now:

Pass 001 materialized the six raw combined moment matrices and showed the real
side13 He + GTO residual MWG representation is available. The final-basis
density-density path now blocks at:

```text
:missing_residual_mwg_density_density_kernel_for_final_basis_projection
```

Current state:

- `route_global_residual_gto_mwg_representation(...)` can produce residual MWG
  centers and widths from final-basis moments.
- `route_global_combined_gto_final_basis_density_density_matrix(...)` accepts
  the existing decomposed WL `V_gg` block and the residual MWG representation.
- It still refuses to fake `V_gR` / `V_RR` with raw GTO density-density.
- The old MWG helper shape is known:
  - `_qwrg_mwg_interaction_components(...)` builds gausslet/residual and
    residual/residual components in the parent-product gausslet space.
  - `_qwrg_final_residual_mwg_component_blocks(...)` takes a final gausslet
    interaction block, a parent-product-to-final contraction matrix, axis
    bundles, expansion, residual centers, and residual widths, then returns
    fixed/residual, residual/residual, and final interaction blocks.

Exact task:

Add the narrow decomposed-WL route-global residual MWG density-density assembly
needed to materialize final-basis:

```text
V_final =
[ V_gg   V_gR
  V_Rg   V_RR ]
```

where:

- `V_gg` is the existing decomposed WL gausslet density-density matrix;
- `V_gR` / `V_Rg` are gausslet/residual-MWG density-density blocks;
- `V_RR` is the residual-MWG/residual-MWG density-density block.

Code surfaces to inspect first:

- `src/cartesian_pair_block_materialization/route_global_combined_gto_density_density.jl`
- `src/cartesian_pair_block_materialization/white_lindsey_density_density.jl`
- `src/cartesian_pair_block_materialization/route_global_one_body_adapter.jl`
- `src/ordinary_qw_operator_assembly.jl`
  - `_qwrg_mwg_interaction_components`
  - `_qwrg_final_residual_mwg_component_blocks`
- `test/ordinary/mwg_residual_component_helper_runtests.jl`

Named surfaces are starting points, not authority. If a name is stale or
inconsistent with the live repo, stop and report.

Implementation guidance:

The likely implementation shape is:

1. Obtain or build the parent-product-to-final gausslet contraction matrix for
   the decomposed WL retained gausslet sector.
   - This must come from the shellification-backed decomposed WL retained-unit
     authority, not from old fixed-block matrices.
   - Reuse the current factorized retained-basis bridge or unit coefficient
     maps if that is the cleanest available source.
   - If this contraction matrix is not cleanly available, stop and report the
     exact blocker rather than inventing a substitute.

2. Use the existing parent axis bundle object and Coulomb expansion.

3. Use the materialized residual MWG centers/widths from
   `route_global_residual_gto_mwg_representation(...)`.

4. Build `V_gR` and `V_RR` using the MWG component kernels, then assemble
   `V_final` with the existing `V_gg` block.

5. Keep the result internal/provider-level. Do not run He + GTO RHF yet unless
   the implementation is complete and a tiny diagnostic solve is explicitly
   labeled non-acceptance.

Important convention:

- Do not use raw GTO/GTO density-density as final operator data.
- Do not use raw GTO rows except as residual-construction inputs.
- Do not use old fixed-block matrices as authority.
- Preserve the corrected IDA convention:
  - `V_gg` is already the final retained WL density-density matrix with raw
    numerator projected before retained weight division.
  - Residual MWG blocks must follow the old MWG density-normalized convention
    through the MWG kernels.
- For future PQS, weight division happens after projection/Lowdin.

Trust boundary:

Do not add:

- He + GTO RHF acceptance
- generalized-overlap final solve
- full-parent CPB fallback
- old fixed-block matrix authority
- flat `_qwrg_diatomic_*` fallback
- ordinary Cartesian IDA fallback
- PQS changes
- side17 long-running gate
- broad metadata tests

Validation:

- focused test/probe only;
- verify final `V_final` shape is `(final_dimension, final_dimension)`;
- verify `V_gg` block is identical to the supplied decomposed WL density-density
  matrix;
- verify `V_final` is symmetric and finite;
- verify `V_gR` and `V_RR` are materialized for real side13 He + GTO;
- verify raw GTO density-density flags remain false;
- run `julia --project=. test/nested/cartesian_combined_gto_density_density_readiness_runtests.jl`;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- exact MWG kernel surfaces used;
- how the parent-product-to-final gausslet contraction matrix was obtained;
- `V_final`, `V_gR`, and `V_RR` shapes;
- side13 He + GTO residual density-density status;
- whether the final-basis density-density matrix is ready for a side13 He + GTO
  RHF pass;
- validation run;
- deletion/shrinkage report.
