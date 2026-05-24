# Cartesian Nested Endcap/Panel Shared-Shell Route

This page records the guarded opt-in endcap/panel shared-shell route for the
bond-aligned diatomic nested source builder. It is an experimental mainline
route, not the default diatomic nested policy and not a broad OPCU import.

## Pseudocode

1. Start from the ordinary bond-aligned diatomic nested source geometry.
   The default shared-shell policy remains:

   ```julia
   shared_shell_layer_policy = :complete_rectangular
   ```

   The endcap/panel route is selected only when the caller explicitly passes:

   ```julia
   shared_shell_layer_policy = :endcap_panel_owned
   ```

2. Identify the one-layer shared shell around the current inner rectangular
   box. The validated route uses a thin endcap-box/perimeter shell, not a
   generic rectangular annulus.

3. Partition that shell into owned units:
   - two bond-axis-normal endcaps
   - four side panels over the inner bond-axis span

   Each support point in the shell must be owned by exactly one unit. The
   coverage audit must report zero duplicate, missing, and outside support
   points.

4. Build product-doside contraction maps for every owned unit.
   The current public frontend controls are:

   ```julia
   shared_shell_endcap_panel_q = 4
   shared_shell_endcap_panel_L = 4
   ```

   The retained counts are:
   - `q * q` per endcap
   - `q * L` per panel

   `:direct_selector` remains a debug/scaffolding contract only and is not the
   intended integration route.

5. Assemble the owned units into a shared-shell layer and run the existing
   nested shell packet machinery. The shell layer must have finite packet
   weights and agree with the support-reference packet oracle at roundoff on
   the bounded validation cases.

6. Convert the nested source to a fixed block. The fixed block may no longer be
   exactly representable as one fully factorized final Cartesian basis, so it
   carries staged owned-unit metadata as an auxiliary sidecar.

7. When by-center nuclear sidecars are requested, choose the fastest valid
   contraction path:
   - `:factorized_final` for exactly factorized final blocks
   - `:product_staged_factorized` for product-owned staged metadata
   - `:staged_factorized` for generic staged metadata
   - `:general_parent_dense` as the dense fallback/oracle

   If exact factorized final-basis extraction succeeds, `:factorized_final` is
   the correct fastest path. The product-staged sidecar is the intended
   nonfactorized fallback for validated product-owned layers.

8. Pass the fixed block into the ordinary Qiu-White operator builder. The
   molecular residual-Gaussian interaction treatment should normally use
   `:mwg`; `:ggt_nearest` remains a fallback/debug route and is still useful
   for reproducing older chemistry rows.

The source frontend accepts `gausslet_backend = :pgdg_localized_experimental`
on the opt-in `:endcap_panel_owned` path. The older/default
`:complete_rectangular` nested source frontend remains numerical-reference-only
until it has a separate analytic-backend validation pass.

9. Treat the route as opt-in until consumer-side chemistry checks justify a
   broader default. It is ready for CR2 and related consumers to compare, but
   not to assume as a settled production policy.

## Code Pointers

- Public frontend keyword plumbing:
  - `src/ordinary_qw_nested_frontends.jl:bond_aligned_diatomic_nested_fixed_source`
  - `src/ordinary_qw_nested_frontends.jl:bond_aligned_diatomic_nested_fixed_block`
- Internal diatomic source policy:
  - `src/cartesian_nested_diatomic.jl:_nested_bond_aligned_diatomic_source`
- Owned-unit construction:
  - `src/cartesian_nested_owned_units.jl:_nested_endcap_panel_owned_units`
  - `src/cartesian_nested_owned_units.jl:_nested_endcap_panel_shell_layer`
- Staged by-center sidecar selection:
  - `src/cartesian_nested_faces.jl:_nested_by_center_sidecar_path`
  - `src/cartesian_nested_faces.jl:_nested_product_staged_by_center_sidecar_cache`
- Product-staged nuclear contraction:
  - `src/ordinary_qw_raw_blocks.jl:_qwrg_bond_aligned_staged_by_center_nuclear_one_body_by_center`

## Validation Boundary

The imported mainline slice has been validated as a narrow internal and
experimental path:

- H2-like source construction produces an endcap/panel fixed block smaller
  than the default complete-rectangular path while preserving clean overlap
  diagnostics.
- Ordinary Qiu-White operator construction succeeds for both `:ggt_nearest`
  and `:mwg` residual interactions on the bounded H2 fixture.
- The old-standard H2 `R = 4.0` HF/ED comparison was reproduced through the
  guarded mainline route.
- Product-staged by-center nuclear sidecars agree with dense fallback to
  roundoff on the bounded staged-sidecar probe and are faster than the dense
  fallback there.

The validation does not yet prove:

- broad Cr2-scale readiness
- arbitrary-molecule placement
- that the endcap/panel route should become the default
- that all future OPCU/high-order construction ideas have been imported

## Diagnostics To Check

For consumer-side trials, record at least:

- fixed-block dimensions for the default and endcap/panel routes
- final operator dimension and residual count
- overlap error
- H/V symmetry errors
- residual owner coverage
- MWG width finiteness
- by-center nuclear sidecar path from `_nested_by_center_sidecar_path`, including
  whether `:factorized_final`, `:product_staged_factorized`, or a dense fallback
  was used
- wall time and allocation for the operator build

If the by-center path reports `:general_parent_dense` on a case that was
expected to use either the factorized-final or staged-product route, that is a
performance-contract issue to investigate before treating the calculation as
representative. Reports should record the actual path instead of inferring it
from the shared-shell policy.

## References

- Default diatomic box policy:
  [Cartesian nested diatomic box policy](cartesian_nested_diatomic_box_policy.md)
- Primitive shell packet construction:
  [Cartesian nested face construction](cartesian_nested_face_construction.md)
- Developer validation memo:
  [High-order endcap/panel mainline validation](../developer/high_order_endcap_panel_mainline_validation_2026-05-16.md)
- Developer chemistry reproduction:
  [High-order endcap/panel H2 chemistry reproduction](../developer/high_order_endcap_panel_h2_chemistry_reproduction_2026-05-16.md)
