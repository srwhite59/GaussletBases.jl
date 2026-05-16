# High-Order Endcap/Panel Mainline Validation

Date: 2026-05-16

This memo records the current mainline status of the limited high-order
endcap/panel import slice. The imported path is internal, experimental, and
disabled by default.

## Current State

The mainline repository now contains the bounded pieces needed to package the
validated endcap/panel shared-shell idea through ordinary/QW operator
construction:

- owned-unit support and coverage audit scaffolding
- product-doside owned-unit contraction maps
- endcap/panel shell-layer assembly and packet preflight
- guarded internal diatomic source integration
- nested molecular QW operator preflight on a small H2-like case

The public/default diatomic nested source path is unchanged:

```julia
shared_shell_layer_policy = :complete_rectangular
```

The experimental path is available only through the internal source builder:

```julia
_nested_bond_aligned_diatomic_source(
    basis,
    bundles;
    shared_shell_layer_policy = :endcap_panel_owned,
    shared_shell_endcap_panel_q = 4,
    shared_shell_endcap_panel_L = 4,
)
```

No public frontend keyword has been added.

## Validated Operator Preflight

The checked preflight uses a small bond-aligned H2-like basis with a nested
fixed block and an H/cc-pVTZ `lmax = 0` molecular Gaussian supplement with
`max_width = 1.0`.

The default source path remains valid:

- fixed block size: `(539, 347)`
- MWG operator dimension: `349`
- residual count: `2`

The guarded endcap/panel source path also reaches the ordinary/QW operator
layer:

- fixed block size: `(539, 313)`
- operator dimension: `315`
- residual count: `2`
- overlap error: about `5.54e-13`
- one-body symmetry error: `0.0`
- interaction symmetry error: `0.0`
- `interaction_treatment = :ggt_nearest` builds
- `interaction_treatment = :mwg` builds
- nearest residual widths are `NaN`
- MWG residual widths are finite and positive
- residual owners cover centers `1` and `2`

Warm-process timing sanity on the endcap/panel fixed block was approximately:

- `:ggt_nearest`: `0.0366 s`, `94 MB`
- `:mwg`: `0.0424 s`, `99 MB`

These timings are a developer preflight sanity check, not a benchmark contract.

## Boundary

This is not a broad OPCU import, FSB/FBU import, or high-order branch merge.
The imported slice is only the guarded endcap/panel shared-shell path needed to
verify that the source can be packaged and consumed by the existing ordinary/QW
operator machinery.

No H2 HF/ED chemistry, BO reference, CR2 workflow, or DMRG validation is claimed
by this memo.

## Reproduction Artifact

The current preflight can be rerun with:

```sh
julia --project=. tmp/work/h2_endcap_panel_qw_operator_preflight.jl
```

The script prints fixed-block dimensions, operator dimensions, residual counts,
overlap and symmetry errors, width status, owner sets, and warm timing/allocation
for the endcap `:ggt_nearest` and `:mwg` operator builds.

## Next Decision

The next step should be chosen explicitly by repo-manager:

- decide whether and how to expose an internal/public frontend control for the
  guarded source policy; or
- run an old-standard-aligned H2 HF/ED chemistry reproduction using this path.

Until then, the route should remain internal and disabled by default.
