# PQS/WL Mapping `s_factor`

Status: implemented expert input and provenance contract under
`HP-PQS-MAP-SFACTOR-FN-01` and `HP-PQS-MAP-SFACTOR-TEST-01`.

## Public Convention

`basis.s_factor` is optional, defaults to `1.0`, and must convert to a finite
positive `Float64`. The canonical driver exposes the same scalar as
`s_factor`. Omission and explicit `1.0` select the standard mapping family;
the producer does not tune the value.

For every center with explicit charge `Z` and common `core_spacing = c`:

```text
mapping_s_standard = sqrt(Z * c)
mapping_s_effective = s_factor * mapping_s_standard
```

`core_spacing` remains the physical near-nucleus target spacing. `s_factor`
changes mapping shape independently; it does not replace `core_spacing`,
radius, padding, reference spacing, or tail spacing.

## One-Center Mapping

An origin-centered atom uses the literal mapping:

```text
AsinhMapping(
    c = core_spacing,
    s = s_factor * sqrt(Z * core_spacing),
    tail_spacing = tail_spacing,
)
```

Equivalently, its asinh range is:

```text
a = sqrt(core_spacing / Z) / s_factor
```

The explicit nuclear charge supplies `Z`; the atom label has no mapping
authority.

## Multicenter Mapping

For each Cartesian axis, the implemented center-list path fits one positive
combined inverse-sqrt mapping. For center `i` it uses:

```text
core_range_i = sqrt(core_spacing / Z_i) / s_factor
target_spacing_i = core_spacing
```

Centers sharing the same coordinate on that axis are combined using the
smallest core range and target spacing. The fit solves for positive amplitudes
`A_i` in the existing local-density condition:

```text
1 / target_spacing_j = 1 / tail_spacing
    + sum_i A_i / sqrt((x_j - x_i)^2 + core_range_i^2)
```

Thus the same scalar factor is applied before each axis's combined mapping is
fit; it is not a post-fit coordinate rescaling or a route-specific switch.

## Provenance

Facade-written base and supplemented artifacts preserve mapping truth in their
existing provenance groups:

```text
s_factor
mapping_s_factor
mapping_s_standard
mapping_s_effective
```

For one center, standard and effective values are scalars. For multicenter
systems they are per-center values in input-center order. Parent metadata also
carries the per-atom standard/effective values used by the center-list fit.
Protected ladder recipe readback preserves `s_factor`; legacy recipes without
it resolve to the unchanged `1.0` default.

Provenance is descriptive only. Construction must use the resolved live input,
not recover mapping controls from an artifact sidecar.

## Source Ownership

- `src/mappings.jl` owns the one-center atomic mapping and positive combined
  inverse-sqrt fit.
- `src/pqs_source_box_route_driver_helpers.jl` owns factor normalization,
  multicenter fit inputs, and parent mapping metadata.
- `src/cartesian_base_hamiltonian.jl` owns the public basis input and ordinary
  artifact provenance.
- `bin/cartesian_ham_builder.jl` owns the normal expert driver input.
- `src/cartesian_protected_ladder_bundle.jl` may only preserve and read recipe
  provenance for this control.

No package export or source owner outside those surfaces is introduced by the
`s_factor` contract.

## Failure Behavior

Nonfinite or nonpositive `s_factor` throws `ArgumentError`. One-center mapping
also rejects invalid `Z`, `core_spacing`, or `tail_spacing`. Multicenter fitting
rejects invalid lengths or scales, target spacing not smaller than tail
spacing, nonfinite amplitudes, or a nonpositive-amplitude solution. These
errors propagate; the producer does not silently restore `1.0`, clamp a fit,
or substitute an element-specific value.

## Non-Goals

This contract does not revive public `d`, `parent_mapping_d`,
`parent_mapping_Z`, or route-specific mapping controls. It does not add
element tables, automatic tuning/recommendations, custom mapping objects,
solver workflow, EGOI, rho0/P0, protected-localized interaction changes,
Residual Gaussian/injection policy changes, or Cr2-specific behavior.
