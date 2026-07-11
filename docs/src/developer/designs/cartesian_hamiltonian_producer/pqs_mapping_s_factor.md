# PQS/WL Mapping `s_factor`

Status: implemented expert input and provenance facility under
`HP-PQS-MAP-SFACTOR-FN-01` and
`HP-PQS-MAP-SFACTOR-TEST-01`.

This is a narrow expert knob for Cartesian/PQS/WL parent mapping shape. It is
not a new default policy, not automatic tuning, and not a statement that small
`core_spacing` is bad. The purpose is to let expert consumers such as CR2 scan
mapping strength independently from the near-core physical scale when the
standard `s = sqrt(Z * core_spacing)` family is too restrictive.

## Public Convention

Add optional positive `s_factor`, default `1.0`.

```text
standard_s = sqrt(Z * core_spacing)
effective_s = s_factor * standard_s
```

`s_factor = 1.0` must preserve current behavior. Omitted `s_factor` and
explicit `s_factor = 1.0` should be numerically identical up to ordinary
floating-point/reconstruction effects.

For one-center White-Lindsey/atom mapping, the convention is literal:

```text
AsinhMapping(c = core_spacing,
             s = s_factor * sqrt(Z * core_spacing),
             tail_spacing = tail_spacing)
```

The durable public knob is only `s_factor`. Public `d`,
`parent_mapping_d`, `parent_mapping_Z`, and route-specific mapping internals
remain unsupported. `core_spacing` remains the near-nucleus physical scale and
box/radius/padding remain separate concepts.

For multicenter PQS mapping, the intended behavior is the analogous
per-center mapping-strength factor applied before the combined inverse-sqrt
mapping is fit. Doer must record exactly how `s_factor` maps into the
combined-invsqrt construction. If that mapping is ambiguous, implement the
one-center path only and report the exact multicenter design question before
touching CR2 production scripts.

## Provenance

Any producer/artifact/manifest path that records mapping controls should
record enough truth to reproduce the mapping:

- `mapping_s_factor`;
- `mapping_s_standard`;
- `mapping_s_effective`;
- `mapping_c` / `mapping_d` / `core_spacing` as already appropriate for that
  path.

For multicenter paths, provenance should also record whether the values are
per-center, per-axis, or fitted combined mapping values. Do not hide
multicenter ambiguity behind a single scalar provenance field.

## Approved Source Surface

- `src/mappings.jl`;
- `src/pqs_source_box_route_driver_helpers.jl`;
- `src/cartesian_base_hamiltonian.jl`;
- `bin/cartesian_ham_builder.jl` only if needed for the normal expert input
  path;
- `src/cartesian_protected_ladder_bundle.jl` only to preserve/read recipe
  provenance.

No package export, solver workflow, EGOI, rho0/P0, protected-localized
interaction convention, or automatic element/default-table behavior is
approved by this lane.

## Guardrails

- `s_factor` must be finite and positive.
- Default behavior must be unchanged for omitted `s_factor`.
- The driver may expose the knob as an expert scalar, but must not tune it.
- Do not add element-table defaults, charge-derived policy beyond the existing
  explicit `Z`, or automatic recommendations.
- Do not reinterpret this as a replacement for `core_spacing`.
- Do not change `Vee`, IDA/MWG, EGOI, rho0/P0, solver workflow, or protected
  artifact conventions.

## Validation

Approved validation for `HP-PQS-MAP-SFACTOR-TEST-01`:

- `git diff --check`;
- package load;
- default H/H2 or small base artifact/readback path unchanged with omitted
  `s_factor`;
- explicit one-center atom with `s_factor != 1` records provenance and
  changes the mapping;
- small multicenter smoke if the multicenter path supports the knob;
- no Cr2 production run required.

Failure rule: if multicenter combined-invsqrt mapping cannot unambiguously
support the same `s_factor` semantics, stop at one-center implementation and
report the exact missing mapping fact. Do not guess the combined mapping rule
and do not patch CR2 scripts around it.
