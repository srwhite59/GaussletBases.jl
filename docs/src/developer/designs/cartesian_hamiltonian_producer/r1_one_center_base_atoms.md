# R1 One-Center Base Atoms

Status: implemented public-facade contract for explicit origin-centered
all-electron one-center base atoms under `HP-R1-ATOM-*`. The general facade
and base artifact contracts remain canonical in
[R1 public base producer](r1_public_base_producer.md).

## Meaning And Interface

A base atom Hamiltonian is unsupplemented, uncorrected, and all-electron. It
uses the shared terminal Cartesian basis, exact assembled one-body operators,
and localized IDA interaction. The exported call remains:

```julia
cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

No atom-specific function, result wrapper, report, or route object is part of
the public interface.

## System Contract

`system` must contain exactly:

```text
atom_symbols, nuclear_charges, atom_locations, nup, ndn
```

For a one-center atom:

- the three center collections are `AbstractVector` values of length one;
- the location is exactly the finite tuple `(0.0, 0.0, 0.0)`;
- the explicit nuclear charge is finite, positive, and integer-valued;
- `nup` and `ndn` are nonnegative integers and are not `Bool`;
- neutrality requires `nup + ndn == round(Int, nuclear_charge)`.

The atom symbol is converted to a string and persisted as a label. It does not
infer charge, electron count, spin, basis, mapping policy, or ECP behavior.

## Basis Contract

The atom `basis` is a plain exact-key `NamedTuple`. It requires:

```text
core_spacing, radius, and at least one of ns or legacy q
```

All spacings and extents are finite and positive. Size inputs are positive
integers and are not `Bool`. Implemented optional fields and defaults are:

| Field | Default or rule |
| --- | --- |
| `parent_axis_family` | `:G10`; other families reject |
| `reference_spacing` | `1.0` |
| `tail_spacing` | `10.0` |
| `nesting` | `:pqs`; accepts `:pqs` or `:wl` |
| `source_span` | `:ordinary`; `:mapped_comx` is PQS-only |
| `s_factor` | `1.0`; finite and positive |
| `coulomb_accuracy` | `:compact`; current source also accepts `:high` |

`nesting`, `source_span`, and `coulomb_accuracy` accept symbols or strings and
normalize to symbols. `parent_axis_family` must be the symbol `:G10`.
`:standard` Coulomb accuracy is not implemented at this baseline. Public
`parent_mapping_Z`, `parent_mapping_d`, `parent_mapping_rule`, backend
controls, and parent axis counts are unsupported.

Durable calls use `ns`. The producer derives:

```text
nesting = :pqs  -> q = ns
nesting = :wl   -> q = ns - 2
```

Generic WL input requires `ns >= 3`. Legacy `q` may be supplied alone, in
which case `ns` is reconstructed; if both are present they must agree with the
selected nesting.

`core_spacing` is the single public near-nucleus physical scale. Deprecated
atom-only `d`, if supplied, must be finite, positive, and exactly equal to
`core_spacing`. It is rejected for diatomics and does not revive public
`parent_mapping_d`.

## Mapping And Parent Extent

For charge `Z`, core spacing `c`, and mapping factor `f`, the one-center
White-Lindsey mapping is:

```text
standard_s  = sqrt(Z * c)
effective_s = f * standard_s
AsinhMapping(c = c, s = effective_s, tail_spacing = tail_spacing)
```

The exact expert-factor semantics and provenance are owned by
[PQS/WL mapping `s_factor`](pqs_mapping_s_factor.md).

`radius` is physical parent-extent authority. The producer finds the odd
mapped axis count that covers that radius, derives the direct-core minimum from
public `ns`, and uses:

```text
direct_core_side = isodd(ns) ? ns : ns + 1
parent_side = max(mapped_count_covering_radius, direct_core_side)
```

Consequently `ns` does not replace radius, driver padding, reference spacing,
or tail spacing.

## Shared Construction

After one-center geometry and shellification normalization, atoms use the same
producer stages as diatomics:

```text
working terminal basis
-> product/moment operators
-> unit-nuclear operators
-> localized IDA interaction
-> CartesianIDAHamiltonian assembly
-> optional artifact
```

Atom-specific Hamiltonian builders, parallel materialization, and atom-only
one-body or interaction orchestration are outside this contract.

## Source Ownership

The one-center atom facade normalization and shared-workflow wiring are owned
only by `src/cartesian_base_hamiltonian.jl`. The existing export/include in
`src/GaussletBases.jl` belongs to the broader R1 base facade and is unchanged
by the atom relaxation. No atom-only source file, committed fixture, or new
test owner is part of this contract.

## Artifact Behavior

`hamfile === nothing` writes nothing. A nonempty path writes the ordinary
version-1 Cartesian IDA artifact and returns the same in-memory Hamiltonian;
an empty path rejects. The ordinary reader continues to consume only the
matrix and physical payload and ignores provenance sidecars.

The existing base sidecars record the explicit system and basis truth,
including `ns`, derived `q`, `core_spacing`, radius, axis counts, mapping kind,
`mapping_d`, `mapping_s_factor`, `mapping_s_standard`,
`mapping_s_effective`, electron counts, and final dimension. Route is
`:one_center_pqs_base` or `:one_center_wl_base` according to `nesting`. These
values are provenance, not algorithmic readback input; no atom-specific schema
is introduced.

## Failure Behavior

Missing or unknown keys, non-vector center collections, non-tuple or
translated locations, invalid charge/electron counts, invalid size or
spacing, inconsistent `ns`/`q`, incompatible nesting/source-span choices,
invalid `s_factor` or Coulomb policy, mismatched legacy `d`, and empty
`hamfile` throw before expensive construction where practical. Numerical and
filesystem failures propagate. The facade never returns a blocker/status or
partial result.

## Non-Goals

This contract does not authorize supplements or corrections, translated or
multicenter broadening, element tables or automatic defaults, ECPs,
pseudopotentials, solver/RHF workflow, public API redesign, new artifact
fields/readers, route diagnostics, or changes to terminal, raw-block,
Residual Gaussian, MWG, or IDA semantics. Supported supplemented atoms are a
separate shared-composition contract.
