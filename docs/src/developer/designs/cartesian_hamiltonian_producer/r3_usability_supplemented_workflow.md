# R3 Usability Supplemented Workflow

Status: implemented internal supported facade for residual-GTO/MWG
supplemented Hamiltonians. The function is module-qualified and non-exported.

## Owned IDs

- `HP-R3U-FILE-01` - implemented source and validation surfaces;
- `HP-R3U-FN-01` - implemented non-exported supplemented facade;
- `HP-R3U-WIRE-01` - implemented same-construction composition;
- `HP-R3U-TEST-01` - implemented standalone H2 facade gate.

The homonuclear z-axis scope and canonical-driver wiring are governed by the
implemented `HP-R3U-ZDI-*` contract in
[R3 homonuclear supplemented workflow](r3_homonuclear_diatomic_supplemented_workflow.md).
One-center supplemented composition reuses this implementation under separate
`HP-COMP-SUPPATOM-*` authority; it does not broaden the R3U molecular scope.

## Interface

The implemented call is:

```julia
cartesian_residual_gto_mwg_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    supplement::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

It is defined in `src/cartesian_base_hamiltonian.jl` and is intentionally not
exported from `GaussletBases`. It returns the existing Hamiltonian directly,
not a wrapper, report, status object, payload, or `(value, metadata)` pair.

## Supported Molecular Input

The R3U molecular path accepts explicit neutral all-electron homonuclear
two-center Cartesian z-axis systems:

- `system` has exactly `atom_symbols`, `nuclear_charges`, `atom_locations`,
  `nup`, and `ndn`;
- center collections are vectors of equal length;
- the two symbols and nuclear charges are equal;
- both centers have zero `x` and `y`, finite distinct `z`, and positive finite
  charges;
- `nup` and `ndn` are nonnegative integers whose sum equals the integer total
  nuclear charge.

There is no element-specific H, Be, or Cr branch. Heteronuclear, charged,
ECP, non-z-axis, translated/rotated general molecular, and solver inputs fail
validation. Cr2 may use the generic internal path, but it has no special
default, fixture, or production claim here.

## Basis Input

For the diatomic path, `basis` requires:

- `core_spacing`;
- `xmax_parallel`;
- `xmax_transverse`;
- at least one of `ns` or legacy-compatible `q`.

The shared base producer owns normalization of `ns`, route-local `q`,
`nesting`, `source_span`, `s_factor`, `coulomb_accuracy`,
`parent_axis_family`, `reference_spacing`, and `tail_spacing`. This facade does
not define a second basis policy. Unknown keys and inconsistent `ns`/`q`
values fail before construction. White-Lindsey minimum-`ns` behavior belongs
to the composition contracts.

## Supplement Input

`supplement` requires:

- `basis_by_center::AbstractVector` of string labels;
- integer `lmax` with `0 <= lmax <= 6`.

Optional fields are:

- `uncontracted::Bool = false`;
- `width_filtering = nothing` or exactly `(; max_width)` with positive finite
  `max_width`;
- `basisfile = nothing` or a trusted local/project path string.

The basis-label count must match the center count. The current homonuclear
path requires identical labels on both centers and uses the existing named
Cartesian supplement loader. Fitted fields, ECP data, element defaults, and
heteronuclear basis-by-center composition are not part of this contract.

## Same-Construction Composition

The facade performs one construction:

```text
validated system, basis, and supplement
-> cartesian_base_working_basis(...; supplemented=true)
-> base products, unit nuclear matrices, Vee, and base Hamiltonian
-> named Gaussian supplement representation
-> owner-local residual Gaussian basis
-> exact augmented one-body and moment matrices
-> residual-containing MWG/IDA interaction
-> CartesianIDAHamiltonian
-> optional existing artifact writer
```

The base Hamiltonian, terminal basis, parent bundles, Gaussian supplement,
residual object, and producer-owned Coulomb expansion therefore share one
construction. The facade must not be replaced by a post-hoc augmentation of
an arbitrary dimension-compatible Hamiltonian.

Residual-basis, exact-operator, and MWG formulas belong to
[Residual Gaussian domain module](residual_gaussian_domain_module.md). The
terminal compatibility functions in
`src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` compose
raw blocks, enforce same-construction checks, delegate numerical physics to
the domain module, and assemble the existing Hamiltonian type.

## Artifact Boundary

If `hamfile === nothing`, no file is written. Otherwise the facade writes the
existing Cartesian IDA Hamiltonian artifact, compact supplement provenance,
and the current Hamiltonian manifest, then still returns the in-memory
Hamiltonian. Empty paths fail; ordinary file-system errors propagate.

Artifact identity, keys, provenance, compatibility, and readback are owned by
[Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md).
This facade does not expose residual transforms, raw stage objects, MWG
descriptors, pair factors, or provenance payloads.

## Validation

The committed focused gate is:

```text
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
```

Its facade section checks the returned type, same-construction H2 endpoint,
artifact/readback matrix parity, compact provenance, Coulomb provenance, and
malformed or unsupported input failures. At the Pass 377 baseline, the H2
facade fixture has base dimension `487`, residual dimension `18`, and final
dimension `505`; the accepted lowest-orbital IDA self-Coulomb reference is
`0.4574161883692301` for the current compact fixture.

The file is a standalone focused gate, not normal `Pkg.test` pressure. Be2 and
Cr2 measurements remain ignored/user-run evidence unless separately promoted.

## Failure And Non-Goals

Invalid keys, dimensions, ownership, centers, nuclear charges, carried Coulomb
expansion, PGDG exponent parity, exact base blocks, residual geometry, or
artifact path fail rather than returning a status payload.

This contract does not authorize:

- public export or broad API redesign;
- general molecular orientation or heteronuclear/ECP support;
- solver, RHF/UHF, EGOI, or screened-Hartree workflow;
- new artifact shapes or residual-basis serialization;
- route reports, persistent factor caches, or exposed construction stages;
- Cr2-specific branches, defaults, committed fixtures, or scientific claims.

Implementation history and superseded R3-A/B/C fixture narratives are indexed
in [R3 compatibility history](r3_residual_gto_mwg_augmentation.md).
