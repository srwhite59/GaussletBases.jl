# Export layer

The package’s current export layer is producer-side only. It now exposes both:

- public in-memory payload builders
- JLD2 writers that delegate to those payload builders
- one public current-model atomic interaction accessor

The current atomic IDA Hamiltonian model is exported honestly as dense or
sliced/block density-density data for downstream solver consumers.

The Cartesian one-basis IDA boundary is now exposed as a minimal Hamiltonian
object and versioned JLD2 artifact:

- `cartesian_base_hamiltonian`
- `CartesianIDAHamiltonian`
- `one_body_hamiltonian`
- `nuclear_repulsion`
- `write_cartesian_ida_hamiltonian`
- `read_cartesian_ida_hamiltonian`

This format stores `K`, separated uncharged `{U_A}`, `Vee`, nuclear charges,
`ncenter x 3` positions, and spin counts. It derives nuclear repulsion on load
and does not store route diagnostics or solver results.

## Cartesian Base Hamiltonian

`cartesian_base_hamiltonian(system; basis, hamfile=nothing)` is the current
public producer for the base Cartesian IDA Hamiltonian. It returns a
`CartesianIDAHamiltonian{Float64}` directly. If `hamfile` is not `nothing`, it
writes the existing Cartesian IDA Hamiltonian artifact and records fixed
producer provenance under `producer_provenance/`.

The public scope is intentionally narrow:

- one positive-charge atom centered at the origin;
- two distinct, equal-label/equal-charge centers on the Cartesian z axis;
- unsupplemented, uncorrected, all-electron localized-IDA Hamiltonians.

H and H2 are examples, not element restrictions. The facade does not support
translated atoms, arbitrary molecular orientation, heteronuclear diatomics,
supplements, corrections, solver controls, or public route-stage selection.

### Expert staged Cartesian construction

`cartesian_base_working_basis` is an expert/unstable staged-construction
boundary, not a second Hamiltonian facade or a stable result schema.
Ordinary users should call `cartesian_base_hamiltonian`. The staged constructor
validates its system and basis descriptions, resolves the Coulomb expansion,
and builds the parent and terminal realization in memory. Its current input, parent,
terminal, provenance, inventory, and due-diligence fields are expert/unstable,
non-schema inspection state rather than compatibility promises; inventory and
due diligence may be `nothing` when no terminal basis is realized.

This constructor does not assemble the complete one-body or nuclear operator,
`Vee`, or `CartesianIDAHamiltonian`; write an artifact; run a solver or
correction; provide a PRF wrapper; or automate a paper workflow. Its
`source_mode_overrides` layout remains research-only, and PRF definitions
remain private and unexported.

One-center H:

```julia
h_system = (;
    atom_symbols = ["H"],
    nuclear_charges = [1.0],
    atom_locations = [(0.0, 0.0, 0.0)],
    nup = 1,
    ndn = 0,
)

h_basis = (;
    q = 5,
    core_spacing = 0.5,
    radius = 4.0,
    reference_spacing = 1.0,
)

h_ham = cartesian_base_hamiltonian(h_system; basis = h_basis)
```

Z-axis H2:

```julia
h2_system = (;
    atom_symbols = ["H", "H"],
    nuclear_charges = [1.0, 1.0],
    atom_locations = [(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)],
    nup = 1,
    ndn = 1,
)

h2_basis = (;
    q = 5,
    core_spacing = 0.5,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
)

h2_ham = cartesian_base_hamiltonian(
    h2_system;
    basis = h2_basis,
    hamfile = "h2_cartesian_ida_hamiltonian.jld2",
)
```

For a one-center atom, `d` is a deprecated compatibility input: omit it, as
above, or set it exactly equal to `core_spacing`. Diatomics reject `d`.
For PQS, `q = ns`; for White-Lindsey, `q = ns - 2`. Supplying both size fields
requires consistency, while `ns` is the natural matched-comparison input.

For the current SlicedMRGUtils / HamIO bridge family, the package also exposes a
thin explicit compatibility adapter:

- the native sliced payload keeps the package’s current `l0_desc_mzigzag`
  within-slice ordering
- the HamV6 compatibility payload/writer reorders within each slice to
  `mzigzag_then_l`

Both payload families now also carry one shared package-owned producer/source
manifest in `meta/manifest/...`, including the atomic charge, public radial
extent recipe, mapping parameters, radial dimension, channel convention, and
producer identity.

The angular research track now also exposes one narrow downstream-facing bridge
path:

- `angular_benchmark_exact_hamv6_payload`
- `write_angular_benchmark_exact_hamv6_jld2`

That bridge is intentionally limited to the benchmark line's exact common
low-`l` reference. The full mixed shell-local angular basis is not yet exported
as HamV6, because the current consumer language still assumes definite
per-orbital `l,m` labels.

For the current mixed-basis angular line, the direct in-memory HF handshake is
instead:

- `build_atomic_injected_angular_hfdmrg_payload`

That payload builder hands dense `H`, `V`, seeds, and occupations directly to
`HFDMRG.solve_hfdmrg(...)` without claiming full mixed-basis HamV6
compatibility.

The current branch-point note for that boundary is:

- `docs/angular_consumer_contract_boundary.md`

The same angular line now also exposes one native fixed-radial sequence export
surface for increasing-`N_sph` ladders:

- `build_atomic_fixed_radial_angular_sequence`
- `atomic_fixed_radial_angular_level_dense_payload`
- `write_atomic_fixed_radial_angular_level_jld2`
- `atomic_fixed_radial_angular_overlap_sidecar_payload`
- `write_atomic_fixed_radial_angular_overlap_sidecar_jld2`

This is a producer-side contract for external continuation studies. It exports:

- one dense native level artifact per `N_sph`
- one adjacent shell-local overlap sidecar per `N_sph[k] -> N_sph[k+1]`
- one full non-adjacent upper-triangle family of direct shell-local overlap
  sidecars inside the same sequence
- stable radial shell ids and shell centers across the sequence
- stable within-shell labels from cached shell-local angular profiles

It does not yet include the later common-target embedding/lift layer, and it
does not claim compatibility with older dense consumer formats.

## Experimental Angular Profile And Sequence Producers

These bindings expose experimental producer-side profile, overlap, and
fixed-radial sequence state. They support external continuation and transfer;
they do not provide restart orchestration or a completed angular application
workflow.

```@docs
ShellLocalAngularProfile
ShellLocalAngularProfileOverlap
shell_local_angular_profile
adjacent_shell_local_angular_profile_overlap
AtomicFixedRadialAngularSequenceLevel
AtomicFixedRadialAngularSequenceOverlapSidecar
AtomicFixedRadialAngularSequence
```

For the narrative explanation of the current producer-side story, see:

- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Example guide](../howto/example_guide.md)

## External Cartesian GTO Transfer

The packet, import, and determinant-preparation bindings used by the
[manual workflow](../manual/external_cartesian_gto_transfer.md) are included in
the curated reference below. The overlap operation is documented separately
with the ordinary Cartesian probe API.
External-GTO fingerprints are strict packet-integrity hashes, not numerical-
equivalence, permutation-, tolerance-, or convention-invariant comparisons.

```@docs
cartesian_base_hamiltonian
cartesian_base_working_basis
atomic_ida_density_interaction_matrix
CartesianIDAHamiltonian
cartesian_residual_gto_mwg_system
ExternalGTOOrbitalSpinBlock
ExternalGTOOrbitalPacket
ExternalGTOOrbitalSpinImport
ExternalGTOOrbitalImportResult
import_external_gto_orbitals
read_external_cartesian_gto_packet
closest_external_gto_determinant
external_gto_ordering_fingerprint
external_gto_overlap_fingerprint
PQSH2PlusRow
PQSH2PlusComparison
pqs_h2plus_comparison
ExactRepresentedHartreeField
FittedReferenceHartreeField
ScreenedHartreeCorrection
screened_hartree_correction
screened_hartree_delta_one_body
screened_hartree_energy_constant
screened_hartree_consistency_error
screened_hartree_field_kind
one_body_hamiltonian
nuclear_repulsion
build_atomic_fixed_radial_angular_sequence
atomic_fixed_radial_angular_level_dense_payload
atomic_fixed_radial_angular_overlap_sidecar_payload
fullida_dense_payload
sliced_ham_payload
atomic_hamv6_payload
angular_benchmark_exact_hamv6_payload
write_atomic_fixed_radial_angular_level_jld2
write_atomic_fixed_radial_angular_overlap_sidecar_jld2
write_cartesian_ida_hamiltonian
read_cartesian_ida_hamiltonian
write_fullida_dense_jld2
write_sliced_ham_jld2
write_atomic_hamv6_jld2
write_angular_benchmark_exact_hamv6_jld2
```
