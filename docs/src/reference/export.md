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

The R1 public scope is intentionally narrow:

- origin-centered hydrogen atom;
- hydrogen molecule on the Cartesian z axis;
- unsupplemented, uncorrected, all-electron localized-IDA Hamiltonians.

It does not support arbitrary molecular orientation, x/y-aligned H2, other
atoms, supplements, corrections, solver controls, or public route-stage
selection.

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
    d = 0.3,
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

For H, `d` is the White-Lindsey atomic mapping parameter and is independent of
`core_spacing` and `reference_spacing`. For H2, `d` is not a public input.

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

For the narrative explanation of the current producer-side story, see:

- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Example guide](../howto/example_guide.md)

```@docs
atomic_ida_density_interaction_matrix
CartesianIDAHamiltonian
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
