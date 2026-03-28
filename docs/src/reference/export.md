# Export layer

The package’s current export layer is producer-side only. It now exposes both:

- public in-memory payload builders
- JLD2 writers that delegate to those payload builders
- one public current-model atomic interaction accessor

The current atomic IDA Hamiltonian model is exported honestly as dense or
sliced/block density-density data for downstream solver consumers.

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
- stable radial shell ids and shell centers across the sequence
- stable within-shell labels from cached shell-local angular profiles

It does not yet include the later common-target embedding/lift layer, and it
does not claim compatibility with older dense consumer formats.

For the narrative explanation of the current producer-side story, see:

- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Example guide](../howto/example_guide.md)

```@docs
atomic_ida_density_interaction_matrix
build_atomic_fixed_radial_angular_sequence
atomic_fixed_radial_angular_level_dense_payload
atomic_fixed_radial_angular_overlap_sidecar_payload
fullida_dense_payload
sliced_ham_payload
atomic_hamv6_payload
angular_benchmark_exact_hamv6_payload
write_atomic_fixed_radial_angular_level_jld2
write_atomic_fixed_radial_angular_overlap_sidecar_jld2
write_fullida_dense_jld2
write_sliced_ham_jld2
write_atomic_hamv6_jld2
write_angular_benchmark_exact_hamv6_jld2
```
