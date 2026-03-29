## Fixed-Radial Legacy `dmrgatom.jl` Export Milestone

This note records the narrow legacy-direct dense export added for the
`dmrgatom.jl` consumer family.

### What Exists Now

The fixed-radial angular sequence line now has a dedicated consumer-specific
writer from `AtomicFixedRadialAngularSequenceLevel`:

- `atomic_fixed_radial_legacy_dmrgatom_payload(level; meta=nothing)`
- `write_atomic_fixed_radial_legacy_dmrgatom_jld2(path, level; meta=nothing)`

The writer targets the minimal legacy direct-dense shape used by
`/Users/srw/Dropbox/GaussletModules/bin/dmrgatom.jl`:

- `meta/Z` as the integer atomic charge expected by the legacy consumer
- `H1`
- `Vee`
- `basis_centers`
- `dims_per_shell`

### Geometry Contract

The exported `basis_centers` are not claimed to be exact physical orbital
centers. They are representative operational geometry for the legacy
within-shell ordering heuristic used by `dmrgatom.jl`.

The current center policy is:

- `representative_dominant_prototype_direction_v1`

For each final shell-local orbital, the exported center is:

- shell radius
- times the sphere-point-set direction associated with that orbital's dominant
  prototype-point coefficient

This is applied to both:

- the exact injected block
- the mixed block

and the rows are emitted in the same final shell-major orbital order already
used by the fixed-radial sequence level.

### Provenance / Audit Fields

The writer records enough information to audit the representative-center
construction:

- center policy name
- within-shell dominant prototype index per orbital
- source angular profile id
- gauge version

These provenance fields are carried under `bridge/`.

### Trust Boundary

This is a narrow legacy-direct export for one specific consumer family only:

- `dmrgatom.jl`

It does not settle broader legacy chemistry-style JLD2 compatibility for:

- `dmrgjld2.jl`
- `dmrgb1.jl`
- chemistry-style `Xmat/Ymat/Zmat` or related bundles

### Practical Consequence

The repo now has an explicit writer for direct legacy `dmrgatom.jl` handoff
from the in-memory fixed-radial angular sequence level, without pretending that
the native `angular_fixed_radial_dense_v1` artifact is itself the legacy file
format.
