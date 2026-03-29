## Legacy-Direct `dmrgatom.jl` Export Plan

This note records the narrow first-pass plan for a legacy-direct dense JLD2
export from the fixed-radial angular sequence line to the legacy
`dmrgatom.jl` consumer family.

### Scope

The target is only:

- `/Users/srw/Dropbox/GaussletModules/bin/dmrgatom.jl`

The intended required payload is the minimal direct-dense shape used by that
consumer:

- `meta/Z`
- `H1`
- `Vee`
- `basis_centers`
- `dims_per_shell`

This note does not broaden the target to:

- `dmrgb1.jl`
- `dmrgjld2.jl`
- chemistry-style `Xmat/Ymat/Zmat`, `Vat`, `Kfinal`, or similar extras
- sampled-variance side products

### Clean Repo Boundary

The clean producer boundary is the in-memory
`AtomicFixedRadialAngularSequenceLevel`, not an on-disk
`angular_fixed_radial_dense_v1` file.

Reason:

- the in-memory level still has the full shell/profile geometry
- the current native level artifact does not carry enough per-orbital geometry
  to reconstruct legacy `basis_centers` cleanly after the fact

So the smallest clean repo-side implementation path is a dedicated payload
builder and writer of the form:

- `atomic_fixed_radial_legacy_dmrgatom_payload(level; meta=nothing)`
- `write_atomic_fixed_radial_legacy_dmrgatom_jld2(path, level; meta=nothing)`

where `level` is an `AtomicFixedRadialAngularSequenceLevel`.

### Straightforward Fields

The following pieces are direct from the current level object:

- `H1` from `level.payload.hamiltonian`
- `Vee` from `level.payload.interaction`
- `dims_per_shell` from `level.shell_dimensions`
- `meta/Z` from the radial/operator source manifest carried by the level
  payload

Under the fixed-radial sequence contract, `dims_per_shell` is the final
shell-major orbital count per radial shell for that level, using the current
profile order inside each shell.

### Geometry and Order

The nontrivial field is `basis_centers`.

`dmrgatom.jl` does not treat `basis_centers` as passive provenance. It uses the
per-orbital coordinates inside each shell to derive a shell-local ordering,
including a north/south pole pinning step from the `z` coordinate. So the
exporter cannot treat `basis_centers` as an arbitrary placeholder.

The raw geometric ingredients do exist at the producer boundary:

- shell radii from `level.shell_centers_r`
- shell-local sphere coordinates from `level.profile.basis.point_set.coordinates`

But there is a real contract question: the final shell-local orbitals are the
profile basis order, not the raw sphere-point order. In particular:

- the exact injected block comes first
- the mixed block follows in deterministic seed order
- those orbitals are not literally the raw point functions

So there is a hidden geometry policy to settle before coding:

- what 3D coordinate should be exported for each final shell-local orbital in
  the level order?

This is the main contract risk for the first implementation.

### Current Recommendation

Proceed with a dedicated writer from `AtomicFixedRadialAngularSequenceLevel`,
but do not implement it as a blind field-mapping pass.

First settle one explicit orbital-center policy for the final shell-local
basis. Examples of the policy level that must be made explicit:

- whether centers are taken as raw sphere-point coordinates in some profile
  order
- whether they are derived from the final basis coefficients against the
  prototype point functions
- or whether another deterministic shell-local center construction is used

Until that policy is fixed, the writer would risk emitting `basis_centers` that
look valid structurally but do not match what `dmrgatom.jl` expects
operationally.

### Trust Boundary

If implemented, this should be documented as:

- a legacy-direct dense export for the `dmrgatom.jl` consumer family only
- built from the fixed-radial angular sequence level object
- separate from the native `angular_fixed_radial_dense_v1` export line

It does not settle the broader legacy chemistry-style JLD2 compatibility
problem.

### Practical Outcome

The first code target should be a writer with an explicit name such as:

- `write_atomic_fixed_radial_legacy_dmrgatom_jld2(...)`

That name keeps this path separate from the native angular sequence export and
makes the consumer-specific trust boundary obvious.
