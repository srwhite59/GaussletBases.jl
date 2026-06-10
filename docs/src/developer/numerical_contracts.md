# Numerical Contracts

This page records a small internal engineering policy that is easy for machine
generated code to miss.

It is intentionally developer-facing rather than part of the main user manual.

## Orthonormal blocks

When a block is constructed to be orthonormal, the intended contract is:

1. build it to be orthonormal in the relevant metric
2. check that the resulting overlap is `I +` small unavoidable Float64 noise
3. then treat the block as orthonormal

Examples include:

- finalized PGDG / COMX-cleaned blocks
- nested fixed blocks
- any other internal block whose construction is explicitly meant to produce an
  orthonormal basis

## What not to do

Do not keep propagating a near-identity overlap matrix as if it were meaningful
mathematical data.

In particular:

- do not store `S = I + ε` by default when `ε` is just Float64 residue
- do not build downstream logic that keeps consulting such an `S`
- do not interpret tiny nonorthogonality at the `1e-12` to `1e-14` level as a
  real structural feature

## What to do instead

Use overlap matrices in this regime only for:

- construction
- validation
- assertions
- or one final cleanup step if the deviation from identity is not small enough

For transfers between final orthonormal working bases:

- use only the cross overlap between the two final bases
- treat self-overlaps of those final bases as diagnostic-only
- do not turn the final working path into a generalized-overlap formulation

After that, the working representation should simply regard the block as having
identity overlap.

This is a coding and design policy, not a user-facing scientific statement.

## Nested fixed-block kinetic

`_NestedFixedBlock3D.kinetic` follows the nested packet contract.

That means:

- it is the kinetic matrix carried by the assembled nested packet
- it is the kinetic payload that downstream nested operator routes should use
- it is not automatically interchangeable with "contract the ordinary parent
  kinetic later and call that the same thing"

For current nested diatomic routes, this distinction is real. A nested
fixed-block kinetic can differ measurably from a later contraction of a
separately assembled ordinary parent one-body path even when both live on the
same final basis dimension.

## One-body reassembly

`assembled_one_body_hamiltonian(...)` reassembles from the operator payload that
was actually stored:

- stored `kinetic_one_body`
- stored `nuclear_one_body_by_center`
- requested `nuclear_charges`

So the meaningful contract is:

- compare reassembled one-body matrices against other operators built from the
  same stored kinetic contract

Do not use `assembled_one_body_hamiltonian(...)` to assert equality between two
routes that already disagree about what the kinetic payload is supposed to be.

In particular, for nested fixed-block routes:

- a by-center payload should be compared against another by-center payload on
  the same nested kinetic contract
- not against a separate total-only route that rebuilds one-body terms through a
  different parent-space contraction path

## WL Gausslet-Only H Atom Acceptance Baselines

The WL/Cartesian gausslet-only hydrogen acceptance check is
`test/nested/cartesian_wl_gausslet_h_atom_acceptance_runtests.jl`. It is a
bounded scientific gate, not an exhaustive local-helper test.

Small fitted-strength fixture:

- one proton at `(0.0, 0.0, 0.0)` with `Z = 1.0`
- `MappedUniformBasisSpec(:G10)` with `count = 7`
- `fit_asinh_mapping_for_strength(s = 0.5, npoints = 7, xmax = 6.0)`
- `reference_spacing = 1.0`
- backend `:pgdg_localized_experimental`
- Coulomb expansion `coulomb_gaussian_expansion(doacc = false)`
- basis dimension / parent support size `343`
- generalized solve against the carried Cartesian overlap

The current small-fixture lowest one-electron energy is approximately
`-0.4706400351534759` Hartree, so the acceptance window is deliberately loose:
`-0.5 < E < -0.45`. This verifies the expected variational side of the
hydrogen ground state while keeping the fixture small.

Coarse/distorted WL-style fixture:

- one proton at `(0.0, 0.0, 0.0)` with `Z = 1.0`
- `MappedUniformBasisSpec(:G10)` with `count = 11`
- `AsinhMapping(c = 0.25, s = 1.0, tail_spacing = 10.0)`
- `reference_spacing = 1.0`
- requested old-style core-spacing parameter `d = 0.25` is represented by the
  repo mapping parameter `c = 0.25`
- requested `xmax = 8.0` is represented by the bounded `count = 11` endpoint,
  giving outer centers at approximately `+/-8.182848593885893`
- backend `:pgdg_localized_experimental`
- Coulomb expansion `coulomb_gaussian_expansion(doacc = false)`
- basis dimension / parent support size `1331`
- generalized solve against the carried Cartesian overlap

The current coarse/distorted lowest one-electron energy is approximately
`-0.4966106635473884` Hartree, with acceptance window `-0.5 < E < -0.49`.
The test reports the solve kind, basis size, spacing/core-spacing parameter,
distortion strength, requested `xmax`, observed energy, and distance from
`-0.5` Hartree for both fixtures.

These checks do not include GTO supplements, PQS retained transforms, CPB/GTO
bundle consumption, route/global refactors, or Hamiltonian assembly beyond the
one-electron `H1 = kinetic + nuclear attraction` path.
