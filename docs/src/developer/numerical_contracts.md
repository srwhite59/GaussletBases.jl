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

## Post-CPB WL Gausslet-Only H Atom Acceptance Baseline

The WL/Cartesian gausslet-only hydrogen acceptance check is
`test/nested/cartesian_wl_gausslet_h_atom_acceptance_runtests.jl`. It is a
bounded scientific gate through the post-CPB White-Lindsey local operator
path, not an exhaustive local-helper test.

Active fixture:

- one proton at `(0.0, 0.0, 0.0)` with `Z = 1.0`
- `q = 5`, `n_s = q`
- standard q-to-core-spacing rule, recorded as
  `:standard_pqs_ns_equals_q`
- core spacing `d = 0.15`
- `MappedUniformBasisSpec(:G10)` parent axes with counts `(15, 15, 15)`
- `white_lindsey_atomic_mapping(Z = 1.0, d = 0.15, tail_spacing = 10.0)`
- `reference_spacing = 1.0`
- backend `:pgdg_localized_experimental`
- Coulomb expansion `coulomb_gaussian_expansion(doacc = false)`
- basis dimension / parent support size `3375`
- generalized solve against the carried Cartesian overlap
- overlap, kinetic, and by-center electron-nuclear blocks materialized by
  `CartesianCPBBlockProviders`

The current post-CPB WL lowest one-electron energy is approximately
`-0.4832079279118124` Hartree. The test verifies the variational side of the
hydrogen ground state, records distance from `-0.5` Hartree, and records
distance from the old direct-route transition baselines:

- small fitted direct route `-0.4706400351534759` Hartree
- coarse/distorted direct route `-0.4966106635473884` Hartree

Those direct-route values are reference baselines for the transition only.
They are not the ongoing route under test.

## Post-CPB WL Gausslet-Only H2+ Acceptance Baseline

The WL/Cartesian gausslet-only H2+ acceptance check is
`test/nested/cartesian_wl_gausslet_h2plus_acceptance_runtests.jl`. It is a
bounded one-electron diatomic gate through the same post-CPB local operator
path, not a route/global or helper-level test.

Fixture:

- two protons on the z axis at `(0.0, 0.0, -1.0)` and `(0.0, 0.0, 1.0)`
- internuclear distance `R = 2.0` bohr
- nuclear charges `(1.0, 1.0)`
- `q = 5`, `n_s = q`
- standard q-to-core-spacing rule, recorded as
  `:standard_pqs_ns_equals_q`
- core spacing `d = 0.15`
- `MappedUniformBasisSpec(:G10)` parent axes with counts `(15, 15, 17)`
- `white_lindsey_atomic_mapping(Z = 1.0, d = 0.15, tail_spacing = 10.0)`
- `reference_spacing = 1.0`
- backend `:pgdg_localized_experimental`
- Coulomb expansion `coulomb_gaussian_expansion(doacc = false)`
- basis dimension / parent support size `3825`
- generalized solve against the carried Cartesian overlap
- overlap, kinetic, and two by-center electron-nuclear blocks materialized by
  `CartesianCPBBlockProviders`

The electronic Hamiltonian is

```text
H_elec = kinetic - 1/r_A - 1/r_B
```

and the Born-Oppenheimer total reported by the test is

```text
E_total = E_electronic + 1/R.
```

At `R = 2.0`, the reference electronic energy is approximately
`-1.1026342144949465` Hartree and the reference total energy is approximately
`-0.6026342144949465` Hartree. The current post-CPB WL fixture reports

- electronic energy `-1.0971828374927926` Hartree
- proton-proton repulsion `0.5` Hartree
- total energy `-0.5971828374927926` Hartree
- total-energy distance from the `R = 2.0` reference of approximately
  `0.005451377002153923` Hartree

The old direct-route H2+ transition baseline was

- electronic energy `-1.0654839328172023` Hartree
- total energy `-0.5654839328172023` Hartree

That direct-route value remains only a transition reference. The active test
asserts the post-CPB path materializes/consumes CPB-local operator blocks and
does not call `ordinary_cartesian_ida_operators` or hand-build the full
Cartesian product matrices directly from 1D factors.

These checks do not include GTO supplements, PQS retained transforms, CPB/GTO
bundle consumption, route/global refactors, or Hamiltonian assembly beyond the
one-electron `H1 = kinetic + nuclear attraction` path.
