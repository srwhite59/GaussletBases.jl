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

- do not store `S = I + eps` by default when `eps` is just Float64 residue
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

## Decomposed WL Gausslet-Only Acceptance Boundary

The intended gausslet-only scientific acceptance path is a true decomposed
White-Lindsey calculation with `q = 5` and `ns = 5`, not a single CPB covering
the full parent product window. The active readiness check is
`test/nested/cartesian_wl_gausslet_h_atom_acceptance_runtests.jl`.

Current status:

- decomposed White-Lindsey overlap and kinetic local pair-block paths exist
- global pilots exist for those supported safe one-body terms
- CPB-local electron-nuclear by-center provider blocks exist
- decomposed White-Lindsey electron-nuclear by-center selector, placement plan,
  and global by-center matrix path do not exist yet
- therefore active H and H2+ scientific acceptance through the decomposed WL
  path is blocked on
  `:missing_decomposed_wl_electron_nuclear_by_center_placement`

The retired transition helpers built one full-parent CPB with role
`:wl_cpb_acceptance_full_parent_window`. That path exercised CPB-local
operators but not White-Lindsey boundary-unit decomposition, so it is not an
active acceptance contract. Its observed values are retained only as transition
notes:

- H full-window CPB energy `-0.4832079279118124` Hartree
- H2+ full-window CPB electronic energy `-1.0971828374927926` Hartree
- H2+ full-window CPB total energy `-0.5971828374927926` Hartree

Older direct-route transition baselines remain historical comparison points
only:

- H small fitted direct route `-0.4706400351534759` Hartree
- H coarse/distorted direct route `-0.4966106635473884` Hartree
- H2+ R = 2.0 direct electronic energy `-1.0654839328172023` Hartree
- H2+ R = 2.0 direct total energy `-0.5654839328172023` Hartree

The next implementation needed before restoring active scientific H/H2+
acceptance is a decomposed WL electron-nuclear by-center path that keeps center
records separate, applies nuclear charges at the acceptance/Hamiltonian
assembly boundary, and places the resulting by-center blocks through real
boundary-unit ranges. Do not reintroduce the full-parent CPB helper or a direct
Cartesian product fallback as the active route.
