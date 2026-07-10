# Atomic HF Reference Packets

Status: implemented internal subsystem under
`HP-PQS-ATOMREF-PACKET-FN-01` and
`HP-PQS-ATOMREF-PACKET-TEST-01`. A fixed-policy determinant-moment polish of
the fitted potential is approved, pending implementation, under
`HP-PQS-ATOMREF-POTMOM-FN-01` and
`HP-PQS-ATOMREF-POTMOM-TEST-01`.

This page is the canonical contract for reusable one-center atomic
Hartree-Fock (HF) reference packets. The registry owns permission, lifecycle,
and file surfaces; this page owns packet semantics, validation, and failure
behavior.

## Purpose And Scope

An atomic HF reference packet preserves one reproducible atom/basis reference
for screened-Hartree and related internal consumers. It combines three
different physical or numerical objects without confusing their roles:

1. a pure-GTO occupied determinant;
2. a near-exact spherical Gaussian fit to that determinant density;
3. a compact Gaussian fit to the density's radial Hartree potential.

The implemented initial scope is:

- one-center closed-shell RHF references;
- Be `Z = 4`, `2e` core `1s^2` and Ne all-electron `10e`;
- cc-pV5Z with `lmax = 1`;
- packet construction, validation, write/readback, and bounded Be/Ne
  consumption checks.

Atom, nuclear charge, electron count, basis, center, occupancy/fill-shell
convention, and construction controls are explicit inputs. Element tables do
not infer electron count, spin, or filled shells.

## Authoritative Roles

The pure-GTO occupied coefficients and occupations define the reference
density:

```text
P0_AA = C_occ * Diagonal(occupations) * C_occ'
q0_AA = diag(P0_AA)
```

The density fit is a compressed representation of that same determinant
density. It defines the reference cloud and its Coulomb self-energy.

The potential fit is only a fast evaluator of the Hartree potential generated
by the fitted density. It does not define `P0`, `q0`, or the self-energy.

Density-fit and potential-fit Gaussian terms are evaluation data. They are
not supplement orbitals, protected basis functions, or occupied directions.
Consumers needing exact determinant representation must protect or import the
packet's original occupied GTO space.

## Packet Contents

The packet records:

- system facts: atom, nuclear charge, electron count, center, explicit
  fill-shell convention, basis name, `lmax`, and construction controls;
- supplement identity: ordered labels, angular powers, centers, primitive
  exponents and coefficients, overlap matrix, and overlap fingerprint;
- HF reference: occupied coefficients in supplement coordinates,
  occupations, orbital energies, density matrix, energy, convergence, and
  iteration diagnostics;
- density fit: spherical Gaussian widths/exponents and signed weights, radial
  target and fit data, charge, self-energy, and fit diagnostics;
- potential fit: Gaussian coefficients/exponents, fixed broad-tail and
  trim/refit facts, radial errors, tail diagnostics, and an optional compact
  determinant-moment polish record;
- provenance: code/input identity, fit tolerances, reference energy where
  supplied, and role-qualified Coulomb-expansion summaries.

The ordered supplement and its fingerprint are part of packet identity.
Consumers must match the packet to the current supplement by explicit order,
overlap, and fingerprint facts. Label-only orbital selection is not a packet
mapping contract.

## Convergence And Occupancy

Only an explicitly converged RHF result may define a packet. Construction must
stop before density or potential fitting when RHF is unconverged. Writing must
reject an in-memory packet marked unconverged, and readback/validation must
report convergence explicitly. Packet consumers must reject an unconverged
packet before using its determinant or fitted fields, including in diagnostic
mode.

There is no `allow_unconverged` packet-construction option.

Occupations and the fill-shell convention are explicit packet facts. The
occupied block must be orthonormal in the packet supplement metric, and its
density trace must reproduce the stated electron count.

## Coulomb Accuracy Roles

Packet Coulomb provenance is role-qualified rather than one
Hamiltonian-wide policy:

- atomic reference RHF uses the high-accuracy Coulomb expansion;
- density target/self-energy evaluation uses the compact packet-reference
  expansion;
- the fitted-potential tail scaffold uses the compact expansion.

The density-fit `J0_G` path must resolve its compact packet-reference expansion
and pass that object explicitly to the mixed-Hartree helper. It must not
inherit a helper default. The potential fit remains an approximation to the
same density-fit Hartree field and must be checked against that field on a
bounded matrix/anchor case.

These role-qualified packet choices do not change the producer-wide
`coulomb_accuracy` policy of a Hamiltonian consuming the packet.

## Determinant-Moment Potential Polish

The packet determinant, density fit, and fitted potential have distinct roles,
but a non-diagnostic screened-Hartree anchor requires them to agree on the
reference Coulomb energy. For one packet at the origin and an identical packet
displaced by radial distance `d`, define:

```text
m_j(d) = Tr(P_HF * V_j(d))
t(d)   = (rho_fit(0) | rho_fit(d))_compact
```

`V_j(d)` is the exact pure-GTO matrix of one retained spherical Gaussian
potential term with the packet potential exponent and unit coefficient.
`P_HF` is the pure-GTO packet determinant. `t(d)` uses the packet density fit
and its explicit role-qualified compact Coulomb expansion. Thus this is a
consistency fit for determinant-defined `P0`, density-fit-defined `E0`, and
potential-fit-defined `J0`; it is not a scalar energy correction.

For two identical placed packets separated by `d`, the screened-Hartree direct
anchor residual is exactly the sum of the two self and two cross moment
residuals:

```text
Tr(P0 * J0_fit) - E0_fit = 2*(m(0)'c - t(0)) + 2*(m(d)'c - t(d))
```

The moment objective therefore repairs the packet representation that enters
the anchor rather than adding a posterior scalar to the Hamiltonian.

Start from the existing retained `33`-term fitted potential with coefficients
`c0`. Preserve every exponent and keep the first `5` broad-tail coefficients
exactly fixed. For coefficients `6:33`, solve once for `delta` on the fixed
distance grid:

```text
D = [0, 0.1, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 8, 10, 20] bohr

min_delta ||W_radial * A_radial * delta||_2^2
        + ||1e4 * (M_free * delta - (t - M*c0))||_2^2
```

Use the existing radial grid, `_radial_fit_weights(...)`, column-scaled
`_solve_weighted_svd(...)`, and fixed SVD `rtol = 1e-12`. The radial target is
zero change from the already-accepted radial fit. Apply the resulting `delta`
once; do not iterate, tune by element, expose controls, or add Gaussian terms.
The fixed policy ID is
`:determinant_densityfit_coulomb_moment_v1`.

The initial scope is the existing one-center spherical closed-shell Be and Ne
packets and identical-packet radial pairs. The multi-distance grid prevents the
Be2 `R = 4` gate from becoming a one-geometry fit. This lane does not establish
a heteronuclear cross-anchor policy. If the retained potential does not have
the approved `33`-term, five-fixed-term shape, or the weighted solve does not
meet rank, radial, tail, and moment checks, stop rather than generalizing or
silently marking the packet polished.

New packets carry one nested `potential_fit` moment-polish record with:

- policy ID;
- retained SVD rank;
- maximum coefficient change;
- maximum absolute moment error.

The policy ID fixes the distance grid, moment weight, SVD tolerance, and five
fixed broad terms. Readback must return the record when present. Older packets
without it remain readable and report the record as unavailable/`nothing`;
absence must never be interpreted as polished. The artifact kind and existing
packet convention ID do not change.

## Validation And Failure Behavior

Required packet checks are:

- explicit RHF convergence;
- write/readback roundtrip;
- occupied-metric orthogonality and density trace;
- density reconstructed from occupied coefficients agrees with the stored
  density;
- current supplement ordering, overlap, and fingerprint match;
- density-fit charge and self-energy errors are reported and within the
  intended fit tolerance;
- potential-fit radial and far-tail diagnostics are reported;
- potential-fit `J0_G` agrees with the density-fit mixed-Hartree path on a
  bounded matrix and energy-anchor check;
- when determinant-moment polishing is present, the term count and exponents
  are unchanged, the first five coefficients are unchanged, the maximum Be/Ne
  moment error is at most `1e-9 Ha`, and radial/tail diagnostics do not regress
  materially;
- the density-fit path demonstrably receives the explicit compact expansion;
- packet outputs contain no unreported fallback reference or inferred
  occupancy convention.

Fail rather than fitting, writing, or consuming when convergence, packet
identity, occupied orthogonality, density trace, or required fit validation is
not satisfied. Do not repair packet mismatch by relabeling orbitals, changing
occupations, treating fit terms as orbitals, or silently choosing another
Coulomb expansion.

## Ownership And Tests

Implementation and test ownership is recorded in the compact registry entry
for `HP-PQS-ATOMREF-PACKET-*`. The bounded tests cover packet roundtrip,
validation, unconverged-reference rejection, Be/Ne consumption, explicit
compact expansion passage, and the vendored basis-data regression.

The packet facility depends on the historical basis collection described in
[Legacy BasisSets provenance](../../../../legacy_basissets_provenance.md).
That page owns snapshot hashes, parser counts, mixed provenance, unresolved
licensing status, and loader policy. Those data facts are not packet
semantics.

Numerical implementation and measurement evidence remains in the manager
running log, especially Passes 325, 330, 331, and 345.

## Related Consumers

- [Screened Hartree correction assembly](screened_hartree_correction_assembly.md)
  consumes represented packet determinants and fitted reference fields.
- [Protected additive atomic reference correction](protected_additive_reference_correction.md)
  protects placed packet occupied spaces while preserving original packet
  blocks for additive molecular `P0`.
- [Screened Hartree residual-density formalism](screened_hartree_residual_density.md)
  owns the physical residual-density decomposition and measurement history.

## Exclusions

This packet contract does not approve:

- production corrected Hamiltonians or corrected protected artifacts;
- public driver defaults or solver workflow;
- exchange, EGOI, or row-gauge rho0/P0 shortcuts;
- Cr/Cr2 production claims;
- fitted density or potential terms as protected orbitals;
- element-table inference of charge, electron count, spin, or occupancy;
- a general basis download, conversion, or licensing subsystem.
- a public potential-polish option, element-specific tuning, heteronuclear
  moment policy, scalar anchor patch, or weakening of the `1e-8 Ha` screened-
  Hartree anchor.
