# Atomic HF Reference Packets

Status: implemented internal subsystem under
`HP-PQS-ATOMREF-PACKET-FN-01` and
`HP-PQS-ATOMREF-PACKET-TEST-01`. The determinant-moment fitted-potential
polish formerly recorded under `HP-PQS-ATOMREF-POTMOM-FN-01` and
`HP-PQS-ATOMREF-POTMOM-TEST-01` is retired and must not be produced or
consumed.

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
  trim/refit facts, radial errors, tail diagnostics, matrix comparisons, and
  reference-energy consistency diagnostics;
- provenance: code/input identity, fit tolerances, reference energy where
  supplied, and role-qualified Coulomb-expansion summaries.

The ordered supplement and its fingerprint are part of packet identity.
Label-only orbital selection is not a packet mapping contract.

## Packet Integrity And Owner-Local Embedding

Packet self-integrity and molecular embedding equivalence are separate gates.
Before any embedding, recompute the stored overlap fingerprint and require

```text
fingerprint(packet.overlap) == packet.overlap_fingerprint
```

exactly. A stale or corrupted stored packet must fail. This check is not
relaxed by an embedding tolerance.

Mapping a packet into a translated or reconstructed owner-local molecular
supplement must still match exactly on atom and basis identity, packet function
count, owner indices, placement center within the existing tolerance, ordered
labels, angular powers, and packet-to-molecular column order. After those
structural checks, compare the mapped owner-local overlap block numerically:

```text
S_packet = packet.overlap
S_block  = S_AA[packet_to_molecular_order, packet_to_molecular_order]

norm(S_block - S_packet, Inf) <= overlap_atol
```

`overlap_atol` defaults to `1e-10` and must be finite and nonnegative. The raw
SHA-256 fingerprint of `S_block` is diagnostic only: translation or equivalent
reconstruction may change final `Float64` bits without changing the physical
overlap matrix. A mapped-block fingerprint mismatch is therefore not a failure
when every structural check and the numerical overlap bound passes.

Embedding returns one nested internal overlap-mapping summary containing:

- stored packet fingerprint;
- recomputed packet fingerprint;
- mapped-block fingerprint;
- whether the mapped-block fingerprint exactly matches the packet fingerprint;
- maximum absolute element difference;
- matrix infinity-norm difference;
- `overlap_atol`.

Do not copy these into flat stage, artifact, or public-result fields. Packet
self-integrity remains a hard prerequisite even when the mapped block is
numerically equivalent.

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

## Ordinary Radial Potential Fit And Consistency

Packet construction uses one untuned sequence:

```text
atomic determinant
-> spherical density fit
-> ordinary radial Hartree-potential fit
-> reported approximation errors
```

The density fit remains the authority for the reference cloud and no-half
Coulomb self-energy `E0_fit`. The potential fit is only a faster approximation
to that cloud's Hartree field `J0_fit`. It starts from the role-qualified
compact 45-term Coulomb scaffold and currently retains 33 Gaussian potential
terms after the ordinary radial trim/refit. A Hamiltonian or packet RHF built
with the separate high-accuracy 135-term policy does not make this fitted
potential exact.

For a represented packet determinant, report

```text
potential_fit_consistency_error = Tr(P0 * J0_fit) - E0_fit
```

along with radial, tail, matrix, and charge/self-energy fit errors. This value
measures consistency between determinant-defined `P0`, density-fit-defined
`E0_fit`, and approximate potential-fit `J0_fit`. It is not required to be
below `1e-8 Ha`, and it must not be hidden by a scalar correction or by fitting
the potential to molecular separations.

For additive packet references, the consumer also reports available
self/cross contributions:

```text
epsilon_aa = Tr(P_a * J_a_fit) - E_aa
epsilon_ab = Tr(P_a * J_b_fit) + Tr(P_b * J_a_fit) - 2 E_ab
epsilon_total = sum_a epsilon_aa + sum_{a<b} epsilon_ab
```

where packet density fits define `E_aa` and `E_ab`. The total must agree with
`Tr(P0 * J0_fit) - E0_fit` up to numerical assembly error.

The retired determinant-moment polish, separation grid, moment weight,
fixed-tail moment solve, and element/molecule-specific coefficient adjustment
must not run. Packet readback must explicitly reject any packet containing
`potential_fit/moment_polish/*`; polished Be, Ne, and Cr packets require
regeneration through the ordinary pipeline. Do not add a compatibility adapter
or silently strip the retired provenance.

## Validation And Failure Behavior

Required packet checks are:

- explicit RHF convergence;
- write/readback roundtrip;
- occupied-metric orthogonality and density trace;
- density reconstructed from occupied coefficients agrees with the stored
  density;
- exact stored-packet overlap-fingerprint integrity;
- exact owner-local atom/basis, function-count, owner, placement, ordered-label,
  angular-power, and column-order mapping;
- mapped owner-local overlap infinity-norm error at most the unchanged
  `overlap_atol = 1e-10`; a numerically equivalent reconstructed block may have
  a different raw fingerprint;
- density-fit charge and self-energy errors are reported and within the
  intended fit tolerance;
- potential-fit radial, far-tail, matrix, and energy-consistency diagnostics
  are reported, including its compact 45-term scaffold and current 33-term
  retained shape;
- potential-fit `J0_G` is compared with the density-fit mixed-Hartree path on a
  bounded matrix case, but its reported energy-consistency error is not an
  acceptance threshold;
- exact/density-fit oracle paths retain strict energy and field identities;
- packets containing retired moment-polish provenance are rejected;
- the density-fit path demonstrably receives the explicit compact expansion;
- packet outputs contain no unreported fallback reference or inferred
  occupancy convention.

Fail rather than fitting, writing, or consuming when convergence, packet
self-integrity, exact owner-local mapping, numerical overlap equivalence,
occupied orthogonality, density trace, finiteness, symmetry, or required
ordinary-fit construction validation is not satisfied. A finite validated
potential fit is not rejected solely because its reported energy-consistency
error exceeds `1e-8 Ha`. Do not repair packet mismatch by relabeling orbitals,
changing occupations, treating fit terms as orbitals, silently choosing
another Coulomb expansion, or reviving determinant-moment polishing.

## Ownership And Tests

Implementation and test ownership is recorded in the compact registry entry
for `HP-PQS-ATOMREF-PACKET-*`. The bounded tests cover packet roundtrip,
validation, unconverged-reference rejection, exact and numerically equivalent
owner-local embedding, structural mapping failures, Be/Ne consumption,
explicit compact expansion passage, and the vendored basis-data regression.

For the embedding-equivalence follow-on, source edits are limited to
`src/cartesian_reference_density/atomic_hf_reference_packets.jl`. The existing
private additive-reference caller may forward the one nested mapping summary
only if directly needed. No other packet, correction, artifact, or workflow
surface is reopened by this amendment.

The moment-polish retirement follow-on must delete
`_POTENTIAL_MOMENT_DISTANCES`, `_determinant_potential_moments`,
`_polish_atomic_reference_potential`, its packet fields, writer/readback
support, and moment-specific tests from the existing packet source/test owner.
`build_atomic_hf_reference_packet(...)` must consume
`fit_atomic_reference_potential(...)` directly, and readback must reject the
retired keys. The matching correction/additive diagnostic changes remain owned
by their canonical contracts; no replacement fit helper or compatibility
shape is approved. The source/test pass should be materially line-negative.

The packet facility depends on the historical basis collection described in
[Legacy BasisSets provenance](../../legacy_basissets_provenance.md).
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
- a replacement potential-polish policy, element-specific tuning,
  separation-trained fit, or scalar anchor patch.
