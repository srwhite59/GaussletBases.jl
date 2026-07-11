# R3 Homonuclear Z-Axis Diatomic Supplemented Workflow

Status: implemented narrow molecular-scope contract for the supplemented
facade and canonical driver. This is generic homonuclear z-axis authority, not
element-specific Cr2 authority or general molecular support.

## Owned IDs

- `HP-R3U-ZDI-FN-01` - implemented producer input scope;
- `HP-R3U-ZDI-WIRE-01` - implemented canonical-driver wiring;
- `HP-R3U-ZDI-TEST-01` - completed H2/Be2 validation contract.

Source implementation commits are `c57e709e7` for the producer scope and
`3a4933812` for the original compact driver wiring. Later driver staging kept
the same system contract while exposing construction timings.

## Supported Systems

The molecular supplemented path supports:

- exactly two centers;
- equal atom symbols and equal finite positive nuclear charges;
- centers on the Cartesian z axis with distinct finite z coordinates;
- explicit nonnegative integer `nup` and `ndn`;
- neutral all-electron count,
  `nup + ndn == round(Int, sum(nuclear_charges))`;
- explicit base and supplement inputs;
- no element-specific defaults or element-specific source branches.

Unsupported inputs fail before expensive construction where practical:

- heteronuclear systems;
- non-z-axis or general translated/rotated molecular geometries;
- charged or ECP systems;
- solver/RHF workflow;
- Cr2-specific defaults, fixtures, or branches.

One-center supplemented atoms are governed by `HP-COMP-SUPPATOM-*`, not this
molecular-scope contract.

## Producer Boundary

`src/cartesian_base_hamiltonian.jl` validates the system through
`_cartesian_r3_diatomic_inputs(...)` and the supplement through
`_cartesian_r3_supplement_inputs(...)`. Required physical inputs remain
explicit: atom symbols, charges, positions, electron counts, near-core
spacing, longitudinal/transverse extents, and supplement basis labels.

The shared base producer owns `ns` normalization, route-local `q`, nesting,
source span, mapping strength, Coulomb accuracy, and other already-approved
basis controls. Supplement labels must match the two centers and remain equal
for the homonuclear named-basis loader. An optional trusted `basisfile` is
supported; no element basis default or pseudopotential behavior is implied.

The numerical composition is the implemented non-exported facade contract in
[R3 usability supplemented workflow](r3_usability_supplemented_workflow.md).
It remains same-construction and has no Cr2-specific path.

## Driver Boundary

`bin/cartesian_ham_builder.jl` supports supplemented mode through the same
producer stages. The current driver may call the approved visible staged
functions rather than one opaque facade call so it can report stage timings;
this does not create a second numerical workflow.

Driver inputs may provide the explicit homonuclear z-axis system, base basis,
supplement labels, optional trusted `basisfile`, and output artifact path. The
driver must not add route diagnostics, report/status payloads, private
Hamiltonian variants, new artifact schemas, public exports, solver behavior,
or Cr2-specific workflow.

## Validation

Durable validation comprises:

- the standalone H2 supplemented facade/artifact gate in
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`;
- canonical-driver H2 and Be2 supplemented artifact/readback checks from the
  accepted implementation passes;
- malformed homonuclear, orientation, charge, and supplement-input rejection;
- optional ignored/user-run Cr2 stress only after H2/Be2 pass.

Cr2 stress is evidence that the generic path can be exercised. It is not a
committed fixture, committed test, special workflow, or production claim.

Artifact behavior is delegated to
[Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md),
and residual numerical behavior is delegated to
[Residual Gaussian domain module](residual_gaussian_domain_module.md).
