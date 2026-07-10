# Screened Hartree Correction Assembly

Status: approved narrow internal source/design authority under
`HP-PQS-SCREEN-HARTREE-CORR-FN-01` and
`HP-PQS-SCREEN-HARTREE-CORR-TEST-01`.

This lane promotes the screened-Hartree residual-density machinery from
ignored probes to a source-backed internal helper that consumers such as CR2
can call naturally. It does not approve public driver support, production
artifacts, solver workflow, or endpoint claims.

## Purpose

Build an internal correction object for the Hartree residual-density
decomposition:

```text
rho = rho0 + delta_rho
```

where the saved/protected HF determinant defines `P0` and `q0`, and the
atomic reference packet density/potential fields evaluate the Galerkin
reference Hartree field.

The returned object should contain:

```text
ScreenedHartreeCorrection:
    delta_one_body      # Delta_J0
    energy_constant     # C
    q0
    P0 diagnostics
    J0_G diagnostics
    E0_G diagnostics
    anchor checks
    packet/provenance summary
```

`Delta_J0 + C` is part of the screened direct electron-electron interaction in
energy accounting, even though it is represented operationally as a one-body
matrix plus scalar constant. It is not a change to the physical
kinetic-plus-nuclear Hamiltonian and not an arbitrary energy offset.

## Inputs

Approved inputs:

- final orthonormal working basis/operators;
- `V_IDA` in the same final basis and site/order convention;
- one or more `AtomicHFReferencePacket` objects placed on molecule centers;
- imported/protected occupied reference coefficients defining `P0` and `q0`.

For molecules, each one-center packet placement must be explicit. The atom,
charge, electron count, center, basis, `lmax`, and fill-shell/spin convention
come from packet/provenance data and caller placement facts, not element-table
inference.

The occupied packet orbitals define `P0/q0`. Fitted density and
fitted-potential terms only evaluate `J0_G` and `E0_G`; they are not
supplement orbitals and not protected basis content.

Every consumed packet's in-memory or readback RHF convergence flag must be
explicitly true. Reject an unconverged packet before using its occupied
coefficients, density fit, or potential fit. Validation metadata is not
permission to continue with an unconverged reference.

## Operation

Build `q0` from the represented reference determinant:

```text
q0 = diag(P0)
```

Build `J0_G` and `E0_G` from the same packet density. Use validated fitted
potential terms for fast `J0_G` where available and where their packet
diagnostics pass. Fall back to the exact density-fit Galerkin path only where
the source pass explicitly supports it and reports the cost.

Return:

```text
Delta_J0 = J0_G - Diagonal(V_IDA * q0)
C        = 0.5 * q0' * V_IDA * q0 - 0.5 * E0_G
```

The helper must report anchor/derivative checks and projection/capture
diagnostics. It must keep the correction in memory unless a later artifact
authority is approved.

## Protected Additive Packet Consumer

`HP-RG-PROTECT-ADDREF-FN-01` approves a narrow molecular consumer of this
same correction core. It supplies separately validated native-`L` occupied
blocks from placed packets and forms:

```text
P0 = sum_a P_a
q0 = diag(P0)
J0 = sum_a J_a
E0 = sum_a E_aa + 2*sum_{a<b} E_ab
```

The packet blocks must remain separate for `P0`; their orthonormalized union is
used only to protect the basis span. An internal additive entry point may
validate each block separately and then delegate to the existing correction
and anchor core. It must not concatenate the blocks and impose global
orthogonality, duplicate the correction formula, rotate `Vee_L`, or write a
corrected artifact. The detailed source and Be2 gates are in
`protected_additive_reference_correction.md`.

## Source Surface

Approved source surface:

- `src/cartesian_reference_density/CartesianReferenceDensity.jl`;
- `src/cartesian_reference_density/screened_hartree_correction.jl`;
- `src/GaussletBases.jl` only for include/qualified access wiring;
- optional narrow use of
  `src/cartesian_reference_density/atomic_hf_reference_packets.jl` for packet
  readback/validation helpers.

Do not place this object in the Hamiltonian artifact writer, solver workflow,
EGOI lane, or external-GTO import lane. Those facilities may supply inputs,
but they do not own this correction object.

## Diagnostics

Required diagnostics:

- packet identity/provenance and placement facts;
- explicit packet RHF convergence status;
- electron count and `q0` charge by packet and total;
- `P0` trace and final-basis representation/capture loss;
- `J0_G` finite/symmetry checks;
- `E0_G` diagnostics and packet self-energy consistency;
- `Delta_J0` finite/symmetry checks;
- direct anchor identity:

  ```text
  E_current_direct[P0] + Tr(P0 * Delta_J0) + C == E_exact_direct[P0]
  ```

- derivative/field anchor check:

  ```text
  F_current_direct[P0] + Delta_J0 == J0_G
  ```

- potential-fit-vs-density-fit agreement when both are evaluated;
- row/sector/locality summaries sufficient for due-diligence review.

## Tests

Approved test surface:

- `test/nested/cartesian_screened_hartree_correction_runtests.jl`

Tests are correctness-only. They should use small Be/Ne-style packets and
bounded constructions.

Required test coverage:

- packet consistency and fingerprint validation;
- rejection of an unconverged packet;
- `q0` charge and `P0` trace;
- finite/symmetric `Delta_J0`;
- direct energy/derivative anchor identities;
- potential-fit agreement with the exact density-fit `J0_G` path on a small
  case;
- rejection or clear failure on mismatched packet/working-basis facts.

No Be2/Cr2 energy assertions, SCF convergence gates, solver tests, or
production endpoint claims are approved.

## Explicit Exclusions

Forbidden in this lane:

- public driver defaults or polished public workflow;
- production artifact schema, writer, or reader changes;
- solver workflow;
- Cr2 production claims;
- exchange correction;
- EGOI changes;
- rho0/P0 row-gauge shortcuts;
- fitted density or fitted-potential terms as protected orbitals;
- Hamiltonian source transforms;
- `Vee` source transforms;
- `C' V C` interaction rotation.

## Decision Rule

If the correction can be assembled from represented packet determinants,
validated packet density/potential fields, and same-basis `V_IDA` with clean
anchor checks, proceed under this authority. If it requires artifact schema,
solver integration, source interaction transforms, exchange correction, or
row-gauge substitutions, stop and request a new design amendment.
