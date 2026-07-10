# Screened Hartree Correction Assembly

Status: implemented internal facility under
`HP-PQS-SCREEN-HARTREE-CORR-FN-01` and
`HP-PQS-SCREEN-HARTREE-CORR-TEST-01`.

This page is the canonical contract for the source-backed in-memory correction
API. The durable physics identity is owned by
[Screened Hartree residual-density formalism](screened_hartree_residual_density.md).
This facility does not provide public driver support, production artifacts,
solver workflow, or endpoint claims.

## Purpose

Build an internal correction object for the Hartree residual-density
decomposition:

```text
rho = rho0 + delta_rho
```

where the saved/protected HF determinant defines `P0` and `q0`, and the
atomic reference packet density/potential fields evaluate the Galerkin
reference Hartree field.

[Atomic HF reference packets](atomic_hf_reference_packets.md) owns packet
identity, determinant/density/potential roles, convergence, and validation.
This page owns only correction assembly from valid represented packet data.

The returned object contains:

```text
ScreenedHartreeCorrection:
    delta_one_body       # Delta_J0
    energy_constant      # C
    q0
    P0
    J0_G
    f_app_direct         # Diagonal(V_IDA * q0)
    Vq0
    E0_G
    q0Vq0
    energy_accounting
    diagnostics
    packet_summary
```

`Delta_J0 + C` is part of the screened direct electron-electron interaction in
energy accounting, even though it is represented operationally as a one-body
matrix plus scalar constant. It is not a change to the physical
kinetic-plus-nuclear Hamiltonian and not an arbitrary energy offset.

## Implemented API And Inputs

The core entry point is:

```text
build_screened_hartree_correction(
    V_IDA, J0_G, E0_G, represented_coefficients, occupations)
```

All matrices and represented coefficients must use one orthonormal final-basis
dimension and ordering. `V_IDA` is the direct IDA/MWG matrix in that basis;
`J0_G` is the Galerkin reference Hartree field in the same basis; and `E0_G`
is the no-half reference self-energy for the same density.

`build_atomic_packet_screened_hartree_correction(...)` is the one-packet
wrapper. It validates a converged packet, evaluates `J0_G` from either its
fitted potential or density fit, takes `E0_G` from the packet density fit, and
delegates to the core entry point. Packet placement is explicit when the
reference is not at its saved center.

`build_additive_screened_hartree_correction(...)` accepts separate represented
coefficient and occupation blocks, sums their density matrices, and preserves
per-block validation. It does not globally orthogonalize physically distinct
packet blocks. Molecular construction of placed fields, cross self-energy,
and native protected-localized coefficients belongs to
[Protected additive atomic reference correction](protected_additive_reference_correction.md).

The occupied packet orbitals define `P0/q0`. The density fit defines the
compressed reference cloud and `E0_G`; fitted-potential terms are only a fast
`J0_G` evaluator. Neither fit is supplement or protected basis content.

Every consumed packet's in-memory or readback RHF convergence flag must be
explicitly true. An unconverged packet is rejected before its occupied
coefficients or fitted fields are used, including in diagnostic-only mode.

## Operation

The core builds `q0` from the represented reference determinant:

```text
q0 = diag(P0)
```

`J0_G`, `E0_G`, and `q0` must refer to the same reference density. The packet
wrapper may use validated fitted-potential terms for fast `J0_G`; the density
fit remains the reference cloud and self-energy authority.

Return:

```text
Delta_J0 = J0_G - Diagonal(V_IDA * q0)
C        = 0.5 * q0' * V_IDA * q0 - 0.5 * E0_G
```

The helper must report anchor/derivative checks and represented-reference
diagnostics. It must keep the correction in memory unless a later artifact
authority is approved.

The energy interpretation depends on the `J0_G` source. For an exact or
density-fit oracle field, retain the strict identity

```text
Tr(P0 * J0_G) = E0_G
```

within the active numerical tolerance. For the ordinary fitted-potential path,
report instead

```text
potential_fit_consistency_error = Tr(P0 * J0_fit) - E0_fit
```

without rejecting an otherwise valid fit solely because this approximation
error exceeds `1e-8 Ha`. The density fit still owns `E0_fit`; the potential fit
does not redefine the reference energy. The derivative identity with respect
to the supplied `J0_G`, representation checks, and finite/symmetry checks
remain strict.

## Protected Additive Packet Consumer

`HP-RG-PROTECT-ADDREF-FN-01` implements a narrow molecular consumer of this
same correction core. It supplies separately validated native-`L` occupied
blocks from placed packets and forms:

```text
P0 = sum_a P_a
q0 = diag(P0)
J0 = sum_a J_a
E0 = sum_a E_aa + 2*sum_{a<b} E_ab
```

The packet blocks remain separate for `P0`; their orthonormalized union is
used only to protect the basis span. The internal additive entry point
validates each block separately and then delegates to the existing correction
and anchor core. It must not concatenate the blocks and impose global
orthogonality, duplicate the correction formula, rotate `Vee_L`, or write a
corrected artifact. The detailed source and Be2 gates are in
`protected_additive_reference_correction.md`.

## Source Surface

Implemented source surface:

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

Implemented diagnostics include:

- correction source and bounded packet identity/provenance summary;
- represented dimension, orbital count, `P0` trace/trace loss, occupation sum,
  and occupied orthogonality error;
- total `q0` charge and its minimum/maximum entries;
- input and resulting symmetry/finite checks for `V_IDA`, `J0_G`, and
  `Delta_J0`;
- `E0_G`, `q0' * V_IDA * q0`, current direct energy, correction expectation,
  and corrected direct energy;
- exact/density-fit direct anchor identity:

  ```text
  E_current_direct[P0] + Tr(P0 * Delta_J0) + C == E_exact_direct[P0]
  ```

- fitted-potential consistency error `Tr(P0 * J0_fit) - E0_fit`, reported as
  an approximation rather than an anchor rejection;
- derivative/field algebra check:

  ```text
  F_current_direct[P0] + Delta_J0 == J0_G
  ```

- optional fitted-potential-versus-density-fit matrix relative-Frobenius and
  maximum-entry differences when both are evaluated;
- for additive references, block count/traces and maximum inter-block occupied
  overlap, plus total and available self/cross fitted-potential consistency
  contributions.

Packet placement, protected-basis capture, row/sector locality, and terminal
due diligence belong to the packet embedding and protected additive-reference
consumer, not this correction object.

## Tests

Implemented test surface:

- `test/nested/cartesian_screened_hartree_correction_runtests.jl`

Tests are correctness-only and use small Be/Ne-style packets and bounded
constructions.

Required test coverage:

- packet consistency and fingerprint validation;
- rejection of an unconverged packet;
- `q0` charge and `P0` trace;
- finite/symmetric `Delta_J0`;
- strict derivative/algebra identities for every source;
- strict direct energy identity for exact/density-fit oracle sources;
- acceptance and reporting of a finite ordinary fitted-potential consistency
  error without a `1e-8 Ha` rejection;
- potential-fit agreement with the exact density-fit `J0_G` path on a small
  case as a reported matrix/energy approximation;
- additive blocks remain separate while their density matrices sum;
- additive total consistency error agrees with its self/cross decomposition;
- packets containing retired `potential_fit/moment_polish/*` provenance are
  rejected;
- rejection or clear failure on mismatched packet/working-basis facts.

No Be2/Cr2 energy assertions, SCF convergence gates, solver tests, or
production endpoint claims are approved.

## Failure Behavior

Dimension mismatch, nonfinite inputs, negative occupations, input asymmetry,
reference trace/orthogonality loss, negative `q0`, unconverged packet
consumption, and failed derivative/algebra checks are hard failures. Exact and
density-fit oracle energy identities also remain strict.

An ordinary fitted-potential result is not rejected solely because
`Tr(P0 * J0_fit) - E0_fit` exceeds `1e-8 Ha`; that value and the associated
radial, tail, and matrix errors must remain visible. `diagnostic_only = true`
does not bypass structural, representation, finiteness, symmetry, convergence,
or algebraic failures and is not required merely to observe a finite
fitted-potential consistency error.

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
