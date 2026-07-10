# Screened Hartree Residual-Density Formalism

Document role: canonical durable physics contract for the screened-Hartree
residual-density construction. The completed measurement lanes are summarized
at the end of this page; they are historical evidence, not active source
authority. The implemented API is owned by
[Screened Hartree correction assembly](screened_hartree_correction_assembly.md),
and reusable reference objects are owned by
[Atomic HF reference packets](atomic_hf_reference_packets.md).

The point nucleus remains an exact/Galerkin one-body operator. IDA/MWG is used
only for the residual fluctuation density.

## Operator Identity

Start from the exact density decomposition

```text
rho_hat(r) = rho0(r) + delta_rho_hat(r)
delta_rho_hat(r) = rho_hat(r) - rho0(r)
```

Then the Hartree Coulomb term satisfies the exact identity

```text
1/2 (rho_hat | rho_hat)
  = (rho_hat | rho0)
  + 1/2 (delta_rho_hat | delta_rho_hat)
  - 1/2 (rho0 | rho0).
```

This identity is true for any chosen `rho0`. The approximation is introduced
only afterward:

```text
treat the rho0 one-body field and rho0 self-energy accurately/Galerkin;
treat only the smaller residual fluctuation term with IDA/MWG.
```

This is the physical target. It is not a row-gauge rho0 diagnostic and not a
correction to bare nuclear attraction alone.

## Reference Density

Choose a pure-GTO atomic closed-shell reference determinant once. For occupied
orbitals

```text
phi_m(r) = sum_a C_GTO[a,m] g_a(r)
```

with occupations `n_m`, define

```text
rho0(r) = sum_m n_m |phi_m(r)|^2.
```

The occupied/reference GTO directions that define this determinant are
protected basis content. They must be represented exactly in the final
gausslet+supplement working span:

```text
phi_m(r) = sum_i C0[i,m] chi_i(r).
```

Then the final-basis reference density matrix and IDA density vector are

```text
P0[i,j] = sum_m n_m C0[i,m] C0[j,m]
q0[i]   = P0[i,i].
```

This is not a projection approximation when the reference directions are
protected exactly. The pure-GTO and gausslet+supplement descriptions are two
coordinate representations of the same real-space determinant.

## Galerkin Reference Field

The nuclear attraction stays exact/Galerkin:

```text
Vnuc_G[i,j] =
    sum_A integral chi_i(r) * (-Z_A / |r - R_A|) * chi_j(r) dr.
```

The exact/Galerkin Hartree potential of the reference density is

```text
v0_G(r) = integral rho0(r') / |r - r'| dr'
```

and its one-body matrix is

```text
J0_G[i,j] = integral chi_i(r) * v0_G(r) * chi_j(r) dr.
```

The reference self-energy is

```text
E0_G = double_integral rho0(r) rho0(r') / |r - r'| dr dr'.
```

Equivalently, if `P0` is represented in the final basis and exact four-center
Galerkin integrals are available,

```text
J0_G[i,j] = sum_k sum_l P0[k,l] (ij|kl)_G
E0_G      = sum_i sum_j P0[i,j] J0_G[i,j].
```

For production-shaped work, `J0_G` and `E0_G` may be evaluated from the pure-GTO
reference density or from a compact fitted atom-centered Gaussian density, but
the fitted density must represent the same `rho0` and its fit error must be
reported or negligible.

## Screened Hartree-Only Functional

Let `P` be the spin-summed density matrix in the final orthonormal working
basis, and let

```text
q(P)[i] = P[i,i].
```

Let `V_IDA` be the current two-index IDA/MWG direct interaction matrix in that
same final basis convention. The screened Hartree-only functional is

```text
E_screen[P] =
    sum_i sum_j P[i,j] * (T[i,j] + Vnuc_G[i,j] + J0_G[i,j])
  + 1/2 sum_i sum_k (q(P)[i] - q0[i]) * V_IDA[i,k] *
                    (q(P)[k] - q0[k])
  - 1/2 E0_G.
```

At `P = P0`, the residual IDA term vanishes and the energy is

```text
sum_i sum_j P0[i,j] * (T[i,j] + Vnuc_G[i,j]) + 1/2 E0_G,
```

the exact Hartree energy of the reference density in the Galerkin one-body
field.

## Correction Relative To Current Direct IDA

Relative to the current direct IDA Hartree model

```text
E_current[P] =
    sum_i sum_j P[i,j] * (T[i,j] + Vnuc_G[i,j])
  + 1/2 q(P)' * V_IDA * q(P),
```

the screened functional adds

```text
Delta_J0[i,j] = J0_G[i,j] - delta_ij * sum_k V_IDA[i,k] q0[k]
```

and the constant

```text
C = 1/2 q0' * V_IDA * q0 - 1/2 E0_G.
```

Thus

```text
E_screen[P] = E_current[P] + Tr(P * Delta_J0) + C.
```

Although `Delta_J0` is represented operationally as a one-particle matrix and
`C` as a scalar constant, both belong to the screened direct electron-electron
interaction in energy accounting. They are the correction produced by
rewriting the Hartree term around `rho0`; they are not a change in the
physical kinetic-plus-nuclear one-body Hamiltonian and not an arbitrary energy
offset. Energy comparisons and error decompositions must group
`Delta_J0 + C` with the Hartree/IDA interaction model.

The nuclear attraction remains Galerkin; no IDA nuclear external-potential
term appears in this formulation.

## Protected Reference Directions

If a GTO direction is needed to define the reference determinant or its
closed-shell density, it must be protected. It is not an optional residual
candidate and must not be discarded by:

- RG selection;
- broad-direction rejection;
- injection cleanup;
- compactness filtering;
- any later convenience truncation.

The final basis may still include the normal RG/MWG, protected-injection, and
localized machinery. The rule is only that the reference determinant remains
exactly representable, so `q0 = diag(P0)` is exact in the final working basis.

## Compressed Density And Potential Roles

The pure-GTO occupied determinant defines `rho0`, `P0`, and `q0`. A validated
atom-centered Gaussian density fit may compress that same `rho0`; the density
fit defines the compressed reference cloud and its `E0_G`. It is not a new
model density or protected orbital content.

A fitted Gaussian potential may accelerate the matrix evaluation of `J0_G`.
It represents the Hartree potential of the density fit only. It does not
redefine `rho0`, `P0`, `q0`, or `E0_G`, and its Gaussian terms are not
supplement or protected orbitals.

The durable convention is:

```text
pure-GTO RHF defines rho0;
protected GTO directions make the same determinant exactly representable;
q0 comes from that final-basis determinant;
the density fit compresses the same rho0 and defines E0_G;
the potential fit is only a validated fast evaluator for J0_G.
```

Do not mix `q0`, `J0_G`, and `E0_G` from different reference densities.
Density- and potential-fit errors must be reported and small enough for the
intended correction accuracy. Packet-specific fit construction, provenance,
and validation belong to
[Atomic HF reference packets](atomic_hf_reference_packets.md).

## Exchange Is Separate

This is a Hartree/direct screened residual-density correction. It is not exact
frozen-core HF.

A true exchange analogue would decompose the one-particle density matrix,

```text
P = P0 + Delta_P,
```

and expand the nonlocal HF exchange functional around `P0`. That is a
separate future branch. It must not be folded into this Hartree lane.

## Historical Measurement Evidence

The measurement IDs below record completed or superseded evidence lanes.
Their detailed numerical records remain in the append-only manager running
log; they do not authorize new probes, source changes, or workflow behavior.

- `HP-PQS-SCREEN-HARTREE-AUDIT-01` is the completed H/Be/Be2 static and
  endpoint measurement. Passes 317-319 established roundoff reference
  representation and anchor consistency, compact/core-local correction modes,
  and basis-dependent endpoint shifts. The evidence supports the construction's
  stability, not a universal endpoint improvement.
- `HP-PQS-SCREEN-HARTREE-NE-AUDIT-01` is historical and operationally
  superseded. Its direct occupied-pair density construction remains an oracle;
  the practical endpoint path was closed through the fitted-cloud and
  packet-driven measurements in Passes 320 and 326.
- `HP-PQS-SCREEN-HARTREE-NE-FITCLOUD-AUDIT-01` is a completed historical
  measurement. Passes 321-326 established that the fitted density is a
  near-exact compression of the determinant density, while `P0/q0` must still
  come from the represented occupied determinant.
- `HP-PQS-SCREEN-HARTREE-POTFIT-AUDIT-01` is a completed historical
  measurement. Passes 327-329 established the fitted potential as a fast
  `J0_G` evaluator only. Its durable packet shape and validation rules moved
  to [Atomic HF reference packets](atomic_hf_reference_packets.md) in
  Passes 330-331.

The implemented correction object is governed separately by
`HP-PQS-SCREEN-HARTREE-CORR-FN-01` and
`HP-PQS-SCREEN-HARTREE-CORR-TEST-01`; see
[Screened Hartree correction assembly](screened_hartree_correction_assembly.md).

## Reusable Atomic HF Reference Packets

The implemented packet facility is governed by
`HP-PQS-ATOMREF-PACKET-FN-01` and `HP-PQS-ATOMREF-PACKET-TEST-01`.
[Atomic HF reference packets](atomic_hf_reference_packets.md) owns packet
roles, contents, identity, Coulomb conventions, validation, and failure
behavior. This formalism consumes the packet's occupied determinant as the
reference density, its density fit as the reference cloud/self-energy, and
its potential fit only as a fast evaluator of that cloud's Hartree field.
