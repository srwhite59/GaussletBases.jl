# Screened Hartree Residual-Density Formalism

Status: approved measurement-only audit authority under
`HP-PQS-SCREEN-HARTREE-AUDIT-01`. This page records the crystallized
Hartree-only screened-density branch. It is not source authority.

This branch supersedes the all-IDA screened-nucleus route as the main
Hartree-only correction target. The point nucleus remains an exact/Galerkin
one-body operator. IDA/MWG is used only for the residual fluctuation density.

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

No `uN_IDA` appears in this formulation.

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

## Relation To Fitted Densities

A fit of the pure-GTO RHF density to 20-30 atom-centered `s` Gaussians may be
used as an evaluation device for `J0_G` and `E0_G`. It is not a separate model
density unless explicitly declared so.

The safe convention is:

```text
pure-GTO RHF defines rho0;
protected GTO directions make the same determinant exactly representable;
q0 comes from that final-basis determinant;
J0_G and E0_G come from the same rho0, directly or through a reported
near-exact fit.
```

Do not mix `q0` from one density with `J0_G` or `E0_G` from another.

## Exchange Is Separate

This is a Hartree/direct screened residual-density correction. It is not exact
frozen-core HF.

A true exchange analogue would decompose the one-particle density matrix,

```text
P = P0 + Delta_P,
```

and expand the nonlocal HF exchange functional around `P0`. That is a
separate future branch. It must not be folded into this Hartree lane.

## Approved Measurement Shape

`HP-PQS-SCREEN-HARTREE-AUDIT-01` approves a protected-GTO Hartree screened
residual-density audit.

Allowed:

- ignored `tmp/work` probes only;
- durable `/Users/srw/dmrgtmp` outputs;
- pure-GTO closed-shell atomic RHF reference, done once;
- H/Be/Be2 first, then Cr atom `[Ar]`-like core;
- explicitly protect the occupied/reference GTO directions in the final
  working span;
- construct `P0_final` and `q0 = diag(P0_final)`;
- construct `J0_G` and `E0_G` from the same pure-GTO reference density;
- form `Delta_J0` and `C`;
- test the direct Hartree energy/derivative identity;
- inspect low-mode behavior, locality, spectra, and orbital expectations.

Forbidden:

- tracked source edits;
- artifacts, schema, or public workflow;
- solver workflow;
- Cr2;
- production corrected Hamiltonian;
- dependence on `uN_IDA`;
- row-gauge rho0 shortcuts;
- discarding reference GTO directions;
- exchange correction;
- EGOI changes;
- Cr2 production claims.

Decision rule: continue only if `P0` is represented exactly or at roundoff,
`J0_G`, `E0_G`, and `q0` refer to the same `rho0`, `Delta_J0` is
core-local/moderate, and H/Be/Be2 low modes remain benign.

Stop if the protected reference determinant cannot be represented exactly,
`J0_G`, `E0_G`, and `q0` come from mismatched densities, or the correction
creates broad or valence-destabilizing modes.
