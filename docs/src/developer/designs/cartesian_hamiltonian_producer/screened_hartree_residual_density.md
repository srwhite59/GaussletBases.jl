# Screened Hartree Residual-Density Formalism

Document role: canonical screened-Hartree residual-density formalism plus a
transitional history of its measurement lanes. The formalism/audit itself is
measurement-only under `HP-PQS-SCREEN-HARTREE-AUDIT-01`; the atomic packet
section later in this page records separately implemented source authority
under `HP-PQS-ATOMREF-PACKET-*` and will be extracted into its own subsystem
contract in a later topology pass. This page does not independently broaden
either authority.

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

## Ne Endpoint Measurement

`HP-PQS-SCREEN-HARTREE-NE-AUDIT-01` approves one narrow endpoint extension of
this measurement branch.

Allowed:

- ignored `tmp/work` probe only;
- durable `/Users/srw/dmrgtmp` output;
- Ne atom only;
- closed-shell RHF;
- cc-pV5Z supplement with `lmax = 1`;
- screen by all electrons: protect the pure-GTO all-electron Ne reference
  determinant `1s^2 2s^2 2p^6` in the final basis;
- use this screened Hartree residual-density formalism:
  - `Vnuc_G` remains Galerkin;
  - `J0_G`, `E0_G`, and `q0` come from the protected GTO determinant;
  - IDA/MWG acts only on `q - q0`;
- compare against radial-gausslet Ne reference
  `E_ref = -128.547098109 Ha`;
- run standard-scaled PQS points, at least `ns = 5` and `ns = 7` if feasible.

Standard core spacings for Ne, `Z = 10`, use
`core_spacing = 1.2 / (Z * (ns - 1))`:

```text
ns = 5: core_spacing = 0.030
ns = 7: core_spacing = 0.020
```

Required reporting:

- dimensions, residual counts, and candidate counts;
- protected determinant representation loss;
- `Tr(P0)`, `q0` charge, and per-orbital projection loss;
- anchor energy and derivative errors;
- uncorrected RHF energy/error;
- screened-Hartree RHF energy/error;
- screened minus uncorrected shift;
- `Delta_J0` eigenvalue range, diagonal range, and low-mode
  locality/sector makeup;
- whether `lmax = 1` residual directions are actually retained.

Forbidden:

- tracked source edits;
- artifacts or public workflow;
- solver or driver integration;
- Cr/Cr2;
- exchange correction;
- EGOI changes;
- rho0/P0 revival;
- mapping default or fitting-policy changes;
- treating Ne as a broad first-row endpoint claim before this bounded
  measurement is reviewed.

## Ne Fitted-Cloud Endpoint Variant

`HP-PQS-SCREEN-HARTREE-NE-FITCLOUD-AUDIT-01` approves a measurement-only
fitted-cloud variant of the Ne endpoint. The exact determinant-density path is
an oracle/validation path. The fitted compact Gaussian cloud is the intended
practical endpoint path.

The fitted cloud is not a tunable physical model and not a rough screening
approximation. It is a near-exact compressed representation of the same
pure-GTO all-electron Ne reference density, used to avoid the expensive oracle
construction from all occupied determinant pair densities.

Allowed:

- ignored `tmp/work` probe only;
- durable `/Users/srw/dmrgtmp` output;
- Ne atom only;
- closed-shell RHF endpoint;
- cc-pV5Z supplement with `lmax = 1`;
- standard-scaled `ns = 5` and `ns = 7` if feasible;
- fit the all-electron Ne reference density, total charge `10e`, to a compact
  atom-centered sum of spherical Gaussian density terms;
- protect/represent exactly the GTO/cloud generators used by the fit, so the
  fitted density is exactly reproducible in the gausslet+supplement working
  space;
- build `J0_G`, `E0_G`, and `q0` from the fitted cloud, not from all occupied
  determinant pair densities;
- compare endpoint RHF energies to `E_ref(Ne) = -128.547098109 Ha`.

Fit standard:

- the fit is a compression of the pure-GTO reference density, not a new model;
- increase the number/flexibility of atom-centered Gaussian density terms until
  the fit reaches about `1e-8` relative error in Coulomb-relevant diagnostics,
  or until the fit is clearly limited by singular/ill-conditioned linear
  algebra;
- endpoint energies may be interpreted only if fit error is well below the
  observed screened-Hartree energy shift and below the target mHa-scale endpoint
  discussion.

Required fit diagnostics:

- number of Gaussian density terms;
- fit rank and condition/singular spectrum;
- max/min Gaussian exponents or widths;
- total fitted charge error;
- radial/enclosed-charge error if available;
- density residual norm or sampled radial residual;
- Coulomb self-energy absolute/relative error versus the exact
  determinant-density oracle where feasible, at least on `ns = 5` or a small
  oracle check;
- `J0_G` matrix/action error versus the exact determinant-density oracle on the
  smallest feasible point;
- representation/protection loss of every cloud generator;
- whether the stop was accuracy-satisfied or singular-math-limited.

Endpoint diagnostics:

- uncorrected RHF energy/error;
- screened-fitted-cloud RHF energy/error;
- screened minus uncorrected shift;
- anchor energy/derivative errors for the fitted cloud;
- `Delta_J0` eigenvalue and diagonal range;
- low-mode locality/sector makeup;
- retained residual counts, including p-channel confirmation.

Forbidden:

- tracked source edits;
- artifacts or public workflow;
- solver or driver integration;
- Cr/Cr2;
- exchange correction;
- EGOI changes;
- rho0/P0 row-gauge shortcuts;
- unreported density fits;
- discarding any protected cloud direction;
- broad first-row claims.

Interpretation rule: if the fitted cloud differs materially from the exact
determinant density, report the difference as fit error before interpreting
energy shifts.

## Fitted-Potential Measurement Variant

`HP-PQS-SCREEN-HARTREE-POTFIT-AUDIT-01` approves a narrow
measurement/prototype amendment for accelerating screened-Hartree endpoint
probes. It does not change the physics object:

- the saved HF determinant still defines `P0` and `q0`;
- the near-exact spherical Gaussian density fit still defines the reference
  cloud and `E0_G`;
- the fitted-potential object is only a fast radial Gaussian representation of
  that same cloud's Hartree potential for building `J0_G`.

Allowed:

- ignored `tmp/work` probes only;
- durable `/Users/srw/dmrgtmp` outputs;
- atomic one-center reference packets first;
- Be, Ne, and Be2 consumption after the one-center packet gates pass;
- optional `potential_fit/*` groups in ignored atomic HF reference density-fit
  prototype packets;
- use of the fitted-potential packet only as a faster `J0_G` builder in
  screened-Hartree probes.

Construction convention:

1. Start from the saved atomic HF determinant and its near-exact spherical
   Gaussian density fit.
2. Evaluate the analytic radial Hartree potential of that density fit:

   ```text
   J0(r) = sum_i w_i * erf(sqrt(beta_i) * r) / r
   ```

   with the `r = 0` limit handled analytically.
3. Use the repo Coulomb Gaussian expansion, scaled by total cloud charge, as a
   fixed long-range tail scaffold.
4. Fit only the short/intermediate residual, or add/refit local terms, so the
   far `Q/r` tail is protected.
5. Fit on a radial grid out to a large radius, for example `100` bohr.
6. Drive the potential-fit error to about `1e-8` in Coulomb-relevant
   diagnostics unless the linear algebra becomes singular or ill-conditioned.

Required packet/probe diagnostics:

- charge and density-fit provenance;
- Coulomb expansion or tail scaffold identity/fingerprint;
- potential-fit term count, exponents, coefficients, and signed/positive
  status;
- radial absolute and relative errors, including near-core and tail bins;
- `r * J_fit(r) -> Q` tail check;
- self-energy or anchor mismatch;
- matrix-level `J0_G` comparison against the existing density-fit exact
  Galerkin path on at least one small case;
- endpoint sensitivity only if the fit error is far below the
  screened-Hartree shift.

Forbidden:

- tracked source edits;
- production artifact schema or reader changes;
- solver or public workflow;
- changes to the saved HF determinant convention;
- treating potential-fit Gaussians as supplement or protected orbitals;
- rho0/P0 row-gauge shortcuts;
- EGOI or exchange changes;
- Cr/Cr2 production claims.

Decision rule: if the fitted potential cannot make matrix/anchor errors
negligible compared with endpoint shifts, stop and report fit limitations. If
it passes, it may be used as the fast `J0_G` builder in the next Be2 or Ne
screened-Hartree measurement, while the density fit remains the reference
cloud.

## Reusable Atomic HF Reference Packets

`HP-PQS-ATOMREF-PACKET-FN-01` and `HP-PQS-ATOMREF-PACKET-TEST-01` govern the
implemented reusable one-center atomic HF reference packet facility. It
replaces ignored ad hoc packets with a standard atom/basis reference object for
screened-Hartree probes.

This packet is analogous to defining a standard atom/basis reference:

- atom and charge;
- electron count and explicit filled-shell/spin convention;
- basis name and `lmax`;
- pure-GTO HF occupied orbitals;
- near-exact spherical Gaussian density fit;
- fast fitted Hartree potential for `J0_G`.

Required conventions:

- HF occupied orbitals define `P0` and `q0`;
- packet construction requires converged RHF and must stop before fitting when
  RHF is unconverged;
- packet writing rejects an in-memory unconverged packet, validation/readback
  report convergence explicitly, and packet-driven screened-Hartree
  consumption rejects an unconverged packet;
- no `allow_unconverged` packet-builder option is approved;
- the density fit defines the reference cloud and self-energy;
- the potential fit is only a fast representation of that same cloud's Hartree
  potential;
- fitted density and fitted-potential Gaussian terms are not supplement
  orbitals and not protected basis content;
- filled-shell and occupancy choices are explicit and must not be inferred
  from element tables.

Approved initial packet contents:

- `system`: atom, charge, electron count, spin/fill-shell convention, basis
  name, `lmax`, and geometry center;
- `supplement_basis`: labels, angular powers, exponents, coefficients,
  overlap, and overlap fingerprint;
- `hf`: occupied coefficients in supplement coordinates, occupations, orbital
  energies, density matrix, and explicit RHF/UHF convention metadata;
- `density_fit`: spherical Gaussian density coefficients/exponents, fit
  diagnostics, charge, and self-Coulomb;
- `potential_fit`: fitted radial `J0(r)` coefficients/exponents, fixed
  broad-tail provenance, refit/trim metadata, and radial/matrix/anchor
  diagnostics;
- `provenance`: code version, input controls, fit tolerances, and reference
  energies where available.

Initial scope:

- one-center atom packets only;
- Be core `Z=4`, `2e`, and Ne all-electron `10e`;
- cc-pV5Z, `lmax = 1`;
- writer/readback plus validation helpers;
- small Be/Ne smoke consumption that builds `P0/q0` and `J0_G` from packet
  contents.

Implemented source surface:

- `src/cartesian_reference_density/CartesianReferenceDensity.jl`;
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl`;
- `src/GaussletBases.jl` only for include/qualified access wiring;
- `src/cartesian_reference_density/screened_hartree_correction.jl` only for
  rejecting an unconverged packet at the consumption boundary;
- optional narrow calls to existing exact Hartree helpers in
  `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` for validation,
  without changing their numerical contract.

The density-fit `J0_G` path must obtain the packet-local compact Coulomb
expansion by its documented role and pass that exact object explicitly to the
mixed-Hartree helper. A helper default must not choose this approximation.

Approved test surface:

- `test/nested/cartesian_atomic_hf_reference_packet_runtests.jl`, limited to
  packet roundtrip, validation helpers, and small Be/Ne consumption smokes with
  temporary packet output and no committed binary fixtures;
- `test/nested/cartesian_screened_hartree_correction_runtests.jl`, only for
  unconverged-packet consumer rejection;
- `test/misc/runtests.jl`, only for the vendored-basis hash/parser regression
  below.

Validation requirements:

- packet readback roundtrip;
- unconverged RHF builder failure, in-memory packet write rejection, explicit
  convergence reporting, and packet-consumer rejection;
- occupied orthogonality and density trace;
- current supplement overlap/fingerprint match;
- density-fit charge and self-energy errors;
- potential-fit radial tail and matrix/anchor checks against exact
  density-fit `J0_G`;
- small Be/Ne smoke consuming packet to build `P0/q0` and `J0_G`;
- explicit compact Coulomb expansion passage for density-fit `J0_G`;
- no committed large binary fixtures.

Vendored-basis reconciliation under `HP-PQS-ATOMREF-PACKET-TEST-01` is limited
to correcting the header of `data/legacy/BasisSets`, correcting
`docs/legacy_basissets_provenance.md`, and adding a cheap regression in
`test/misc/runtests.jl`. The normalized scientific body must remain unchanged.
Its SHA-256 is
`b83883f4589234dd512eb760c95280854a2f42d007ab6e3533abda39a2829051`;
the full current file at `6cfe15d0c` is
`b40763c17cdbd8825cecf05aedbf34b35084fec7a3375661a65b49e749d33251`
and contains `60` basis blocks. The regression checks the normalized-body hash
and parser counts for H cc-pVTZ (`6/8` shells/rows), H cc-pVQZ (`10/12`), Be
cc-pV5Z (`21/42`), Ne cc-pV5Z (`21/54`), and Cr cc-pV5Z (`32/434`). Mixed
historical provenance and unresolved license/redistribution status must be
reported honestly; do not claim every custom block came from Basis Set
Exchange.

`HP-RG-PROTECT-ADDREF-FN-01` is the narrow molecular follow-on. It protects
the orthonormalized union of placed packet occupied spaces but forms `P0` from
the original per-packet blocks, sums placed fitted-potential fields, and adds
all packet density-fit cross energies before calling the existing native-`L`
screened-Hartree correction core. See
`protected_additive_reference_correction.md`. This does not broaden the packet
artifact itself or change fitted density/potential terms into orbitals.

Explicit exclusions:

- screened-Hartree production Hamiltonian assembly;
- artifact workflow integration beyond the packet itself;
- public driver defaults;
- solver workflow;
- EGOI, exchange, rho0/P0 row-gauge shortcuts, Cr/Cr2 production claims;
- treating fitted density or potential terms as protected GTOs.

Evidence: ignored packets and probes show this shape works. Ne packet-driven
endpoints removed hard-coded occupied labels, and refined potential fitting
reduced `J0_G` construction from about `118` s exact to about `0.96` s with
`33` terms and sub-microhartree anchor error.
