# Invariants

These guardrails apply to the Cartesian/PQS base Hamiltonian producer.

## Architecture

- Geometry-specific code may produce terminal regions differently for
  one-center and bond-aligned diatomic systems.
- After terminal support, retained, and transform records exist, terminal basis
  realization, one-body assembly, IDA assembly, and
  `CartesianIDAHamiltonian` construction are shared.
- Terminal basis realization must not dispatch on atom count, system
  classification, route kind, bond axis, or terminal role names. It may dispatch
  on terminal lowering/transform kind.
- Do not create recursive stage embedding, route-status payload chains, flat
  field clouds, or report mirrors for construction data.
- Do not carry numerical matrices, factor tensors, terminal bases, or raw pair
  data through metadata.

## Pair And Box Terminology

Use these terms consistently:

| Term | Meaning | Production judgment |
| --- | --- | --- |
| Full parent 3D space | Entire Cartesian parent grid/basis. | Avoid dense 3D work. |
| Local product box / CPB | Full rectangular local box chosen to expose product structure. | Preferred intermediate. |
| Axis pair-factor terms | One-dimensional tensors such as `raw_x[k,i,j]`. | Preferred fast data. |
| Local product-box pair action | Pair contraction through 1D factors on local boxes. | Preferred fast kernel. |
| Terminal block pair | Pair of realized final-basis blocks used by K/U/V assembly. | Preferred final-basis operation. |
| Route pair inventory | Metadata enumeration of unit pairs and helper labels. | Planning only; not numerical authority. |
| Unstructured shell 3D contraction | Generic shell coefficient contraction without 1D simplification. | Avoid as routine production. |
| Global parent pair matrix | Dense pair object on the entire parent. | Forbidden production shape. |

"Full box" often means the fast local product-box route. "Full parent" is the
dangerous dense shape.

## Basis

- Direct core/slab/boundary blocks use implicit identity maps; do not allocate
  dense identity matrices.
- Parent gausslet rows are orthonormal to machine precision. Terminal regions
  own disjoint parent rows, so block-local terminal basis supports are
  structurally orthogonal across blocks.
- Terminal basis blocks are represented on disjoint owned terminal regions.
  `support_indices` and `support_states` are authoritative owned rows, not
  post-projection enlarged supports.
- PQS shell realization uses the full source box only to generate boundary
  product-mode columns. Those columns are restricted to the shell-owned support
  rows before shell-local Gram construction and symmetric Lowdin.
- PQS shell blocks must not project against previous terminal blocks, recurse
  through previous blocks, or grow support onto previous terminal regions.
- No global Lowdin is permitted.
- No global support matrix or global dense coefficient matrix is permitted as
  production authority.
- Final-basis self-overlaps are construction checks only, not downstream working
  data.
- Cross-block overlap is zero by construction because owned parent-row supports
  are disjoint. It is not a physical residual to compute, minimize, or repair.
  A nonzero structural overlap means duplicated support rows, incorrect row
  restriction, wrong support ownership, or an indexing error.
- Block-sparse terminal-basis representation does not imply block-diagonal
  operators. Cross-block kinetic, nuclear-attraction, and localized IDA matrix
  elements may be nonzero and remain assembled over terminal block pairs.

## Parent-Axis Numerical Ownership

- Once the Cartesian parent lattice and axis bundle are realized, that
  construction is the sole authority for reusable parent-only one-dimensional
  numerical data: overlap, kinetic, coordinate, second moment, integral
  weights, Gaussian factor terms, raw pair-factor terms, and exponent ordering.
- Parent-supplement cross tables are not parent data. They additionally depend
  on the validated supplement, Gaussian expansion, and physical centers, so
  they are construction-local augmentation work data derived from the
  authoritative parent-axis source.
- Downstream code may materialize rectangular parent-by-supplement cross
  matrices and project them through terminal blocks. It must not rebuild the
  one-dimensional cross tables independently per terminal block or per
  operator in production.
- Parent-supplement cross tables are numerical work data, not metadata, report
  fields, route-stage fields, artifacts, public API, or global mutable caches.
- A shared upstream source owner for supplement-dependent cross tables requires
  a separate docs-only amendment when a second production consumer or measured
  parent-side source cost justifies it.

## R3 Residual Selection

- Residual-GTO candidates are partitioned by physical owner center before
  residual-content selection.
- For owner `a`, the residual Gram
  `M_a = S_AaAa - X_a' X_a` is interpreted as the residual occupation spectrum
  of unit-occupied projected owner candidates. Its eigenvalues control
  retention; they are not merely numerical-rank diagnostics.
- The residual occupation cutoff `eta_RG = 1.0e-8` is separate from numerical
  negative-eigenvalue and stabilization tolerances. The negative-eigenvalue
  tolerances remain `tau_neg_abs = 1.0e-12` and `tau_neg_rel = 1.0e-12`.
- Owner-local modes below `eta_RG` are discarded. Eigenvalue flooring must not
  be used to retain nearly absent residual modes.
- Retained owner-local residual sectors are orthonormalized locally, then
  concatenated and merged by one final symmetric Lowdin over the inter-owner
  overlap matrix.
- Final merge thresholds are `tau_merge_abs = 1.0e-12` and
  `tau_merge_rel = 1.0e-12`. Any merge eigenvalue below
  `-max(tau_merge_abs, tau_merge_rel * max(lambda_max(S_merge), 1.0))` is a
  construction error, and any merge eigenvalue at or below that positive
  threshold is a near-singular merge error. Final `G' S R` and
  `R' S R - I` errors must be below `1.0e-10`.
- Global raw-candidate symmetric Lowdin and global raw-column pivoted-Cholesky
  selection are not the R3 residual algorithm.
- MWG centers and widths are computed from the final merged residual functions.
  Residual integral weights may be near zero or sign-changing and are not
  base-PQS IDA weights.
- Width filtering changes the candidate span; it is not a conditioning repair
  for owner-local residual selection.

## Lowdin

Use symmetric Lowdin through the matrix inverse square root:

```julia
transform = inv(sqrt(Symmetric(overlap)))
cleaned = coefficients * transform
```

Do not hand-roll one-sided eigendecomposition transforms for Lowdin.
For the R3 final inter-owner merge, a stable symmetric eigensystem
implementation is acceptable when it implements the same symmetric inverse
square root. It must diagnose near-singular merge spectra rather than floor
eigenvalues to preserve low-occupation directions.

## Gauge And Weights

- Slice A sign-canonicalizes completed block columns so localized final
  integrals/IDA weights are positive.
- Gausslet final integrals and final IDA weights must be finite, positive, and
  not near zero.
- Residual-Gaussian near-zero weights are a different non-base/supplement lane;
  they are not part of this base producer.
- Slice C recomputes final localized IDA weights from the supplied
  sign-canonicalized basis and parent axis weights. Negative recomputed weights
  mean the supplied basis and weights are inconsistent; Slice C must not repeat
  sign flips.

## Operators

- Kinetic energy is assembled from final-basis block products:

  ```text
  K = Tx*S y*Sz + Sx*Ty*Sz + Sx*Sy*Tz
  ```

- Unit nuclear attraction matrices are uncharged:

  ```text
  U_A = -1/r_A
  H1 = K + sum_A Z_A * U_A
  ```

- Nuclear charges are not applied to `U_A` before constructing
  `CartesianIDAHamiltonian`.
- Gaussian-expanded nuclear attraction and localized IDA use term-first
  contraction with the Gaussian expansion index as the short inner reduction.
- At most one terminal support-pair workspace should be live. Tile or stream
  any simultaneous local workspace above 64 MiB.

## Materialization

- PQS materialization receives `terminal_basis_realization` directly.
- Requested base PQS materialization returns `CartesianIDAHamiltonian{Float64}`
  directly.
- Missing terminal basis for requested PQS Hamiltonian construction is an error,
  not a blocked-result wrapper.
- Artifact writing uses the existing `write_cartesian_ida_hamiltonian` shape.
- R1 public facade artifact writes may add only the `HP-R1-ART-01`
  `producer_provenance/` keys in the final Hamiltonian file. Those keys are
  consumer provenance, not staged algorithm inputs.
- R3 supplemented Hamiltonian artifact writes may add only the `HP-R3-ART-01`
  `supplement_provenance/` keys in the final Hamiltonian file. Those keys are
  consumer provenance, not staged algorithm inputs, and do not change the
  `CartesianIDAHamiltonian` matrix artifact shape.
