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

## Lowdin

Use symmetric Lowdin through the matrix inverse square root:

```julia
transform = inv(sqrt(Symmetric(overlap)))
cleaned = coefficients * transform
```

Do not hand-roll one-sided eigendecomposition transforms for Lowdin.

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
