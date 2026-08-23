# PQS Residual-GTO Same-Construction Working Basis

## Status And Authority

`HP-REP-PQS-RG-WORKING-FN-01` and
`HP-REP-PQS-RG-WORKING-TEST-01` own the implemented narrow expert construction
and its completed bounded validation. Commit `346589a6d` is the accepted source
implementation; both records now grant maintenance only.

This contract closes one representation boundary. It does not change residual
selection, one-body or interaction physics, external-packet semantics, or the
opaque `CartesianIDAHamiltonian` contract.

## Consumer Need

Before this boundary was implemented,
`cartesian_residual_gto_mwg_hamiltonian(...)` constructed a PQS terminal basis
`G`, a Gaussian supplement `A`, retained residuals `R`, exact augmented
one-body matrices, and the MWG interaction in one pass, then returned only
`CartesianIDAHamiltonian`. The implemented system constructor retains the
compact same-construction data needed to overlap that final basis with later
GTO probes.

The external importer correctly requires

```text
C_B = S_BX * C_X,
S_BX = <B|X>,
```

for an orthonormal final basis `B` and explicit external GTO basis `X`.
Reconstructing a representation from an opaque Hamiltonian is impossible and
is not authorized.

## Public Boundary

The package exports exactly one expert root constructor:

```julia
cartesian_residual_gto_mwg_system(
    system::NamedTuple;
    basis::NamedTuple,
    supplement::NamedTuple,
)
```

The call returns one concrete in-memory result from the same construction as
the Hamiltonian. The concrete result type remains unexported. Its supported
consumer surface is:

```julia
result.hamiltonian
gto_overlap_matrix(result, probes; block_indices = nothing)
import_external_gto_orbitals(result, packet)
```

`result.hamiltonian` is the existing `CartesianIDAHamiltonian{Float64}`. The
result retains only the terminal realization, factorized parent-axis donor,
explicit supplement representation, and existing residual object needed by
the overlap action. It must not retain the full base stage, operator stages,
or a dense parent-by-terminal coefficient matrix.

The existing non-exported
`cartesian_residual_gto_mwg_hamiltonian(...; hamfile)` call keeps its exact
signature, direct Hamiltonian return, artifact behavior, and numerical output.
Both public and compatibility calls must share one private construction rather
than rebuild independently.

No `basis_representation(result)` method is required or approved. This object
supports exact GTO cross overlap and external-orbital import; it is not a new
general Cartesian representation family.

## Cross-Overlap Contract

For terminal basis `G`, explicit supplement `A`, and retained residuals

```text
R = G*T_G + A*T_A,
B = [G,R],
```

reuse the existing factorized parent/GTO and contracted-GTO overlap kernels to
form

```text
S_GX = <G|X>,
S_AX = <A|X>,
S_RX = T_G' * S_GX + T_A' * S_AX,
S_BX = [S_GX; S_RX].
```

Rows are returned in native `[G,R]` Hamiltonian order. Probe columns retain
the packet/probe AO order. Existing `block_indices` ordering and bounds
semantics remain unchanged.

The terminal projection must use the support-blockwise terminal realization
and factorized parent data. The unavoidable final-by-probe output and a
parent-by-probe analytic temporary are allowed. A dense parent-by-terminal
map, dense final self-overlap, generalized final metric, or assembled
`[I T_G; 0 T_A]` raw-to-final matrix is forbidden.

Before returning or applying an overlap, validate:

- Hamiltonian dimension equals `nG + nR`;
- the residual base dimension equals `nG`;
- supplement orbital count equals the residual candidate count;
- `T_G` and `T_A` have exact `(nG,nR)` and `(nA,nR)` shapes;
- the public constructor accepts no separately supplied terminal, supplement,
  transform, representation, or Hamiltonian component; all retained pieces are
  captured inside the one producer call;
- requested probes and all computed cross-overlap entries are finite.

The result stores no external packet, imported orbitals, density, screening
field, solver state, or consumer policy.

## Implementation Boundary

Approved source owners are limited to:

- `src/GaussletBases.jl` for the one export;
- `src/cartesian_base_hamiltonian.jl` for the compact result and shared
  same-construction assembly;
- `src/cartesian_gto_probes.jl` for the exact cross-overlap dispatch;
- `docs/src/reference/export.md` for the exported function's doc entry.

Reuse the existing terminal projection, supplement-overlap, residual, and
external-import owners. Add no source file, module, exported type, accessor,
adapter, provider hierarchy, metadata key, status symbol, or artifact field.

Source budget, including the export and docstring:

- preferred: at most 110 added source lines;
- hard stop: at most 150 added source lines.

These are stop-and-report bounds, not reasons to obscure validation or merge
unrelated responsibilities. If coherent implementation exceeds the hard
bound, make no source commit and return the exact missing reusable operation.

## Accepted Implementation Validation

The original implementation validation used:

- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` for exact
  same-construction, formula, import, artifact, and provenance evidence;
- `test/docs/runtests.jl` for the single export/reference parity entry.

That supplemented H2 gate caught the information-loss and transpose/order
defects that the older Hamiltonian tests could not detect:

1. Construct the result through the new public function and preserve exact
   Hamiltonian matrices, dimensions, due diligence, and direct-facade behavior.
2. Compare `gto_overlap_matrix(result, probes)` with a direct assembly of
   `S_GX`, `S_AX`, and `T_G'*S_GX + T_A'*S_AX`.
3. Verify native row order, `block_indices`, finiteness, and dimensions.
4. Build a bounded valid `ExternalGTOOrbitalPacket`; verify the unchanged
   importer returns exactly `S_BX*C_X` and its existing capture diagnostics.
5. Reject malformed probe identity, dimensions, nonfinite data, and any
   dimension/identity inconsistency at the private assembly boundary before
   import.
6. Inspect and report terminal due diligence for the endpoint.

Accepted implementation evidence from commit `346589a6d` is:

- the H2 basis has `487` terminal plus `18` residual functions, for native
  dimension `505`;
- the system and direct facades produce exactly equal Hamiltonian matrices,
  and artifact readback deltas remain zero;
- direct assembly and the implemented `S_GX`/`S_AX`/`S_RX` formula agree
  exactly, as do imported coefficients and `S_BX*C_X`;
- the imported occupied capture is within `1e-8` of unity;
- the complete existing H2 owner passed `464/464`, and the supplemented
  facade/import owner passed `92/92`;
- terminal due diligence remains `9x9x15` parent axes, approximately
  `+/-4.871`, `+/-4.871`, and `+/-6.110` bounds, `4`-bohr padding, retained
  rows `[275,114,98]`, and only the pre-existing shape, axis-display, and
  large-identity warnings.

## Public Cartesian CI Extraction

`HP-REP-PQS-RG-WORKING-CI-FN-01` and
`HP-REP-PQS-RG-WORKING-CI-TEST-01` own the implemented public-CI extraction.
Commit `50b5ee9a1` moved the supported contract out of the unwired,
mixed-purpose nested suite without changing production behavior.

Maintain exactly one narrow public owner:

```text
test/driver_public/cartesian_residual_gto_mwg_system_runtests.jl
```

Keep that owner included once in the existing `:cartesian` group in
`test/runtests.jl`. The existing Julia `1.10` Supported-floor row owns the
remote regression. Do not add a CI row, alter the PQS or screening paper
groups, or wire the complete R3A file or any other `test/nested` suite.

The new owner exercises only documented public behavior:

1. Construct the bounded H2 system with the accepted `q=5`, `0.5` core
   spacing, `6.0/4.0` extents, and contracted H/cc-pVTZ `lmax=1` supplement.
2. Treat the concrete result as opaque except for `result.hamiltonian`.
   Validate native dimension `487 + 18 = 505`, finite symmetric matrices,
   the `1/1` particle sector, nuclei, and the existing
   `0.4574161883692301` self-Coulomb fingerprint.
3. Build a normalized one-orbital external GTO probe entirely through public
   representation and packet constructors. Exercise full and indexed
   `gto_overlap_matrix`, either restricted or alpha/beta import, source
   orthogonality, capture, and stale packet fingerprint rejection. Do not call
   `_cartesian_supplement_cross_overlap`, inspect private result fields, or
   reconstruct the residual transform formula in this owner.
4. Reject a small representative set of malformed public `system`, `basis`,
   and `supplement` inputs. Do not reproduce every historical blocked-path
   assertion.

The old R3A owner must lose the public-facade assertions duplicated by this
gate. Preserve there only unique private exact-formula, assembly-boundary,
artifact/readback, provenance, and due-diligence evidence until that suite has
its separate maintenance-versus-quarantine review. The existing
`HP-REP-PQS-RG-WORKING-TEST-01` continues to own those private oracle and docs
checks; the new CI records own only the public extraction and runner wiring.

The accepted owner is `68` lines and the runner adds one line; `75` obsolete
public-facade lines were removed from the old owner, for net `-5` tracked test
lines. Future maintenance must not use this record to expand the fixture,
restore duplicate nested assertions, or add a helper. The original
preferred/hard owner limits were `100/140` lines and remain upper bounds, not
targets.

Accepted evidence from commit `50b5ee9a1` is:

- the focused public owner passed `21/21` and the complete `:cartesian` group
  passed `253/253`;
- the retained nested owners passed `464/464` plus `64/64` private checks;
- remote CI run `32613923516` passed Supported floor, PQS paper, and Screening
  paper without changing either paper gate;
- terminal due diligence remains `9x9x15` parent axes, `4`-bohr padding, `487`
  terminal functions, `18` residuals, and final dimension `505`.

With repository acceptance complete, the separate REQ-101 consumer may rerun
only its frozen `R=2.35`, `ns=5`, PQS early case. It must use the unchanged
external
determinant and supplement, the public packet/importer, source metric
`C_X' S_XX C_X`, Euclidean final capture, the existing `1e-8` capture bounds,
and worst occupied loss at most `1e-4`. This is external acceptance evidence,
not a committed C2 fixture or a repo-owned paper campaign.

## Exclusions And Failure Rule

Do not change `CartesianIDAHamiltonian`, residual selection or cutoff, exact
one-body assembly, MWG/IDA interaction, external packet/importer behavior,
screened-Hartree algebra, public basis inputs, artifacts, drivers, solvers, or
defaults. Add no PySCF dependency or parser, C2/REQ-101 branch, molecular
screening automation, generalized overlap, final-final transfer, Hamiltonian
transform, interaction transform, persistent sidecar, or dense representation.

If the exact cross overlap cannot be produced from the same construction with
the existing factorized and supplement kernels, or if direct-facade parity,
the H2 gate, metric identities, or resource boundary fails, make no source
commit. Report the exact missing seam; do not reconstruct from the Hamiltonian,
fall back to a private consumer map, or add partial scaffolding.
