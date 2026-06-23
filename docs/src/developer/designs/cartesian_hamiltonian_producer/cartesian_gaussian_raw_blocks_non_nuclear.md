# Cartesian Gaussian Raw Blocks - Non-Nuclear Slice

Status: approved narrow source authority for non-nuclear Cartesian Gaussian
raw-block reuse. This extends the neutral `CartesianGaussianRawBlocks` owner
under new `HP-CGRB-NN-*` IDs. It must not be implemented under
`HP-CGRB-FN-02`, which remains nuclear-only family reuse.

## Decision

Post-`HP-CGRB-FN-02` Cr2 q4 profiling shows that uncharged nuclear raw blocks
are no longer the dominant exact-operator cost. The remaining measured costs
are non-nuclear mixed/self Gaussian raw-block construction and final-basis
`G-G` product matrices:

- full Cr2 q4 exact augmented-operator wrapper:
  `8.3579s / 19107.314 MiB`;
- neutral nuclear raw blocks after family reuse: `0.6364s / 19.828 MiB`;
- residual setup mixed overlap `X`: `3.0855s / 10990.106 MiB`;
- Qiu-White self moment `A-A`: `2.1954s / 8453.301 MiB`;
- `G-G` product matrices: `2.8728s / 7024.456 MiB`;
- Qiu-White cross moment `G-A`: `0.4495s / 2489.566 MiB`.

The approved next source authority is therefore a narrow non-nuclear raw-block
slice in the existing neutral owner:

```text
src/cartesian_gaussian_raw_blocks/
  non_nuclear_blocks.jl
```

This is not a new route, cache framework, public API, or Residual Gaussian
algorithm change. It is a shared numerical kernel owner for exact
parent-supplement and supplement-supplement non-nuclear Gaussian raw blocks.

`G-G` product-matrix allocation is a measured remaining cost but is explicitly
deferred by this amendment. Optimizing final-basis `G-G` product matrices
needs a separate source-surface decision.

## Approved IDs

- `HP-CGRB-NN-FILE-01` - non-nuclear raw-block file under the existing neutral
  Cartesian Gaussian raw-block module.
- `HP-CGRB-NN-FN-01` - exact non-nuclear Gaussian `G-A` and `A-A` raw-block
  construction.
- `HP-CGRB-NN-WIRE-01` - behavior-preserving rewiring of Residual Gaussian and
  Qiu-White consumers to the neutral non-nuclear kernel.
- `HP-CGRB-NN-TEST-01` - focused parity and endpoint validation for the
  extraction.

## Scope

`HP-CGRB-NN-FN-01` approves only exact Cartesian Gaussian raw blocks for:

- overlap;
- kinetic;
- coordinate moments `x`, `y`, `z`;
- second moments `x^2`, `y^2`, `z^2`.

Approved block families:

- parent-supplement `G-A`;
- supplement-supplement `A-A`.

Approved construction details:

- analytic one-dimensional table construction;
- unique supplement axis-family reuse;
- canonical `A-A` family-pair table keys and orientation handling;
- upper-triangular `A-A` assembly and mirroring;
- function-local scratch/workspace reuse;
- coupled product-axis contraction preserving the existing Qiu-White numerical
  convention;
- reuse of a once-built overlap `G-A` block for both residual setup mixed
  overlap `X = G' S A` and exact augmented-operator assembly when both are
  constructed in the same call.

The kernel may return a compact fixed-field internal result containing only
the approved raw matrices. It must not be a status object, report payload,
route stage, persistent cache, metadata carrier, or artifact data.

## Not Approved

This amendment does not approve:

- nuclear raw-block changes;
- final-basis `G-G` product-matrix optimization;
- terminal projection;
- residual Gaussian selection, residual orientation, or augmented-operator
  transform changes;
- Qiu-White semantic changes or Qiu-White route objects;
- pair factors or matched-width Gaussian interaction;
- parent construction or parent-stage fields;
- persistent caches or broad provider bundles;
- metadata, report, status, or payload fields;
- artifact schema changes;
- public API or exports;
- Cr2 facade support, full Cr2 Hamiltonian construction, or Cr2 artifact
  workflow.

## Source Surfaces

Approved owner file:

```text
src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl
```

Allowed module plumbing:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
```

Only the include needed to load `non_nuclear_blocks.jl` is approved there.
Root include changes in `src/GaussletBases.jl` are not expected because the
neutral module already exists; if implementation discovers a real root include
or order issue, stop and request an amended source surface.

Allowed caller rewiring surfaces:

- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, only
  to replace current non-nuclear Qiu-White donor calls and to share the
  once-built overlap block between residual setup and exact operators when
  local wiring permits;
- `src/ordinary_qw_raw_blocks.jl`, only to replace duplicated non-nuclear
  Gaussian raw-block loops with the neutral kernel;
- `src/ordinary_qw_operator_assembly.jl`, only where Qiu-White consumers need
  adjustment to consume the neutral output.

No source edit outside these surfaces is approved.

## Implementation Sequence

1. Extract current Qiu-White non-nuclear `G-A`/`A-A` behavior into the neutral
   owner without changing numerical conventions.
2. Rewire Residual Gaussian exact-operator construction and residual mixed
   overlap setup to consume the neutral output where they currently use the
   Qiu-White donor organization.
3. Rewire Qiu-White consumers to the neutral output.
4. Delete duplicate route-local non-nuclear loops once parity is established
   and no live caller remains.
5. Optimize allocation inside the neutral owner only after extraction parity.

Do not add a persistent raw-block bundle/cache object to avoid one duplicate
overlap build. Construction-local fixed-field return values are acceptable
only when consumed immediately by the approved callers.

## Validation

`HP-CGRB-NN-TEST-01` approves the following validation only:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged as ignored performance/usability
  validation if needed;
- Cr2 q4 non-nuclear `G-A` and `A-A` raw blocks for overlap, kinetic,
  coordinate moments, and second moments match the current implementation at
  roundoff, as ignored measurement only;
- residual setup mixed overlap `X` matches the current construction at
  roundoff;
- one small Qiu-White non-nuclear parity fixture.

Approved committed standalone parity file, if no existing test can host it
cleanly:

```text
test/nested/cartesian_gaussian_raw_blocks_non_nuclear_runtests.jl
```

That test file, if added, must stay standalone and must not be added to
`test/runtests.jl` without a later amendment. It may validate only the neutral
non-nuclear raw-block contract and small Qiu-White parity. It must not assert
route status fields, report mirrors, payload fields, Cr2 workflow, artifacts,
public API, or Residual Gaussian internals.

## Relation To Nuclear Raw Blocks

`HP-CGRB-FN-01` and `HP-CGRB-FN-02` remain nuclear-only. Non-nuclear work must
use the `HP-CGRB-NN-*` IDs above. The shared module may reuse local utility
patterns across files, but nuclear and non-nuclear behavior must remain
separately reviewable.

## Relation To Residual Gaussian

`CartesianResidualGaussians` still owns residual selection, exact augmented
operator transformation, moment-matched Gaussian descriptors, and residual
interaction assembly. It may consume neutral non-nuclear raw blocks as exact
operator inputs, but it does not own the raw analytic Gaussian non-nuclear
block construction.

## Relation To Qiu-White

Qiu-White route code remains a consumer and parity reference during extraction.
Qiu-White route objects, route status fields, and route-specific staging must
not move into `CartesianGaussianRawBlocks`.
