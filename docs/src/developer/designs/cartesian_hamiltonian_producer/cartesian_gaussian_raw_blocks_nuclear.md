# Cartesian Gaussian Raw Blocks - Nuclear Slice

Status: approved narrow source authority for a neutral Cartesian Gaussian
raw-block owner, nuclear slice only. This amendment does not approve source
work outside the listed files and caller rewiring surfaces.

## Decision

Exact uncharged by-center Cartesian Gaussian nuclear block construction is a
shared numerical operation. Residual Gaussian augmentation and Qiu-White route
code both need the same parent-supplement and supplement-supplement nuclear
blocks, and Cr2 profiling shows this construction is the dominant allocation
cost in the current exact-operator path.

The operation should no longer be owned by Residual Gaussian helper code or by
Qiu-White route code. It moves to a neutral raw-block owner:

```text
src/cartesian_gaussian_raw_blocks/
  CartesianGaussianRawBlocks.jl
  nuclear_blocks.jl
```

This owner is numerical kernel code, not a route stage, cache framework,
artifact workflow, or public API.

## Approved IDs

- `HP-CGRB-FILE-01` - neutral Cartesian Gaussian raw-block module files and
  root include plumbing.
- `HP-CGRB-FN-01` - exact uncharged by-center Gaussian nuclear `G-A` and `A-A`
  raw-block construction.
- `HP-CGRB-FN-02` - one-dimensional supplement axis-family reuse and
  term-first nuclear table organization inside the neutral raw-block owner.
- `HP-CGRB-WIRE-01` - behavior-preserving rewiring from Residual Gaussian and
  Qiu-White callers to the neutral nuclear kernel, with duplicate route-local
  loops deleted after parity.
- `HP-CGRB-TEST-01` - focused parity and endpoint validation for the extraction.
- `HP-CGAI-FN-01` - optional low-level Cartesian Gaussian axis helper. It is
  allowed only as support for `HP-CGRB-FN-02`, not as the accepted performance
  endpoint.

## Scope

`HP-CGRB-FN-01` approves only:

- exact parent-supplement nuclear attraction blocks, `G-A`;
- exact supplement-supplement nuclear attraction blocks, `A-A`;
- by-center uncharged nuclear matrices, `U_A = -1/r_A`;
- analytic one-dimensional nuclear factor construction;
- reuse across unique nuclear coordinates;
- upper-triangular `A-A` assembly and mirroring;
- function-local scratch and workspace reuse;
- term-first contraction over the Gaussian expansion.

The raw-block owner returns uncharged by-center matrices. Physical nuclear
charges are applied by Hamiltonian or one-body assembly consumers, not inside
this kernel.

Later `HP-PQS-COULOMB-ACCURACY-FN-01` authority narrowly permits
`_factor_axis_integral` to replace its cancellation-prone weighted-center
variance with the algebraically equivalent pairwise weighted-distance
identity. That amendment changes neither this owner's physical scope nor its
polynomial moment, normalization, or term-first contraction conventions.

## Not Approved

This nuclear-slice amendment does not approve:

- overlap, kinetic, coordinate moments, or second moments;
- pair factors or matched-width Gaussian interaction;
- terminal projection;
- residual Gaussian selection or augmented-operator transforms;
- Qiu-White route objects;
- parent construction;
- persistent caches or raw-block bundles;
- metadata, report, status, or payload fields;
- artifact schema changes;
- public API or exports;
- Cr2 facade support or Cr2 artifact workflow.

Non-nuclear overlap/kinetic/moment raw-block work, if needed, must use the
separate `HP-CGRB-NN-*` authority in
`cartesian_gaussian_raw_blocks_non_nuclear.md`. It must not be implemented
under `HP-CGRB-FN-02`.

## Source Surfaces

Approved owner files:

```text
src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

`src/GaussletBases.jl` may add only the internal include needed to load
`cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl`. Include-order
changes are allowed only when required to make existing Gaussian integral
dependencies available to the neutral module and to make the neutral module
available to approved consumers. No public export is approved.

Allowed caller rewiring surfaces:

- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`, only
  to replace the existing residual-GTO nuclear raw-block loop with the neutral
  kernel and delete the replaced helper once no live caller remains;
- `src/ordinary_qw_raw_blocks.jl`, only to replace duplicated Qiu-White
  Gaussian nuclear raw-block loops with the neutral kernel;
- `src/ordinary_qw_operator_assembly.jl`, only where a Qiu-White nuclear
  consumer must be adjusted to call the neutral kernel output.

If implementation discovers a required source edit outside these files, stop
and request a docs-only source-surface amendment before coding that edit.

## Implementation Sequence

1. Extract the existing nuclear `G-A` and `A-A` behavior into the neutral owner
   without changing numerical conventions.
2. Rewire Residual Gaussian and Qiu-White callers to consume the neutral
   nuclear kernel.
3. Delete duplicate route-local nuclear loops once parity is established and no
   live caller remains.
4. Optimize allocation inside the neutral owner only after extraction parity.

Optimization must remain behavior-preserving. Do not add persistent cache
objects, status fields, provider payloads, route objects, or artifact data to
justify reuse.

## Approved Family-Reuse Optimization

Cr2 q4 profiling after `47d9b2a3` shows that the remaining neutral nuclear
raw-block allocation is one-dimensional analytic axis-integral table work, not
Residual Gaussian logic, Qiu-White route logic, or wrapper/result allocation:

- full neutral nuclear raw blocks: about `13.880s / 44552.840 MiB`;
- `G-A` primitive axis tables only: about `3.188s / 10787.069 MiB`;
- `G-A` 3D assembly from cached factors: about `0.529s / 20.625 MiB`;
- streamed upper-triangle `A-A` fill: about `9.890s / 33071.848 MiB`;
- scalar one-dimensional integral calls: `51,539,760`;
- parity against the current neutral function: `G-A = 0.0`, `A-A = 0.0`.

The failed in-place-table WIP reduced Cr2 q4 nuclear raw-block allocation only
from about `44552.840 MiB` to about `43471.448 MiB`. Therefore
`HP-CGAI-FN-01` is superseded as a performance endpoint. It may remain an
optional small helper, but the approved optimization target is
`HP-CGRB-FN-02`.

`HP-CGRB-FN-02` approves reorganizing nuclear raw-block construction inside:

```text
src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl
```

The kernel must identify unique supplement one-dimensional axis families
independent of flattened 3D orbital labels, then reuse their factor tables
across every 3D orbital or orbital pair that references the same family.

Approved implementation concepts:

- identify unique supplement axis families from primitive exponents, axis
  center, Cartesian power, and primitive normalization prefactors;
- build integer maps
  `orbital_axis_family[orbital, axis] -> family_id`;
- build unique `G-A` table keys
  `(axis, supplement_axis_family, nuclear_axis_coordinate)`;
- build unique `A-A` table keys
  `(axis, canonical(left_family, right_family), nuclear_axis_coordinate)`;
- keep transpose/orientation flags where canonical `A-A` family-pair order is
  reversed;
- fill each required one-dimensional table at most once per Coulomb Gaussian
  term;
- reuse those tables across all flattened 3D orbitals and orbital pairs that
  reference the same axis families;
- preserve the coupled primitive-pair contraction
  `sum_pq c_p c_q Ix[p,q] Iy[p,q] Iz[p,q]`;
- use only function-local workspaces and integer lookup plans.

Independent contraction of x/y/z axis tables into separate scalar contractions
is forbidden. The correct contracted 3D matrix element keeps the coupled
primitive-pair form above and is generally not the product of independently
contracted one-dimensional shell integrals.

Optional low-level support:

- `src/cartesian_gaussian_axis_integrals.jl` may add a specialized
  nonallocating nuclear-factor scalar integral for the `:factor` term if
  needed;
- the same file may keep/add a tiny in-place table-fill helper if it supports
  the family-reuse kernel;
- neither helper is sufficient as the optimization endpoint.

The optional in-place table-fill helper shape remains:

```julia
_cartesian_gaussian_axis_integral_table!(
    destination,
    left_exponents,
    left_centers,
    left_powers,
    left_prefactors,
    right_exponents,
    right_centers,
    right_powers,
    right_prefactors,
    term;
    factor_exponent = 0.0,
    factor_center = 0.0,
)
```

Its required contract is:

- fill an already allocated destination matrix;
- allocate no result matrix;
- return the same numerical values as
  `_cartesian_gaussian_axis_integral_table(...)`;
- preserve existing scalar `_cartesian_gaussian_axis_integral(...)` behavior.

Approved source surfaces for this follow-on optimization:

- `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`, for the
  one-dimensional family inventory, integer maps, unique table inventories,
  term-first table fill, and coupled primitive-pair accumulation;
- `src/cartesian_gaussian_axis_integrals.jl`, only for the optional
  specialized nonallocating nuclear-factor scalar integral or tiny table-fill
  helper described above.

No other source file is approved. If implementation needs another source
owner, stop for a docs-only amendment.

Implementation sequence for the later source pass:

1. Add a count-only/structural check, in ignored code if useful, reporting
   unique axis-family counts, unique `G-A` and `A-A` table counts, and current
   versus unique table-fill/scalar-call counts.
2. Implement a function-local axis-family inventory and integer
   orbital-to-family maps.
3. Fill unique one-dimensional tables term-first and reuse them across all 3D
   orbital/orbital-pair accumulations.
4. Add a specialized nonallocating nuclear-factor scalar integral only if the
   family-reuse kernel still needs it for material allocation reduction.
5. Preserve H2, Be2, Qiu-White, and Cr2 parity.
6. Remeasure Cr2 q4 nuclear raw-block allocation.

Acceptance criteria for the source pass:

- H2 Residual Gaussian endpoint is unchanged;
- Be2 facade/readback remains unchanged;
- Qiu-White nuclear parity remains unchanged;
- Cr2 q4 raw nuclear `G-A`/`A-A` parity remains unchanged;
- Cr2 q4 nuclear raw-block allocation drops substantially from about
  `44552.840 MiB`;
- report unique family counts, unique table counts, scalar-call count
  reduction, time, and allocation;
- no persistent object or broadened raw-block framework is introduced;
- no independent axis contraction is introduced;
- if low-level helper work is included, exact table/scalar values match the old
  helper at roundoff;
- no metadata, status, cache, route object, payload, report, artifact, public
  API, or persistent raw-block bundle is added.

Rationale: the Cr2 q4 fixture has only a modest number of real
one-dimensional Gaussian families, but the flattened 3D orbital-pair loop
rebuilds equivalent axis tables many times. The current `A-A` path already
fills preallocated tables in place, so an in-place result-table helper alone
cannot solve the measured `A-A` allocation. The primary optimization must
recover the shell/axis structure, reuse one-dimensional family tables, and
avoid heap allocation in genuinely distinct scalar nuclear factor evaluations.

## Deprecated Optimization Target

`HP-CGAI-FN-01` is no longer the accepted performance endpoint for Cr2 nuclear
raw-block allocation. It remains allowed only as an optional helper surface:

- `src/cartesian_gaussian_axis_integrals.jl`, for a specialized nonallocating
  nuclear-factor scalar integral or tiny table-fill helper if needed by
  `HP-CGRB-FN-02`;
- `src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`, only as the consumer
  of that helper inside the neutral nuclear kernel.

Do not proceed under the assumption that eliminating result-matrix allocation
solves the Cr2 q4 nuclear profile.

## Validation

`HP-CGRB-TEST-01` approves the following validation only:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian endpoint unchanged, as an ignored performance/usability
  measurement if it is too slow for committed validation;
- Cr2 q4 exact nuclear `G-A` and `A-A` blocks match the current implementation
  at roundoff, as an ignored measurement only;
- one small Qiu-White nuclear parity fixture.

Approved committed parity gate, if no existing test can host it cleanly:

```text
test/nested/cartesian_gaussian_raw_blocks_nuclear_runtests.jl
```

That test file, if added, must stay standalone and must not be added to
`test/runtests.jl` without a later amendment. It may validate only the neutral
nuclear raw-block contract and small Qiu-White parity. It must not assert route
status fields, report mirrors, payload fields, Cr2 workflow, artifacts, public
API, or Residual Gaussian internals.

## Relation To Residual Gaussian

`CartesianResidualGaussians` still owns residual basis selection, exact
augmented operator transformation, moment-matched Gaussian descriptors, and
residual interaction assembly. It may consume neutral nuclear raw blocks as
inputs to `transform_augmented_operator(...)`, but it does not own raw analytic
Gaussian nuclear formula construction.

## Relation To Qiu-White

Qiu-White route code remains a consumer and oracle source for parity during the
extraction. Qiu-White route objects, status fields, and route-specific staging
must not move into `CartesianGaussianRawBlocks`.
