# Raw Product Source And Retained Transform Policy

This note records the intended next abstraction for high-order Cartesian
construction. It is a design contract, not a claim that all routes already use
this machinery.

## Core Idea

Most retained Cartesian units should be understood as:

```text
raw product source space -> retained-space transform -> retained unit
```

The raw product source space is a small local product block built from the
parent gausslet axes. The retained-space transform then performs the local
selection, projection, cleanup, or identity map needed for the construction.

Operator blocks should be assembled in the raw product source spaces first,
then finished by small dense transforms:

```text
O_retained = T_left' * O_raw_product * T_right
```

This is the preferred framing for ordinary product slabs, midpoint caps,
rectangular shells, and projected q-shells. It should reduce the number of
special mixed-block cases because each side presents the same shape: a raw
product source plus a retained transform.

## Examples

A midpoint slab `q x q x 1` is the simple case. Its raw source is the product
space

```text
X_q * Y_q * fixed_z
```

and its retained transform is usually identity or close to identity. This is a
special case of a product-source retained unit, not a separate Hamiltonian
idea.

A projected q-shell, `PQS(q, L)`, starts from the full local product block

```text
X_q * Y_q * Z_L
```

The shell is not defined by subtracting a contracted inner block. The clean
definition is:

1. build the full local product transform;
2. project it onto the raw boundary rows of the shell;
3. keep the boundary-projected product-mode span;
4. apply full-rank symmetric Lowdin cleanup.

For a rectangular shell this gives the counting picture
`q*q*L - (q-2)*(q-2)*(L-2)`, but that count is not the construction rule. The
construction rule is boundary projection of the full product block, followed by
cleanup. This avoids any dependence on how an interior cube would have been
contracted.

## Operator Consequence

Mixed blocks such as PQS/product should not be treated as strange special
cases. Both sides should expose:

- raw product source metadata;
- a retained transform;
- support and boundary ownership diagnostics;
- metric/operator capability labels.

Then overlap, position, kinetic, `x2`, and Gaussian local one-body blocks can
be formed in the raw product source spaces and transformed down. For
factorized two-electron or MWG-style terms, the same philosophy may apply, but
the contraction must go through the factor terms carefully. It must not assume
positive quadrature weights for functions that are not quadrature-carrying.

## Raw Source Pair Operator Packets

A real retained basis has several retained units: atom cores, midpoint slabs,
PQS shells, contact caps, mismatch shells, supplement blocks, and fallback
pieces. The final operator matrix is assembled from blocks between every pair
of retained units.

The missing intermediate object is therefore a pair packet for raw product
source spaces:

```text
source_i, source_j
-> RawProductSourcePairOperatorPacket3D
-> T_i' * O_raw(i,j) * T_j
-> retained operator block
```

Only the upper triangle needs to be built for symmetric operators. The packet
should be lazy or cacheable so large constructions do not materialize all raw
pairs too early.

A useful vocabulary is:

```text
RawProductSource3D
    parent axes
    local box
    axis transforms
    source dimension
    support map / boundary map
    raw source quadrature weights

RetainedTransform3D
    source id
    retained dimension
    transform T
    selection / projection / cleanup diagnostics
    quadrature weight role for retained columns

RawProductSourcePairOperatorPacket3D
    left source id
    right source id
    operator kind
    raw overlap / position / kinetic / gaussian factors / pair factors
    backend provenance
    supported term list
    symmetry status
```

This is the structure that makes a PQS/product mixed block ordinary. The PQS
side has a full `q x q x L` product source and a boundary-projection/Lowdin
retained transform. The product or midpoint-slab side has its own product
source and usually a simple retained transform. The operator block is built
between the two raw sources, then transformed on both sides.

## Weight Contract

Retained units need an explicit weight role.

The safest positive quadrature object is the raw product source. Raw
gausslet-derived product sources should have positive source-point weights in
non-pathological cases:

- product slabs;
- midpoint caps treated as product slabs;
- atom-local cubes and shells;
- rectangular product shells;
- projected q-shells.

For these raw source spaces, negative or zero source weights in normal geometry
are a construction problem, conditioning problem, or gauge problem.

Retained-column weights need a stricter label. If the retained transform is
identity-like, retained weights may remain positive. If the retained transform
includes projection and Lowdin cleanup, as in PQS, induced per-column weights
can become small or signed depending on the gauge and conditioning. The
intended PQS construction should induce positive retained weights in ordinary
cases, but code must not silently rely on that fact. IDA/MWG-style paths that
divide by retained weights should either operate at the raw source level or
require an explicit retained-column `:positive_required` role plus a finite
positive check.

Supplement or final-basis functions are different. Angular GTO supplements
such as p and d functions, residual directions, and other non-quadrature
additions do not carry positive quadrature weights. Their integral or weight
can be zero, signed, or not meaningful as a quadrature mass. Nothing in the
construction should divide by those weights, and IDA-style positive-weight
paths must be forbidden for those functions.

A useful payload field would be:

```text
quadrature_weight_role =
    :raw_source_positive
    :positive_required
    :not_quadrature_weight
    :debug_reference_only
```

The intended meaning is:

- `:raw_source_positive`: finite positive weights are attached to raw product
  source points and can be used by raw-source quadrature contractions.
- `:positive_required`: finite positive retained-column weights are part of
  the construction contract before IDA/MWG-like use.
- `:not_quadrature_weight`: weights are not valid quadrature masses; reject
  any path that divides by them.
- `:debug_reference_only`: the payload is a reference/debug object and should
  not be promoted to production quadrature behavior without a separate review.

## Fit With The Current Framework

The existing framework direction is:

```text
Cartesian parent
-> shell/region descriptor
-> contraction rule
-> resolved payload
-> metric/operator helper
```

This note refines what an execution-ready payload should contain. A useful
payload should not merely say "PQS" or "product". It should expose the raw
product source and retained transform that let operator code build the raw
block first and finish with small linear algebra.

For pairwise operator construction, the next useful object is a
`RawProductSourcePairOperatorPacket3D`-style packet. That packet should be the
place where raw product-source operator factors, backend provenance, symmetry,
and supported-term lists are recorded before retained transforms are applied.

Current private PQS sidecar fixtures are still low-order/reference-scoped.
They prove payload discovery and low-order checks, but they do not yet install
this full raw-product-source contract into production metric, QW, Hamiltonian,
or CR2 routes.

## Private Metadata Checkpoint

The current private implementation now has the first raw-source planning
scaffold:

- raw product source metadata;
- retained transform metadata;
- an upper-triangular raw source pair plan;
- a resolved raw source pair audit;
- a private raw overlap packet;
- a private toy raw `:axis_index_x` packet for the identity product/slab
  fixture;
- private physical raw `:position_x`, `:position_y`, and `:position_z` packets
  for product/slab fixtures, consuming explicit 1D axis metric data with
  non-integer positions;
- private retained low-order blocks for materialized product/slab transforms;
- a nonidentity materialized product/slab retained transform;
- product/product cross-pair physical position packets for same-support
  identity/nonidentity fixtures;
- product/product cross-pair physical position packets for distinct shifted
  supports with nontrivial explicit axis metrics;
- private factorized product/doside raw packets for `:overlap` and physical
  `:position_x`, `:position_y`, and `:position_z`, built from explicit 1D axis
  factors over staged local intervals;
- private retained-block comparisons against the existing product-staged metric
  block kernel for overlap and physical position terms.

The retained-block calculations are intentionally tiny product/slab fixtures.
The baseline commits for the raw-packet line include `8ba7ab4 Rename toy
product packet axis index term`, `d0fb995 Add private physical product position
packets`, `3ef0fa1 Generalize private physical product packets`, `9589f9a Test
nonidentity product retained position blocks`, `26f88d6 Add private product
cross position packets`, and `c6d4f41 Test product cross packets on distinct
supports`. The factorized packet checkpoint continues through `107e68c Add
private factorized product overlap packets`, `3c1eb95 Add private factorized
product position packets`, and `1c28f8f Compare factorized product packets to
staged metrics`.

```text
source dimension 4 -> retained dimension 4
retained transform kind = :product_axis_transform
T' * I_raw * T == I_4
T' * axis_index_x_raw * T == diag(1, 1, 2, 2)
T' * position_x_raw * T == diag(0.25, 0.25, 1.75, 1.75)
T' * position_y_raw * T == diag(-0.5, 0.5, -0.5, 0.5)
T' * position_z_raw * T == diag(3.25, 3.25, 3.25, 3.25)
```

The `:axis_index_x` packet is a toy axis-index diagnostic for the fixture. It
is not a physical `position_x` operator and should not be extended as if it
were one.

The private physical position packets are separate. They build raw matrices
over raw product source support rows from explicit 1D axis metric data and do
not build dense full-parent 3D matrices. Retained blocks are checked only
where the retained transforms are materialized. The current private checks
cover self-pair retained overlap/position, nonidentity retained transforms,
same-support product/product cross pairs, and distinct-support product/product
cross pairs. All of these remain fixture/reference checks, not production
metric execution.

The private factorized product/doside packets are the next validation layer.
They build materialized raw support-row matrices for these small fixtures, but
the entries are assembled from explicit 1D overlap/position factors restricted
to staged local intervals. Retained blocks are then checked with the same
`T_left' * O_raw * T_right` algebra and compared against
`_fill_product_staged_metric_blocks!(...)`, the existing product-staged
retained metric block kernel. This proves the low-order product/product
factorization for the private fixtures; it is still not a large-region fast
operator assembly path.

The follow-up retained-block checkpoint removes that raw support-row matrix
from the private product/doside comparison path. The helper
`CartesianContractedParentMetrics._product_doside_retained_low_order_block(...)`
directly returns retained product/product blocks from staged 1D factors for
`:overlap` and physical `:position_x`, `:position_y`, and `:position_z`. It
matches both `_fill_product_staged_metric_blocks!(...)` and the private
raw-packet retained path on the focused fixtures. Guard tests pin that
unsupported terms reject, non-`:product_doside` units reject, and axis metadata
or axis metric mismatches reject. This checkpoint is recorded by
`bde3dac Add private product retained block helper` and
`0470f27 Test product retained block helper guards`.

The product low-order adoption checkpoint then uses these private helpers only
inside the existing product/doside low-order metric plumbing. As of
`250b3d7 Adopt product retained block fill helper`,
`_fill_product_staged_metric_blocks!(...)` delegates product/doside
`:overlap` and physical `:position_x`, `:position_y`, and `:position_z`
retained block fills to `_product_doside_retained_low_order_block(...)`. As of
`57738aa Factor product linear vector helper`,
`_staged_unit_linear_vectors(...)` delegates product/doside retained weights
and first moments to `_product_doside_retained_linear_vectors(...)`.
Support-dense and generic fallback paths remain unchanged. Pair blocks and
one-unit linear vectors remain distinct concepts; weights and first moments
are not pair packets.

The helper seam is still private low-order infrastructure. Its adoption only
removes duplicate product/doside internals in the existing metric path. It does
not change public/default routes, QW/Hamiltonian builders, backend/default
policy, PGDG/quadrature policy, IDA or positive-weight assumptions, CR2 paths,
or science validation. It also does not add PQS/product mixed blocks,
support-dense mixed packet generalization, all-pairs production matrices,
kinetic, nuclear/local, Gaussian, MWG, or interaction terms. GTO supplement
functions remain outside retained-column positive quadrature assumptions.

The first private non-low-order operator checkpoint is kinetic. As of
`04b7095 Add private product kinetic block helper`,
`_product_doside_retained_separable_sum_block(left_unit, right_unit,
axis_factor_terms)` builds retained product/doside blocks from explicit
separable axis-factor triples, and
`_product_doside_retained_kinetic_block(left_unit, right_unit, axis_ops)` is a
thin wrapper for

```text
(Kx, Sy, Sz) + (Sx, Ky, Sz) + (Sx, Sy, Kz)
```

The focused tests compare product/doside self and cross-unit kinetic blocks
against support-local retained references formed from `_staged_unit_entries(...)`
and `_contract_pair_block(...)`. This is private signed operator-block
infrastructure only. Kinetic is not a weight object, not a positive quadrature
claim, and not an IDA division path.

The real-block shadow checkpoint is recorded by
`ecad56a Test product kinetic helper on real blocks`. The existing q4
bond-aligned endcap-panel shared-shell source-policy fixture exposes six
`:product_doside` units through its staged sidecar. The test compares
`_product_doside_retained_kinetic_block(...)` against the already-built real
route matrix block
`fixed_block.kinetic[left.column_range, right.column_range]` for all 21
upper-triangular product/product pairs. Cross-pair mirror blocks are checked
against the transposed kinetic block, and the observed maximum block error was
about `8.9e-16`. This remains private shadow validation only, not construction
adoption.

This is a contract checkpoint, not production metric execution. The PQS raw
overlap packet exists only as raw packet/reference plumbing. PQS retained
transforms are not applied, and the PQS Lowdin cleanup matrix is not treated
as the full raw-to-retained transform. PQS/product mixed retained blocks remain
unsupported. Support-dense mixed packets, production product/product
operator assembly, all-pairs matrix construction, kinetic adoption outside the
private product/doside helper, nuclear/local, Gaussian, MWG, interaction,
QW/Hamiltonian, public/default, backend/default, PGDG/quadrature, CR2, and
science/energy behavior remain unchanged and outside this private fixture line.

Before any PQS retained-block execution, the code needs an explicit way to
represent or resolve the full factored PQS transform:

```text
raw_product_modes
-> raw_boundary_projection
-> full_rank_symmetric_lowdin_cleanup
-> retained_columns
```

General physical raw source pair operator packets beyond these private
product/slab fixtures, kinetic adoption into existing real route consumers,
nuclear, Gaussian/local, interaction/MWG, and mixed product/PQS or
product/support-dense pairs need separate design. The current checkpoint
changes no all-pairs matrix construction, QW or Hamiltonian path,
public/default route, backend/default policy, PGDG or quadrature policy, CR2
path, or science status.

## Non-Goals

This policy does not change:

- public/default construction routes;
- QW or Hamiltonian kernels;
- backend defaults;
- PGDG or numerical-reference policy;
- supplement handling;
- CR2 or science validation status.

The next design choice is either to build a private product-only kinetic shadow
matrix/packet from product/product blocks, or to stop before broader adoption.
Full kinetic adoption requires support-dense and mixed fallback blocks, so it
should not be inferred from the product/product shadow result. Local/Gaussian
one-body terms, PQS/product mixed blocks, support-dense mixed packet
generalization, and all-pairs production assembly remain separate design
questions. Any adoption pass should remain private until it has explicit
operator checks, signed-block diagnostics, and a clear production-readiness
decision.
