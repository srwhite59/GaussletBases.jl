# Projected q-Shell Policy

## Purpose

Projected q-Shell (PQS) is the intended high-order shell abstraction for
regular atom-local and molecular/exterior shell construction. It replaces the
earlier preferred mental model of stitching separate endcap, panel, annulus,
or leftover slab pieces in the regular middle of the construction.

This page records the policy direction only. It does not make PQS the default
builder route, does not change any Hamiltonian or QW construction path, and
does not claim production science validation.

## Definition

PQS now has two distinct stages that must not be conflated.

The raw product-box reference stage is:

```text
orthogonal 1D COMX/source transforms
-> orthogonal q x q x L product-box modes
-> boundary COMX-product mode selection
-> raw-box retained operator reference
```

At this stage, no shell-row projection is applied and no Lowdin cleanup is
used. The retained functions live in the `q x q x L` product-box mode span,
not in a shell-row support subspace. Operator references are obtained by
building the product-box operator from 1D factors and selecting the same
boundary COMX-product columns on both sides:

```text
O_boundary = P_boundary' * O_product_box * P_boundary
```

The shell-realization stage is separate:

```text
selected product-box modes
-> project to shell rows
-> full-rank symmetric Lowdin cleanup
-> isometric shell-realized representation
```

Lowdin belongs to this later projection-plus-realization stage. It should not
be pulled into the raw product-box operator reference.

For a `q x q x L` local block:

1. Build or use the full local product/block transform for the whole block.
2. Select the boundary COMX-product modes: product columns whose local mode
   index is first or last on at least one axis.
3. Build raw-box retained operator references by boundary-column selection
   when staying in product-box mode space.
4. Only when a shell-supported representation is needed, project those
   selected product modes onto the raw outer-shell coordinates and apply
   full-rank symmetric Lowdin cleanup.

The shell must not depend on how the interior block was contracted. It must
not be defined by subtracting a previously locked contracted inner cube, and
it must not depend on locked prior spans. Interior contraction is a separate
choice; the selected PQS modes are boundary COMX-product modes from the full
local block transform. Shell-row projection plus Lowdin is the later
realization map, not the raw-box operator contract.

The subtraction formula is a counting picture, not the construction contract:

```text
dim PQS(q, L) = q^2 L - (q - 2)^2 (L - 2)
```

For `q = 5`, this gives `16L + 18`. The cubic case is `L = q`, so
`PQS(5, 5)` has dimension `5^3 - 3^3 = 98`.

## Policy

- Cubic atom shells are `PQS(q, q)`.
- Rectangular molecular/exterior shells are `PQS(q, L)`.
- Regular shells should be represented as PQS boundary traces, not as
  separate endcap/panel stitching.
- PQS must be independent of the contraction chosen for the interior block.
- Fit-to-box adjustments happen only at the outermost parent boundary.
- Generous boxes may omit awkward outer padding, but diagnostics must say
  exactly what was omitted.
- Regular mid-region construction should not leave leftover slabs.
- Regular construction should not create far split-box shells.
- PGDG analytic integrals remain mandatory unless an explicit
  `:numerical_reference` or debug path is requested.

The existing q-row and endcap/panel work in mainline is transitional
infrastructure and validation scaffolding. It remains useful for diagnostics,
receipt wrapping, sidecar plumbing, and CR2 handoffs, but it is not the final
preferred shell abstraction.

## External First-Gate Evidence

High-order ran a q=5 unified block-shell C2 `R = 4.7` preflight. The unified
shell row matched the q5 annulus dimension, beat the q4 and q5 annulus rows on
RHF/EGOI for that point, had clean support accounting, avoided dense fallback,
used the full analytic PGDG-style path, had exact H/V symmetry, and had clean
by-center/total H agreement.

The reported rows were:

| route | dimension | RHF total | EGOI relative Frobenius |
|---|---:|---:|---:|
| q4 four-panel/direct5 | `1482` | `-74.881526646276` | `0.00135160` |
| q5 true annulus | `1554` | `-74.888408985519` | `0.00142188` |
| q5 unified block shell/PQS candidate | `1554` | `-74.891723662832` | `0.000866239` |

This is positive first-gate evidence for PQS as the intended abstraction. It
deserves an accepted-R ladder follow-up before any default-policy decision.
It is not enough by itself to make PQS the default route.

## Mainline Implementation Checkpoint

As of 2026-05-28, mainline has a private opt-in PQS construction smoke. The
implementation trail is:

- `06483f8`: added the private projected q-shell local layer helper
- `bd66e51`: exposed PQS as a realization candidate in metadata
- `e4fbcca`: added opt-in PQS source and fixed-block construction
- `42103e3`: added the pure nested PGDG QW smoke for the opt-in PQS fixed
  block

The helper follows the intended boundary-mode contract: build the full local
COMX transforms, select boundary COMX-product modes, project those modes onto
raw boundary rows, and apply full-rank symmetric Lowdin cleanup. It does not
subtract a contracted inner cube and does not project against locked prior
spans.

The path remains private and explicit opt-in. Existing defaults continue to
use the endcap/panel route where they did before. The first small fixture has
full-parent dimension `735`; that is construction and QW-smoke evidence only,
not a compression-quality, energy, CR2, or production science claim.

A later private contract-correction checkpoint adds
`_pqs_raw_product_box_reference_block(...)` as a raw product-box self-block
reference path. It supports only overlap, `position_x`, `position_y`,
`position_z`, `x2_x`, `x2_y`, `x2_z`, and kinetic. It builds operator blocks
in the selected product-box mode span, uses boundary-column selection as the
oracle, records numerical 1D/product/selected-overlap checks, and explicitly
reports that shell projection is postponed. It does not use the stored
shell-row support coefficients, does not use Lowdin at the raw-box stage, does
not divide by retained PQS weights, and does not change packet construction,
QW/Hamiltonian, public/default, local/ECP/Gaussian/MWG, CR2, or IDA/MWG
behavior.

A follow-up private helper, `_pqs_product_box_realization_plan(...)`, now
builds both setup pieces together while keeping them separate: a raw
product-box plan with 1D operator factors and boundary selector, and a
shell-realization plan with shell-row projection, Lowdin cleanup, and isometry
diagnostics. Operators still start from the raw product-box 1D factor form;
the shell-realization map is not used while building raw-box operators.

The raw product-box side is now its own private plan boundary:
`_pqs_raw_product_box_plan(...)`. PQS self references, PQS/product mixed
source-box references, and the PQS/GTO cross-overlap shadow consume that raw
plan directly. `_pqs_shell_realization_plan(...)` remains the separate
projection-plus-Lowdin shell-row realization object, and
`_pqs_product_box_realization_plan(...)` is only the compatibility/diagnostic
wrapper that returns both plans together. This keeps the code aligned with the
intended staging: build source-box operators from 1D factors first, then apply
any shell realization later and explicitly.

The PQS raw plan can now wrap the shared private
`_cartesian_raw_product_box_plan(...)` helper when that shared source-box plan
is supplied. That helper owns intervals, total source-mode dimensions, z-fast
mode ordering, 1D axis transforms, axis-local coefficients, and retained-rule
free diagnostics. Descriptor-only PQS remains a structural fallback and reports
that the shared plan is unavailable because descriptors do not carry axis
bundles. Focused PQS/product and PQS/GTO source-box shadows now prove
shared-backed raw-plan consumption where available, still without packet
adoption, fixed-block construction changes, QW/Hamiltonian changes, IDA/MWG
changes, shell projection, Lowdin in raw-box operators, or public/default route
changes.

The raw-plan-first cleanup has now removed the remaining duplicate
descriptor-specific raw-box operator plumbing. Descriptor-level PQS raw-box
reference is only a convenience adapter that creates a raw plan and delegates.
The PQS/product shadow layout consumes the raw plan directly, while the
descriptor method only checks consistency and delegates. The PQS/GTO
source-box shadow likewise consumes the raw plan directly, with descriptor
input kept as an adapter only. The final-basis GTO handoff remains a separate
authoritative path.

The raw-source policy now records the corresponding private migration plan:
raw product-box records own intervals/source dimensions/1D transforms,
retained-rule records own boundary selection or shell realization, and pair
operator plans own cross-source operator factors. PQS should continue to move
through those boundaries rather than adding more descriptor-specific mixed
block paths.

The first PQS/product source-box mixed-block helper then uses this raw plan
with product/doside retained metadata. It is validated against explicit
source-box references for overlap, position, `x2`, and kinetic, including a
non-identity product-axis retained-transform check. This records that
PQS/product mixed source-box blocks handle nontrivial product retained
transforms, still private/reference-only and without adopting packet
construction or any QW/Hamiltonian route.

That mixed-block helper is already the direct source-box path: it does not
materialize a dense 3D raw source-box pair matrix. It builds 1D cross-axis
factors, selects PQS boundary modes, applies product/doside retained mode
metadata, and assembles retained blocks directly. The multi-term private
helpers `_pqs_product_source_box_reference_blocks_from_pair_plan(...)` and
`_pqs_product_source_box_reference_blocks(...)` reuse one PQS/product pair plan
across overlap, position, `x2`, and kinetic requests; `_pqs_product_source_box_shadow_blocks(...)`
uses that path for its PQS/product component blocks. This remains
private/shadow-only source-box infrastructure, with no shell projection,
Lowdin, retained PQS weight division, packet adoption, QW/Hamiltonian route, or
public/default behavior change.

`_pqs_product_source_box_shadow_blocks(...)` is the private two-block
layout/reference consumer for this path. It places one mode-selected PQS
source-box unit and one product/doside retained unit into a shadow layout and
fills PQS/PQS, PQS/product, product/PQS by transpose for symmetric real terms,
and product/product blocks. Its mixed PQS/product component blocks now use the
multi-term pair-plan reuse path. Covered terms are `:overlap`,
`:position_x/y/z`, `:x2_x/y/z`, and `:kinetic`. The tests include a
rectangular PQS source box and a non-identity product/doside transform. This
remains private source-box shadow evidence only: no shell-row projection, no
Lowdin, no `support_coefficient_matrix` PQS oracle, no retained PQS weight
division, no packet adoption, and no QW/Hamiltonian, public/default, CR2,
local/ECP/Gaussian/MWG/interaction, or IDA/MWG behavior change.

A private GTO cross-overlap shadow now extends the same source-box boundary.
`_pqs_source_box_gto_cross_overlap_shadow(...)` uses
`_pqs_source_box_gto_axis_projection(...)` to project existing 1D
Cartesian/GTO primitive-axis overlap tables through the PQS source-box axis
coefficients, then assembles overlaps for selected boundary COMX-product
modes. The rectangular `5 x 5 x 7` test compares this result against an
explicit dense source-box reference:

```text
full source-box coefficients over product-box parent rows
-> parent-row/GTO overlap
-> transpose(C_source_box) * S_parent_gto
-> boundary-column selection
```

This is source-box shadow evidence only. The current final-basis GTO handoff
remains unchanged and authoritative for handoff use. The helper does not use
shell-row projection, Lowdin, `support_coefficient_matrix`, retained PQS
weights, IDA division, packet construction, QW/Hamiltonian, public/default,
CR2, local/ECP/Gaussian/MWG/interaction, or IDA/MWG behavior.

The next implementation, if scoped, should consume the shell-realization plan
through a separate helper:

```text
selected product-box modes
-> project to shell rows
-> Lowdin cleanup
-> isometric shell-realized representation
```

That helper should not be merged into the raw product-box operator reference.

The pure nested QW smoke uses the existing receipt path with
`gausslet_backend = :pgdg_localized_experimental`,
`interaction_treatment = :ggt_nearest`, and `nuclear_term_storage =
:total_only`. It reports clean source/sidecar agreement, finite symmetric
overlap/one-body/interaction matrices, zero residuals, and no warning-level
numerical-quadrature logs for that path.

The missing next design problem is product-staged PQS sidecar and performance
work. PQS should not be adopted for by-center, supplement, or
performance-sensitive routes until that sidecar/performance contract exists.

## Fixed-Only Be2 Broad-Parent Diagnostic

A private fixed-only Be2 PQS probe was run from the ignored repo script
`tmp/work/be2_pqs_fixed_only_mvp_probe.jl`, using the existing output directory
`tmp/work/be2_pqs_fixed_only_mvp_outputs/`. It must be classified as a
broad-parent-boundary diagnostic and construction-contract failure relative to
the intended q=5-local PQS MVP, not as a passed PQS MVP gate. It is not a
supplement-coupled route, not HF or energy validation, not CR2 science
validation, and not public/default route adoption.

The target was Be2/cc-pVDZ at `R = 5.0`, all-electron `Z = 4`, `q = 5`, with
the CR2 parent grid `(13, 13, 25)` and parent dimension `4225`. The probe
confirmed the z-fast parent ordering and parent grid coordinates/weights
against the repo-built parent basis. The fixed dimension was `3409`, which is
near-parent behavior for a parent dimension of `4225`. The resulting
capture/H1 numbers are therefore not meaningful evidence that the intended q=5
PQS compression works.

The QW construction used the PGDG backend
`gausslet_backend = :pgdg_localized_experimental`, did not use numerical
fallback, and reported `dense_parent_matrix_used = false`. The fixed-only
capture/H1 readback was:

| spin | fixed capture | worst orbital | max fixed H1 delta |
|---|---:|---|---:|
| alpha | `0.999997644256` | `alpha_mo_4` | `0.136217 mHa` |
| beta | `0.999997644256` | `beta_mo_4` | `0.136217 mHa` |

Those numbers are expected to look good for a near-full-parent fixed object
and should not be compared as a successful q=5 PQS MVP gate. The timing
checkpoint was `88.27 s` for fixed-block build, `7.57 s` for parent QW
construction, and `16.60 s` for PQS fixed QW construction.

The construction audit found that the current shared PQS regions used broad
region dimensions rather than the selected recipe q. The two regular shared
regions were constructed as:

- `PQS(13, 23)`, retained count `1346`;
- `PQS(11, 21)`, retained count `1002`.

The intended q=5-local rectangular shell counts for the same bond-axis lengths
would be much smaller:

- `PQS(5, 23) = 386`;
- `PQS(5, 21) = 354`.

The immediate bug is that `raw_q` currently comes from the full shared-region
transverse size, not from the selected recipe q. The next fix is to make
shared-shell PQS construction q-local, or to explicitly label broad
parent-boundary PQS as reference-only and keep it out of MVP/pass gates. Do
not run supplement-coupled PQS until this construction contract is fixed.

## Corrected Be2 Fixed-Only Interpretation

After the q-local shared-shell correction (`7fc94b7`), the fixed-only Be2 PQS
probe was rerun on the corrected CR2 target: Be2/cc-pVDZ at `R = 5.0`,
all-electron `Z = 4`, `d = 0.15`, parent axes `(15, 15, 27)`, parent
dimension `6075`, and `q = 5`. The run remained private/report-only, used the
PGDG backend, did not use numerical-reference fallback, and reported
`dense_parent_matrix_used = false`.

The corrected q-local PQS source diagnostics are the important construction
checkpoint:

- shared physical boxes remain broad: `(15, 15, 25)`, `(13, 13, 23)`,
  `(11, 11, 21)`;
- shared source-mode dimensions are q-local/adaptive: `(5, 5, 5)`,
  `(5, 5, 5)`, `(5, 5, 6)`;
- shared `pqs_retained_count` values are `98`, `98`, and `114`;
- broad physical PQS forms such as `PQS(15,25)`, `PQS(13,23)`, and
  `PQS(11,21)` are not normal MVP PQS routes.

The source-mode dimensions above are total side lengths. The private
`_nested_diatomic_source_box_dimension_plan(...)` helper now records that
contract explicitly: angular-spacing policy chooses the source-box dimensions,
then COMX/source transforms are deterministic for those dimensions, and the
PQS retained rule acts afterward. Any lower-level selector counts are
diagnostic details, not the PQS-facing side lengths.

The strict q=5 PQS fixed-only result was:

| route | fixed dim | fixed capture | max fixed H1 delta |
|---|---:|---:|---:|
| q-1 panel fixed-only | `1461` | `0.998655961965` | `4.987085 mHa` |
| strict PQS q=5 fixed-only | `1483` | `0.999918551122` | `0.125908 mHa` |
| q5 q-row/endcap-panel fixed-only | `1623` | `0.999955518077` | `0.211455 mHa` |

The interpretation is therefore two-sided:

- strict PQS q=5 passes the compact fixed-only comparison against the q-1
  endcap/panel route;
- strict PQS q=5 remains below the richer q5 endcap/panel fixed-capture
  baseline;
- the q5 endcap/panel route should not be the only fixed-space comparator,
  because it retains a richer shared boundary object than strict q=5 PQS.

A diagnostic-only broader PQS variant using
`shared_shell_angular_resolution_scale = 1.1` reached fixed dimension `1549`,
fixed capture `0.999967462999`, and max fixed H1 delta `0.126753 mHa`. It
passes the q5 endcap/panel fixed threshold, but it is not a promoted policy:
it changes shared source-mode dimensions to `(6, 6, 5)`, `(5, 5, 6)`, and
`(5, 5, 7)`.

The supplement-coupled strict PQS q=5 probe was then run on the same corrected
Be2 target: Be2/cc-pVDZ at `R = 5.0`, all-electron `Z = 4`, `d = 0.15`,
physical box `16 x 16 x 21` bohr, parent axes `(15, 15, 27)`, parent
dimension `6075`, and strict `q = 5`. It used the PGDG backend, did not use
numerical-reference fallback, and reported `dense_parent_matrix_used = false`.
The source diagnostics stayed q-local/adaptive: shared source-mode dimensions
`(5, 5, 5)`, `(5, 5, 5)`, `(5, 5, 6)` with shared `pqs_retained_count` values
`98`, `98`, and `114`.

The fair supplement-coupled comparison is against the compact q-row route with
`5x5` endcaps and `4xL` panels (`protected_atom_side_count = 5`, `q_min = 4`,
`shared_q = 4`, `shared_order = 4`):

| route | fixed/final/residual dims | fixed capture | final capture | max fixed H1 delta | max final H1 delta |
|---|---:|---:|---:|---:|---:|
| fair q-row q-1 panel | `1461 / 1479 / 18` | `0.998655961965` | `0.999980227448` | `4.987085 mHa` | `2.333419 mHa` |
| strict PQS q=5 | `1483 / 1501 / 18` | `0.999918551122` | `0.999989138049` | `0.125908 mHa` | `2.339571 mHa` |

Strict PQS q=5 clearly beats the fair fixed q-row baseline. In the
supplement-coupled final basis it also beats the fair final capture, while its
max final H1 delta is slightly worse by about `0.006152 mHa`. Relative to the
richer q5 q-row/endcap-panel supplement baseline (`1623 / 1641 / 18`, final
capture `0.999989372730`, max final H1 delta `1.874443 mHa`), strict PQS q=5
remains slightly behind in both final capture and final H1. This is therefore a
competitive private MVP result, not a production/default route adoption.

The receipt-audit compatibility checkpoint is recorded by
`a47d2bf Add PQS-compatible hybrid overlap fallback`. Diatomic hybrid overlap
sidecars now support a dense exact fallback for fixed columns that do not have
product/factorized parent data. PQS fixed columns are therefore not marked as
product-factorized. Product/q-row fixed columns still keep the existing
factorized sidecar path and its reconstruction checks.

The dense fallback is an exact parent-to-supplement overlap audit. Its
coefficient factor is a block/support-sparse parent-to-fixed coefficient map
materialized densely for exact parent-to-supplement overlap audit, then
contracted as `parent_to_fixed_coefficients' * overlap_ga`. It is not a dense
parent-parent operator and must not be confused with
`dense_parent_matrix_used` in QW construction.

The strict PQS q5 supplement probe now reports
`final_receipt_audit_available = true`. Its hybrid overlap sidecar kind is
`:dense_bond_aligned_diatomic_mixed_raw`, and the receipt-built final
operators match the direct supplement builder output with max matrix error
`0.0`. This is an audit/representation compatibility fix only: it does not
change supplement construction, QW/Hamiltonian assembly, public/default route
adoption, CR2 artifacts, optimized PQS/product kernels,
local/ECP/Gaussian/MWG handling, or numerical fallback policy.

Interaction and weight diagnostics must keep two contracts separate. Active
fixed/final operators use the intended gausslet IDA/MWG electron-electron
route through the existing nested pair-sum interaction data. Retained PQS
per-column weights in private executable payloads are different: they are
debug/reference metadata, not generic positive quadrature masses and not IDA
division weights. PQS diagnostics should therefore report retained-column
weight role `:debug_reference_only`, `ida_weight_division_allowed = false`, and
no claimed quadrature-weight semantics. Those flags do not disable or weaken
the active IDA/MWG interaction path; they only prevent PQS retained-transform
weights from being mistaken for the positive weights owned by the pair-sum
construction.

The private density-density SCF smoke then used the existing matrix-level
`_closed_shell_density_density_scf` helper on the final
`OrdinaryCartesianOperators3D` data: overlap, H1, and `interaction_matrix`.
Both routes used `interaction_treatment = :mwg`, the PGDG backend, no
numerical-reference fallback, finite/symmetric H and V matrices, and finite
positive residual metadata.

| route | final dim | residuals | converged | iterations | SCF energy |
|---|---:|---:|---:|---:|---:|
| strict PQS q5 final | `1501` | `18` | `true` | `44` | `-32.338778224290 Eh` |
| fair q-row q-1 panel final | `1479` | `18` | `true` | `44` | `-32.342180822067 Eh` |

The strict PQS q5 occupied orbital energies in this smoke were
`(-4.7317411478, -4.7317331742, -0.3803502823, -0.2518522168)`. Strict PQS q5
and the fair q-row final route both converge in this density-density IDA/MWG
SCF smoke. The strict PQS q5 final energy is higher than the fair q-row final
by `+0.003402597777 Eh`. IDA/MWG is the intended gausslet electron-electron
Hamiltonian form; the repo does not normally use a Galerkin four-index Coulomb
tensor with gausslet bases. The caution is narrower: lower IDA/MWG energy
should not be reported as automatically better basis quality in the usual
basis-set variational sense.

This proves private HF-smoke compatibility of the strict PQS q5 final
operators for the intended gausslet RHF Hamiltonian form. It is not production
HF validation, does not add a public HF API, and does not change
QW/Hamiltonian construction, IDA/MWG semantics, CR2 artifacts, or
public/default routes.

The external conventional-Coulomb RHF reference ladder for all-electron Be2 at
`R = 5.0` bohr is:

| PySCF basis | RHF total energy | delta vs 5Z |
|---|---:|---:|
| cc-pVDZ | `-29.136053493620 Eh` | `+1.622 mHa` |
| cc-pVTZ | `-29.137357901570 Eh` | `+0.317 mHa` |
| cc-pVQZ | `-29.137579392916 Eh` | `+0.096 mHa` |
| cc-pV5Z | `-29.137675195692 Eh` | `0.000 mHa` |

`aug-cc-pVQZ` lowers `cc-pVQZ` by only `0.018 mHa`. These PySCF totals include
nuclear repulsion. For Be2 with `Z = 4` and `R = 5.0` bohr, the nuclear
repulsion is `3.2 Eh`, so the PySCF cc-pV5Z electronic-only energy is
approximately `-32.337675195692 Eh`. Repo private SCF smoke energies appear to
be electronic energies unless nuclear repulsion is added separately. Compare
total-to-total or electronic-to-electronic only, and keep the Hamiltonian model
difference explicit: PySCF uses the conventional Gaussian four-index Coulomb
Hamiltonian, while the repo smoke uses the intended gausslet IDA/MWG
Hamiltonian.

The near-term scientific step should be a Be2 gausslet-RHF ladder under the
intended IDA/MWG Hamiltonian, tracking q/source-mode/basis-size convergence.
The PySCF RHF ladder should be used as an external conventional-Coulomb
reference, not as a direct same-model variational target. `cc-pVDZ` remains
useful as a small plumbing/capture target; `cc-pVQZ` or `cc-pV5Z` is the more
appropriate conventional-Coulomb reference for a nearly converged Be2 RHF
number.

## Sidecar Prototype Checkpoint

As of 2026-05-28, mainline also has a private descriptor/prototype line for
future PQS sidecar work. This line is internal scaffolding only; it is not
installed into fixed blocks, contracted-parent metrics, QW builders, or public
routes.

The staged descriptor now records the data needed to reproduce the PQS shell
metric without reinterpreting the construction:

- boundary COMX-product mode indices selected from the full local block;
- raw-boundary support rows;
- axis-local COMX transforms and axis intervals;
- the stored full-rank symmetric Lowdin cleanup transform.

Two private metric prototypes have been checked:

- a support-local descriptor prototype that reproduces overlap and weights and
  acts as the overlap-invariant debug oracle;
- a slab/product prototype that decomposes the raw boundary into six
  rectangular pieces without building a support-local boundary overlap matrix
  or any dense full-parent matrix.

The overlap contract was then tightened: PQS overlap is an orthonormality
invariant, not an operator-contraction target. The slab/product helper returns
identity overlap after the PQS cleanup contract has been established, while the
support-local prototype records the debug overlap error. The nontrivial
slab/product checks now cover weights and first moments (`x`, `y`, `z`) against
the existing support-reference packet.

The first benchmark, run before the overlap-invariant correction, was
direction-setting rather than publication-quality. On small q=5 fixtures it
found:

| fixture | support/retained | support-local prototype | slab/product prototype |
|---|---:|---:|---:|
| `PQS(5,5)` | `98 / 98` | `0.000084 s / 0.471 MiB` | `0.000334 s / 0.510 MiB` |
| `PQS(5,7)` | `130 / 130` | `0.000093 s / 0.847 MiB` | `0.000521 s / 0.798 MiB` |
| `PQS(5,13)` | `226 / 226` | `0.000348 s / 2.349 MiB` | `0.001611 s / 1.995 MiB` |

The interpretation is deliberately conservative: the slab/product helper is
correct and allocation-promising for longer rectangular shells, but it is not
yet a speed win on small q=5 fixtures. Before sidecar adoption, the next
technical work should optimize the product path with caching, scratch reuse,
and boundary-mode grouping. This checkpoint does not justify immediate
by-center, supplement, QW, Hamiltonian, CR2, or science-route adoption.

## Private Sidecar-Fixture Checkpoint

The private checkpoint first added a sidecar-shaped executable resolved-payload
fixture for PQS low-order metric checks. Descriptor-only PQS remains
unsupported/prototype: without an explicit `column_range` it still reports
missing installed sidecar payload fields and is not consumed by production
metric paths or public/default routes.

The fixture combines the existing PQS descriptor with an explicit
`column_range`, raw-boundary support coefficients, boundary COMX mode data,
and the stored full-rank symmetric Lowdin cleanup transform. It is fixture-only
and does not make PQS production-supported.

A follow-up checkpoint now adds a private single-PQS-layer `_NestedFixedBlock3D`
fixture that stores `_CartesianProjectedQShellSidecarFixture3D` in a
fixed-block sidecar slot, but only exposes it behind a PQS-specific private
accessor. Ordinary by-center, QW, and Hamiltonian consumers do not consume this
sidecar, and `_nested_staged_by_center_sidecar` remains incompatible/loud for
the PQS fixture. This proves fixed-block sidecar attachment/discovery only; the
mixed q4 recipe fixed block is not covered by this fixture.

Validated low-order checks are deliberately narrow:

- PQS self blocks use identity self-overlap only as the post-cleanup
  orthonormality invariant, and validate weights, position matrices, and first
  moments against the existing reference path.
- PQS/support-dense mixed blocks use the support-local reference path. They do
  not use the identity shortcut for mixed overlap.
- PQS/product optimized or adopted metric blocks remain explicitly unsupported.
  Private reference helpers now exist for overlap, `position_x/y/z`, and
  kinetic; they are separately named reference/debug paths, not silent
  fallbacks from the optimized path and not packet adoption.

This is low-order metric/reference readiness only. The PQS/product kinetic
helper is a private signed-operator reference using
`(K,S,S) + (S,K,S) + (S,S,K)`; it does not imply readiness for `x2`,
nuclear/local one-body terms, Gaussian or pair terms, interactions,
QW/Hamiltonian construction, H1, energy, CR2 validation, or any science route.
No default builders, public APIs, backend policy, PGDG/quadrature policy, or
QW/Hamiltonian paths changed.

If scoped later, the next implementation should still stay in a tiny private
fixture lane unless a separate design explicitly covers mixed q4 recipe sidecar
coverage and non-low-order metric terms.

## Consequences For Mainline

The next source-construction design should treat PQS as the regular shell
target:

- atom-local q shells should move toward `PQS(q, q)`;
- molecular/exterior q shells should move toward `PQS(q, L)`;
- existing endcap/panel helpers should be treated as transitional support and
  audit scaffolding unless explicitly retained for a boundary-only special
  case;
- diagnostics should distinguish construction by raw-boundary projection from
  any counting-only annulus formula;
- any route that trims or omits outer support for fit-to-box reasons must
  report that as an outer-boundary adjustment, not as a regular shell policy.

Before adoption, repo-side work still needs a dedicated implementation and
validation plan. That plan should compare PQS against the transitional
endcap/panel and annulus rows without changing public defaults, backend
policy, Hamiltonian kernels, or quadrature behavior.
