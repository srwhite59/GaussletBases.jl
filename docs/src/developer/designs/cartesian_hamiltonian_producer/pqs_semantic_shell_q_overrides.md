# Semantic Per-Shell PQS Source-q Overrides

Status: implemented internal opt-in facility with completed bounded H2 and
padded Be2 validation.
Authority is recorded under `HP-PQS-SHELLQ-OVERRIDE-FN-01` and
`HP-PQS-SHELLQ-OVERRIDE-TEST-01`.

This page is the canonical contract for bounded source-mode refinement of
selected complete shells in bond-aligned homonuclear diatomic PQS
construction. The registry owns lifecycle, permission, and exact source/test
ceilings. This page owns input semantics, construction order, failure
behavior, validation, and exclusions.

## Purpose

The global route `q` controls the ordinary producer construction. A
measurement may need to ask a narrower question:

> What changes if one semantic complete shell receives a larger PQS source
> span while parent geometry, shell support, global `ns/q`, and every other
> terminal region remain fixed?

This lane permits that question only inside the private numerical-complete
additive-reference composition with `source_span = :ordinary`. It is a general
producer diagnostic, not a Cr-specific branch, production policy, or automatic
refinement rule.

Changing a shell source span changes its retained dimension and therefore the
final basis and Hamiltonian. It does not change which parent points the shell
owns.

## Private Input

The private input is `source_mode_overrides`, for example:

```julia
source_mode_overrides = [
    (;
        role = :atom_local_shell,
        shell_index = 7,
        owner = :all,
        source_q = 8,
    ),
]
```

Omitting the keyword and supplying an empty collection are equivalent and
must reproduce ordinary construction exactly.

Each record has exactly four fields:

| Field | Contract |
| --- | --- |
| `role` | `:atom_local_shell` or `:shared_molecular_shell` |
| `shell_index` | positive native semantic shell index |
| `owner` | exactly `:all` |
| `source_q` | non-Boolean integer strictly larger than the route `q` |

Unknown fields and alternate spellings are errors. In particular, callers may
not provide `q`, lowering/retained `q`, `L`, a full source shape, a terminal
region key, an order index, or a route-local ordinal.

Duplicate `(role, shell_index, owner)` records are errors even when their
`source_q` values agree. There is no last-wins rule. Record order has no
semantic effect; normalized application order must be deterministic.

`source_q` is a shell-local source-span size. It is not public/route `q` and
does not alter `ns`, the route's PQS `q = ns`, direct-core side, parent axes,
shellification, or lowering policy. Mapped-COMX source spans are rejected in
this first lane because their later carried-axis enrichment has not yet been
proved compositional with an overridden shell shape.

## Semantic Matching

Overrides are matched only after shellification has produced terminal regions.
The identity is:

```text
region.role
raw_region.shell_index
```

and the region must also be a `:complete_shell` selected for
`:pqs_filled_source_cpb` lowering.

`shell_index` is the positive native counter within the named semantic role.
It is not `terminal_region_22`, `order_index`, array position, or another
unstable key.

For `role = :atom_local_shell`, `owner = :all` means the symmetric pair of
owner-local shells at that semantic index. Both owner regions must exist and
must be refined together. A one-sided match, unequal source request, or
geometry with no such pair is an error.

For `role = :shared_molecular_shell`, the semantic index identifies the
corresponding shared complete shell. The supported homonuclear z-axis path has
one shared region for that index. Zero or multiple matches are errors.

No owner-specific or asymmetric request is accepted in this lane.

## Construction Order

The implementation must preserve this order:

```text
resolve parent and route q
build parent axes
shellify terminal support
build the ordinary PQS lowering plan
validate and match semantic overrides
rewrite only matched complete-shell source dimensions
freeze lowering inventory and retained/support/transform records
realize terminal basis and operators
```

The existing post-shellification terminal low-order route is the owner. An
override must not be applied in the early region-contract builder, where the
shared-shell angular selector lacks parent/bundle context, or after retained
and support records have already frozen dimensions.

The normalized match map must rewrite both `available_contracts` and
`selected_contracts` from the same source dimensions, as the current aspect-
shell enrichment does. Inventory and retained-unit construction then consume
the rewritten plan.

### Atom-local complete shell

For a matched atom-local shell:

```text
source_mode_shape = (source_q, source_q, source_q)
```

The existing PQS boundary-product retained rule is recomputed from this
authoritative shape. Global route `q`, parent axes, support boxes, shell
ownership, and region order remain unchanged.

The lowering contract's existing `q` field remains `route_q`. The atom-local
rewrite sets `source_mode_shape` and recomputes the retained count from that
shape; it does not fabricate shared-shell `raw_q`, `raw_L`, or `L` selector
facts.

### Shared molecular complete shell

For a matched shared shell, `source_q` is the transverse source size. The
existing angular-band dimension selector must be rerun with that value and the
same parent, atomic locations, outer/inner boxes, bond axis, support count, and
retention policy. It selects the longitudinal dimension `L`:

```text
source_mode_shape = (source_q, source_q, L)
```

The caller cannot choose `L`. Physical aspect estimates remain advisory
due-diligence comparisons and do not replace the angular-band selector.

The argument split is exact:

```text
lowering-contract q                 = route_q
complete-shell retention resolver  = source_q
angular selector nside             = source_q
angular selector selected_q        = source_q
```

The rewritten contract records the selector's existing
`source_mode_dims/raw_source_dims`, `selected_q`, `raw_q`, `raw_L`, `L`, and
axis-retained-count fields. `selected_q` is `source_q`; `raw_q/raw_L` are the
selector results; `L = raw_L`; and `pqs_retained_count` is recomputed from the
authoritative returned shape. None of these fields redefines route `q`.

## One Authoritative Shape

The rewritten lowering contract owns the shape. Existing downstream consumers
must derive all affected dimensions from it:

- lowering-contract inventory;
- retained-unit and support records;
- retained-unit transform contracts;
- raw-product retained rules and multilayer source plans;
- terminal realization and final column ranges;
- source-mode provenance already available in memory;
- terminal inventory and due-diligence rows;
- exact one-body, IDA, residual, MWG, packet-reference, and correction
  construction rebuilt for the resulting basis.

No downstream stage may independently reconstruct the route-default shape or
patch a second copy. The implementation adds no persistent override summary,
new status/warning symbol, artifact field, or alternate shell builder. The
caller-owned private recipe may retain its input, but producer results gain no
new field.

Due-diligence acceptance identifies rows by their actual terminal `role` and
semantic `shell_index`. The coarser `owner_contact_shared` classification is a
reported classification, not an override key.

If the existing authoritative `source_mode_shape` does not reach any required
consumer, stop and report that exact missing seam. Do not broaden the approved
source surface to repair several downstream copies.

## Approved Source Boundary

Implementation is limited to:

- `src/pqs_source_box_route_driver_helpers.jl` for strict normalization,
  semantic matching, and lowering-contract shape enrichment;
- `src/cartesian_base_hamiltonian.jl` for a private keyword on
  `cartesian_base_working_basis(...)` and narrow internal-stage forwarding;
- `src/cartesian_protected_ladder_bundle.jl` for private forwarding from the
  `_plb_build_numerical_complete_additive_reference_member(...)` call site
  through an explicit `_plb_build_inputs(...)` keyword;
- `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` for bounded
  correctness coverage.

The private keyword is not a basis `NamedTuple` key and is not forwarded by
the exported base facade or canonical driver. No new module, file, struct,
persistent result shape, public input, artifact field, or committed large
fixture is approved.

Existing multilayer realization files are intentionally outside the source
ceiling. Their current authoritative-shape path must work unchanged or the
implementation stops.

`_plb_build_inputs(...)` must default to an empty override and must not inspect
the recipe implicitly. Only
`_plb_build_numerical_complete_additive_reference_member(...)` may extract the
private recipe value and pass it explicitly. Compact/protected,
numerical-complete without additive references, protected additive, ladder,
and artifact-producing entry points reject a nonempty override. A nonempty
override with empty packet placements is also an error; it must not fall back
to an ordinary numerical-complete member.

The approved path is in-memory only. Protected member/ladder artifact writers,
manifests, recipe persistence, and readback do not consume or preserve this
input. The permitted call site may pass its caller-owned value through the
explicit `_plb_build_inputs(...)` keyword, but no writer or returned producer
object gains a field.

## Validation Contract

### Ordinary parity

- Omitted and explicitly empty overrides produce identical parent, terminal
  basis, `H1`, `Vee`, residual, packet, correction, and due-diligence results.
- PQS construction outside the private numerical-complete additive path is
  unchanged.

### Bounded construction gates

- Exercise one atom-local and one shared-shell refinement on bounded H2
  constructions.
- Exercise a physically padded Be2 numerical-complete additive-reference
  construction after the H2 gates pass.
- Confirm that exactly the selected semantic shell rows change
  `source_mode_shape`; atom-local `owner = :all` changes both paired rows.
- Confirm parent axes, physical/index shell boxes, support ownership, region
  order, direct cores, angular-extension slabs, thin slabs, and all unmatched
  regions are unchanged.
- Confirm retained counts agree with the active boundary formula for each new
  authoritative shape.
- Rebuild and report final dimension, numerical residual count, packet
  occupied capture, `J0`, `E0`, correction anchors, source shapes, and every
  terminal due-diligence warning.
- Require finite/symmetric `H1` and `Vee`, valid residual metrics, and the
  existing strict packet-capture checks.

No HF or energy endpoint is approved in the source pass.

### Rejection gates

Tests must reject:

- duplicate or unmatched semantic targets;
- zero/negative/noninteger/Boolean shell indices or source sizes;
- `source_q <= route_q`;
- partial-owner or asymmetric atom-local requests;
- unsupported roles or extra fields;
- terminal-region/order/ordinal keys;
- direct `L`, full-shape, or lowering-q requests;
- one-center and White-Lindsey construction;
- mapped-COMX source spans and any artifact-producing use.

## Consumer Measurement

With the source implementation and H2/Be2 gates accepted, an ignored CR2
measurement may refine shell 7 and shell 8 independently at
`source_q = 8`, then their pair only if the individual measurements justify
it.

The baseline and variants must use the same commit, system geometry, parent
controls, supplement, Coulomb policy, packet, and numerical-complete settings;
only `source_mode_overrides` may differ.

Every variant imports the same explicit external-GTO/PySCF occupied
determinant independently:

```text
C_variant = <variant | G_external> C_external
```

Compare occupied capture and density trace in each variant. Do not add or
materialize a dense baseline-to-variant final-basis overlap; at dimensions near
`7053`, that object is hundreds of megabytes and is unnecessary for this
question. An occupied-only baseline/variant comparison may remain an ignored
CR2-local diagnostic.

Before HF, report the imported-state `H1`, valence `Vee`,
core/core-valence/valence decomposition, residual weight, capture, dimensions,
source shapes, and due-diligence warnings. This remains measurement evidence,
not a Cr2 production claim. These IDs do not themselves authorize HF; any
later run must already have separate consumer-side authority and must wait for
review of the static comparison.

## Failure Behavior

Stop if:

- a semantic target cannot be matched uniquely and symmetrically;
- the angular selector cannot produce a valid shared-shell shape;
- a matched override changes parent geometry, support, ownership, or an
  unmatched region;
- downstream records disagree on shape or retained count;
- a requested span exceeds the shell's realizable support/rank or makes the
  shell projection numerically invalid;
- implementation requires edits outside the approved source surfaces;
- H2 or padded Be2 develops invalid metrics, poor packet capture, nonfinite or
  asymmetric operators, or a bad low mode.

Do not weaken capture/rank tolerances, floor residual eigenvalues, add occupied
injection, or create compatibility metadata to force the study through.

## Explicit Non-Goals

This authority does not approve:

- a production/default `q` change or automatic source refinement;
- public facade, basis, canonical-driver, or solver controls;
- parent-grid, shell-support, shellification, direct-core, slab, or
  aspect-policy replacement;
- White-Lindsey or one-center behavior;
- mapped-COMX source-span behavior or a second COMX realization path;
- arbitrary owner-specific asymmetry;
- residual cutoff, MWG, screened-Hartree, EGOI, or artifact changes;
- shell-projection, Lowdin, sign, rank, or overlap-policy changes;
- `C' V C`, source-Hamiltonian transforms, or interaction rotation;
- a dense final-final overlap or durable transfer API;
- Cr2-specific source behavior, production energies, or paper claims.
