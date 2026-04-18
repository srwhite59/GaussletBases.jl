# Cartesian Nested Decomposition Plan

## Purpose

This note records the first structural decomposition plan for the two largest
mixed-responsibility files in the repo:

- `src/cartesian_nested_faces.jl`
- `src/ordinary_qiu_white_rg.jl`

The goal is to reduce code-shape bloat and improve maintainability without
deleting real capability or silently changing numerical behavior.

This is a planning note, not a claim that the decomposition has been completed.

## Scope And Non-Goals

This plan is about:

- ownership boundaries
- file/module split
- hot-path visibility
- lower-risk refactor chunking

This plan is not about:

- public API redesign beyond narrow seam cleanup
- deleting experimental chain/square capability
- removing reference/legacy behavior unless it is separately justified
- changing packet math, geometry policy, or operator formulas

## Audit Summary

### `src/cartesian_nested_faces.jl`

This file currently combines at least seven distinct responsibilities:

1. shared nested data structures and coefficient-map utilities
2. packet kernels and compact term-sum contraction paths
3. generic shell / sequence assembly
4. one-center atomic public fixed-block helpers
5. bond-aligned diatomic geometry and source assembly
6. experimental chain and square-lattice source assembly
7. timing, trace, and reporting helpers

The most important live seams already exist in the code:

- timing/reporting around `nested_fixed_block_timing_report` and the diatomic
  doside trace helpers
- packet construction around the compact-only packet builders and weighted
  term-sum helpers
- one-center atomic helpers around `_build_one_center_atomic_shell_sequence`
  and the public atomic fixed-block / diagnostics surface
- diatomic source logic around `_nested_bond_aligned_diatomic_sequence_for_box`,
  `_nested_bond_aligned_diatomic_split_geometry`, and
  `_nested_bond_aligned_diatomic_source`
- experimental chain/square source logic around
  `_nested_bond_aligned_homonuclear_chain_source` and
  `_nested_axis_aligned_homonuclear_square_lattice_source`

### `src/ordinary_qiu_white_rg.jl`

This file currently combines at least six distinct responsibilities:

1. public QW basis containers and basis constructors
2. public nested diatomic/chain/square front ends and geometry diagnostics
3. QW timing/report helpers
4. residual-space diagnostics, keep policy, and stabilization logic
5. raw atomic / diatomic block construction and supplement plumbing
6. final operator assembly for ordinary, nested-fixed-block, and
   supplement-bearing routes

The main live seams are:

- basis construction at the top of the file
- nested-source/fixed-block front ends in the `1054:1631` region
- residual analysis in the `1712:3819` region
- ordinary and nested operator assembly in the `4375:end` region

## Current Role Classification

### Current Production Path

The current production or production-adjacent path is:

- compact-only nested packet construction
- one-center atomic fixed-block construction
- bond-aligned diatomic nested fixed-source/fixed-block/diagnostics
- ordinary QW operator assembly on:
  - direct bond-aligned ordinary bases
  - one-center atomic nested fixed blocks
  - bond-aligned diatomic nested fixed blocks
- exact Cartesian overlap/projector/transfer support that consumes these
  payloads elsewhere

### Current Experimental But Still Live

The current experimental but still live path is:

- bond-aligned homonuclear chain basis, nested geometry, and consumer payloads
- axis-aligned homonuclear square-lattice basis, nested geometry, and consumer
  payloads
- molecular supplement / residual-space analysis on QW routes
- diatomic doside / odd-retain diagnostics and geometry export support

These are not dead code. They are real repo capability, but they should not
obscure the compact production path.

### Legacy / Internal / Reference-Only

The following are still live enough to preserve, but they are not the primary
production story:

- one-center atomic legacy-profile working-box helpers
- compatibility timing-report labels and timing collector scaffolding
- doside trace writers and debug-geometry helpers
- residual keep-policy compatibility aliases and narrow legacy interaction
  treatment diagnostics

The problem is not that these exist. The problem is that they currently live in
the same files as the hot production kernels and the primary public front ends.

## Suspected Bloat Vs Justified Complexity

### Suspected Bloat

- Public front ends, experimental front ends, and low-level kernels are mixed
  in the same files.
- Timing/reporting helpers are interleaved with hot kernel code.
- Chain/square experimental paths live beside production diatomic and atomic
  paths.
- Tests reach into internal nested structs directly, which makes file-local
  cleanup harder than it should be.
- `src/ordinary_qiu_white_rg.jl` interleaves residual diagnostics with raw
  operator construction and final operator assembly.

### Justified Complexity

- The repo really does support two distinct nested packet kernels:
  `:factorized_direct` and `:support_reference`.
- The one-center atomic and bond-aligned diatomic geometry policies are not the
  same problem and should not be forced into one policy file.
- Chain and square-lattice experimental geometries are genuinely different
  split problems and deserve their own space.
- Residual-space analysis is numerically nontrivial and deserves a dedicated
  implementation area rather than being flattened into operator wrappers.

## Performance Pressure Points From Code Shape

The present code shape is likely hurting performance in a few concrete ways:

- Hot packet builders still branch on runtime `Symbol` mode selectors such as
  `packet_kernel`, `interaction_treatment`, and related front-door mode values.
  That makes it harder to isolate specialization-friendly hot paths.
- `_CartesianCoefficientMap` is a `Union` over dense and sparse matrix storage.
  That is sometimes justified, but keeping both paths interleaved in the same
  large file makes inference and allocation auditing harder.
- `_CartesianNestedShellPacket3D` and `_NestedFixedBlock3D` still carry
  optional fields for compact and noncompact term storage, even though compact
  is now the real production path.
- Source/bundle reuse seams are harder to see than they should be. When the
  front-end path is hidden inside a mixed file, callers naturally rediscover
  basis-level rebuilding instead of reusing source-level objects.
- Timing/reporting code living next to kernel code makes it harder to profile
  or simplify the hot path without touching unrelated behavior.

This note does not claim these issues are the only bottlenecks. It records the
places where structure is now plausibly also a performance tax.

## Proposed Decomposition

### `src/cartesian_nested_faces.jl`

Proposed destination files:

1. `src/cartesian_nested_common.jl`
   - shared nested types
   - coefficient-map storage helpers
   - axis-bundle helpers
   - shared interval/box geometry utilities
   - contract-audit helpers

2. `src/cartesian_nested_reporting.jl`
   - `NestedFixedBlockBuildTimingSummary`
   - timing collector/report helpers
   - diatomic doside trace structs and writers
   - narrow text-report helpers

3. `src/cartesian_nested_packets.jl`
   - compact packet structs
   - factorized-direct and support-reference packet kernels
   - weighted term-sum helpers
   - `_nested_fixed_block(...)` adapters

4. `src/cartesian_nested_sequences.jl`
   - generic doside / face / edge / shell / sequence construction
   - complete-shell retention contracts
   - `_nested_complete_rectangular_shell(...)`
   - `_nested_shell_sequence_from_core_block(...)`

5. `src/cartesian_nested_atomic.jl`
   - one-center atomic shell-sequence builders
   - one-center atomic fixed-block front ends
   - one-center atomic structure diagnostics/report
   - legacy-profile atomic helpers

6. `src/cartesian_nested_diatomic.jl`
   - diatomic split geometry
   - diatomic nonuniform core logic
   - diatomic source builder
   - diatomic-specific diagnostic/trace glue

7. `src/cartesian_nested_experimental_geometries.jl`
   - homonuclear chain source types and builders
   - square-lattice source types and builders
   - experimental geometry diagnostics/report helpers

Ownership rule:

- generic shell/packet math belongs in the common/packet/sequence files
- public atomic/diatomic/experimental front ends belong in their dedicated path
  files
- reporting belongs outside the kernel files unless the report is intrinsic to
  a path object

### `src/ordinary_qiu_white_rg.jl`

Proposed destination files:

1. `src/ordinary_qw_types_and_bases.jl`
   - `OrdinaryCartesianOperators3D`
   - `QWRGResidualSpaceDiagnostics`
   - QW basis container structs
   - basis coordinate helpers and public basis constructors
   - basic chain/square basis diagnostics

2. `src/ordinary_qw_nested_frontends.jl`
   - `_resolved_nested_term_coefficients`
   - public diatomic fixed-source/fixed-block/diagnostics
   - public chain/square nested diagnostics wrappers
   - nested-source/fixed-block front-door helpers

3. `src/ordinary_qw_residuals.jl`
   - residual keep policy
   - residual stabilization
   - residual-space diagnostics and moment extraction
   - `diagnose_qwrg_residual_space`

4. `src/ordinary_qw_raw_blocks.jl`
   - atomic axis/block helpers
   - diatomic shell/raw-block helpers
   - supplement primitive/block helpers
   - 1D/3D raw matrix assembly used by operator construction

5. `src/ordinary_qw_operators.jl`
   - `assembled_one_body_hamiltonian`
   - `ordinary_cartesian_vee_expectation`
   - `ordinary_cartesian_1s2_check`
   - ordinary/direct operator assembly
   - nested fixed-block operator assembly
   - supplement-bearing diatomic operator assembly

6. `src/ordinary_qw_experimental_paths.jl`
   - experimental chain/square consumer path wrappers
   - narrow report/timing helpers that are path-specific rather than general

Ownership rule:

- basis/front-door construction should be separate from raw block math
- residual-space analysis should not sit inside the same file segment as final
  operator assembly
- experimental chain/square path wrappers should not live in the same file area
  as the core diatomic and atomic production front ends

## Dependency Shape After The Split

### Nested Cartesian Core

The intended dependency order is:

1. `cartesian_nested_common.jl`
2. `cartesian_nested_reporting.jl`
3. `cartesian_nested_packets.jl`
4. `cartesian_nested_sequences.jl`
5. `cartesian_nested_atomic.jl`
6. `cartesian_nested_diatomic.jl`
7. `cartesian_nested_experimental_geometries.jl`

`bond_aligned_diatomic_geometry.jl` and
`bond_aligned_diatomic_geometry_export.jl` should continue to depend only on
the diatomic source/geometry types, not on atomic or experimental chain/square
helpers.

### Ordinary / QW Layer

The intended dependency order is:

1. `ordinary_qw_types_and_bases.jl`
2. `ordinary_qw_nested_frontends.jl`
3. `ordinary_qw_residuals.jl`
4. `ordinary_qw_raw_blocks.jl`
5. `ordinary_qw_operators.jl`
6. `ordinary_qw_experimental_paths.jl`

The operator file should depend on the residual/block files, not the other way
around.

## Highest-Risk Coupling Points

These are the places where careless file motion is most likely to cause drift:

- `_NestedFixedBlock3D` is consumed by:
  - ordinary/QW operator assembly
  - cartesian basis representation
  - geometry export
  - many tests
- `_CartesianNestedBondAlignedDiatomicSource3D` is consumed by:
  - diatomic geometry export
  - diatomic diagnostics
  - source reuse tests
- the current test suite reaches into internal nested structs directly, so
  type/file moves must preserve include order and names until tests are
  deliberately narrowed
- chain/square experimental exports rely on path structs assembled in
  `ordinary_qiu_white_rg.jl` but on source types assembled in
  `cartesian_nested_faces.jl`
- ordinary operator assembly depends on nested fixed-block fields that are
  intentionally compact-only now; that contract must not drift during motion

## Pure File Motion Vs Logic Refactor

### Mostly Pure File Motion

- timing/reporting extraction from `cartesian_nested_faces.jl`
- one-center atomic helper extraction from `cartesian_nested_faces.jl`
- chain/square experimental geometry extraction from both source files
- basis/type/front-end extraction from `ordinary_qiu_white_rg.jl`

### Small Logic Refactor Required

- packet-kernel extraction, because packet setup and packet assembly are
  currently interleaved with shared sequence code
- ordinary operator split, because final operator assembly currently sits next
  to the raw-block helpers it consumes
- residual-analysis split, because keep-policy/stabilization helpers are
  interwoven with the routes that use them

## Repo-Doer Chunk Plan

### Chunk 1: Extract nested timing and reporting helpers

- Files to edit:
  - `src/cartesian_nested_faces.jl`
  - new `src/cartesian_nested_reporting.jl`
  - `src/GaussletBases.jl`
- Preserve:
  - timing labels
  - `nested_fixed_block_timing_report(...)`
  - diatomic doside trace text format
- Minimal tests/smokes:
  - `GAUSSLETBASES_TEST_GROUPS=nested julia --project=. test/runtests.jl`
    for:
    - `One-center atomic fixed-block timing surface`
  - `GAUSSLETBASES_TEST_GROUPS=diatomic julia --project=. test/runtests.jl`
    for:
    - `Bond-aligned diatomic doside / COMX trace diagnostics`
    - `Bond-aligned diatomic plane projection export`
- Do not touch:
  - packet math
  - atomic/diatomic geometry policy

### Chunk 2: Isolate one-center atomic nested core

- Files to edit:
  - `src/cartesian_nested_faces.jl`
  - new `src/cartesian_nested_atomic.jl`
  - `src/GaussletBases.jl`
- Preserve:
  - public atomic shell-sequence/fixed-block names
  - legacy-profile behavior
  - compact-only packet contract
- Minimal tests/smokes:
  - `One-center atomic full-parent nested contract`
  - `One-center atomic legacy-profile nested contract`
  - `One-center atomic compact fixed-block term storage`
  - `One-center atomic factorized direct packet kernel`
- Do not touch:
  - diatomic source builders
  - chain/square experimental sources

### Chunk 3: Isolate diatomic geometry and source logic

- Files to edit:
  - `src/cartesian_nested_faces.jl`
  - new `src/cartesian_nested_diatomic.jl`
  - `src/GaussletBases.jl`
  - only include-order users if required
- Preserve:
  - one-build source reuse contract
  - current split/no-split behavior
  - compact packet path
  - geometry export compatibility
- Minimal tests/smokes:
  - `Bond-aligned diatomic nested source reuse path`
  - `Bond-aligned diatomic split geometry`
  - `Bond-aligned diatomic nested fixed block`
  - `Bond-aligned diatomic compact nested fixed-block contract`
  - `Bond-aligned diatomic raw source geometry and 3d export`
  - `Bond-aligned diatomic shared-shell odd-retain experiment`
- Do not touch:
  - chain/square split policy
  - operator formulas in `ordinary_qiu_white_rg.jl`

### Chunk 4: Isolate experimental chain/square nested producers

- Files to edit:
  - `src/cartesian_nested_faces.jl`
  - `src/ordinary_qiu_white_rg.jl`
  - new `src/cartesian_nested_experimental_geometries.jl`
  - new `src/ordinary_qw_experimental_paths.jl`
  - `src/GaussletBases.jl`
- Preserve:
  - chain/square exported names
  - dense export payloads
  - current experimental status and behavior
- Minimal tests/smokes:
  - `Bond-aligned homonuclear chain nested geometry diagnostics`
  - `Experimental bond-aligned homonuclear chain nested QW consumer path`
  - `Axis-aligned homonuclear square-lattice nested geometry diagnostics`
  - `Experimental axis-aligned homonuclear square-lattice nested QW consumer path`
- Do not touch:
  - production atomic and diatomic path semantics

### Chunk 5: Isolate QW basis/front-end layer

- Files to edit:
  - `src/ordinary_qiu_white_rg.jl`
  - new `src/ordinary_qw_types_and_bases.jl`
  - new `src/ordinary_qw_nested_frontends.jl`
  - `src/GaussletBases.jl`
- Preserve:
  - public basis constructors
  - public diatomic source/fixed-block/diagnostics names
  - one-build diatomic source API
- Minimal tests/smokes:
  - `Bond-aligned diatomic QW reference path`
  - `Bond-aligned diatomic nested source reuse path`
  - `Bond-aligned homonuclear chain ordinary QW reference path`
  - `Axis-aligned homonuclear square-lattice ordinary QW reference path`
- Do not touch:
  - raw block formulas
  - residual keep-policy semantics

### Chunk 6: Isolate residual diagnostics and stabilization

- Files to edit:
  - `src/ordinary_qiu_white_rg.jl`
  - new `src/ordinary_qw_residuals.jl`
  - `src/GaussletBases.jl`
- Preserve:
  - `QWRGResidualSpaceDiagnostics`
  - stabilization behavior
  - public `diagnose_qwrg_residual_space`
- Minimal tests/smokes:
  - `QW residual-space keep policy is near-null-only and stabilized`
  - `One-center atomic legacy-profile residual completion contract`
  - `Atomic residual keep policy rejects relative_case_scale on public QW routes`
  - `One-center atomic ns=9 legacy-profile residual stabilization closes center-extraction failure`
- Do not touch:
  - front-end basis/source APIs
  - ordinary operator naming surface

### Chunk 7: Isolate raw block builders and final operator assembly

- Files to edit:
  - `src/ordinary_qiu_white_rg.jl`
  - new `src/ordinary_qw_raw_blocks.jl`
  - new `src/ordinary_qw_operators.jl`
  - `src/GaussletBases.jl`
- Preserve:
  - `OrdinaryCartesianOperators3D`
  - `assembled_one_body_hamiltonian`
  - compact nested fixed-block contract
  - current per-center nuclear-term behavior
- Minimal tests/smokes:
  - `Atomic lmax=0 supplement uses the explicit 3D shell route in QW consumers`
  - `Active atomic lmax=1 supplement is explicit and physical in QW routes`
  - `Qiu-White residual Gaussian reference path`
  - `Per-center nuclear one-body reassembly on diatomic routes`
  - `Bond-aligned diatomic nested QW consumer path`
  - `Bond-aligned diatomic molecular supplement ordinary QW path`
- Do not touch:
  - nested geometry policy
  - experimental chain/square split policy

## Recommended Execution Order

Use this order:

1. reporting/timing extraction
2. one-center atomic extraction
3. diatomic extraction
4. experimental chain/square extraction
5. QW basis/front-end extraction
6. residual extraction
7. raw block/operator extraction

That order keeps the compact production paths visible early and postpones the
highest-coupling operator refactor until the front-end seams are cleaner.

## Review Discipline

For every chunk:

- prefer pure file motion first
- preserve names and behavior before narrowing anything
- keep compact production paths easy to find
- do not silently widen public APIs
- do not claim performance wins without measurement

The problem is not solved when the file count increases. The problem is solved
only when ownership boundaries become clear enough that later optimization and
maintenance stop fighting the file shape.
