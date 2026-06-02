# Be Atom Supplement Source Facts Checkpoint - 2026-06-02

This checkpoint records the private Be atom supplement target facts gathered
for the source-box component lane. It is a developer checkpoint only. It does
not make any construction path public and does not adopt QW, Hamiltonian, or
packet routing for the source-box route.

## Target

- Atom: Be.
- Nuclear charge: `Z = 4`.
- Electron model: all electron, no ECP.
- Geometry: one center at the origin, while preserving one-center
  by-center/component bookkeeping.
- q ladder: `q = 4, 5, 6`, with `q = 7` shape-only.
- Initial uncontracted core side rule:
  - odd q: `core_side = q`;
  - even q: `core_side = q + 1`.
- Distortion rule:
  - `dZ = 1.2 / (core_side - 3)`;
  - `d = dZ / Z`;
  - `s = white_lindsey_atomic_mapping(...).s`.
- Parent target radius: `10.0`.
- Reference spacing: `1.0`.
- Tail spacing: `10.0`.

The accepted fixed-only dimensions are:

| q | core side | parent side | fixed dimension |
|---:|---:|---:|---:|
| 4 | 5 | 15 | 615 |
| 5 | 5 | 15 | 615 |
| 6 | 7 | 23 | 2087 |
| 7 | 7 | 23 | 2087 shape-only |

## Supplement Data

The Be supplement target uses `legacy_atomic_gaussian_supplement("Be",
"cc-pV5Z"; lmax = 5, basisfile = "/Users/srw/BasisSets")`.

The explicit `basisfile` is a machine-local data source. The repo-vendored
legacy basis data is not sufficient for this Be `cc-pV5Z` target, so no
tracked test should require this file. The loaded supplement contains 126
Cartesian shell orbitals.

## Source Facts

The private probe uses existing source/facts front doors and direct residual
inputs only:

- existing `cartesian_qw_operator_build_source(fixed_block, supplement; ...)`
  for source metadata/front-door validation;
- direct fixed/GTO overlap and supplement/supplement overlap construction;
- existing near-null residual-space construction.

The build-source front door resolves the gausslet backend to
`:pgdg_localized_experimental`, with source input kind
`:atomic_nested_fixed_block_input`.

Accepted constructed rows:

| q | fixed/GTO overlap | supplement overlap | residual kept | residual discarded | raw-to-final shape | final/source dim | final overlap error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | `(615, 126)` | `(126, 126)` | 125 | 1 | `(741, 740)` | 740 | `4.905e-10` |
| 5 | `(615, 126)` | `(126, 126)` | 125 | 1 | `(741, 740)` | 740 | `4.905e-10` |
| 6 | `(2087, 126)` | `(126, 126)` | 122 | 4 | `(2213, 2209)` | 2209 | `9.179e-10` |

The q7 row is shape-only:

| q | expected fixed/GTO overlap | expected supplement overlap | source/residual construction |
|---:|---:|---:|---|
| 7 | `(2087, 126)` | `(126, 126)` | not run |

## Component-Route Diagnostics

The same ignored private probe now records component-route diagnostics for the
Be atom supplement target. These diagnostics assemble private matrices for
source/final-space checks only. They do not adopt a Hamiltonian, packet, QW,
or public/default route.

The retained component inventory is:

- a fixed carried component from the nested all-electron Be core;
- a raw GTO supplement component from the `cc-pV5Z` Cartesian shell data;
- a near-null residual final supplement component after overlap residual-space
  filtering;
- constructed q4/q5/q6 rows and q7 shape-only metadata.

Accepted component layout rows:

| q | status | fixed raw/final range | raw GTO range | residual final range | raw-to-final shape | final dim |
|---:|---|---|---|---|---:|---:|
| 4 | constructed | `1:615` | `616:741` | `616:740` | `(741, 740)` | 740 |
| 5 | constructed | `1:615` | `616:741` | `616:740` | `(741, 740)` | 740 |
| 6 | constructed | `1:2087` | `2088:2213` | `2088:2209` | `(2213, 2209)` | 2209 |
| 7 | shape-only | not constructed | expected 126 GTOs | not constructed | not constructed | not constructed |

### Overlap Component Consumer

The overlap-only component consumer builds the raw fixed/GTO/GTO component
overlap, contracts it with `raw_to_final`, and compares the result to the
source final overlap. The source-final comparison error is exactly `0.0` for
q4/q5/q6.

| q | raw overlap shape | raw-to-final shape | final overlap shape | source-final comparison error | max final overlap error |
|---:|---:|---:|---:|---:|---:|
| 4 | `(741, 741)` | `(741, 740)` | `(740, 740)` | `0.0` | `4.905e-10` |
| 5 | `(741, 741)` | `(741, 740)` | `(740, 740)` | `0.0` | `4.905e-10` |
| 6 | `(2213, 2213)` | `(2213, 2209)` | `(2209, 2209)` | `0.0` | `9.179e-10` |

### Kinetic Component Consumer

The kinetic component consumer assembles the fixed/fixed, fixed/GTO, and
GTO/GTO kinetic blocks, contracts through `raw_to_final`, and compares against
the existing final one-body block authority for the kinetic-only input.

| q | final kinetic shape | authority | authority error |
|---:|---:|---|---:|
| 4 | `(740, 740)` | `_qwrg_final_one_body_matrix_from_blocks` | `3.893e-10` |
| 5 | `(740, 740)` | `_qwrg_final_one_body_matrix_from_blocks` | `3.893e-10` |
| 6 | `(2209, 2209)` | `_qwrg_final_one_body_matrix_from_blocks` | `9.131e-10` |

### By-Center Nuclear Component Consumer

The by-center nuclear-attraction component consumer uses the physical
electron-nuclear sign convention: `-Z` times the positive one-center Gaussian
factor expansion, with `Z = 4` for Be, so the explicit scale is `-4`.

The source front door used by the probe still requests total-only nuclear term
storage. The direct component diagnostic keeps the one-center by-center
nuclear attraction matrix separate and compares it to
`_qwrg_final_nuclear_one_body_by_center`.

| q | final nuclear shape | scale | authority error |
|---:|---:|---:|---:|
| 4 | `(740, 740)` | `-4` | `1.356e-9` |
| 5 | `(740, 740)` | `-4` | `1.356e-9` |
| 6 | `(2209, 2209)` | `-4` | `1.104e-9` |

### Combined One-Body Diagnostic

The combined one-body diagnostic forms
`H1 = kinetic + one-center nuclear attraction` in the fixed + residual final
space. This is private diagnostic output only. It is not a Hamiltonian route
adoption and is not an SCF target.

The comparison authority is the total-only
`_qwrg_final_one_body_matrix_from_blocks` path. The by-center nuclear
component remains separately reported, so the authority storage caveat is:
total-only one-body authority compared while the by-center Be nuclear
component is retained separately in diagnostics.

| q | H1 shape | finite | symmetry error | authority error | kinetic inf norm | nuclear inf norm | H1 inf norm |
|---:|---:|---|---:|---:|---:|---:|---:|
| 4 | `(740, 740)` | true | `0.0` | `1.625e-9` | `329.741` | `63.888` | `290.423` |
| 5 | `(740, 740)` | true | `0.0` | `1.625e-9` | `329.741` | `63.888` | `290.423` |
| 6 | `(2209, 2209)` | true | `0.0` | `8.669e-9` | `850.102` | `127.278` | `763.346` |

Observed H1 assembly and authority-comparison timings from the private probe
were:

| q | H1 assembly seconds | H1 assembly bytes | authority seconds | authority bytes |
|---:|---:|---:|---:|---:|
| 4 | `0.004287041` | `17563968` | `0.025175918` | `40143328` |
| 5 | `0.004464307` | `17563968` | `0.068722806` | `40143328` |
| 6 | `0.094632635` | `156434752` | `0.177662453` | `374143456` |

### Final-Residual MWG Interaction Component Diagnostic

The final-residual MWG interaction diagnostic follows the existing ordinary
MWG authority shape: raw GTO supplement rows are residual-construction inputs,
not final interaction rows. The fixed/fixed block is taken from the existing
fixed IDA path, and fixed/residual plus residual/residual blocks are built
from effective residual Gaussian centers and widths through
`_qwrg_mwg_interaction_components`.

The comparison authority is the existing final interaction matrix from
`ordinary_cartesian_qiu_white_operators(fixed_block, supplement;
interaction_treatment = :mwg)`. The component assembly matches that authority
exactly on the private q4/q5/q6 probe.

The component extraction is also available as the private helper
`_qwrg_final_residual_mwg_component_blocks(...)`. It returns the fixed/fixed,
fixed/residual, residual/residual, and assembled final interaction blocks plus
diagnostics. The older full-matrix helper,
`_qwrg_fixed_block_interaction_matrix_mwg(...)`, delegates to this component
helper and preserves its existing final-matrix authority behavior.

| q | fixed dim | raw GTO dim | residual count | final dim | raw-to-final shape | fixed/fixed | fixed/residual | residual/residual | authority error |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 615 | 126 | 125 | 740 | `(741, 740)` | `(615, 615)` | `(615, 125)` | `(125, 125)` | `0.0` |
| 5 | 615 | 126 | 125 | 740 | `(741, 740)` | `(615, 615)` | `(615, 125)` | `(125, 125)` | `0.0` |
| 6 | 2087 | 126 | 122 | 2209 | `(2213, 2209)` | `(2087, 2087)` | `(2087, 122)` | `(122, 122)` | `0.0` |

The q7 row remains shape-only: expected fixed dimension `2087`, expected raw
GTO dimension `126`, and expected raw dimension `2213`.

Residual center and width checks passed for q4/q5/q6:

- residual centers are finite;
- residual widths are finite;
- residual widths are positive;
- final interaction symmetry error is `0.0`.

Owner-index caveat: the current one-center atomic authority reports owner
index `0` for residual rows. The private component diagnostic may separately
record logical one-center owner `1` metadata for reporting, but that is
report-only and must not be generalized to multi-center owner semantics.

The diagnostic explicitly does not build raw GTO/GTO or fixed/raw-GTO MWG
interaction blocks. Raw GTO rows remain residual-construction inputs only.

Durable tracked ordinary coverage now lives in
`test/ordinary/mwg_residual_component_helper_runtests.jl`. That test uses a
small in-repo one-center fixture, checks the helper output against both the
delegating full-matrix helper and the existing ordinary final MWG authority,
and asserts the same private-route boundaries: caller-owned residual metadata
without generalized owner semantics, no raw GTO/GTO or fixed/raw-GTO MWG
blocks, raw GTO rows as residual-construction inputs only, and no public,
default, packet, fixed-block, QW/Hamiltonian, Be2/CR2, ECP, SCF/HF, or
MWG/IDA semantic adoption.

Observed q6 final-residual MWG interaction timing/allocation rows were:

| q | phase | elapsed seconds | allocated bytes |
|---:|---|---:|---:|
| 6 | `existing_final_mwg_authority` | `9.003502754` | `9226055184` |
| 6 | `raw_cartesian_gto_blocks` | `4.056610367` | `5980174752` |
| 6 | `residual_space` | `0.498724446` | `582154288` |
| 6 | `residual_mwg_provenance` | `1.289357810` | `505031696` |
| 6 | `mwg_component_extraction` | `0.827999912` | `1089336624` |
| 6 | `final_residual_mwg_assembly` | `0.054681351` | `117117168` |

These timings are private diagnostic scale signals, not production readiness
or public-route performance claims.

## Timing And Allocation Caution

The expensive diagnostic phase is raw Cartesian/GTO block construction. The
private probe timing rows observed:

| q | phase | elapsed seconds | allocated bytes |
|---:|---|---:|---:|
| 4 | `raw_cartesian_gto_blocks` | 2.723 | 3977100480 |
| 5 | `raw_cartesian_gto_blocks` | 2.637 | 3977099520 |
| 6 | `raw_cartesian_gto_blocks` | 3.852 | 5988878752 |

The q6 residual-space phase also allocated about `573696048` bytes. These
numbers are acceptable for a private source/facts probe, but they are not a
production-route performance claim.

## Non-Goals

This checkpoint does not:

- build QW operators;
- build a Hamiltonian or SCF target;
- build raw GTO/GTO or fixed/raw-GTO MWG interaction blocks;
- adopt packet/fixed-block/QW/Hamiltonian routing for source-box work;
- change public/default route behavior;
- introduce ECP behavior;
- move to Be2;
- make a CR2 or science-quality claim;
- change MWG/IDA semantics;
- assign retained PQS/source-box columns positive quadrature weights;
- use retained-weight or IDA division.

The final-residual MWG component diagnostic is now mechanically available for
one-center Be q4/q5/q6 as a private authority-matched probe. It remains
separate from the IDA source-box electron-electron lane and does not promote
MWG supplement/residual coupling to a public/default, Hamiltonian, QW, packet,
fixed-block, or multi-center production route.
