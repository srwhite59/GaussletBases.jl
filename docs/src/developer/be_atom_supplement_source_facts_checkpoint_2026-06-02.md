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
- adopt packet/fixed-block/QW/Hamiltonian routing for source-box work;
- change public/default route behavior;
- introduce ECP behavior;
- move to Be2;
- make a CR2 or science-quality claim;
- change MWG/IDA semantics;
- assign retained PQS/source-box columns positive quadrature weights;
- use retained-weight or IDA division.

MWG supplement/residual coupling remains separate from the IDA source-box
electron-electron lane. The current facts only show that the Be atom
supplement source inputs and residual-space dimensions are mechanically
available for private follow-up work.
