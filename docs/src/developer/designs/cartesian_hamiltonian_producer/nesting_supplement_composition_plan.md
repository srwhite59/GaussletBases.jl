# Nesting/Supplement Composition Plan

This note freezes the target shape for the canonical Cartesian Hamiltonian
producer after `nesting = :wl` became a real construction-family input.

It is a planning and authority-boundary amendment only. It approves no source
work by itself.

## Target Contract

The long-term supported producer shape is three successive user choices that
compose:

```text
geometry:   atom | z-axis diatomic
nesting:    :pqs | :wl
supplement: off | on
```

The driver may expose these as visible public construction choices, but it must
not encode the current partial implementation as the permanent contract. The
implementation target is a common staged construction:

```text
public system / basis / optional supplement contract
-> geometry normalization
-> nesting-specific parent and terminal-basis realization
-> common base product, unit-nuclear, IDA, and Hamiltonian assembly
-> optional residual-GTO/MWG augmentation
-> existing Hamiltonian artifact
```

`nesting` may affect parent construction, shellification, retained rules, and
terminal-basis realization. After it produces a
`CartesianTerminalBasisRealization`, downstream base operators, residual
Gaussian augmentation, IDA/MWG interaction, and artifact writing should not
branch on whether the terminal basis came from PQS or White-Lindsey.

`supplement` may affect supplement loading, residual Gaussian basis selection,
exact augmented operator transforms, MWG descriptors, residual-containing IDA
blocks, and supplemented artifact provenance. It must not resurrect route
reports, route-stage controls, pair/assembly public stages, or driver-specific
helper choreography.

## Current Matrix

Current status is intentionally explicit. Unsupported cells must fail clearly;
they must not be hidden by driver defaults or mislabeled artifact provenance.

| Geometry | Supplement | `nesting = :pqs` | `nesting = :wl` |
| --- | --- | --- | --- |
| atom | off | implemented for explicit origin-centered all-electron base atoms, with H as the committed endpoint | implemented for one-center base atoms through the WL terminal-basis seam |
| atom | on | not approved / not wired | not approved / not wired |
| z-axis diatomic | off | H2 base works; broader generic base diatomic support remains limited | blocked: no native WL diatomic terminal records |
| z-axis diatomic | on | supported for explicit homonuclear z-axis diatomics through the residual-GTO/MWG path | blocked first by missing WL diatomic base terminal records |

## Common Boundary Rules

- Atoms and diatomics must share the same producer workflow after
  geometry/shellification normalization. Atom-only Hamiltonian builders,
  atom-only materialization paths, and atom-specific artifact shapes are
  forbidden.
- PQS and White-Lindsey must converge to the same terminal-basis boundary:
  a `CartesianTerminalBasisRealization` with disjoint owned terminal supports.
- Residual-GTO/MWG supplementation consumes the terminal basis, public system
  facts, parent bundles, supplement specification, and same-construction base
  operators. It must not consume route skeletons, route reports, or old WL
  H1/H1+J materialization objects.
- Artifact provenance records the public choices (`geometry`, `nesting`, and
  supplement state) truthfully. It must not infer route identity from helper
  file names or preserve PQS labels for WL artifacts.
- The canonical driver remains compact and copyable. It can print public
  contracts and coarse physics-stage timings, but it must not grow stop-after
  controls, raw-provider switches, route diagnostics, allocation probes, or
  status/report payloads.

## Dependency Order

### 1. White-Lindsey Diatomic Base

Candidate goal: make `nesting = :wl`, `supplement = off` work for z-axis
diatomic base artifacts by producing native WL terminal records and the common
`CartesianTerminalBasisRealization`.

This should extend the WL terminal-basis seam from one-center atoms to
two-center z-axis diatomics. It should not adapt the old WL H1/H1+J
materialization path, change the canonical driver contract, or add route
diagnostics.

### 2. Supplemented Atoms

Candidate goal: make `geometry = atom`, `supplement = on` work through the same
Residual Gaussian path used by supplemented diatomics.

One-center residual selection is the one-owner case of the same owner-local
residual Gaussian algorithm. It must not introduce a separate atom supplement
algorithm, atom-specific Hamiltonian builder, or atom-only artifact schema.

### 3. Supplemented White-Lindsey

Candidate goal: after WL atom and WL diatomic base terminal bases are real,
allow residual-GTO/MWG supplementation to consume WL terminal bases through the
same Residual Gaussian/raw-block/IDA boundary used for PQS.

This lane should not start by branching the driver on WL supplemented cases.
It should first prove that the RG augmentation boundary is genuinely
nesting-neutral once supplied with a valid terminal basis and same-construction
base operators.

## Candidate IDs

These IDs are placeholders for later docs-only amendments. They do not
authorize implementation until promoted in `registry.md` with exact files,
functions, validation, forbidden surfaces, and line budgets.

- `HP-COMP-WLDIAT-FN-01` / `HP-COMP-WLDIAT-TEST-01`: WL z-axis diatomic base
  terminal-basis and artifact path.
- `HP-COMP-SUPPATOM-FN-01` / `HP-COMP-SUPPATOM-TEST-01`: supplemented
  one-center atom path through common Residual Gaussian augmentation.
- `HP-COMP-SUPPWL-FN-01` / `HP-COMP-SUPPWL-TEST-01`: supplemented
  White-Lindsey path through the common RG boundary after WL base terminal
  bases exist.

## Deferred

- translated atoms;
- non-z-axis diatomics and rotations;
- heteronuclear supplemented workflow beyond explicit future approval;
- ECP, EGOI, solver/RHF, or CR2-specific branches;
- public export/API redesign;
- new Hamiltonian wrapper or artifact matrix format;
- route diagnostics, private stop-after stages, report/status payloads, and
  broad driver switches.
