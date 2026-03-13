# Gausslets Public Repo Context

This note is for a future Codex session working in `codexhome/repositories/Gausslets`.

## The Immediate Situation

The new `Gausslets` repository currently contains only `Gausslets.jl`. That is the right starting point, but it is not yet the full picture of the code and papers that matter.

The purpose of this packet is to supply the missing background so the new repo can be shaped with the older gausslet line and the newer radial-gausslet line in view.

## The Main Layers To Keep Straight

### 1. Core gausslet line

This is the line represented by:

- White 2017
- White and Lindsey 2023
- `Gausslets.jl`

This is the main conceptual and reusable numerical layer.

### 2. Radial half-line / radial-gausslet line

This is the line represented by:

- the current `Radial_Gausslets` draft
- `RadialGGrid.jl`
- boundary-gausslet cache generation
- radial one-electron operators
- radial Coulomb multipole tables

This is the main follow-on line that should inform any public radial-gausslet extraction.

### 3. Ylm / atom-Hamiltonian producer line

This is the line represented by:

- `DiagYlm.jl`
- `atombasisYlmopt.jl`
- the support reports copied here

This is the current partial-IDA / HamV6 producer path used to build atom Hamiltonians in a form later consumed by DMRG code.

### 4. Full-IDA shell-local angular line

This is more recent and not yet cleanly part of the public split. The important producer is:

- `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/Radial/testatom/sphgatomps.jl`

That line is highly relevant to current research, but it is not the cleanest first public-facing extraction target.

## Current Code Locations

### Public-repo staging area

- `/Users/srw/Dropbox/codexhome/repositories/Gausslets`

### Active cleaned producer-side radial tree

- `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/radial`

### Shared general-purpose module area

- `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules`

### DMRG consumer side

- `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/slicedmrgutils`

## What This Suggests For Public Repo Planning

The public repo work should probably start from the cleanest and most reusable pieces:

- core gausslet functionality
- the radial half-line basis construction
- clean operator-building interfaces

It should not start by trying to publish every current atom or DMRG path all at once.

## Why The Code Snapshots Here Matter

- `RadialGGrid.jl` is the most important radial code file because it is the central half-line basis and operator builder.
- `DiagYlm.jl` shows the angular-coupling layer used in the current radial/atomic pipeline.
- `atombasisYlmopt.jl` shows how the radial and angular pieces are assembled into a consumable Hamiltonian artifact.

Together they give a future Codex a practical map of the current software structure, even before looking at the whole `work/radial` tree.
