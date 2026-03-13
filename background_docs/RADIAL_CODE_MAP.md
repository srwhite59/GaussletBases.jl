# Radial Code Map

This note is a compact map of the current radial-side code that is most useful for understanding how a future public gausslet / radial-gausslet repository should be organized.

## 1. `RadialGGrid.jl`

Current active path:

- `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/radial/RadialGGrid.jl`

Real source file at present:

- `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/RadialGGrid.jl`

### Role

`RadialGGrid.jl` is the core half-line radial module. It is the most important radial code file in the current tree.

It provides:

- the boundary-gausslet basis on `r >= 0`
- the mapped quadrature grid
- radial one-electron matrices
- radial Coulomb multipole tables `Veel[L]`
- X-diagonalization / center ordering support
- the main `RadialData` payload used by higher-level drivers

### Important exported entry points

- `RadialData`
- `radial_moments`
- `xdiag_finalize`
- `make_erf_mapped_grid`
- `build_radial_and_one_electron`
- `build_Veel_L`
- `build_Veel_L_cheb`
- `build_Veel_tables`
- `radial_states`

### Why it matters for a public repo

If a public radial-gausslet package is made, `RadialGGrid.jl` is one of the strongest candidates for the architectural center of that repo.

## 2. `DiagYlm.jl`

Current active path:

- `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/radial/DiagYlm.jl`

### Role

This is the angular-coupling layer for the Ylm-based atom-Hamiltonian pipeline.

It provides:

- angular Coulomb kernel construction
- J/K contraction machinery
- Ylm indexing helpers

This is important because it shows how the radial Coulomb multipoles from `RadialGGrid.jl` are turned into usable angular-coupled operators.

## 3. `atombasisYlmopt.jl`

Current active path:

- `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/radial/atombasisYlmopt.jl`

### Role

This is the main partial-IDA / HamV6 producer path.

It:

- builds the radial basis through `RadialGGrid`
- builds angular kernels through `DiagYlm`
- assembles one- and two-electron components
- writes a `HamV6` JLD2 file for later DMRG consumption

This file is especially useful as an integration example.

## 4. Other helpful files in the active radial tree

### `make_boundary_gausslets.jl`

Builds or refreshes the boundary basis cache used by `RadialGGrid.jl`.

### `solveHatom.jl`

Useful for hydrogenic checks and basic sanity validation.

### `HFatom.jl`

Current closed-shell RHF driver on the radial-plus-Ylm basis, though some older script structure is still mixed in.

### `He1p.jl`

Small-atom / one-particle testing and diagnostics.

## 5. Reports copied into this packet

### `Report10.12.25.md`

Best consolidated overview of the VeeYlm / Heatom / HFatom pipeline.

### `ReportHamconstruct.md`

Best concise explanation of the radial-plus-Ylm operator assembly and the two-index integral-diagonal approximation logic.

### `ReportgirdboundaryHatomv2.md`

Best implementation-and-validation note for boundary gausslets, mapped grids, and hydrogen tests.

## 6. Important code that is not copied here

The current full-IDA shell-local angular producer is still in the older test tree:

- `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/Radial/testatom/sphgatomps.jl`

It is relevant to current research and DMRG work, but it is not yet part of the cleaned public-extraction story.

## 7. Public-extraction guidance

The likely first-pass public pieces are:

- core gausslet functionality
- radial boundary-basis machinery
- clean radial operator construction

The more experimental atom-Hamiltonian, DMRG, and full-IDA shell-local code should probably be treated as later-stage or separate extraction work.
