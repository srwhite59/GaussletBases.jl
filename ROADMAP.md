# Roadmap

This roadmap is a short note about the **current pressure points** in the repo.

It is not a schedule. It is the answer to:

- what matters most next, given the current code shape

## Current center of gravity

The repo is no longer mainly an early radial-only or 1D-only project. The live
pressure is now spread across four connected fronts:

1. radial / atomic workflows that are already mature enough to support real use
2. exact Cartesian overlap / projector / transfer as workflow primitives
3. one-center nested Cartesian and bond-aligned diatomic workflow support
4. frozen-core, exactification, and larger-Z carryover questions

## Highest-value next work

### 1. Keep overlap / projector / transfer as first-class workflow primitives

The exact Cartesian transfer surface is now real enough that it should be
treated as infrastructure, not as an isolated diagnostic.

The main pressure here is:

- keep `cross_overlap`, `basis_projector`, and `transfer_orbitals` stable
- keep bundle/sidecar/export contracts aligned with those transfer primitives
- use them as the normal basis-to-basis handoff language for current Cartesian
  workflows

### 2. Deepen the one-center Cartesian line

The one-center nested Cartesian route is now real, compact-only on production
paths, and worth treating as a genuine workflow surface.

The main next questions are:

- frozen-core and exactification contracts
- larger-Z carryover
- stable diagnostics/reporting around nested fixed blocks
- clearer workflow boundaries between compact production paths and archived
  reference logic

### 3. Mature the bond-aligned diatomic workflow

The diatomic route is also now real code with supported source reuse,
diagnostics, geometry payloads, and Qiu-White consumers.

The main next questions are:

- counterpoise / ghost-basis driver support
- larger-Z diatomic carryover
- continued geometry-policy maturation and diagnostics
- supported artifact/reporting patterns for real workflow use

### 4. Carry compact-only cleanup into broader workflow contracts

The repo has now done substantial compact-only cleanup on nested production
paths. The next value is not more optional debug branching; it is keeping the
mainline small and stable.

That means:

- keeping compact summed packet contracts as the production default
- resisting reintroduction of broad slow debug-era branches
- continuing to simplify public production seams where old compatibility no
  longer has a real user

### 5. Keep larger-Z / chromium-facing work connected to the mainline

The repo now has enough one-center and diatomic infrastructure that larger-Z
and chromium-facing experiments should reuse the same:

- transfer primitives
- compact nested packet contracts
- one-build source/fixed-block patterns
- frozen-core / exactification language

This is more valuable than opening entirely separate workflow dialects.

## Important but secondary work

### Radial / atomic exactification and export

The mature radial/atomic line remains important, but the roadmap should no
longer read as if it is the only active frontier.

The next meaningful radial/atomic work is still:

- exactification of interaction structure
- better solver-facing export contracts
- narrower, cleaner atomic workflow documentation

### Experimental chain / square-lattice promotion

These routes should only be promoted further once the same compact, transfer,
and workflow-contract questions are answered cleanly there too.

For now they should remain explicitly experimental.

### Angular research

The angular line remains scientifically important, but it should continue to be
advanced as an experimental research track rather than competing for the main
public workflow story.

## Not the current roadmap center

This roadmap is **not** centered on:

- reviving the old 1D COMX-cleaned hybrid route
- treating flat note-history files as the main authority
- building a broad solver package before the current workflow primitives settle
- pretending the repo is still only a radial or 1D project

Those older framings no longer match the live code.
