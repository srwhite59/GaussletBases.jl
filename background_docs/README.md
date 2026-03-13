## Gausslets Repo Background Packet

This folder is an internal background packet for future Codex work on the public-facing `Gausslets` repository.

It is not a release directory. Some files here are draft manuscripts, internal reports, or code snapshots copied for orientation.

### Suggested Reading Order

1. `White_2017_Hybrid_gridbasis_set_discretizations.pdf`
2. `White_Lindsey_2023_Nested_Gausslet_Basis_Sets.pdf`
3. `Radial_Gausslets_draft.pdf`
4. `GAUSSLETS_REPO_CONTEXT.md`
5. `RADIAL_CODE_MAP.md`
6. The copied radial support reports

### What Is In This Packet

- Two canonical gausslet papers:
  - White 2017
  - White and Lindsey 2023
- Current radial-gausslet manuscript draft:
  - PDF
  - TeX source snapshot
- Internal notes:
  - `GAUSSLETS_REPO_CONTEXT.md`
  - `RADIAL_CODE_MAP.md`
- Radial support reports:
  - `Report10.12.25.md`
  - `ReportHamconstruct.md`
  - `ReportgirdboundaryHatomv2.md`
- Code snapshots:
  - `code_snapshots/RadialGGrid.jl`
  - `code_snapshots/DiagYlm.jl`
  - `code_snapshots/atombasisYlmopt.jl`

### Why These Files

The goal is to give a future Codex enough context to answer:

- what the original gausslet line established
- what the current radial-gausslet paper is trying to say
- how the active radial code is organized now
- which parts are likely candidates for public extraction

### Important Scope Note

The copied radial code snapshots mainly represent the current partial-IDA / Ylm producer line.

The active full-IDA shell-local angular producer still lives outside this packet, in the older tree:

- `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/Radial/testatom/sphgatomps.jl`

That file is referenced in the notes, but not copied here because it is still tied to the older experimental tree rather than the cleaned `work/radial` location.
