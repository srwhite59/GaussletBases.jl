# Legacy BasisSets Provenance

GaussletBases vendors the historical legacy-style basis collection at:

- `data/legacy/BasisSets`

This is compatibility data for the existing `ReadBasis.jl`-style parser and
the Cartesian/PQS supplement path. It is not a general basis-library service.

## Snapshot Identity

The scientific body was copied on `2026-07-08` from:

- `/Users/srw/Dropbox/GaussletModules/BasisSets`

The source file's observed modification time was
`2026-04-17 22:19:26 -0700`. Its raw SHA-256 was:

```text
9c5c0e96917a88b3ccdf713696437b0094e9da60834175ff22e682dc8b90b737
```

The only normalization applied to the historical scientific body was removal
of trailing ASCII spaces and tabs from each line. Line order, blank lines, LF
line endings, and the final newline were preserved. A repo provenance header
was then prepended; no basis exponent, contraction coefficient, shell header,
basis label, or block order was changed.

The normalized scientific body begins at the first `#BASIS SET:` line. Its
SHA-256 is:

```text
b83883f4589234dd512eb760c95280854a2f42d007ab6e3533abda39a2829051
```

At the `6cfe15d0c` reconciliation baseline, the complete vendored file had
SHA-256:

```text
b40763c17cdbd8825cecf05aedbf34b35084fec7a3375661a65b49e749d33251
```

After the approved mixed-provenance header-only correction, the complete
vendored file has SHA-256:

```text
6796d34f2c813e5b627b03d003ddd29e0701cfd18e1b4c2ee7d263ce671dbede
```

The normalized scientific-body hash remains unchanged.

## Contents

The file contains `60` blocks identified by `#BASIS SET:` headers. It includes
standard named bases and custom/historical blocks. The collection has mixed
historical provenance: some standard blocks may derive from Basis Set Exchange
or published basis data, but available evidence does not establish that every
custom block came from Basis Set Exchange. Do not make that blanket claim.

License and redistribution status for the complete mixed collection is
unresolved in the current repo record. Treat it as unresolved unless verified
source/license evidence is added.

## Regression Contract

`HP-PQS-ATOMREF-PACKET-TEST-01` approves a cheap regression in
`test/misc/runtests.jl`. It computes the normalized-body SHA-256 from the first
`#BASIS SET:` line after stripping trailing ASCII spaces/tabs, and checks these
parser results:

| Block | Shells | Primitive rows |
| --- | ---: | ---: |
| H cc-pVTZ | 6 | 8 |
| H cc-pVQZ | 10 | 12 |
| Be cc-pV5Z | 21 | 42 |
| Ne cc-pV5Z | 21 | 54 |
| Cr cc-pV5Z | 32 | 434 |

The regression protects data identity and representative parser coverage. It
does not turn all 60 historical blocks into scientifically validated producer
fixtures.

## Loader Policy

The loader search order is:

1. explicit `basisfile = ...`;
2. `GAUSSLETBASES_BASISSETS_PATH`;
3. vendored repo copy `data/legacy/BasisSets`;
4. legacy fallback `~/BasisSets`.

Users needing a different collection should provide `basisfile` or
`GAUSSLETBASES_BASISSETS_PATH`. This authority does not add a basis download,
conversion, or license-resolution subsystem.
