# Cartesian/PQS Operational Facts

This is an operational memory sheet for agents. It is not design authority and
does not approve new source, artifact, solver, or physics-default behavior.
Use it to avoid repeating known bad probe choices.

## Standard Core Spacing

For the current standard scaled one-center/atom-local ladder, use

```text
core_spacing = 1.2 / (Z * (ns - 1))
```

unless the user or a committed driver input explicitly says otherwise.

Examples:

```text
H  ns=5  -> 0.300000
Be ns=5  -> 0.075000
Be ns=7  -> 0.050000
Ne ns=5  -> 0.030000
Ne ns=7  -> 0.020000
Cr ns=7  -> 0.008333
```

Do not call a fixed-`core_spacing` `ns` ladder a standard scaled ladder. It is
a fixed-spacing refinement check.

The mapping-strength convention is:

```text
mapping_s_standard = sqrt(Z * core_spacing)
mapping_s_effective = s_factor * mapping_s_standard
```

with `s_factor = 1.0` unless an expert scan explicitly sets it.

## Dimer Boxes And Padding

For endpoint curves, do not reuse old small-system helpers with hard-coded
`xmax_parallel = 6.0` and `xmax_transverse = 4.0`. Those are algebra-smoke
fixtures only.

The canonical driver convention is padding beyond each nucleus:

```text
xmax_parallel = R / 2 + padding
xmax_transverse = padding
```

The current driver default is `padding = 10.0` bohr. Treat this as a visible
starting point, not a universal adequacy proof. Helium can sometimes tolerate
padding/radius around `7` bohr, but other atoms and molecules should generally
start at `10` bohr or larger. A padding/radius around `20` bohr is the
conservative safe choice when endpoint accuracy or curve continuity matters.
Diffuse atoms, especially Be and Be-like endpoint work, can need these larger
boxes. For Be2 curves, scan or justify the padding instead of assuming 10 bohr
is enough.

For an `R` ladder, either:

- use fixed parent axis counts / the same parent where possible; or
- choose the box from the largest `R` plus a real per-atom padding and verify
  that no parent-count or retained-basis discontinuity contaminates the curve.

Before interpreting any endpoint curve, report:

```text
requested padding
nominal atom padding = xmax_parallel - R/2
actual parent physical bounds
actual padding beyond each nucleus
parent axis counts
final dimension
retained residual counts by channel
shell/slab topology changes
```

Reject or quarantine curve points where the parent counts, final dimension, or
retained channel identities change unless the pass is explicitly studying that
change.

## Driver Versus Ignored Helpers

Use the canonical driver or driver-style basis construction for endpoint and
curve work. The driver prints terminal due diligence, including:

- normalized system and geometry;
- spacing and box controls;
- parent physical bounds and axis counts;
- nucleus snap errors;
- shell/slab rows and retained counts.

Ignored `tmp/work` helpers are measurement history, not reusable defaults. If a
helper returns a basis with fixed `xmax_parallel` independent of `R`, assume it
is unsafe for a molecular curve until proven otherwise.

## Repo Due-Diligence Contract

The formal repo contract is `HP-DRV-SHELLDD-FN-01` /
`HP-DRV-SHELLDD-TEST-01`, recorded in:

```text
docs/src/developer/designs/cartesian_hamiltonian_producer/terminal_shellification_due_diligence.md
```

That contract requires Cartesian/PQS producer/driver workflows to expose a
bounded terminal due-diligence report before consumers interpret energies,
residual behavior, injection behavior, or high-cost production artifacts. The
required report includes:

- normalized system and geometry facts;
- validated atom locations, bond axis/length, and snapped nuclear indices;
- parent physical bounds, parent axis counts, center/spacing summaries, and
  gausslet/IDA weight statistics;
- dimension and compression accounting;
- shell-by-shell terminal order/key, role, region kind, shell index,
  owner/contact/shared classification, index and physical boxes, physical side
  lengths/aspect ratios, source-mode shape, expected aspect-balanced source
  shape, source-mode count, retained count, final column range,
  lowering/retained/realization rules, slab metadata, and advisory warning
  flags.

Warning flags are advisory by default. The due-diligence report is a review
surface, not permission to change shellification, source-mode selection,
artifacts, public inputs, numerical construction, solver workflow, or Cr2
workflow.

## Capture Due Diligence

Before blaming screened-Hartree, EGOI, IDA/MWG, or SCF for a molecular endpoint
problem, check representation quality:

1. Supplement -> parent lattice capture. This says whether the parent grid has
   enough resolution for the GTO supplement at all.
2. Supplement -> final gausslet-only contracted lattice capture. This says
   whether contraction/shellification kept supplement-relevant content before
   residuals, injection, EGOI, or HF enter.
3. Protected/reference determinant representation loss, if a packet or
   screened-Hartree reference density is involved.

Parent bad means the base lattice settings are suspect. Parent good but final
bad means contraction/shellification is losing the space. Both good pushes the
debug target toward residual selection, injection, H1/Vee transforms, IDA/MWG,
or correction bookkeeping.

## Screened-Hartree Accounting

For screened-Hartree residual-density work, `Delta_J0` and the scalar constant
`C` are represented operationally as a one-particle matrix plus a scalar, but
they are part of the screened direct electron-electron/Hartree model. Do not
count them as a physical kinetic/nuclear one-body change or as an arbitrary
energy offset.

The fitted density and fitted-potential Gaussians in atomic reference packets
are evaluation devices for the reference density field and self-energy. They
are not protected orbitals and not supplement functions.
