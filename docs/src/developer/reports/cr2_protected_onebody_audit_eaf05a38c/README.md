# Cr2 Protected Fixed-Sector One-Body Audit

Status: measurement evidence for `HP-RG-PROTECT-ONEBODY-AUDIT-01`.

Raw output:

```text
/Users/srw/dmrgtmp/cr2_protected_onebody_audit_eaf05a38c/
```

## Purpose

This audit tested whether the source-backed staged protected-original geometry
can receive exact one-body operators in the injected fixed sector without
changing production Hamiltonian construction.

The geometry is:

```text
M = [G, R_compact]
Z = [Z_protected, Z_broad]
F = [Z, M Qperp]
```

The probe formed in-memory:

```text
F' K F
F' U_A F by nuclear center
F' H1 F
```

No production source, artifact schema, IDA/MWG interaction, public driver,
or Cr2 HF workflow was changed by the audit.

## Geometry Match

The source-backed geometry exactly matched the staged-filter target:

| quantity | value |
| --- | ---: |
| source commit | `eaf05a38c` |
| base dimension | 6915 |
| compact residual dimension | 30 |
| protected originals | 30 |
| broad retained originals | 87 |
| `Z` dimension | 117 |
| `F` dimension | 6945 |
| `B_min` | 0.993465824505872 |
| `B < 0.99` | 0 |
| source `B_min` delta | 0.0 |
| source fake-trace delta | 0.0 |

Fixed-sector orthogonality diagnostics:

```text
F' S F - I block estimate = 1.1641532182693481e-9
Z' S M Qperp              = 9.872889823506115e-16
protected span min sv     = 0.9999999999999996
protected identity max    = 5.10702591327572e-15
```

## One-Body Diagnostics

Symmetry samples:

```text
K sample symmetry       = 1.9679191609611735e-10
max U_A sample symmetry = 3.7439938302408216e-12
H1 sample symmetry      = 6.59525767332525e-11
```

Per-center uncharged nuclear diagnostics:

| center | charge | trace | charged trace | symmetry sample | low min |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 24.0 | -21314.69035137109 | -511552.5684329062 | 3.7439938302408216e-12 | -310.3713251562225 |
| 2 | 24.0 | -21314.69018160351 | -511552.56435848423 | 1.9693136010801027e-12 | -310.37119254003386 |

Trace diagnostics:

```text
trace K_FF        = 4.0271122908679575e7
trace U_total_FF  = -1.0231051327913904e6
trace H1_FF       = 3.9248017775889225e7
```

Converged low `H1_FF` Ritz values:

```text
-295.5569201561
-295.5569039890
-79.8649001575
-79.8648899852
-79.5459243157
-79.5453658961
-79.2688084735
-79.2687478275
```

The matching compact-main `H1_M` low values also converged, with minimum
`-295.5569819352`.

Lowest `H1_FF` modes are dominated by protected originals, not broad injected
directions:

| mode | eigenvalue | protected weight | broad-Z weight | `Qperp` weight |
| ---: | ---: | ---: | ---: | ---: |
| 1 | -295.5569201561 | 0.9997640435 | 0.0001804204 | 0.0000555362 |
| 2 | -295.5569039890 | 0.9997641639 | 0.0001612109 | 0.0000746252 |
| 3 | -79.8649001575 | 0.9679441832 | 0.0225794079 | 0.0094764089 |
| 4 | -79.8648899852 | 0.9688254563 | 0.0224326869 | 0.0087418568 |

The audit did not find an obvious low-`H1` broad injected mode before any
IDA/MWG interaction design.

## Caveat

The `K_FF` low Ritz values and `M Qperp` `H1` subblock low Ritz values hit
Krylov convergence limits. They are diagnostic only. The converged `H1_FF`
and `H1_M` checks are the meaningful positive one-body evidence.

## Design Consequence

The measurement supports a narrow source lane for protected-original exact
one-body transformation. It does not support IDA/MWG interaction transforms,
artifact writing, public driver/API changes, Cr2 HF, residual default changes,
or treating rejected broad directions as MWG residuals.
