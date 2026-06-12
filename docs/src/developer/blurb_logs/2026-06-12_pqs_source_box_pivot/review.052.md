Review 052: accepted; the pre-final density seam is the right bridge into RHF.

The pass made the density-interaction convention concrete without promoting the
bad alternatives:

```text
accepted gauge:      localized pre-final positive weights
final orbital map:   c_prefinal = combined_lowdin_cleanup * c_final
H1 J:                0.6397851751855723
fixed oracle J:      0.6397857768997106
oracle delta:       -6.0e-7
```

The important point is that the final orthonormal basis remains the one-body
solve basis, while the electron-electron diagnostic is evaluated after mapping
the occupied orbital into the pre-final density gauge. The signed final-weight
division and raw no-division variants remain rejected diagnostic values, not
alternate contracts.

The helper returns a report-shaped object, which is acceptable for this narrow
seam, but do not let it grow into a broad route result. The next step should be
an exploratory complete PQS He RHF probe that consumes this object directly and
keeps the final/pre-final density transformation explicit.

-- repo-manager@macmini
