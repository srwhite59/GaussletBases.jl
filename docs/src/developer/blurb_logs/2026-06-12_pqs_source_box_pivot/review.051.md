Review 051: accepted; density must live in the localized pre-final gauge.

The probe resolved the major convention fork:

```text
pre-final weights: all positive
final weights:     signed
```

and the pre-final positive-weight contraction matches the same-geometry
fixed-block J oracle:

```text
pre-final J:         0.6397851751855723
fixed-block J:       0.6397857768997106
delta:              -6.0e-7
exact 1s reference:  0.625
```

Rejected conventions:

```text
final signed-weight division: 18.132786333403647
raw projected no-division:     0.56242162015939745
stored normalized terms:       0.05282601375665979
```

The next blocker is no longer "which weights"; it is how a final-basis orbital
coefficient vector consumes a density interaction that is defined in the
localized pre-final density gauge. The likely map is through the combined
Lowdin cleanup carried by the final-basis object, but it needs to be made
explicit and tested on the H1 self-Coulomb diagnostic before RHF.

-- repo-manager@macmini
