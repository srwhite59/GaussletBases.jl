Review 050: accepted as a blocker; do not accept signed final-gauge division.

The density/J probe correctly refused to retain a production helper. The
candidate convention:

```text
project raw numerator into the complete final basis
divide by final_ida_weights[i] * final_ida_weights[j]
```

is not physically usable for the current complete PQS final gauge:

```text
J = 18.132786333403647
reference 1s J = 0.625
fixed-block oracle J = 0.6397857768997106
```

The important clue is that the complete final gauge has signed final function
integrals. The old nested fixed-block density builders require positive
retained weights and normalize coefficient functions by those weights before
contracting raw pair terms. That is a localized density gauge assumption, not
an arbitrary orthonormal final-gauge invariant operation.

The raw projected numerator with no final-weight division gives:

```text
J = 0.56242162015939745
```

which is much more reasonable, but adopting that blindly would be another
unreviewed convention. The next pass should audit the density gauge explicitly:
compare pre-final/localized complete core-shell weights and pair contraction
against the same-geometry fixed-block oracle. Do not proceed to RHF yet.

-- repo-manager@macmini
