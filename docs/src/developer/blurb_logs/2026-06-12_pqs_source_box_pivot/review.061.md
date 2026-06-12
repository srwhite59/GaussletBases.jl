Review 061: accepted; side13 multi-layer PQS is ready for an RHF probe.

The WL-style pre-RHF checks are coherent:

```text
final dimension: 1549
Z=2 H1:         -1.9755618232013417
J:               1.2169264388860319
WL side13 H1:   -1.9748150892830352
WL side13 J:     1.2158294767735702
```

The density interaction used the intended pre-final positive-weight gauge:

```text
support weights:      all positive
pre-final weights:    all positive
pair matrix:          finite, symmetric
signed final weights: not used
raw no-division:      not used
fixed-block pairs:    not authority
```

That is enough to run RHF as the next probe. This still should not become an
acceptance gate; it is a WL-aligned physics probe for the multi-layer PQS seam.
Keep the spacing/Z rule provisional and avoid turning this into a q-only ladder.

-- repo-manager@macmini
