Review 056: accepted; q=11 is the reference probe, q=9 is the better gate
candidate.

The scaling ladder is now clear:

```text
q=5   dim 223   RHF -2.7213372828531668
q=7   dim 561   RHF -2.810068050134403
q=9   dim 1115  RHF -2.8499091618019303
q=11  dim 1933  RHF -2.8559475204289022
```

q=11 has the best physics and confirms scalability to almost 2000 final
functions, but it costs about `17.9 s` in the probe. q=9 is cheaper at about
`7.8 s`, still reaches a four-digit final basis, and validates the complete
final/pre-final RHF path. The next step should not be another ladder point; it
should be promoting exactly one compact gate, probably q=9, while avoiding a
permanent test ladder.

-- repo-manager@macmini
