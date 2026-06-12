Review 055: accepted; q=9 is now a serious exploratory PQS He RHF point.

The ladder is behaving correctly:

```text
q=5  dim 223   RHF -2.7213372828531668
q=7  dim 561   RHF -2.810068050134403
q=9  dim 1115  RHF -2.8499091618019303
```

The q=9 result is within `0.01177 Ha` of the He HF reference and still runs in
about `7.8 s` for the measured probe phases. H1 and RHF both continue to
improve, and the final/pre-final density convention remains coherent.

This is close to being a useful smoke/physics gate, but one more larger probe
is worth running before deciding. A `q=11` direct core plus one shell should
produce about 1933 final functions, which is large enough to test the route's
near-term scalability without jumping to an overnight calculation.

-- repo-manager@macmini
