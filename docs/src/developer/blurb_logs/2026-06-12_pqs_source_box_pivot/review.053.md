Review 053: accepted as an exploratory compact-box He RHF smoke result.

The complete core/shell PQS route now reaches closed-shell RHF with the intended
split:

```text
one-electron basis: final orthonormal basis
density gauge:      localized pre-final positive-weight gauge
orbital map:        c_prefinal = combined_lowdin_cleanup * c_final
RHF total:          -2.7213372828531668
He HF error:        +0.14034271275907217
iterations:         9
```

This is a route-coherence result, not an accepted physics baseline. The compact
`5^3` direct core plus one shell has a Z=2 H1 error of about `0.152 Ha`, so the
RHF error is plausibly fixture quality rather than an immediate convention
failure. The next useful question is whether the same route scales to a larger
one-shell fixture and improves H1/J/RHF without changing the production path.

-- repo-manager@macmini
