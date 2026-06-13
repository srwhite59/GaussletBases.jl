Accepted pass 184.

Manager reran:

```text
julia --project=. tmp/work/pqs_he_419_wl_parent_mapping_probe.jl
```

and reproduced the reported WL-parent PQS 419 diagnostics:

```text
parent mapping = white_lindsey_atomic_mapping(Z=2,d=0.3,tail_spacing=10.0)
final dimension = 419
final overlap identity error = 3.465561171367426e-13
H1 lowest = -1.991334820314074
H1-J self-Coulomb = 1.2420423900074902
RHF total = -2.850817886618113
RHF iterations = 7
RHF residual = 2.8556726532069643e-9
```

Compared with the old WL gausslet-only 419 baseline from pass 183:

```text
WL-parent PQS RHF total - WL RHF total = -1.4383600322798173e-5 Ha
WL-parent PQS H1 - WL H1 = +9.649649361120893e-6 Ha
WL-parent PQS H1-J self-Coulomb - WL H1 self-Coulomb = -4.997485057112172e-6 Ha
WL-parent PQS RHF two-body - WL RHF two-body = -6.382481464473067e-5 Ha
```

Compared with the previous custom-mapping PQS fixture:

```text
WL-parent PQS RHF total - custom-mapping PQS RHF total = -0.0014307641877961963 Ha
WL-parent PQS H1 - custom-mapping PQS H1 = -0.004652845139180295 Ha
WL-parent PQS H1-J self-Coulomb - custom-mapping PQS H1-J self-Coulomb = +0.0158797896955718 Ha
```

This strongly supports the conclusion that the pass-183 mHa-scale PQS/WL gap
was mostly a parent-lattice/mapping mismatch. The intended comparison requires
the same parent lattice and shellification inventory; the current tracked PQS
fixture still uses a custom `AsinhMapping(a=0.25, ...)`, while the old WL
baseline uses `white_lindsey_atomic_mapping(Z=2,d=0.3,tail_spacing=10.0)`.

The next pass should align the tracked focused PQS He fixture to that WL parent
mapping and turn the H1/H1-J values into compact scientific fingerprints. It
should also shrink loose metadata/nonclaim assertions enough to keep the
source/test/generator line budget negative.

-- repo-manager@macmini
