Accepted pass 185.

The tracked focused PQS He fixture now uses the same parent mapping as the old
WL q=5/n_s=5 gausslet-only baseline:

```julia
white_lindsey_atomic_mapping(Z = 2.0, d = 0.3, tail_spacing = 10.0)
```

Manager inspected the diff and reran:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
```

which passed:

```text
PQS fixed-q complete core-shell inventory gate | 32 pass, 32 total, 6.4s
```

The test now asserts the live scalar fingerprints:

```text
final dimension = 419
H1 lowest = -1.991334820314074
H1-J self-Coulomb = 1.2420423900074902
```

and deleted the old custom mapping setup plus loose nonclaim/finite/positive
checks. The source/test line budget was:

```text
6 added, 28 deleted, net -22
```

This is the right tracked fixture correction. The next small improvement is to
make the WL comparison visible in the tracked gate without growing the suite:
add compact WL baseline delta assertions in this focused test and delete the
stale integration-only fixed-block term-storage test, whose remaining contract
is already covered elsewhere.

-- repo-manager@macmini
