Pass 163 review - accepted

Reviewed commit:

```text
51aace02 Add CR2 read-only Hamiltonian inspection view
```

Verdict: accepted.

The pass did the right small thing. It made the existing diatomic Hamiltonian
consumer readiness tuple say exactly what CR2 can inspect today, without adding
a new payload layer, solver path, export path, public API, or test file.

Important resulting facts:

- `cr2_read_only_inspector_ready = true`
- `cr2_solver_ready = false`
- `cr2_export_ready = false`
- `cr2_handoff_blocker = :missing_cr2_solver_handoff_format`
- `two_body_representation_kind = :pre_final_density_interaction`
- `density_gauge = :pre_final_localized_positive_weight`
- `raw_pair_factor_convention = :raw_numerator`

The main readiness blocker correctly remains:

```text
:missing_hfdmrg_density_density_contract
```

This preserves the boundary: the current Be2 PQS route can be inspected by a
downstream consumer, but it is not yet a CR2 solver handoff, HamV6 export,
HFDMRG density-density object, dense `Vee`, public API, artifact, H1/J product,
or RHF product.

Line budget check:

```text
src/test added:   23
src/test deleted: 31
net:              -8
```

Manager validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

passed 58/58 in 47.6s.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed.

```text
git diff --check HEAD~1..HEAD
```

passed.

Next direction:

Do not add another GaussletBases payload immediately. The next useful step is
to prepare a CR2-facing read-only handoff blurb so the downstream CR2 agent can
say what concrete solver/handoff format it wants for Be2 WL/PQS comparison and
future Qiu-White atom-local-HF corrections.

-- repo-manager@macmini
