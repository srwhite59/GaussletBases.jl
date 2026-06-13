Pass 181 - add H1-J/density-interaction smoke on the corrected 419-dimensional fixed-q PQS He basis.

Current accepted state:

- Pass 179 fixed route-owned fixed-q PQS source modes.
- Pass 180 restored H1 on the corrected one-center He basis.
- The live atom fixture now has:

```text
final dimension = 419
H1 lowest energy = -1.9866819751748936
H1 matrix finite/symmetric
```

This pass should add the next diagnostic seam: H1-J / pre-final density interaction on the same 419 basis. Do not go to RHF.

Physics target:

```text
one-center He q=5/n_s=5 fixed-q PQS basis
H1-J / density-interaction diagnostic only
no RHF/SCF/DIIS
no DMRG/HFDMRG
no AHGBS residual/supplement layer
no Be2/Cr2 artifact generation or comparison
no exports/artifacts/public API
```

Implementation scope:

- Prefer editing only `test/nested/pqs_direct_retained_final_h1_runtests.jl`.
- You probably should not need source changes.
- If the existing H1-J helper fails on the 419 basis, write `.agent_handoffs/ATTENTION.md` with the exact blocker rather than patching broad source code.

Use the existing route-owned/internal density input convention:

```julia
provenance =
    GaussletBases.CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(
        fixture.plan.bundles;
        expected_term_count = length(fixture.expansion.coefficients),
    )

h1_j_payload = GaussletBases.pqs_multilayer_complete_core_shell_h1_j_payload(
    fixture.plan;
    final_basis = fixture.final_basis,
    h1_payload,
    axis_weights = provenance.axis_weights,
    raw_pair_factor_terms = provenance.raw_axis_pair_factor_terms,
    coulomb_expansion = fixture.expansion,
    metadata = (; fixture = :pqs_fixed_q_he_h1j_gate),
)
```

Adjust names to match local style. Do not use retained diagnostic weights, density-normalized pair terms, signed-final-weight division, raw-no-division, or old fixed-block oracle data.

The focused test should assert only the live H1-J diagnostic contract:

```text
h1_j_payload.status == :materialized_pqs_multilayer_complete_core_shell_h1_j_payload
h1_j_payload.summary.final_dimension == 419
h1_j_payload.summary.h1_energy matches the H1 payload energy
h1_j_payload.summary.h1_energy_reconstruction_error is small
density_interaction.status == :materialized_pqs_complete_core_shell_pre_final_density_interaction
density_gauge == :pre_final_localized_positive_weight
pre_final_pair_matrix_finite == true
pre_final_weights_all_positive == true
self_coulomb is finite and positive
RHF/GTO/driver/export/artifact flags remain false
```

Do not assert a tight reviewed self-Coulomb or total HF energy yet. Report the observed self-Coulomb and H1-J summary values in the response.

Line-budget rule:

```text
git diff --numstat -- src test
sum(deleted) > sum(added)
```

Try to keep the H1-J assertions compact rather than adding another scaffold. Good deletion/shrink options:

1. First, trim low-value repeated nonclaim assertions from `test/nested/pqs_direct_retained_final_h1_runtests.jl` if they do not protect this seam.

2. If that is not enough, `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl` is now a legitimate deletion candidate:
   - Be2/Cr2 work is paused.
   - The previous Be2/WL/PQS artifact was explicitly marked not same-basis comparison-ready.
   - `rg` shows this file is not included by `test/nested/runtests.jl`.
   - It is a route/artifact fingerprint scaffold, not the current atom-first physics target.

Before deleting the Be2 file, verify no live `src`/`test` caller or default-runner include remains. If it is only standalone scaffold plus historical handoff references, delete it outright rather than moving it.

Validation:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

Do not run a broad suite. If H1-J makes the focused test much slower, report the elapsed time.

Reporting requirements:

- exact files changed/deleted;
- observed final dimension, H1 energy, self-Coulomb, density gauge, and density-interaction status;
- confirmation of the raw provenance source used for `axis_weights` and `raw_pair_factor_terms`;
- source/test line additions, deletions, and net;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Write the result to `.agent_handoffs/response.181.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.181.md
```

-- repo-manager@macmini
