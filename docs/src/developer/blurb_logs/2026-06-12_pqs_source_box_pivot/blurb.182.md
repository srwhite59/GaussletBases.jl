Pass 182 - local ignored RHF diagnostic probe on the corrected 419-dimensional fixed-q PQS He H1-J payload.

Current accepted state:

```text
fixed-q one-center He PQS final dimension = 419
H1 energy = -1.9866819751748936
H1-J self-Coulomb = 1.2261626003119184
density gauge = pre_final_localized_positive_weight
```

This pass is a local diagnostic probe only. The purpose is to get the first
RHF/HF scalar from the corrected 419-dimensional PQS Hamiltonian path so we can
compare directionally against the old WL Fig.8-style He n_s=5 result later.

Do not edit tracked `src/`, `test/`, docs, generator, or artifacts except for
the required response handoff copy. Do not commit. If you unexpectedly need a
tracked source/test change, stop with `.agent_handoffs/ATTENTION.md`.

Allowed:

- create an ignored `tmp/work/...` probe script and, if useful, an ignored
  summary file;
- include `test/nested/pqs_direct_retained_final_h1_runtests.jl` to reuse
  `_pqs_h1_complete_fixture()`;
- build the same H1 and H1-J payload as pass 181;
- use the existing private RHF diagnostic helpers as a Hamiltonian validator.

Forbidden:

- no route-driver RHF wiring;
- no tracked test;
- no source changes;
- no Be2/Cr2 artifact/generator work;
- no HFDMRG/DMRG/CR2 run;
- no AHGBS residual/supplement layer;
- no public API/export/artifact behavior.

Suggested probe shape:

1. Include the focused H1 test file to get the fixture helper.

2. Build:

   ```text
   fixture
   H1 payload
   density inputs from `_pqs_source_box_route_driver_complete_core_shell_density_inputs(...)`
     or from `_pqs_source_box_ida_factor_provenance(...)` wrapped in the existing expected shape
   H1-J payload
   RHF input contract with electron_count = 2 and fixture_role = :route_smoke
   RHF SCF payload
   ```

3. Run at least one conservative control:

   ```text
   mixing_kind = :fock_diis
   max_iterations = 50 or 75
   residual_atol = 1e-8
   density_atol = 1e-8
   energy_atol = 1e-10
   ```

   If that is slow or unstable, also try fixed-point only if cheap. Do not tune
   endlessly; this is a probe, not an SCF development pass.

Report:

```text
input contract status/blocker
SCF status/blocker
converged?
iterations
final total energy
one-body energy
two-body energy
density trace
commutator/residual if available
idempotency error if available
H1 energy and H1-J self-Coulomb from the input path
elapsed time
```

Interpretation constraints:

- This RHF result is private diagnostic/probe behavior only.
- Do not call it a physics endpoint yet.
- Do not compare it as accepted against WL unless you explicitly state the WL reference is old Fig.8-style supplemented WL and this PQS value is gausslet-only/no AHGBS residual.
- If SCF does not converge, report the residual/blocker. Do not patch SCF controls in this pass.

Line-budget rule:

- Since this is a no-tracked-source/test pass, `git diff --numstat -- src test` should be empty.
- If any tracked `src`/`test` edit appears, it must be net-negative or you must stop with `ATTENTION.md`.

Validation:

```text
julia --project=. tmp/work/<your_probe>.jl
git status --short --branch
git diff --numstat -- src test
```

Write the result to `.agent_handoffs/response.182.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.182.md
```

-- repo-manager@macmini
