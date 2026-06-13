Pass 156 - retire RHF seam-test pressure from the default suite

Role: repo-doer@macmini

Task type: deletion/shrinkage implementation.

Purpose:

Remove low-value private RHF seam-test pressure from the default nested suite.
RHF is a private Hamiltonian validator only; HFDMRG is the serious HF/DMRG
package. The repo should not keep four standing RHF seam tests in the default
runner just to preserve private staged payload vocabulary.

Governing policy:

- `AGENTS.md` test scope/deletion policy: tests are code with carrying cost;
  development scaffolding tests should be deleted or quarantined once
  superseded.
- `docs/src/developer/pqs_source_box_fixture_policy.md`: compact H1/J/RHF
  route facts are route-smoke/convention diagnostics, not physics endpoints.
- `docs/src/developer/pqs_source_box_operator_framework.md`: PQS remains
  source-box-first; RHF/SCF is not the product path.

Current state:

`test/nested/runtests.jl` includes four RHF files by default:

```text
pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl
pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl
pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl
pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

Pass 155 found that the first three are pure seam/scaffolding tests and the
fourth should shrink to one compact private-Hamiltonian validator if RHF remains
tested at all.

Exact task:

1. Remove all four RHF includes from `test/nested/runtests.jl`.

2. Delete these three test files:

```text
test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl
test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl
test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl
```

3. Shrink `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
   to one compact private RHF validator. Keep only the smallest synthetic check
   that proves the live convention:

   - closed-shell spin-summed density has trace 2 for the two-electron toy
     fixture;
   - the tiny synthetic Hamiltonian converges;
   - trace/idempotency residuals are small;
   - the final one-step payload is recomputed for the returned final density;
   - the result is private/diagnostic, not public/export/artifact/route
     production.

   Delete broad control-payload inventory, DIIS option inventory, blocked-path
   seam tests, repeated status vocabulary, and repeated nonclaim flag checks.

4. Do not delete or edit `src/pqs_multilayer_complete_core_shell_rhf.jl` in this
   pass. Source deletion is a separate decision; for now RHF remains a private
   Hamiltonian validator.

5. Do not touch the Be2 Ham payload fingerprint test in this pass. That file is
   active route-spine pressure and will be shrunk separately.

Trust boundary:

- Tests only.
- No source code changes.
- No route-driver wiring.
- No H1/J, dense `Vee`, Hamiltonian consumer, export/HamV6, CR2/HFDMRG, WL,
  RHF feature, SCF control, DIIS tuning, fixture promotion, or physics endpoint
  work.

Decision rules:

- If shrinking the SCF test exposes a real source bug, stop and report the bug
  rather than patching source in this deletion pass.
- If the compact validator cannot be kept without importing helpers from one of
  the deleted files, inline only the minimal synthetic fixture setup in the SCF
  file.
- If the remaining compact validator is still mostly status/blocker vocabulary,
  shrink further; the useful test is convention/numerical behavior, not payload
  names.

Validation:

Run:

```text
julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
rg -n "pqs_multilayer_complete_core_shell_rhf_(input_contract|initial_density|one_step|scf)_runtests" test/nested/runtests.jl test/runtests.jl
git diff --check
```

The `rg` command should return no matches in runner files after deletion. Do
not run the full nested default suite unless you first explain why it is needed;
this is a deletion/shrinkage pass, not a broad integration gate.

Report back:

- files deleted;
- runner includes removed;
- remaining RHF validator line count/check count;
- validation commands and results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Commit the change if validation passes, with a message like:

```text
Retire RHF seam test pressure
```

-- repo-manager@macmini
