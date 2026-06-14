# Pass 235 manager review - accepted

Accepted.

This pass paid down the pass-234 line-budget exception by deleting old
route-shadow density-density fixture pressure. It removed private
PQS/PQS/product density-density producer/consumer helpers and the corresponding
slow integration-test sections, while preserving compact lower-level math
convention checks.

Key checks:

- Deleted helper names have no remaining callers in `src`, `test`, or `bin`.
- Lower-level density-normalized/raw-weighted pair convention checks remain.
- Nuclear charge/sign convention checks remain.
- Boundary retained count checks remain.
- Independent H2 PQS support-plan work, fake-PQS guard fields, and fake-PQS
  golden regression were untouched.

Manager validation:

- Reviewed the diff and response.
- Ran the deleted-name caller search; no matches.
- Ran `git diff --check`.
- Ran package load:
  `julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'`
  with `elapsed_s=0.650992375`.
- Did not rerun the focused integration check. Doer ran it successfully:
  `3820/3820` passed, reported test time `6m04.6s` and elapsed
  `845.430618416s`.

Line-budget check:

- Scoped `src + test + bin` diff is `19` added / `1587` deleted, net `-1568`.
- This more than pays down the pass-234 `+108` exception.

Guardrail:

- The deleted path was route-shadow fixture pressure, not active independent H2
  PQS route authority. Do not reintroduce it to satisfy future metadata tests.

Next:

- Return to independent H2 PQS recovery with a no-edit audit of retained-rule
  and source-plan authority for `atom_contact_core`, `shared_shell_1`, and
  `shared_shell_2`.

-- repo-manager@macmini
