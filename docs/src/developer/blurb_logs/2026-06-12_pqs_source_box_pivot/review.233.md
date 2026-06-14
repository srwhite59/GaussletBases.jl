# Pass 233 manager review - accepted

Accepted.

The audit found a plausible route-owned support-region path and, importantly,
kept the support counts classified as target constants until generated. The
recommended next seam is narrow enough for implementation: use shellification
geometry to derive primitive support regions, group the two atom-local cores and
midpoint slab into `:atom_contact_core`, and order shared molecular shells
outside-in.

Manager read:

- `atom_contact_core = 275` is derivable geometrically as two `5^3` atom cores
  plus one `5*5*1` contact/midpoint slab.
- Shared shell counts are also derivable geometrically:
  inner `7*7*13 - 5*5*11 = 362`, outer `9*9*15 - 7*7*13 = 578`.
- The target vocabulary wants outside-in order:
  `:shared_shell_1 => 578`, `:shared_shell_2 => 362`.
- Existing shellification/route-core objects own primitive support. The
  three-unit H2 physical support vocabulary should be a compact private
  route-owned grouping object.

Guardrail:

- Do not use `bond_aligned_diatomic_nested_fixed_source(...)`, fake-PQS
  `source.sequence.coefficient_matrix`, or any WL/QW retained-transform data as
  independent support authority.

Next:

- Pass 234 should implement only the support-region materializer/fingerprint.
  It should keep source-plan, retained rules, final basis, H1, H1-J, RHF,
  supplements, CR2, export, and public API blocked.

-- repo-manager@macmini
