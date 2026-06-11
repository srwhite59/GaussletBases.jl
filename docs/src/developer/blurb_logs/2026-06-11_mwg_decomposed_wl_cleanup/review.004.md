Review result:

Accepted as diagnostic-only. The probe used the intended final-basis path:
combined one-electron matrices, residual MWG density-density blocks, ordinary
symmetric RHF, and no raw GTO density-density final operator or fallback route.

Scientific interpretation:

The side13 He + GTO result improves substantially over side13 gausslet-only
RHF:

- side13 gausslet-only RHF total: `-2.8364979997009137`
- side13 He + GTO RHF total: `-2.864394277839975`
- improvement: `0.02789627813906126`

It is not an acceptance result because it is below the He HF reference
`-2.861679995612234` by about `2.7 mHa`. That violates the intended
variational sanity check for an HF acceptance fixture unless the MWG/IDA
approximation is explicitly being judged against a different internal
Hamiltonian reference.

Corrections made:

None. This pass was intentionally probe-only and made no tracked source or test
changes.

Validation reviewed:

- `julia --project=. tmp/work/side13_he_gto_final_basis_rhf_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/side13_he_gto_final_basis_rhf_summary.txt`

Deletion/shrinkage review:

No old surface became obsolete. No test was added. The probe exposed the next
scientific/convention question, so the right next action is an audit, not an
acceptance test or more feature expansion.

Commit/push:

Pending manager commit/push of the tracked response/review log.

Next target:

Audit the final-basis density-density and RHF convention before making any
He + GTO acceptance claim.
