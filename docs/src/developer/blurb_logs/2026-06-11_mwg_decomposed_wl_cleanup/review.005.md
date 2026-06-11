Review result:

Accepted as diagnostic-only. The pass did not find a clear implementation bug
in the side13 He + GTO final-basis RHF path, and it correctly avoided adding an
acceptance test or changing production code.

Scientific interpretation:

The audit confirms that the probe-local RHF formula follows the same
closed-shell two-index density-density convention as the gausslet-only He
probe. It also confirms that the final-basis one-electron path is ordinary
symmetric, with no generalized final solve and no final self-overlap used as
working metric data.

The result remains diagnostic because the side13 He + GTO total
`-2.864394277839975` is below the He HF reference
`-2.861679995612234` by about `2.7 mHa`. The audit's most plausible
interpretation is that the residual MWG/IDA Hamiltonian is approximate and not
necessarily variational relative to exact all-electron HF, but that still needs
a matched old-QW/MWG oracle or a located White-Lindsey Fig. 8 source table.

Important convention note:

The artifact records that residual MWG interaction code density-normalizes
effective-Gaussian pair factors inside `_qwrg_mwg_interaction_components`.
That is the old MWG convention; it is not identical to the retained
weight-boundary convention used by the decomposed WL gausslet block. This is
not automatically a bug, but it is exactly why the next comparison should be a
matched old-QW/MWG oracle rather than an acceptance baseline.

Validation reviewed:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/side13_he_gto_final_basis_rhf_convention_audit.txt`

Deletion/shrinkage review:

No old surface became obsolete. No test was added. The audit correctly
preserved the boundary between probe evidence and acceptance tests.

Commit/push:

Pending manager commit/push of the tracked response/review log.

Next target:

Run a matched old nested/QW MWG oracle comparison for the side13 He + GTO
diagnostic, and separately frame the White-Lindsey Fig. 8 `ns = 5`, `d = 0.3`,
447-basis-function S-only reference point if the source or a reproducible
fixture is available.
