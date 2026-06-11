Purpose:

Audit why the side13 He + GTO final-basis RHF diagnostic falls below the He HF
reference before any acceptance claim is made.

Why now:

The first side13 He + GTO final-basis RHF probe used the intended final-basis
route and improved over gausslet-only, but the total energy was
`-2.864394277839975`, about `2.7 mHa` below the He HF reference
`-2.861679995612234`. That may indicate a final-basis density-density or RHF
convention issue, or it may mean the MWG/IDA approximate Hamiltonian should not
be judged directly against the all-electron HF reference. Do not turn this into
an acceptance test until that distinction is understood.

Current state:

- side13 gausslet-only He RHF total: `-2.8364979997009137`
- side13 He + GTO final-basis RHF diagnostic total: `-2.864394277839975`
- He HF reference: `-2.861679995612234`
- final overlap identity error in the diagnostic: about `2.5e-11`
- final solve kind: ordinary symmetric
- raw GTO density-density is not accepted as final operator data
- no full-parent CPB, direct Cartesian product, ordinary Cartesian IDA, PQS, or
  generalized final solve is allowed
- White-Lindsey Fig. 8 has relevant comparison data for this question. The
  manager's working recollection is: for `ns = 5`, Fig. 8 shows an error of a
  few tenths of a mH using `d = 0.3` and 447 basis functions, with a very
  accurate S-only GTO supplement; the error goes up for smaller `d`. Locate the
  repo-local paper/source/reference first if available and verify these details
  against it. If the data is not present in the repo or known local references,
  report that explicitly and treat the manager note as a lead, not a citation.

Exact task:

Run a focused audit of the side13 He + GTO final-basis density-density and RHF
convention. Use `tmp/work` probes and existing code inspection first. Do not
add an acceptance test in this pass.

Inspect and report on these surfaces:

- `tmp/work/side13_he_gto_final_basis_rhf_probe.jl`
- `tmp/work/side13_he_gto_final_basis_rhf_summary.txt`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_density_density.jl`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_final_basis.jl`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_matrix_assembly.jl`
- `test/nested/cartesian_wl_gausslet_he_atom_acceptance_runtests.jl`
- old MWG references such as `_qwrg_mwg_interaction_components` and
  `_qwrg_final_residual_mwg_component_blocks`

Audit questions:

1. Confirm the RHF formula in the probe matches the gausslet-only He RHF
   convention:
   - closed-shell one spatial orbital density trace convention;
   - Fock form `h + 2J - K` for the two-index density-density approximation;
   - energy decomposition and electron count;
   - no accidental factor-of-two or spin-counting mismatch.

2. Confirm the final-basis one-electron transform is the active final-basis
   ordinary-Hermitian path:
   - no generalized final solve;
   - no use of final self-overlap as working data;
   - final overlap is identity within the existing tolerance.

3. Audit final-basis density-density assembly:
   - `V_gg` comes from the corrected decomposed WL raw-numerator then
     retained-weight-division convention;
   - `V_gR` and `V_RR` use residual MWG effective centers/widths;
   - raw GTO/GTO density-density is not accepted as final data;
   - signs and symmetry are correct;
   - the residual transform/order matches the final-basis projection.

4. Compare diagnostic components:
   - gausslet-only side13 RHF one-electron and electron-electron pieces;
   - He + GTO final-basis one-electron and electron-electron pieces;
   - one-electron-only lowest orbital behavior if useful;
   - a simple closed-shell one-orbital scalar check if the probe has a local
     RHF helper.

5. Decide whether the below-HF result is:
   - a code/convention bug to fix now;
   - a known MWG/IDA approximation effect that should be documented and judged
     against an internal old-MWG oracle instead of the exact HF reference;
   - or unresolved, in which case no acceptance test should be added.

6. Compare against White-Lindsey Fig. 8 if the source data is available:
   - identify the exact paper/source file or documented local reference used;
   - verify the `ns = 5`, `d = 0.3`, 447-basis-function, S-only GTO supplement
     point and the reported few-tenths-of-a-mH error if present;
   - check the reported trend that the error worsens for smaller `d`;
   - record whether the current side13 `d = 0.1`, `s = 1.0`, 520-dimensional
     diagnostic is consistent with that historical benchmark or is using a
     materially different fixture;
   - do not treat an uncited remembered value as authority.

Trust boundary:

Do not add a new test by default. Do not add public driver/export/artifact paths.
Do not add a raw GTO density-density fallback. Do not use full-parent CPB, flat
`_qwrg_diatomic_*` route authority, direct Cartesian product assembly, ordinary
Cartesian IDA, PQS, or a generalized final-basis solve.

If a clear bug is found:

- fix it narrowly;
- rerun the side13 He + GTO probe;
- run the most local relevant existing test;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- report the before/after energy and the exact convention corrected.

If no clear bug is found:

- leave production code unchanged;
- write the diagnostic artifact under `tmp/work`;
- report the most likely interpretation and the next needed oracle/comparison;
- do not add an acceptance test.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- artifact paths;
- whether a bug was found;
- whether any source/test/doc files changed;
- side13 He + GTO energy after any fix, or unchanged diagnostic energy;
- whether the result can become acceptance, must stay diagnostic, or needs an
  old-MWG/internal-Hamiltonian oracle;
- validation run;
- deletion/shrinkage report.
