Implemented the Be atom S+P old nested/QW oracle readiness probe.

Artifact:

- `tmp/work/be_atom_sp_qw_oracle_probe.jl`
- `tmp/work/be_atom_sp_qw_oracle_summary.txt`

Result:

- status: `blocked`
- blocker: `:missing_required_be_ccpv5z_basisfile`
- required basisfile: `/Users/srw/BasisSets`
- `isfile("/Users/srw/BasisSets") == false`

I did not substitute another basis source. The blurb explicitly required:

```julia
legacy_atomic_gaussian_supplement(
    "Be",
    "cc-pV5Z";
    lmax = 1,
    basisfile = "/Users/srw/BasisSets",
)
```

and said to stop/report if that source is missing. Because the trusted old nested/QW oracle could not be built first, I did not audit the newer decomposed/final-basis Be S+P path in this pass.

What did not run:

- old nested/QW Be S+P fixed-block construction
- old nested/QW Be S+P MWG operators
- four-electron closed-shell RHF
- decomposed/final-basis Be S+P one-body or electron-electron audit

Validation:

- Ran `julia --project=. tmp/work/be_atom_sp_qw_oracle_probe.jl`
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Ran `git diff --check`

Deletion/shrinkage report:

- No production code, tests, docs, metadata, or compatibility paths became obsolete.
- Nothing was deleted or simplified because the pass stopped at the required machine-local basisfile availability blocker.
- No new tracked test was added. The probe remains under `tmp/work`.
- Remaining stale/duplicate surfaces to retire next: none identified. The next step is to restore/provide `/Users/srw/BasisSets` or revise the blurb to authorize a specific alternate Be `cc-pV5Z` source.

-- repo-doer@macmini
