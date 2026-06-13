Pass 100 response

Moved the stale ignored H1/J probe artifacts out of active `tmp/work/`
circulation.

Files moved:

- `tmp/work/pqs_direct_retained_final_h1_timing_probe.jl`
- `tmp/work/pqs_direct_retained_final_h1_probe_summary.txt`
- `tmp/work/pqs_complete_core_shell_h1_probe.jl`
- `tmp/work/pqs_complete_core_shell_h1_probe_summary.txt`
- `tmp/work/pqs_complete_core_shell_final_density_j_probe.jl`
- `tmp/work/pqs_complete_core_shell_final_density_j_probe_summary.txt`

Archive directory used:

- `tmp/work/archive_stale_pqs_h1j_2026-06-12/`

Validation:

- Verified each source file existed before moving.
- Ran `git check-ignore -v` on each source file before moving; all were ignored by `.gitignore:12:tmp/`.
- Created `tmp/work/archive_stale_pqs_h1j_2026-06-12/`.
- Moved the files with `mv`, not `rm`.
- Ran `git check-ignore -v` on each archived path; all remain ignored by `.gitignore:12:tmp/`.
- Did not edit source, tests, docs, tracked files, or run Julia.

Archived paths now present:

- `tmp/work/archive_stale_pqs_h1j_2026-06-12/pqs_direct_retained_final_h1_timing_probe.jl`
- `tmp/work/archive_stale_pqs_h1j_2026-06-12/pqs_direct_retained_final_h1_probe_summary.txt`
- `tmp/work/archive_stale_pqs_h1j_2026-06-12/pqs_complete_core_shell_h1_probe.jl`
- `tmp/work/archive_stale_pqs_h1j_2026-06-12/pqs_complete_core_shell_h1_probe_summary.txt`
- `tmp/work/archive_stale_pqs_h1j_2026-06-12/pqs_complete_core_shell_final_density_j_probe.jl`
- `tmp/work/archive_stale_pqs_h1j_2026-06-12/pqs_complete_core_shell_final_density_j_probe_summary.txt`

Git status after move, before writing this response:

```text
## main...origin/main
```

Deletion/shrinkage report:

- deleted: none; files were moved, not removed.
- simplified: active `tmp/work/` no longer presents these stale H1/J timing,
  H1, and rejected final-density convention probes as current validation
  pressure.
- quarantined: six ignored local artifacts under
  `tmp/work/archive_stale_pqs_h1j_2026-06-12/`.
- not deleted because: baton requested `mv`, not `rm`, preserving local history
  while removing it from active circulation.
- exact remaining caller/blocker: no code caller; remaining blocker for tracked
  H1/J test shrinkage is still the absence of a compact tracked driver H1/J
  smoke over the accepted one-center source-box diagnostic path.

-- repo-doer@macmini
