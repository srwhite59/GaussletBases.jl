Purpose:
Move stale ignored H1/J probe artifacts out of active `tmp/work/` circulation.

Why now:
Pass 099 identified local ignored development probes whose questions are now
answered by tracked tests or the driver-owned private H1/J diagnostic path.
They should not keep acting like active validation pressure.

Governing framework:
Use `AGENTS.md` deletion/shrinkage policy and
`docs/src/developer/pqs_source_box_operator_framework.md`.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
already-approved commands. If approval would be required, write `ATTENTION.md`
with the exact command, reason, and blocker, then stop.

Exact task:
Move the following ignored local files into an ignored archive directory:

`tmp/work/archive_stale_pqs_h1j_2026-06-12/`

Files to move:

- `tmp/work/pqs_direct_retained_final_h1_timing_probe.jl`
- `tmp/work/pqs_direct_retained_final_h1_probe_summary.txt`
- `tmp/work/pqs_complete_core_shell_h1_probe.jl`
- `tmp/work/pqs_complete_core_shell_h1_probe_summary.txt`
- `tmp/work/pqs_complete_core_shell_final_density_j_probe.jl`
- `tmp/work/pqs_complete_core_shell_final_density_j_probe_summary.txt`

Use `mv`, not `rm`.

Before moving, verify each file exists and is ignored by git. If any file is not
ignored, stop and report instead of moving it.

Do not edit source, tests, docs, or tracked files.
Do not delete tracked files.
Do not add tests.
Do not run Julia.
Do not touch RHF/SCF/Fock, GTO, IDA/MWG, exports, artifacts, fixture promotion,
or production route behavior.

Validation:

- `git check-ignore -v <each file before move>`
- `git status --short --branch` after move

Report back:

- files moved;
- archive directory used;
- confirmation moved files remain ignored;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.100.md`, continue polling for
`blurb.101.md`, `STOP.md`, or `ATTENTION.md`.
