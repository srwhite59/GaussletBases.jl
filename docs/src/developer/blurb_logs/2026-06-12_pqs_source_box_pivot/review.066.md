Review 066: accepted.

The driver-spine audit gives the right next seam. The successful PQS probes
should not become a separate route family; they should feed the standard stages.

Most important conclusion:

```text
cartesian_assembly is the first missing seam.
```

The route has module-owned pieces for source planning, retained one-body blocks,
multi-layer shell source plans, final-basis realization, final IDA weights, and
pre-final density interaction. But probes still locally rebuild the complete
core/shell final basis from a multi-layer plan by slicing support overlap blocks
and calling `pqs_complete_core_shell_final_basis(...)`.

That repeated probe-local assembly is now the best next shrink target. A narrow
helper that consumes `pqs_multilayer_shell_source_plan(...)` and publishes the
complete core/shell final-basis payload would reduce private probe code without
adding driver wiring, H1, RHF, density, exports, or fixture acceptance.

-- repo-manager@macmini
