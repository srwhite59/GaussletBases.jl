Review 031

Accepted as an incomplete validation checkpoint.

The response followed the updated 10-minute slow-validation policy. It tried
the requested slow nested harness, handled the Juliaup lockfile issue by using
the concrete Julia binary, then stopped the still-running validation after
about 11.5 minutes. The edited pass-029 integration file was not reached, so
there is no evidence of a regression in that file and no full slow-harness
validation claim.

Recorded status:

- pass-029 integration cleanup remains committed but not fully slow-harness
  validated;
- pass-030 route-driver stack-overflow fix remains focused-test validated;
- slow nested harness is too broad for unattended baton validation under the
  current 10-minute silence/wait policy.

Deletion/shrinkage:

- No files changed in pass 031 except the tracked response log.
- No tests were added.
- No old surface became obsolete from validation alone.

Next:

- Continue with the near-term PQS cleanup plan: status/ownership wording,
  final-basis realization module boundary, then direct retained-boundary
  overlap/kinetic kernels.

-- repo-manager@macmini
