Review 040: accepted, and the baton loop should pause here.

The physical one-center H1 probe did what it was meant to do: it separated
route correctness from physical basis adequacy.

Route correctness result:

```text
source dims/count: 5 x 5 x 5 / 125
boundary/final retained count: 98 / 98
direct retained overlap/kinetic/nuclear flags: all true
active raw source one-body materialization: false
final overlap identity error: 2.24e-14
final Hamiltonian vs shell-support oracle: 1.33e-15
H1 eigenvalue vs shell-support oracle: 3.75e-16
ordinary symmetric eigensolve: true
generalized overlap solve: false
```

Physical interpretation:

```text
Z = 1 H reference: -0.5
boundary-shell final-basis H1: -0.08171962129085239
error: +0.41828037870914758
```

That large H error is expected for the current fixture because it is only the
PQS boundary-shell source-mode final basis, not a complete one-center
parent/final basis. This should not become a physical H acceptance baseline.

Design fork:

- Continue with final IDA-weight machinery on the boundary-shell route as a
  convention/kernel target, knowing it is not a physical acceptance basis.
- Or first define a complete one-center PQS/final basis route so H/He+ H1 is
  physically meaningful before adding final IDA/density-density.

This is a manager/user decision point. Do not continue the unattended baton loop
into IDA/RHF until that direction is chosen.

-- repo-manager@macmini
