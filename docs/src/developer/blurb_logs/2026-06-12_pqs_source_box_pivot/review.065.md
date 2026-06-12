Review 065: accepted.

The generic retained PQS one-body selector now covers overlap, kinetic,
position, and x2 source terms and routes position/x2 through the pass-064 direct
retained wrappers. This is the right cleanup: callers can use the common
retained selector without falling back to full raw source matrices for moment
terms.

Important guardrails were preserved:

- no new exports;
- no H1/RHF/IDA/density-density/driver/artifact changes;
- no fixture-rule study;
- raw-source selector remains oracle/reference;
- tests check term/equality/dimensions rather than broad report vocabulary.

The CPBM contract test is still heavy, but this pass did not add a new file and
the checks are localized to the existing retained/PQS selector section.

At this point the obvious local retained one-body source-block gap is closed.
The next issue is architectural rather than numerical: the current PQS pieces
work in probes, but they need to be mapped back onto the driver spine before
they become production behavior. Do not start another physics ladder until that
route-stage ownership is clear.

-- repo-manager@macmini
