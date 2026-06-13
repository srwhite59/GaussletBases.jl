Pass 101 review

Accepted timing/design probe.

The result is decisive enough for planning:

- Accepted one-center source-box H1/J dry-run, cold fresh process:
  `120.102067` seconds.
- Same route shape, second run in the same Julia process: `0.210283` seconds.
- The accepted dry-run materializes the private H1/J diagnostic with:
  - final dimension `223`
  - H1 energy `-5.6629907690725245`
  - self-Coulomb `1.8691288063594704`
  - density gauge `pre_final_localized_positive_weight`
  - RHF/GTO/export/artifact flags false.

Conclusion:
Do not add a normal cold-process tracked driver H1/J smoke now. The blocker is
first-call/compilation latency, not repeated route runtime or fixture geometry.
Smaller lattice is not the right lever for the accepted route contract.

The local timing probe and summary are ignored under `tmp/` and can stay as
temporary review evidence.

Next pass should move to RHF contract design only, with no implementation.

-- repo-manager@macmini
