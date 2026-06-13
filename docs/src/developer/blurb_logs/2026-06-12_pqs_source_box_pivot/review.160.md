Pass 160 review - diatomic consumer duplication audit

Accepted. The audit stayed read-only and identified a concrete shrink pass.

Key finding:

The private Hamiltonian consumer contract added in pass 159 correctly replaced
the active blocker, but it repeats many fields already owned by the handoff:
final dimension, one-body status/reference, density gauge, raw pair convention,
support count, pair shapes, nuclear/electron/spin metadata, and downstream false
flags. Tests currently assert some of those copied fields through the consumer,
which would make the duplication harder to remove later.

Decision:

Pass 161 should be a line-negative source/test shrink:

- thin the consumer payload to source-handoff reference plus compact
  readiness/summary/metadata;
- remove copied scalar fields from the consumer struct and builder;
- consolidate repeated downstream nonclaim flags into one small helper or local
  bundle;
- keep readiness blocked on `:missing_hfdmrg_density_density_contract`;
- trim Be2 test assertions that only preserve copied consumer fields.

Line-budget rule:

Pass 161 must be net-line-negative across tracked `src/` + `test/` files:

```text
git diff --numstat -- src test
```

with total deleted lines greater than total added lines. The target is at least
70 deletions and no more than 30 additions, preferably a net reduction of 65+
lines. If that cannot be achieved safely, doer must write `ATTENTION.md` and
stop rather than making a line-positive or line-neutral source/test pass.

Deletion/shrinkage accounting:

- deleted: none by pass 160
- simplified: pass 161 shrink surface identified
- quarantined: none
- not deleted because: audit only
- exact remaining caller/blocker: active readiness remains blocked on
  `:missing_hfdmrg_density_density_contract`; consumer/readiness duplication
  should be shrunk before any HFDMRG, CR2, HamV6, dense `Vee`, or export format
  work

-- repo-manager@macmini
