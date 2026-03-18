# Dense Hamiltonian Export: `fullida_dense_v1`

This note records the first Hamiltonian export layer for downstream solver
codes.

## 1. Why export is the next step

The package now has coherent in-memory Hamiltonian objects, especially on the
atomic IDA line.

The next high-value step is therefore not more in-package solver work.

It is a clean producer/export layer for external consumers such as:

- HFDMRG-style workflows
- DMRG consumers
- later correlation-method consumers

That keeps the solver layers outside this package while making the Hamiltonian
data easy to reuse.

## 2. Why dense `fullida_dense_v1` comes first

The first export target should be the dense bridge format already used in the
existing ecosystem:

- format note:
  `~/Dropbox/codexhome/docs/Format.FullIDADenseBridge.v1.md`
- existing consumer:
  `work/slicedmrgutils/runs/dmrg_fullida_dense.jl`

That format is the highest-value first target because:

- it is already defined
- it is already consumed
- it matches the current dense in-memory atomic IDA objects well

## 3. Why sliced/block export comes second

The sliced/block export should come later.

That format is structurally richer, but it is also a larger design step. The
first milestone should be one clean dense bridge producer.

## 4. What the current model means physically

The current export must reflect the present model honestly.

For the current atomic and ordinary IDA lines, the interaction is:

- density-density
- two-index IDA
- not a full four-index Coulomb tensor

So the export should write the actual current Hamiltonian model, not a
fictitious more general one.

## 5. First implementation scope

The first implementation scope is:

- `AtomicIDAOperators`
- dense `H1`
- dense `Vee`
- `fullida_dense_v1` bridge metadata

The ordinary Cartesian line can come second once it is clearly in the same
solver-facing state. The dense bridge format has no overlap channel, so the
ordinary line should not be exported casually before its basis side is settled
enough for downstream consumers.
