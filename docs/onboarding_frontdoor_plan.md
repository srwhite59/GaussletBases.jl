# Front-Door Onboarding Plan

This note records the current onboarding pass.

## 1. What is wrong with the present front door

The package documentation is now organized, but the first reader experience is
still weaker than it should be.

The main problems are:

- the README begins too quickly in package-structure language instead of first
  explaining what gausslets are
- the first radial quickstart introduces setup knobs before the reader knows
  which ones matter for a first calculation
- the first README example stops too early instead of carrying through to one
  useful physical result

So the issue is no longer information architecture. It is onboarding and
teaching.

## 2. What belongs in the README

The README should do four things:

1. explain what gausslets are in plain scientific language
2. explain why they are interesting
3. show what this package can currently do with them
4. carry one first radial example through to a hydrogen energy

It should not try to teach every setup knob.

## 3. What belongs in the radial quickstart

The radial quickstart should:

- teach the basic basis -> diagnostics -> quadrature -> operators -> hydrogen
  workflow
- explain the first appearance of `s`, `c`, and `AsinhMapping(c = ..., s = ...)`
- make clear that `c = ...` and `s = ...` are ordinary Julia keyword
  arguments to the mapping constructor

It can introduce more detail than the README, but it should still teach the
workflow before it teaches the tuning knobs.

## 4. What belongs only in the deeper setup note

The deeper setup note is the right place for:

- `reference_spacing`
- `tails`
- `odd_even_kmax`
- `xgaussians`
- fuller discussion of how `s` and `c` may need to change together

Those details matter, but they should not dominate the first page a new reader
opens.
