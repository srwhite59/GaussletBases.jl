# Documentation Presentation Plan

This note records the current presentation and navigation refinement pass.

## 1. Current state

The docs site is now structurally sound:

- a Documenter site exists
- curated manual pages exist
- a curated API reference exists
- docs CI exists
- a first doctest slice exists

So the next problem is no longer documentation structure or missing
infrastructure. It is presentation.

## 2. Goal of this pass

The next improvement should make the site feel more like a standard
user-facing Julia package docs site:

- a small number of clear top-level sections
- a left sidebar that reads as a package manual, not as a research tree
- user-oriented prose first
- developer and history material visibly demoted

This is closer in spirit to package docs like KrylovKit than to a notebook of
accumulated project notes.

## 3. Intended top-level shape

The target top-level navigation is:

- Home
- Manual
- Reference
- Developer Notes

The exact labels can vary, but the important point is the hierarchy:

- user-facing reading and usage material first
- reference second
- development and history material last

## 4. What belongs under each bucket

### Manual

The user-facing manual should contain:

- first radial workflow
- recommended atomic setup
- current atomic branch
- current ordinary branch
- example guide

These pages are about using the package and understanding its current
scientific lines. They should not read as developer notes.

### Reference

The reference section should remain the current curated API reference built
from real docstrings.

### Developer Notes

Developer Notes should contain:

- architecture and current direction
- supporting note maps
- transition and hardening notes
- links back out to the flat supporting-note tree

This keeps the development history available without letting it visually
compete with the user manual.

## 5. Scope discipline

This pass should not become another broad content rewrite. The content is
already strong enough. The goal is:

- simpler sidebar organization
- small wording changes where needed
- clearer demotion of developer/history material

That should be enough to make the site feel much closer to a standard
user-facing Julia package docs site.
