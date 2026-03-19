# Documentation Navigation Feel Note

The main remaining documentation problem is no longer missing content or
missing structure. The current issue is navigation feel.

The docs site already has the right pieces:

- a landing-page README
- a rendered Documenter site
- a user manual
- a curated API reference
- a demoted developer/supporting-note section

But the site can still feel too much like a visible directory tree if the
sidebar foregrounds every child page equally.

The next polish step is therefore to make the site feel more like a small set
of primary user documents:

- Home
- Manual
- Reference
- Developer Notes

Those should read as the main clickable documents of the site, not just as
section labels wrapped around a file tree.

The intended reading experience is closer to the feel of packages like
KrylovKit:

- a small visible top-level surface
- user-oriented prose first
- a manual that does the routing work
- a clearly separate reference section
- developer/history notes that remain available without competing with the
  user manual

This should be done with minimal content churn:

- strengthen the landing pages
- simplify the visible sidebar
- keep the underlying source-file organization if it is still useful
- avoid another broad rewrite of the scientific pages
