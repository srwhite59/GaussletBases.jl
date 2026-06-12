Review 069: accepted.

The tracked H1 gate now reuses `GaussletBases._pqs_multilayer_support_product_matrix`
instead of carrying a duplicate local rectangular product helper. This is the
right kind of small cleanup: no production code changed, no tests were added,
and the H1 gate stayed focused.

The remaining local helpers in that test are kinetic and nuclear support
operator assembly. Kinetic is probably straightforward, but nuclear carries
by-center/sign/charge/Gaussian-factor-source convention risk. The next step
should be an audit of that seam before moving it into production.

-- repo-manager@macmini
