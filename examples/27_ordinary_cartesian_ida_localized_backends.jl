using LinearAlgebra
using GaussletBases

function truncate_expansion(expansion::CoulombGaussianExpansion, nterms::Int)
    return CoulombGaussianExpansion(
        expansion.coefficients[1:nterms],
        expansion.exponents[1:nterms];
        del = expansion.del,
        s = expansion.s,
        c = expansion.c,
        maxu = expansion.maxu,
    )
end

Z = 2.0
npoints = 5
rmax = 6.0
expansion = truncate_expansion(coulomb_gaussian_expansion(doacc = false), 3)

mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = npoints, xmax = rmax)
basis = build_basis(MappedUniformBasisSpec(:G10;
    count = npoints,
    mapping = mapping,
))

reference_1d = mapped_ordinary_one_body_operators(
    basis;
    exponents = expansion.exponents,
    backend = :numerical_reference,
)
proxy_1d = mapped_ordinary_one_body_operators(
    basis;
    exponents = expansion.exponents,
    backend = :pgdg_experimental,
)
localized_1d = mapped_ordinary_one_body_operators(
    basis;
    exponents = expansion.exponents,
    backend = :pgdg_localized_experimental,
)

reference = ordinary_cartesian_ida_operators(
    basis;
    expansion = expansion,
    Z = Z,
    backend = :numerical_reference,
)
localized = ordinary_cartesian_ida_operators(
    basis;
    expansion = expansion,
    Z = Z,
    backend = :pgdg_localized_experimental,
)

println("Ordinary Cartesian IDA localized-backend comparison")
println("  basis: ", basis)
println("  mapping: ", mapping)
println("  expansion terms: ", length(expansion))
println("  reference 1D overlap error: ", norm(reference_1d.overlap - I, Inf))
println("  pre-COMX proxy 1D overlap error: ", norm(proxy_1d.overlap - I, Inf))
println("  localized proxy 1D overlap error: ", norm(localized_1d.overlap - I, Inf))
println("  reference 3D overlap error: ", norm(reference.overlap_3d - I, Inf))
println("  localized 3D overlap error: ", norm(localized.overlap_3d - I, Inf))
println("  ||H1_localized - H1_reference||_inf: ", norm(localized.one_body_hamiltonian - reference.one_body_hamiltonian, Inf))
println("  ||Vee_localized - Vee_reference||_inf: ", norm(localized.interaction_matrix - reference.interaction_matrix, Inf))
println("  localized operators: ", localized)
