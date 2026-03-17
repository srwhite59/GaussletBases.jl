using LinearAlgebra
using GaussletBases

function localized_reference_1d(
    basis::MappedUniformBasis,
    exponents::AbstractVector{<:Real},
)
    representation = basis_representation(basis; operators = (:overlap, :position, :kinetic))
    transform, _ = GaussletBases._cleanup_comx_transform(
        representation.basis_matrices.overlap,
        representation.basis_matrices.position,
        integral_weights(basis),
    )
    overlap = Matrix(transpose(transform) * representation.basis_matrices.overlap * transform)
    kinetic = Matrix(transpose(transform) * representation.basis_matrices.kinetic * transform)
    gaussian_factors = Matrix{Float64}[
        Matrix(transpose(transform) * gaussian_factor_matrix(basis; exponent = exponent, center = 0.0) * transform)
        for exponent in exponents
    ]
    return overlap, kinetic, gaussian_factors
end

function one_body_hamiltonian_3d(
    overlap_1d::AbstractMatrix{<:Real},
    kinetic_1d::AbstractMatrix{<:Real},
    gaussian_factors::AbstractVector{<:AbstractMatrix{<:Real}},
    expansion::CoulombGaussianExpansion;
    Z::Real,
)
    overlap = Matrix{Float64}(overlap_1d)
    kinetic = Matrix{Float64}(kinetic_1d)
    overlap_3d = kron(overlap, kron(overlap, overlap))
    hamiltonian =
        kron(kinetic, kron(overlap, overlap)) +
        kron(overlap, kron(kinetic, overlap)) +
        kron(overlap, kron(overlap, kinetic))

    for term in eachindex(expansion.coefficients)
        factor = Matrix{Float64}(gaussian_factors[term])
        hamiltonian .-= expansion.coefficients[term] .* kron(factor, kron(factor, factor))
    end

    return overlap_3d, hamiltonian
end

full_expansion = coulomb_gaussian_expansion(doacc = false)
expansion = CoulombGaussianExpansion(
    full_expansion.coefficients[1:3],
    full_expansion.exponents[1:3];
    del = full_expansion.del,
    s = full_expansion.s,
    c = full_expansion.c,
    maxu = full_expansion.maxu,
)

basis = build_basis(MappedUniformBasisSpec(:G10;
    count = 5,
    mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0),
    reference_spacing = 1.0,
))

prototype = GaussletBases.mapped_pgdg_logfit_prototype(basis)
localized_old = mapped_pgdg_localized(prototype)
localized_new = mapped_ordinary_one_body_operators(
    basis;
    exponents = expansion.exponents,
    center = 0.0,
    backend = :pgdg_localized_experimental,
)
localized_oracle = GaussletBases._mapped_ordinary_localized_oracle_operators(
    basis;
    exponents = expansion.exponents,
    center = 0.0,
)
reference_overlap, reference_kinetic, reference_factors = localized_reference_1d(
    basis,
    expansion.exponents,
)
old_overlap = overlap_matrix(localized_old)
old_kinetic = kinetic_matrix(localized_old)
old_factors = gaussian_factor_matrices(
    localized_old;
    exponents = expansion.exponents,
    center = 0.0,
)

_, h1_reference = one_body_hamiltonian_3d(
    reference_overlap,
    reference_kinetic,
    reference_factors,
    expansion;
    Z = 2.0,
)
_, h1_old = one_body_hamiltonian_3d(
    old_overlap,
    old_kinetic,
    old_factors,
    expansion;
    Z = 2.0,
)
_, h1_new = one_body_hamiltonian_3d(
    localized_new.overlap,
    localized_new.kinetic,
    localized_new.gaussian_factors,
    expansion;
    Z = 2.0,
)
_, h1_oracle = one_body_hamiltonian_3d(
    localized_oracle.overlap,
    localized_oracle.kinetic,
    localized_oracle.gaussian_factors,
    expansion;
    Z = 2.0,
)

vee_reference = ordinary_cartesian_ida_operators(
    basis;
    expansion = expansion,
    Z = 2.0,
    backend = :numerical_reference,
)
vee_corrected = ordinary_cartesian_ida_operators(
    basis;
    expansion = expansion,
    Z = 2.0,
    backend = :pgdg_localized_experimental,
)

old_factor_diff = maximum(
    norm(old_factors[index] - reference_factors[index], Inf) for index in eachindex(old_factors)
)
new_factor_diff = maximum(
    norm(localized_new.gaussian_factors[index] - reference_factors[index], Inf)
    for index in eachindex(localized_new.gaussian_factors)
)
oracle_factor_diff = maximum(
    norm(localized_oracle.gaussian_factors[index] - reference_factors[index], Inf)
    for index in eachindex(localized_oracle.gaussian_factors)
)

println("Ordinary PGDG one-body fidelity comparison")
println("  basis: ", basis)
println("  backend: ", localized_new)
println("  expansion terms: ", length(expansion))
println("  1D overlap error (old pure analytic): ", norm(old_overlap - reference_overlap, Inf))
println("  1D overlap error (new pure analytic): ", norm(localized_new.overlap - reference_overlap, Inf))
println("  1D kinetic error (old pure analytic): ", norm(old_kinetic - reference_kinetic, Inf))
println("  1D kinetic error (new pure analytic): ", norm(localized_new.kinetic - reference_kinetic, Inf))
println("  1D kinetic error (oracle): ", norm(localized_oracle.kinetic - reference_kinetic, Inf))
println("  max Gaussian-factor error (old pure analytic): ", old_factor_diff)
println("  max Gaussian-factor error (new pure analytic): ", new_factor_diff)
println("  max Gaussian-factor error (oracle): ", oracle_factor_diff)
println("  3D H1 error (old pure analytic): ", norm(h1_old - h1_reference, Inf))
println("  3D H1 error (new pure analytic): ", norm(h1_new - h1_reference, Inf))
println("  3D H1 error (oracle): ", norm(h1_oracle - h1_reference, Inf))
println("  3D overlap error (corrected localized): ", norm(vee_corrected.overlap_3d - vee_reference.overlap_3d, Inf))
println("  3D Vee error (corrected localized): ", norm(vee_corrected.interaction_matrix - vee_reference.interaction_matrix, Inf))
