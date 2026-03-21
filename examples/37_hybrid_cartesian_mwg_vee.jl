using GaussletBases

function run_case(basis_name::String; count::Int = 11, s::Float64 = 0.6, xmax::Float64 = 6.0)
    source_basis = build_basis(MappedUniformBasisSpec(:G10;
        count = count,
        mapping = fit_asinh_mapping_for_strength(s = s, npoints = count, xmax = xmax),
        reference_spacing = 1.0,
    ))
    legacy = legacy_s_gaussian_data("He", basis_name)
    hybrid_basis = hybrid_mapped_ordinary_basis(
        source_basis;
        core_gaussians = legacy.gaussians,
        backend = :pgdg_localized_experimental,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    pure_operators = ordinary_cartesian_ida_operators(
        source_basis;
        expansion = expansion,
        Z = 2.0,
        backend = :pgdg_localized_experimental,
    )
    combined_operators = ordinary_cartesian_ida_operators(
        hybrid_basis;
        expansion = expansion,
        Z = 2.0,
        interaction_treatment = :combined_basis,
    )
    nearest_operators = ordinary_cartesian_ida_operators(
        hybrid_basis;
        expansion = expansion,
        Z = 2.0,
        interaction_treatment = :residual_gaussian_nearest,
    )
    mwg_operators = ordinary_cartesian_ida_operators(
        hybrid_basis;
        expansion = expansion,
        Z = 2.0,
        interaction_treatment = :residual_gaussian_mwg,
    )

    pure_check = GaussletBases.ordinary_cartesian_1s2_check(pure_operators)
    combined_check = GaussletBases.ordinary_cartesian_1s2_check(combined_operators)
    nearest_check = GaussletBases.ordinary_cartesian_1s2_check(nearest_operators)
    mwg_check = GaussletBases.ordinary_cartesian_1s2_check(mwg_operators)
    mwg_data = GaussletBases._hybrid_residual_gaussian_mwg_data(hybrid_basis)

    return (
        source_basis = source_basis,
        legacy = legacy,
        hybrid_basis = hybrid_basis,
        pure = pure_check,
        combined = combined_check,
        nearest = nearest_check,
        mwg = mwg_check,
        mwg_centers = mwg_data.residual_centers,
        mwg_widths = mwg_data.residual_widths,
    )
end

function print_case(basis_name::String)
    result = run_case(basis_name)
    reference_value = 1.25
    sample_count = min(3, length(result.mwg_widths))

    println("He ", basis_name, " s-primitives")
    println("  source basis: ", result.source_basis)
    println("  hybrid basis: ", result.hybrid_basis)
    println("  exponents used: ", result.legacy.primitive_exponents)
    println("  widths used: ", result.legacy.widths)
    println("  representative MWG residual centers: ", result.mwg_centers[1:sample_count])
    println("  representative MWG residual widths: ", result.mwg_widths[1:sample_count])
    println("  pure ordinary E1: ", result.pure.orbital_energy)
    println("  pure ordinary <Vee>: ", result.pure.vee_expectation)
    println("  hybrid combined E1: ", result.combined.orbital_energy)
    println("  hybrid combined <Vee>: ", result.combined.vee_expectation)
    println("  hybrid nearest residual <Vee>: ", result.nearest.vee_expectation)
    println("  hybrid MWG residual <Vee>: ", result.mwg.vee_expectation)
    println("  combined - 1.25: ", result.combined.vee_expectation - reference_value)
    println("  nearest - 1.25: ", result.nearest.vee_expectation - reference_value)
    println("  MWG - 1.25: ", result.mwg.vee_expectation - reference_value)
    println()
end

println("Hybrid ordinary Cartesian MWG residual-interaction check")
println("  count = 11")
println("  s = 0.6")
println("  xmax = 6.0")
println("  backend = :pgdg_localized_experimental")
println("  comparing :combined_basis, :residual_gaussian_nearest, and :residual_gaussian_mwg")
println()

for basis_name in ("cc-pVTZ", "cc-pVQZ")
    print_case(basis_name)
end
