# Legacy/internal experimental example.
# The old 1D COMX-cleaned hybrid route is not part of the supported public
# GaussletBases workflow. This file is kept only for surrogate regression and
# historical comparison.

using GaussletBases

function report_case(count::Int, s::Float64)
    source_basis = build_basis(MappedUniformBasisSpec(:G10;
        count = count,
        mapping = fit_asinh_mapping_for_strength(s = s, npoints = count, xmax = 6.0),
    ))
    core_gaussians = [
        Gaussian(center = 0.0, width = 0.2),
        Gaussian(center = 0.0, width = 0.6),
    ]
    expansion = coulomb_gaussian_expansion(doacc = false)
    pure_ops = ordinary_cartesian_ida_operators(
        source_basis;
        expansion = expansion,
        Z = 2.0,
        backend = :pgdg_localized_experimental,
    )
    hybrid_basis = GaussletBases.hybrid_mapped_ordinary_basis(
        source_basis;
        core_gaussians = core_gaussians,
        backend = :pgdg_localized_experimental,
    )
    combined_ops = ordinary_cartesian_ida_operators(
        hybrid_basis;
        expansion = expansion,
        Z = 2.0,
    )
    residual_ops = ordinary_cartesian_ida_operators(
        hybrid_basis;
        expansion = expansion,
        Z = 2.0,
        interaction_treatment = :residual_gaussian_nearest,
    )

    pure_check = GaussletBases.ordinary_cartesian_1s2_check(pure_ops)
    combined_check = GaussletBases.ordinary_cartesian_1s2_check(combined_ops)
    residual_check = GaussletBases.ordinary_cartesian_1s2_check(residual_ops)
    reference_value = 1.25
    center_values = sort(Float64.(centers(source_basis)))
    midpoint = argmin(abs.(center_values))
    distances = abs.(center_values .- center_values[midpoint])
    filter!(>(1.0e-12), distances)
    dx_core = minimum(distances)

    println("count = ", count, ", s = ", s, ", dx_core = ", dx_core)
    println("  source basis: ", source_basis)
    println("  hybrid basis: ", hybrid_basis)
    println("  core Gaussian widths: ", [gaussian.width for gaussian in core_gaussians])
    println("  pure ordinary lowest one-body orbital energy: ", pure_check.orbital_energy)
    println("  pure ordinary 1s^2 IDA Vee: ", pure_check.vee_expectation)
    println("  combined hybrid lowest one-body orbital energy: ", combined_check.orbital_energy)
    println("  combined hybrid 1s^2 IDA Vee: ", combined_check.vee_expectation)
    println("  residual-Gaussian 1s^2 IDA Vee: ", residual_check.vee_expectation)
    println("  hydrogenic reference (5/8)Z: ", reference_value)
    println("  combined difference: ", combined_check.vee_expectation - reference_value)
    println("  residual difference: ", residual_check.vee_expectation - reference_value)
    println()
end

println("Hybrid ordinary Cartesian residual-Gaussian IDA check")
println("  mild mapping family with paper-like core spacing")
println("  backend: :pgdg_localized_experimental")
println("  interaction treatments: :combined_basis vs :residual_gaussian_nearest")
println()

for (count, s) in ((9, 0.8), (11, 0.5), (11, 0.6))
    report_case(count, s)
end
