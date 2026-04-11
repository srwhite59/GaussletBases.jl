using GaussletBases

basis_name = "cc-pVTZ"
count = 9
s = 0.8
xmax = 6.0
z_value = 2.0
reference_value = 1.25

source_basis = build_basis(MappedUniformBasisSpec(:G10;
    count = count,
    mapping = fit_asinh_mapping_for_strength(s = s, npoints = count, xmax = xmax),
))
legacy = legacy_s_gaussian_data("He", basis_name)
expansion = coulomb_gaussian_expansion(doacc = false)
run_experimental_mwg = get(ENV, "GAUSSLETBASES_RUN_EXPERIMENTAL_MWG", "0") == "1"
qiu_white_nearest = ordinary_cartesian_qiu_white_operators(
    source_basis,
    legacy;
    expansion = expansion,
    Z = z_value,
    interaction_treatment = :ggt_nearest,
)
nearest_check = GaussletBases.ordinary_cartesian_1s2_check(qiu_white_nearest)
mwg_result, mwg_error = if run_experimental_mwg
    try
        qiu_white_mwg = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = expansion,
            Z = z_value,
            interaction_treatment = :mwg,
        )
        (
            (
                operators = qiu_white_mwg,
                check = GaussletBases.ordinary_cartesian_1s2_check(qiu_white_mwg),
            ),
            nothing,
        )
    catch error
        (nothing, error)
    end
else
    (nothing, nothing)
end

println("Qiu-White residual-Gaussian reference check")
println("  basis: ", source_basis)
println("  He basis supplement: ", legacy.basis_name)
println("  primitive exponents: ", legacy.primitive_exponents)
println("  primitive widths: ", legacy.primitive_widths)
println("  contraction matrix size: ", size(legacy.contraction_matrix))
println("  contracted widths: ", legacy.widths)
println("  count = ", count, ", s = ", s, ", xmax = ", xmax)
println("  note: this is a slow full-expansion light reference run, not a quick smoke check.")
println("  hydrogenic 1s^2 target: ", reference_value)
println()

println("Paper-faithful Qiu-White nearest / GGT path")
println("  E1: ", nearest_check.orbital_energy)
println("  <Vee>: ", nearest_check.vee_expectation)
println("  difference from 1.25: ", nearest_check.vee_expectation - reference_value)
println()

println("Paper-faithful Qiu-White MWG path")
if mwg_result === nothing
    if mwg_error === nothing
        println("  skipped by default")
        println("  set ENV[\"GAUSSLETBASES_RUN_EXPERIMENTAL_MWG\"] = \"1\" to run the still-experimental MWG branch")
    else
        println("  failed: ", sprint(showerror, mwg_error))
    end
else
    sample_count = min(3, mwg_result.operators.residual_count)
    println("  E1: ", mwg_result.check.orbital_energy)
    println("  <Vee>: ", mwg_result.check.vee_expectation)
    println("  difference from 1.25: ", mwg_result.check.vee_expectation - reference_value)
    println("  representative residual centers:")
    println(mwg_result.operators.residual_centers[1:sample_count, :])
    println("  representative residual widths:")
    println(mwg_result.operators.residual_widths[1:sample_count, :])
end
