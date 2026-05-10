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
run_nearest = get(ENV, "GAUSSLETBASES_RUN_GGT_NEAREST", "0") == "1"
qiu_white_mwg = ordinary_cartesian_qiu_white_operators(
    source_basis,
    legacy;
    expansion = expansion,
    Z = z_value,
)
qiu_white_mwg_check = GaussletBases.ordinary_cartesian_1s2_check(qiu_white_mwg)
nearest_result, nearest_error = if run_nearest
    try
        qiu_white_nearest = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = expansion,
            Z = z_value,
            interaction_treatment = :ggt_nearest,
        )
        (
            (
                operators = qiu_white_nearest,
                check = GaussletBases.ordinary_cartesian_1s2_check(qiu_white_nearest),
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

println("Paper-faithful Qiu-White MWG path")
sample_count = min(3, qiu_white_mwg.residual_count)
println("  E1: ", qiu_white_mwg_check.orbital_energy)
println("  <Vee>: ", qiu_white_mwg_check.vee_expectation)
println("  difference from 1.25: ", qiu_white_mwg_check.vee_expectation - reference_value)
println("  representative residual centers:")
println(qiu_white_mwg.residual_centers[1:sample_count, :])
println("  representative residual widths:")
println(qiu_white_mwg.residual_widths[1:sample_count, :])
println()

println("Paper-faithful Qiu-White nearest / GGT fallback path")
if nearest_result === nothing
    if nearest_error === nothing
        println("  skipped by default")
        println("  set ENV[\"GAUSSLETBASES_RUN_GGT_NEAREST\"] = \"1\" to run the fallback/debug nearest-GGT branch")
    else
        println("  failed: ", sprint(showerror, nearest_error))
    end
else
    println("  E1: ", nearest_result.check.orbital_energy)
    println("  <Vee>: ", nearest_result.check.vee_expectation)
    println("  difference from 1.25: ", nearest_result.check.vee_expectation - reference_value)
end
