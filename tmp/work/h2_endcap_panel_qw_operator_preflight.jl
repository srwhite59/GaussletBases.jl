using GaussletBases
using LinearAlgebra
using Printf

set_timing!(false)

function truncated_expansion(expansion::CoulombGaussianExpansion, nterms::Int)
    return CoulombGaussianExpansion(
        expansion.coefficients[1:nterms],
        expansion.exponents[1:nterms];
        del = expansion.del,
        s = expansion.s,
        c = expansion.c,
        maxu = expansion.maxu,
    )
end

function symmetry_error(matrix::AbstractMatrix)
    return norm(matrix - transpose(matrix), Inf)
end

function width_status(widths::AbstractMatrix)
    all(isnan, widths) && return "all_nan"
    all(isfinite, widths) && all(>(0.0), widths) && return "finite_positive"
    return "mixed_or_invalid"
end

function owner_set(operators)
    return sort!(collect(Set(operators.residual_nucleus_indices)))
end

function build_source_and_fixed(basis, bundles, expansion, policy::Symbol)
    source = GaussletBases._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = :z,
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        shared_shell_layer_policy = policy,
        shared_shell_endcap_panel_q = 4,
        shared_shell_endcap_panel_L = 4,
    )
    return source, GaussletBases._nested_fixed_block(source)
end

function build_operators(fixed_block, supplement, expansion, interaction_treatment::Symbol)
    return ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        expansion,
        interaction_treatment,
    )
end

function operator_summary(operators)
    ndimension = operators.gausslet_count + operators.residual_count
    return (
        dimension = ndimension,
        residual_count = operators.residual_count,
        overlap_error = norm(operators.overlap - I, Inf),
        h_symmetry_error = symmetry_error(operators.one_body_hamiltonian),
        v_symmetry_error = symmetry_error(operators.interaction_matrix),
        width_status = width_status(operators.residual_widths),
        owner_set = owner_set(operators),
    )
end

function print_source_summary(label::AbstractString, source, fixed_block)
    shared_columns = [size(layer.coefficient_matrix, 2) for layer in source.shared_shell_layers]
    layer_types = [nameof(typeof(layer)) for layer in source.shared_shell_layers]
    println(label)
    println("  shared_layer_types = ", layer_types)
    println("  shared_layer_columns = ", shared_columns)
    println("  fixed_block_size = ", size(fixed_block.coefficient_matrix))
    @printf("  fixed_overlap_error = %.6e\n", norm(fixed_block.overlap - I, Inf))
end

function print_operator_summary(label::AbstractString, summary)
    println(label)
    println("  dimension = ", summary.dimension)
    println("  residual_count = ", summary.residual_count)
    @printf("  overlap_error = %.6e\n", summary.overlap_error)
    @printf("  h_symmetry_error = %.6e\n", summary.h_symmetry_error)
    @printf("  v_symmetry_error = %.6e\n", summary.v_symmetry_error)
    println("  residual_width_status = ", summary.width_status)
    println("  residual_owner_set = ", summary.owner_set)
end

full_expansion = coulomb_gaussian_expansion(doacc = false)
expansion = truncated_expansion(full_expansion, 3)
basis = bond_aligned_homonuclear_qw_basis(
    family = :G10,
    bond_length = 1.4,
    core_spacing = 0.7,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
    bond_axis = :z,
)
bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
    "H",
    "cc-pVTZ",
    basis.nuclei;
    lmax = 0,
    max_width = 1.0,
)

default_source, default_fixed = build_source_and_fixed(
    basis,
    bundles,
    expansion,
    :complete_rectangular,
)
endcap_source, endcap_fixed = build_source_and_fixed(
    basis,
    bundles,
    expansion,
    :endcap_panel_owned,
)

default_mwg = build_operators(default_fixed, supplement, expansion, :mwg)

# Warm both endcap operator paths before the timing sanity measurement.
build_operators(endcap_fixed, supplement, expansion, :ggt_nearest)
build_operators(endcap_fixed, supplement, expansion, :mwg)
GC.gc()

nearest_timing = @timed build_operators(endcap_fixed, supplement, expansion, :ggt_nearest)
mwg_timing = @timed build_operators(endcap_fixed, supplement, expansion, :mwg)
endcap_nearest = nearest_timing.value
endcap_mwg = mwg_timing.value

println("H2 endcap/panel QW operator preflight")
println("  expansion_terms = ", length(expansion.coefficients))
println("  public_default_policy = :complete_rectangular")
println("  experimental_policy = :endcap_panel_owned")
println("  q = 4")
println("  L = 4")

print_source_summary("default source", default_source, default_fixed)
print_source_summary("endcap/panel source", endcap_source, endcap_fixed)

print_operator_summary("default MWG operators", operator_summary(default_mwg))
print_operator_summary("endcap nearest operators", operator_summary(endcap_nearest))
print_operator_summary("endcap MWG operators", operator_summary(endcap_mwg))

@printf("endcap nearest warm_time_seconds = %.6f\n", nearest_timing.time)
println("endcap nearest warm_alloc_bytes = ", nearest_timing.bytes)
@printf("endcap MWG warm_time_seconds = %.6f\n", mwg_timing.time)
println("endcap MWG warm_alloc_bytes = ", mwg_timing.bytes)

println("chemistry_validation = not_run")
