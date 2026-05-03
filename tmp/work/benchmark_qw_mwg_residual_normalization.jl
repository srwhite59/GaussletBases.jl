using LinearAlgebra
using Printf

using GaussletBases

function _benchmark_expansion(nterms::Int)
    expansion = coulomb_gaussian_expansion(doacc = false)
    return CoulombGaussianExpansion(
        expansion.coefficients[1:nterms],
        expansion.exponents[1:nterms];
        del = expansion.del,
        s = expansion.s,
        c = expansion.c,
        maxu = expansion.maxu,
    )
end

function _nested_atomic_mwg_fixture()
    count = 9
    a = 0.25
    xmax = 8.0
    tail_spacing = 10.0
    endpoint = (count - 1) / 2
    s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
    source_basis = build_basis(MappedUniformBasisSpec(:G10;
        count = count,
        mapping = AsinhMapping(a = a, s = s, tail_spacing = tail_spacing),
        reference_spacing = 1.0,
    ))
    expansion = _benchmark_expansion(3)
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        source_basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    interval = 2:(length(source_basis) - 1)
    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        interval,
        interval,
        interval;
        retain_xy = (3, 2),
        retain_xz = (3, 2),
        retain_yz = (3, 2),
        term_coefficients = term_coefficients,
    )
    fixed_block = GaussletBases._nested_fixed_block(shell, bundle)
    supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
    return fixed_block, supplement, expansion
end

function _block_norms(nearest, mwg)
    fixed_range = 1:mwg.gausslet_count
    residual_range = (mwg.gausslet_count + 1):(mwg.gausslet_count + mwg.residual_count)
    fixed_delta = norm(
        mwg.interaction_matrix[fixed_range, fixed_range] -
        nearest.interaction_matrix[fixed_range, fixed_range],
        Inf,
    )
    nearest_fixed_residual = norm(nearest.interaction_matrix[fixed_range, residual_range], Inf)
    mwg_fixed_residual = norm(mwg.interaction_matrix[fixed_range, residual_range], Inf)
    nearest_residual_residual = norm(nearest.interaction_matrix[residual_range, residual_range], Inf)
    mwg_residual_residual = norm(mwg.interaction_matrix[residual_range, residual_range], Inf)
    return (
        fixed_delta = fixed_delta,
        nearest_fixed_residual = nearest_fixed_residual,
        mwg_fixed_residual = mwg_fixed_residual,
        fixed_residual_ratio = mwg_fixed_residual / nearest_fixed_residual,
        nearest_residual_residual = nearest_residual_residual,
        mwg_residual_residual = mwg_residual_residual,
        residual_residual_ratio = mwg_residual_residual / nearest_residual_residual,
    )
end

function main()
    set_timing!(false)
    set_timing_live!(false)
    fixed_block, supplement, expansion = _nested_atomic_mwg_fixture()

    GC.gc()
    nearest_timed = @timed ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement;
        expansion = expansion,
        Z = 2.0,
        interaction_treatment = :ggt_nearest,
    )
    GC.gc()
    mwg_timed = @timed ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement;
        expansion = expansion,
        Z = 2.0,
        interaction_treatment = :mwg,
    )
    nearest = nearest_timed.value
    mwg = mwg_timed.value
    norms = _block_norms(nearest, mwg)

    println("case=nested_atomic_he_cc_pvtz_count9_lmax0_nterms3")
    println("timing_kind=fresh_process")
    @printf(
        "nearest_s=%.6f nearest_alloc_mb=%.3f nearest_gc_s=%.6f\n",
        nearest_timed.time,
        nearest_timed.bytes / 1024^2,
        nearest_timed.gctime,
    )
    @printf(
        "mwg_s=%.6f mwg_alloc_mb=%.3f mwg_gc_s=%.6f\n",
        mwg_timed.time,
        mwg_timed.bytes / 1024^2,
        mwg_timed.gctime,
    )
    println("fixed_count=$(mwg.gausslet_count)")
    println("residual_count=$(mwg.residual_count)")
    println("interaction_dim=$(size(mwg.interaction_matrix, 1))")
    @printf("fixed_fixed_delta_inf=%.6e\n", norms.fixed_delta)
    @printf(
        "fixed_residual_nearest_inf=%.6e fixed_residual_mwg_inf=%.6e ratio=%.6f\n",
        norms.nearest_fixed_residual,
        norms.mwg_fixed_residual,
        norms.fixed_residual_ratio,
    )
    @printf(
        "residual_residual_nearest_inf=%.6e residual_residual_mwg_inf=%.6e ratio=%.6f\n",
        norms.nearest_residual_residual,
        norms.mwg_residual_residual,
        norms.residual_residual_ratio,
    )
end

main()
