using LinearAlgebra
using Printf

using GaussletBases

function high_l_basis_text()
    return "#BASIS SET: He repo-spdfghi\n" *
           "He    S\n" *
           "      1.2000000              1.0000000\n" *
           "He    P\n" *
           "      1.0500000              1.0000000\n" *
           "He    D\n" *
           "      0.9000000              1.0000000\n" *
           "He    F\n" *
           "      0.7500000              1.0000000\n" *
           "He    G\n" *
           "      0.6000000              1.0000000\n" *
           "He    H\n" *
           "      0.4500000              1.0000000\n" *
           "He    I\n" *
           "      0.3000000              1.0000000\n" *
           "END\n"
end

function write_high_l_basisfile()
    path, io = mktemp()
    write(io, high_l_basis_text())
    close(io)
    return path
end

function probe_case(path::AbstractString, lmax::Int)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.8, npoints = 5, xmax = 4.0),
        reference_spacing = 1.0,
    ))
    supplement = legacy_atomic_gaussian_supplement(
        "He",
        "repo-spdfghi";
        lmax,
        basisfile = path,
    )
    representation = basis_representation(supplement)
    overlap = gto_overlap_matrix(basis, supplement)
    occupancy = gto_occupancy_matrix(basis, supplement)
    return (
        lmax = lmax,
        orbital_count = length(representation.orbitals),
        overlap_size = size(overlap),
        occupancy_size = size(occupancy),
        overlap_checksum = norm(overlap),
        occupancy_trace = tr(occupancy),
    )
end

function qw_smoke_case(path::AbstractString, lmax::Int)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 3,
        mapping = fit_asinh_mapping_for_strength(s = 0.8, npoints = 3, xmax = 3.0),
        reference_spacing = 1.0,
    ))
    supplement = legacy_atomic_gaussian_supplement(
        "He",
        "repo-spdfghi";
        lmax,
        basisfile = path,
    )
    full_expansion = coulomb_gaussian_expansion(doacc = false)
    expansion = CoulombGaussianExpansion(
        full_expansion.coefficients[1:1],
        full_expansion.exponents[1:1];
        del = full_expansion.del,
        s = full_expansion.s,
        c = full_expansion.c,
        maxu = full_expansion.maxu,
    )
    operators = ordinary_cartesian_qiu_white_operators(
        basis,
        supplement;
        expansion,
        Z = 2.0,
        interaction_treatment = :ggt_nearest,
        residual_keep_policy = :near_null_only,
    )
    check = GaussletBases.ordinary_cartesian_1s2_check(operators)
    return (
        lmax = lmax,
        dimension = size(operators.overlap, 1),
        residual_count = operators.residual_count,
        overlap_error = check.overlap_error,
        orbital_energy = check.orbital_energy,
        vee = check.vee_expectation,
    )
end

function run_timed(label::AbstractString, build)
    GC.gc()
    timed = @timed build()
    result = timed.value
    println("case=$label")
    println("timing_kind=fresh_process")
    @printf("time_s=%.6f alloc_mb=%.3f gc_s=%.6f\n", timed.time, timed.bytes / 1024^2, timed.gctime)
    for field in propertynames(result)
        println("$(field)=$(getproperty(result, field))")
    end
    return nothing
end

function main()
    path = write_high_l_basisfile()
    probe_case(path, 3)
    qw_smoke_case(path, 2)
    println("warmup=probe_lmax3_qw_lmax2")
    run_timed("atomic_probe_lmax4", () -> probe_case(path, 4))
    run_timed("atomic_probe_lmax6", () -> probe_case(path, 6))
    run_timed("atomic_qw_smoke_lmax6", () -> qw_smoke_case(path, 6))
    return nothing
end

main()
