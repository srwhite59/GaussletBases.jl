using LinearAlgebra
using Printf

using GaussletBases

function high_l_basis_text()
    return "#BASIS SET: H repo-sgi\n" *
           "H    S\n" *
           "      2.0000000              1.0000000\n" *
           "      0.1250000              0.2500000\n" *
           "H    G\n" *
           "      2.0000000              1.0000000\n" *
           "      0.1250000              0.2500000\n" *
           "H    I\n" *
           "      2.0000000              1.0000000\n" *
           "      0.1250000              0.2500000\n" *
           "END\n"
end

function write_high_l_basisfile()
    path, io = mktemp()
    write(io, high_l_basis_text())
    close(io)
    return path
end

function truncated_expansion(nterms::Int)
    full = coulomb_gaussian_expansion(doacc = false)
    return CoulombGaussianExpansion(
        full.coefficients[1:nterms],
        full.exponents[1:nterms];
        del = full.del,
        s = full.s,
        c = full.c,
        maxu = full.maxu,
    )
end

function h2_basis()
    return bond_aligned_homonuclear_qw_basis(
        bond_length = 2.0,
        core_spacing = 1.0,
        xmax_parallel = 0.5,
        xmax_transverse = 0.5,
        bond_axis = :z,
    )
end

function high_l_supplement(path::AbstractString, basis, lmax::Int)
    return legacy_bond_aligned_diatomic_gaussian_supplement(
        "H",
        "repo-sgi",
        basis.nuclei;
        lmax,
        basisfile = path,
        max_width = 1.0,
    )
end

function owner_locality_summary(operators, basis)
    distance(a, b) = sqrt(sum((a[axis] - b[axis])^2 for axis in 1:3))
    if operators.residual_count == 0
        return (
            max_owner_distance = 0.0,
            min_midpoint_distance = Inf,
            owner_indices = Int[],
            owner_assignment_matches_nearest = true,
            finite_widths = all(isfinite, operators.residual_widths),
            positive_widths = all(>(0.0), vec(operators.residual_widths)),
        )
    end
    owner_distances = Float64[]
    midpoint_distances = Float64[]
    nearest_owner_matches = Bool[]
    for (index, owner) in pairs(operators.residual_nucleus_indices)
        center = Tuple(operators.residual_centers[index, :])
        push!(owner_distances, distance(center, basis.nuclei[owner]))
        push!(midpoint_distances, distance(center, (0.0, 0.0, 0.0)))
        other = owner == 1 ? 2 : 1
        push!(nearest_owner_matches, distance(center, basis.nuclei[owner]) < distance(center, basis.nuclei[other]))
    end
    return (
        max_owner_distance = maximum(owner_distances),
        min_midpoint_distance = minimum(midpoint_distances),
        owner_indices = sort(unique(operators.residual_nucleus_indices)),
        owner_assignment_matches_nearest = all(nearest_owner_matches),
        finite_widths = all(isfinite, operators.residual_widths),
        positive_widths = all(>(0.0), vec(operators.residual_widths)),
    )
end

function qiu_white_case(path::AbstractString, lmax::Int, treatment::Symbol)
    basis = h2_basis()
    supplement = high_l_supplement(path, basis, lmax)
    supplement3d = GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(supplement)
    expansion = truncated_expansion(1)
    operators = ordinary_cartesian_qiu_white_operators(
        basis,
        supplement;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        expansion,
        interaction_treatment = treatment,
    )
    locality = owner_locality_summary(operators, basis)
    return (
        lmax = lmax,
        interaction_treatment = operators.interaction_treatment,
        supplement_orbital_count = length(supplement3d.orbitals),
        dimension = size(operators.overlap, 1),
        residual_count = operators.residual_count,
        residual_owner_indices = locality.owner_indices,
        max_owner_distance = locality.max_owner_distance,
        min_midpoint_distance = locality.min_midpoint_distance,
        owner_assignment_matches_nearest = locality.owner_assignment_matches_nearest,
        finite_widths = locality.finite_widths,
        positive_widths = locality.positive_widths,
        overlap_error = norm(operators.overlap - I, Inf),
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
    qiu_white_case(path, 4, :ggt_nearest)
    println("warmup=diatomic_qw_lmax4_ggt_nearest")
    for lmax in (4, 6), treatment in (:ggt_nearest, :mwg)
        run_timed("diatomic_qw_lmax$(lmax)_$(treatment)", () -> qiu_white_case(path, lmax, treatment))
    end
    return nothing
end

main()
