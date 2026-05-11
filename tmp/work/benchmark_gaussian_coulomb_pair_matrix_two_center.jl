using LinearAlgebra
using Printf

using GaussletBases

const C_CC_PVTZ_LEGACY_BASIS = """
#BASIS SET: C cc-pVTZ

C S
      8.236000e+03           5.310000e-04
      1.235000e+03           4.108000e-03
      2.808000e+02           2.108700e-02
      7.927000e+01           8.185300e-02
      2.559000e+01           2.348170e-01
      8.997000e+00           4.344010e-01
      3.319000e+00           3.461290e-01
      9.059000e-01           3.937800e-02
      3.643000e-01          -8.983000e-03
      1.285000e-01           2.385000e-03
C S
      9.059000e-01           1.000000e+00
C S
      8.236000e+03          -1.130000e-04
      1.235000e+03          -8.780000e-04
      2.808000e+02          -4.540000e-03
      7.927000e+01          -1.813300e-02
      2.559000e+01          -5.576000e-02
      8.997000e+00          -1.268950e-01
      3.319000e+00          -1.703520e-01
      9.059000e-01           1.403820e-01
      3.643000e-01           5.986840e-01
      1.285000e-01           3.953890e-01
C S
      1.285000e-01           1.000000e+00
C P
      3.827000e-01           1.000000e+00
C P
      1.871000e+01           1.403100e-02
      4.133000e+00           8.686600e-02
      1.200000e+00           2.902160e-01
      3.827000e-01           5.010080e-01
      1.209000e-01           3.434060e-01
C P
      1.209000e-01           1.000000e+00
C D
      1.097000e+00           1.000000e+00
C D
      3.180000e-01           1.000000e+00
C F
      7.610000e-01           1.000000e+00
END
"""

function truncated_expansion(nterms::Int)
    full = coulomb_gaussian_expansion(doacc = false)
    count = min(nterms, length(full.coefficients))
    return CoulombGaussianExpansion(
        full.coefficients[1:count],
        full.exponents[1:count];
        del = full.del,
        s = full.s,
        c = full.c,
        maxu = full.maxu,
    )
end

function c2_cc_pvtz_orbitals()
    path, io = mktemp()
    try
        write(io, C_CC_PVTZ_LEGACY_BASIS)
        close(io)
        nuclei = [(0.0, 0.0, -2.35), (0.0, 0.0, 2.35)]
        supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            "C",
            "cc-pVTZ",
            nuclei;
            lmax = 6,
            basisfile = path,
        )
        return basis_representation(supplement).orbitals
    finally
        isopen(io) && close(io)
        rm(path; force = true)
    end
end

function compressed_structure_counts(orbitals)
    compact_pairs, _ = GaussletBases._gaussian_coulomb_compact_pair_index(length(orbitals))
    compact_terms = Vector{Vector{GaussletBases._GaussianCoulombPairTerm3D}}(
        undef,
        length(compact_pairs),
    )
    for (compact_index, (p, q)) in pairs(compact_pairs)
        compact_terms[compact_index] = GaussletBases._gaussian_coulomb_pair_terms(
            orbitals[p],
            orbitals[q],
        )
    end
    compact_coefficients, term_descriptors =
        GaussletBases._gaussian_coulomb_global_term_coefficients(compact_terms)
    axis_terms, _ = GaussletBases._gaussian_coulomb_axis_term_indices(term_descriptors)
    return (
        compact_pair_count = length(compact_pairs),
        unique_term_count = length(term_descriptors),
        coefficient_entry_count = sum(length, compact_coefficients),
        unique_axis_term_count = length(axis_terms),
    )
end

function main()
    nterms = parse(Int, get(ENV, "GAUSSLETBASES_TWO_CENTER_GTO_BENCH_TERMS", "45"))
    orbitals = c2_cc_pvtz_orbitals()
    expansion = truncated_expansion(nterms)
    counts = compressed_structure_counts(orbitals)

    gaussian_coulomb_pair_matrix(orbitals; expansion, max_orbitals = nothing)
    GC.gc()
    timed = @timed gaussian_coulomb_pair_matrix(orbitals; expansion, max_orbitals = nothing)
    matrix = timed.value
    println("case=c2_cc_pvtz_70_orbitals")
    println("timing_kind=fresh_process_post_warmup")
    println("threads=", Threads.nthreads())
    println("orbitals=", length(orbitals))
    println("pair_matrix_shape=", size(matrix))
    println("expansion_terms=", length(expansion.coefficients))
    println("compact_pair_count=", counts.compact_pair_count)
    println("unique_term_count=", counts.unique_term_count)
    println("coefficient_entry_count=", counts.coefficient_entry_count)
    println("unique_axis_term_count=", counts.unique_axis_term_count)
    @printf("time_s=%.6f\n", timed.time)
    @printf("allocated_mb=%.3f\n", timed.bytes / 1024^2)
    @printf("gc_time_s=%.6f\n", timed.gctime)
    @printf("checksum_abs=%.17e\n", sum(abs, matrix))
    println("finite=", all(isfinite, matrix))
    println("symmetric=", isapprox(matrix, transpose(matrix); atol = 1.0e-10, rtol = 1.0e-10))
    return nothing
end

main()
