using GaussletBases
using LinearAlgebra
using Printf

function _timed(f, label)
    GC.gc()
    timing = @timed f()
    @printf("%-28s %9.5f s %9.2f MiB\n", label, timing.time, timing.bytes / 2.0^20)
    return timing.value, timing
end

function _cr_like_supplement()
    path, io = mktemp()
    try
        write(
            io,
            "#BASIS SET: Cr repo-cr-sp\n" *
            "Cr    S\n" *
            "      4.0000000              1.0000000\n" *
            "Cr    P\n" *
            "      2.5000000              1.0000000\n" *
            "END\n",
        )
        close(io)
        return legacy_atomic_gaussian_supplement("Cr", "repo-cr-sp"; lmax = 1, basisfile = path)
    finally
        isopen(io) && close(io)
        rm(path; force = true)
    end
end

expansion = coulomb_gaussian_expansion(doacc = false)
basis = build_basis(MappedUniformBasisSpec(:G10;
    count = 5,
    mapping = white_lindsey_atomic_mapping(Z = 24.0, d = 0.05, tail_spacing = 10.0),
    reference_spacing = 1.0,
))
supplement = _cr_like_supplement()

pgdg_fixed, pgdg_fixed_timing = _timed("pgdg fixed block") do
    one_center_atomic_full_parent_fixed_block(basis; expansion, nside = 3)
end
numerical_fixed, numerical_fixed_timing = _timed("numerical fixed block") do
    one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 3,
        gausslet_backend = :numerical_reference,
    )
end
pgdg_ops, pgdg_ops_timing = _timed("pgdg atomic QW MWG") do
    ordinary_cartesian_qiu_white_operators(
        pgdg_fixed,
        supplement;
        expansion,
        Z = 24.0,
        interaction_treatment = :mwg,
        gausslet_backend = :pgdg_localized_experimental,
    )
end
auto_ops, auto_ops_timing = _timed("auto atomic QW MWG") do
    ordinary_cartesian_qiu_white_operators(
        pgdg_fixed,
        supplement;
        expansion,
        Z = 24.0,
        interaction_treatment = :mwg,
    )
end
reference_ops, reference_ops_timing = _timed("numerical QW MWG") do
    ordinary_cartesian_qiu_white_operators(
        numerical_fixed,
        supplement;
        expansion,
        Z = 24.0,
        interaction_treatment = :mwg,
        gausslet_backend = :numerical_reference,
    )
end

function _symmetry_error(matrix)
    return norm(matrix - transpose(matrix), Inf)
end

@printf("\nCr-like atomic QW PGDG backend summary\n")
@printf("fixed backend default       %s\n", pgdg_fixed.gausslet_backend)
@printf("fixed backend reference     %s\n", numerical_fixed.gausslet_backend)
@printf("operator backend explicit   %s\n", pgdg_ops.gausslet_backend)
@printf("operator backend auto       %s\n", auto_ops.gausslet_backend)
@printf("operator backend reference  %s\n", reference_ops.gausslet_backend)
@printf("fixed dimension pgdg/ref    %d / %d\n", size(pgdg_fixed.overlap, 1), size(numerical_fixed.overlap, 1))
@printf("operator dimension pgdg/ref %d / %d\n", size(pgdg_ops.overlap, 1), size(reference_ops.overlap, 1))
@printf("residual count pgdg/ref     %d / %d\n", pgdg_ops.residual_count, reference_ops.residual_count)
@printf("pgdg overlap error          %.6e\n", norm(pgdg_ops.overlap - I, Inf))
@printf("pgdg H symmetry             %.6e\n", _symmetry_error(pgdg_ops.one_body_hamiltonian))
@printf("pgdg V symmetry             %.6e\n", _symmetry_error(pgdg_ops.interaction_matrix))
@printf("pgdg finite centers         %s\n", all(isfinite, pgdg_ops.residual_centers))
@printf("pgdg finite positive widths %s\n", all(isfinite, pgdg_ops.residual_widths) && all(>(0.0), vec(pgdg_ops.residual_widths)))
@printf("overlap diff vs numerical   %.6e\n", norm(pgdg_ops.overlap - reference_ops.overlap, Inf))
@printf("one-body diff vs numerical  %.6e\n", norm(pgdg_ops.one_body_hamiltonian - reference_ops.one_body_hamiltonian, Inf))
