@testset "Mapped ordinary Cartesian 1D working representation uses localized Gaussian contract" begin
    mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
    basis_a = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 5,
            mapping = mapping,
            reference_spacing = 1.0,
        ),
    )
    basis_b = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 7,
            mapping = mapping,
            reference_spacing = 1.0,
        ),
    )

    rep_a = GaussletBases._mapped_ordinary_working_basis_representation(basis_a)
    rep_b = GaussletBases._mapped_ordinary_working_basis_representation(basis_b)
    S_AA = cross_overlap(rep_a, rep_a)
    S_BB = cross_overlap(rep_b, rep_b)
    S_AB = cross_overlap(rep_a, rep_b)
    I_A = Matrix{Float64}(I, size(S_AA, 1), size(S_AA, 2))
    I_B = Matrix{Float64}(I, size(S_BB, 1), size(S_BB, 2))

    @test all(primitive -> primitive isa Gaussian, primitives(primitive_set(rep_a)))
    @test all(primitive -> primitive isa Gaussian, primitives(primitive_set(rep_b)))
    @test norm(S_AA - I_A, Inf) < 1.0e-12
    @test norm(S_BB - I_B, Inf) < 1.0e-12
    @test maximum(svdvals(S_AB)) <= 1.0 + 1.0e-10

    fixture = _atomic_hybrid_he_same_parent_stress_fixture()
    @test all(
        primitive -> primitive isa Gaussian,
        primitives(primitive_set(fixture.source_rep.axis_representations.x)),
    )
    @test all(
        primitive -> primitive isa Gaussian,
        primitives(primitive_set(fixture.target_rep.axis_representations.x)),
    )
end
