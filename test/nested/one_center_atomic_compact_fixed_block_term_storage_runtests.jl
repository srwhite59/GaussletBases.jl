@testset "One-center atomic compact fixed-block term storage" begin
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)

    compact_full = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    compact_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion,
        working_box = 2:12,
        nside = 5,
    )

    @test !hasproperty(compact_full, :term_storage)
    @test !hasproperty(compact_full, :gaussian_terms)
    @test !hasproperty(compact_full, :pair_terms)
    @test !isnothing(compact_full.gaussian_sum)
    @test !isnothing(compact_full.pair_sum)

    @test !hasproperty(compact_legacy, :term_storage)
    @test !hasproperty(compact_legacy, :gaussian_terms)
    @test !hasproperty(compact_legacy, :pair_terms)
    @test !isnothing(compact_legacy.gaussian_sum)
    @test !isnothing(compact_legacy.pair_sum)

    full_contract = GaussletBases._nested_glass_box_contract(
        one_center_atomic_nested_structure_diagnostics(compact_full; nside = 5),
    )
    legacy_contract = GaussletBases._nested_glass_box_contract(
        one_center_atomic_nested_structure_diagnostics(compact_legacy; nside = 5),
    )
    @test full_contract.fixed_dimension == size(compact_full.overlap, 1)
    @test legacy_contract.fixed_dimension == size(compact_legacy.overlap, 1)
    @test full_contract.leaf_count === nothing
    @test legacy_contract.leaf_count === nothing

end
