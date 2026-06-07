# Integration/slow test. Do not include in default nested runner.

@testset "Bond-aligned diatomic atom-growth anatomy policy" begin
    parent29 = (1:29, 1:29, 1:29)
    recipe29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_recipe(
        parent29;
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    anatomy29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(recipe29)

    @test recipe29 isa GaussletBases._BondAlignedDiatomicAtomGrowthRecipe3D
    @test recipe29.contact_cap_policy == :single_shared_contact_cap
    @test recipe29.mismatch_absorption_policy == :outermost_shared_molecular_shell
    @test anatomy29.atom_side_count_ladder == [5, 7, 9]
    @test anatomy29.final_atom_side_count == 9
    @test anatomy29.contact_gap_count == 1
    @test anatomy29.contact_policy == :single_shared_contact_cap
    @test anatomy29.left_atom_box == (11:19, 11:19, 6:14)
    @test anatomy29.right_atom_box == (11:19, 11:19, 16:24)
    @test anatomy29.contact_box == (11:19, 11:19, 15:15)
    @test anatomy29.inner_atom_contact_box == (11:19, 11:19, 6:24)
    @test anatomy29.outer_regular_start_box == (6:24, 6:24, 1:29)
    @test anatomy29.regular_shared_shell_count == 5
    @test anatomy29.outer_mismatch_low_counts == (5, 5, 0)
    @test anatomy29.outer_mismatch_high_counts == (5, 5, 0)
    @test anatomy29.support_coverage.status == :full_parent_covered
    @test anatomy29.support_coverage.expected_support_count == 29^3
    @test anatomy29.support_coverage.atom_contact_support_count == 9 * 9 * 19
    @test anatomy29.support_coverage.shared_molecular_support_count == 29^3 - 9 * 9 * 19
    @test anatomy29.support_coverage.covered_support_count == 29^3
    @test anatomy29.support_coverage.duplicate_count == 0
    @test anatomy29.support_coverage.missing_count == 0
    @test anatomy29.support_coverage.outside_count == 0
    @test anatomy29.support_coverage.coverage_ok

    even_parent = (1:20, 1:20, 1:20)
    even_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        even_parent;
        bond_axis = :z,
        atom_axis_indices = (7, 14),
        protected_atom_side_count = 4,
    )

    @test even_anatomy.atom_side_count_ladder == [4, 6]
    @test even_anatomy.final_atom_side_count == 6
    @test even_anatomy.contact_gap_count == 0
    @test even_anatomy.contact_policy == :touching_atom_boxes
    @test isnothing(even_anatomy.contact_box)
    @test even_anatomy.left_atom_box == (8:13, 8:13, 5:10)
    @test even_anatomy.right_atom_box == (8:13, 8:13, 11:16)
    @test even_anatomy.left_atom_box[1:2] == even_anatomy.right_atom_box[1:2]
    @test 7 - first(even_anatomy.left_atom_box[3]) == 2
    @test last(even_anatomy.left_atom_box[3]) - 7 == 3
    @test 14 - first(even_anatomy.right_atom_box[3]) == 3
    @test last(even_anatomy.right_atom_box[3]) - 14 == 2
    @test even_anatomy.support_coverage.status == :full_parent_covered
    @test even_anatomy.support_coverage.coverage_ok

    asymmetric_parent = (1:11, 1:13, 1:25)
    asymmetric = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        asymmetric_parent;
        bond_axis = :z,
        atom_axis_indices = (8, 17),
        protected_atom_side_count = 5,
    )

    @test asymmetric.contact_policy == :touching_atom_boxes
    @test asymmetric.inner_atom_contact_box == (2:10, 3:11, 4:21)
    @test asymmetric.outer_regular_start_box == (1:11, 2:12, 3:22)
    @test asymmetric.regular_shared_shell_count == 1
    @test asymmetric.outer_mismatch_low_counts == (0, 1, 2)
    @test asymmetric.outer_mismatch_high_counts == (0, 1, 3)
    @test asymmetric.recipe.mismatch_absorption_policy == :outermost_shared_molecular_shell
    @test asymmetric.support_coverage.expected_support_count == 11 * 13 * 25
    @test asymmetric.support_coverage.covered_support_count == 11 * 13 * 25
    @test asymmetric.support_coverage.coverage_ok

    expansion = coulomb_gaussian_expansion(doacc = false)
    real_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 3.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    real_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(real_basis, expansion)
    real_recipe = GaussletBases._nested_bond_aligned_diatomic_atom_growth_recipe(
        real_basis,
        real_bundles;
        protected_atom_side_count = 3,
    )
    real_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(real_recipe)
    @test real_recipe.parent_box == (1:7, 1:7, 1:13)
    @test real_recipe.atom_axis_indices == (5, 9)
    @test real_anatomy.contact_policy == :single_shared_contact_cap
    @test real_anatomy.support_coverage.coverage_ok

    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        (1:7, 1:7, 1:13);
        bond_axis = :z,
        atom_axis_indices = (3, 11),
        protected_atom_side_count = 3,
    )
end
