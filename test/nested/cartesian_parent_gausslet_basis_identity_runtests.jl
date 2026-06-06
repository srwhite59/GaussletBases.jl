@testset "Cartesian parent gausslet basis identity" begin
    CP = GaussletBases.CartesianParentGaussletBases
    atomic_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 7,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    atomic_parent = CP.CartesianParentGaussletBasis3D(atomic_axis)
    construction_allocated = @allocated CP.CartesianParentGaussletBasis3D(atomic_axis)

    diatomic_basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4;
        core_spacing = 0.5,
        xmax_parallel = 3.0,
        xmax_transverse = 2.0,
    )
    chain_basis = bond_aligned_homonuclear_chain_qw_basis(
        natoms = 3;
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_parallel = 2.5,
        xmax_transverse = 2.0,
    )
    square_basis = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 2;
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_in_plane = 2.5,
        xmax_transverse = 2.0,
    )
    diatomic_parent = CP.CartesianParentGaussletBasis3D(diatomic_basis)
    chain_parent = CP.CartesianParentGaussletBasis3D(chain_basis)
    square_parent = CP.CartesianParentGaussletBasis3D(square_basis)

    function _check_parent_index_contract(parent)
        dims = CP.parent_axis_counts(parent)
        states = (
            (1, 1, 1),
            (min(2, dims[1]), min(3, dims[2]), min(4, dims[3])),
            dims,
        )
        axes = CP.parent_axes(parent)
        for state in states
            flat = CP.parent_flat_index(parent, state...)
            @test flat == GaussletBases._cartesian_flat_index(state..., dims)
            @test CP.parent_unflat_index(parent, flat) ==
                GaussletBases._cartesian_unflat_index(flat, dims)
            @test CP.parent_center(parent, state) == (
                centers(axes.x)[state[1]],
                centers(axes.y)[state[2]],
                centers(axes.z)[state[3]],
            )
        end
    end

    @test CP.cartesian_parent_gausslet_basis(atomic_parent) === atomic_parent
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(atomic_axis)) ==
        CP.parent_axis_counts(atomic_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(diatomic_basis)) ==
        CP.parent_axis_counts(diatomic_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(chain_basis)) ==
        CP.parent_axis_counts(chain_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(square_basis)) ==
        CP.parent_axis_counts(square_parent)

    @test CP.parent_axes(atomic_parent).x === atomic_axis
    @test CP.axis_basis(atomic_parent, :y) === atomic_axis
    @test CP.parent_box(atomic_parent) == (1:7, 1:7, 1:7)
    @test CP.parent_axis_counts(atomic_parent) == (7, 7, 7)
    @test CP.parent_dimension(atomic_parent) == 7^3
    @test atomic_parent.axis_sharing == :shared_xyz
    @test atomic_parent.metadata.basis_family == :mapped_uniform_same_axis
    @test construction_allocated < 50_000
    @test fieldnames(typeof(atomic_parent)) == (:axes, :parent_box, :axis_sharing, :metadata)
    @test !hasproperty(atomic_parent, :gausslet_backend)
    @test !hasproperty(atomic_parent, :backend)
    @test !hasproperty(atomic_parent.metadata, :basis_centers)
    @test !hasproperty(atomic_parent.metadata, :parent_centers)

    @test CP.parent_axes(diatomic_parent).x === diatomic_basis.basis_x
    @test CP.parent_axes(diatomic_parent).z === diatomic_basis.basis_z
    @test diatomic_parent.axis_sharing == :shared_xy
    @test diatomic_parent.metadata.basis_family == :bond_aligned_diatomic
    @test diatomic_parent.metadata.bond_axis == diatomic_basis.bond_axis
    @test diatomic_parent.metadata.nuclei == diatomic_basis.nuclei
    @test diatomic_parent.metadata.nuclear_charges == diatomic_basis.nuclear_charges
    @test CP.parent_box(diatomic_parent) == (
        1:length(diatomic_basis.basis_x),
        1:length(diatomic_basis.basis_y),
        1:length(diatomic_basis.basis_z),
    )

    @test chain_parent.axis_sharing == :shared_xy
    @test chain_parent.metadata.basis_family == :bond_aligned_homonuclear_chain
    @test chain_parent.metadata.chain_axis == chain_basis.chain_axis
    @test chain_parent.metadata.chain_coordinates == chain_basis.chain_coordinates
    @test chain_parent.metadata.nuclei == chain_basis.nuclei

    @test square_parent.axis_sharing == :shared_xy
    @test square_parent.metadata.basis_family == :axis_aligned_homonuclear_square_lattice
    @test square_parent.metadata.lattice_size == square_basis.lattice_size
    @test square_parent.metadata.x_coordinates == square_basis.x_coordinates
    @test square_parent.metadata.y_coordinates == square_basis.y_coordinates

    for parent in (atomic_parent, diatomic_parent, chain_parent, square_parent)
        _check_parent_index_contract(parent)
        @test !hasproperty(parent, :gausslet_backend)
        @test !hasproperty(parent, :backend)
    end

    @test_throws ArgumentError CP.axis_basis(atomic_parent, :q)
    @test_throws ArgumentError CP.parent_flat_index(atomic_parent, 0, 1, 1)
    @test_throws ArgumentError CP.parent_unflat_index(atomic_parent, CP.parent_dimension(atomic_parent) + 1)
    @test_throws ArgumentError CP.cartesian_parent_gausslet_basis((not_a_parent = true,))
end
