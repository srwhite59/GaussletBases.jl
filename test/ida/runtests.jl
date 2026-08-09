include(joinpath(@__DIR__, "cartesian_ida_hamiltonian_runtests.jl"))

@testset "PQS source-box IDA PGDG factors" begin
    function _tiny_pgdg_source_box_ida_basis(count::Int)
        xmax = 3.0
        tail_spacing = 10.0
        a = 0.25
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        return build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = AsinhMapping(; a, s, tail_spacing),
            reference_spacing = 1.0,
        ))
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        _tiny_pgdg_source_box_ida_basis(7);
        exponents = expansion.exponents,
        backend = :pgdg_localized_experimental,
        refinement_levels = 0,
    )
    expected_term_count = length(expansion.coefficients)
    metrics_module = GaussletBases.CartesianContractedParentMetrics
    axis_factors = metrics_module._pqs_source_box_ida_axis_factors(
        bundle;
        axis = :x,
        expected_term_count,
    )

    stored_weight_divided_terms =
        axis_factors.density_normalized_pair_factor_terms
    raw_terms = axis_factors.raw_pair_factor_terms
    weights = axis_factors.source_weights
    n = length(weights)

    @test axis_factors.axis == :x
    @test axis_factors.term_count == expected_term_count
    @test axis_factors.factor_dimensions == (n, n)
    @test size(stored_weight_divided_terms) == (expected_term_count, n, n)
    @test size(raw_terms) == (expected_term_count, n, n)
    @test stored_weight_divided_terms === bundle.pgdg_intermediate.pair_factor_terms
    @test raw_terms === bundle.pgdg_intermediate.pair_factor_terms_raw
    @test weights == bundle.pgdg_intermediate.weights
    @test axis_factors.centers == bundle.pgdg_intermediate.centers
    @test all(isfinite, weights)
    @test all(weight -> abs(weight) > 1.0e-12, weights)

    reconstructed = similar(stored_weight_divided_terms)
    weight_outer = weights * transpose(weights)
    for term in 1:expected_term_count
        reconstructed[term, :, :] .= raw_terms[term, :, :] ./ weight_outer
    end
    # This is the historical PGDG storage relation. Retained IDA routes must
    # transform raw numerators and divide only by final retained density weights.
    @test reconstructed ≈ stored_weight_divided_terms atol = 1.0e-12 rtol = 1.0e-12
    @test_throws ArgumentError metrics_module._pqs_source_box_ida_axis_factors(
        bundle;
        axis = :bad,
        expected_term_count,
    )
    @test_throws DimensionMismatch metrics_module._pqs_source_box_ida_axis_factors(
        bundle;
        axis = :x,
        expected_term_count = expected_term_count + 1,
    )
end

@testset "Atomic IDA ingredients" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    nchannels = length(channels)
    radial_dim = length(rb)

    @test ida isa AtomicIDAOperators
    @test length(orbitals(ida)) == nchannels * radial_dim
    @test size(radial_multipole(ida, 1)) == (radial_dim, radial_dim)
    @test size(gaunt_tensor(ida, 1)) == (nchannels, nchannels, 3)
    @test size(angular_kernel(ida, 1)) == (nchannels, nchannels, nchannels, nchannels)
    @test channel_overlap(ida, YlmChannel(1, 0)) ≈ channel_overlap(atom, YlmChannel(1, 0)) atol = 1.0e-12 rtol = 1.0e-12
    @test channel_hamiltonian(ida, YlmChannel(2, 1)) ≈ channel_hamiltonian(atom, YlmChannel(2, 1)) atol = 1.0e-12 rtol = 1.0e-12
    @test channel_range(ida, YlmChannel(1, -1)) == channel_range(atom, YlmChannel(1, -1))

    q0 = angular_kernel(ida, 0)
    for alpha in 1:nchannels, alphap in 1:nchannels, beta in 1:nchannels, betap in 1:nchannels
        expected = (alpha == alphap && beta == betap) ? 1.0 : 0.0
        @test q0[alpha, alphap, beta, betap] ≈ expected atol = 1.0e-12 rtol = 1.0e-12
    end

    @test gaunt_coefficient(ida, 0, 0, YlmChannel(0, 0), YlmChannel(0, 0)) ≈ inv(sqrt(4 * pi)) atol = 1.0e-12 rtol = 1.0e-12
    @test gaunt_coefficient(ida, 1, 0, YlmChannel(0, 0), YlmChannel(0, 0)) == 0.0
    @test gaunt_coefficient(ida, 1, 1, YlmChannel(1, 0), YlmChannel(1, 0)) == 0.0
    @test abs(gaunt_coefficient(ida, 1, 0, YlmChannel(1, 0), YlmChannel(0, 0))) > 1.0e-12

    orbital_data = orbitals(ida)
    @test orbital_data[1].index == 1
    @test orbital_data[1].channel == YlmChannel(0, 0)
    @test orbital_data[1].radial_index == 1
    @test orbital_data[radial_dim + 1].channel == YlmChannel(1, -1)
    @test orbital_data[radial_dim + 1].radial_index == 1
end

@testset "Atomic full-IDA dense export" begin
    _, _, _, channels, _, ida = _quick_radial_atomic_fixture()
    nchannels = length(channels)
    radial_dim = size(ida.radial_operators.overlap, 1)
    norbitals = length(orbitals(ida))
    shell_centers_r = Float64[Float64(value) for value in ida.radial_operators.shell_centers_r]
    expected_channel_l = Int[channel.l for channel in channels]
    expected_channel_m = Int[channel.m for channel in channels]
    package_version = string(Base.pkgversion(GaussletBases))
    perm = GaussletBases._atomic_shell_major_permutation(ida)
    expected_h1 = Matrix{Float64}(ida.one_body.hamiltonian[perm, perm])
    expected_vee = GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[perm])
    expected_dims = fill(nchannels, radial_dim)
    expected_orders = collect(1:radial_dim)
    payload_data = fullida_dense_payload(
        ida;
        nelec = 2,
        meta = (example = "test_atomic_fullida_dense_export",),
    )

    mktempdir() do dir
        path = joinpath(dir, "atomic_fullida_dense_test.jld2")
        @test write_fullida_dense_jld2(
            path,
            ida;
            nelec = 2,
            meta = (example = "test_atomic_fullida_dense_export",),
        ) == path

        jldopen(path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "H1" in top_keys
            @test "Vee" in top_keys
            @test "dims_per_shell" in top_keys
            @test "orders" in top_keys
            @test "basis_centers" in top_keys
            @test "bridge" in top_keys
            @test "meta" in top_keys

            @test size(file["H1"]) == (norbitals, norbitals)
            @test size(file["Vee"]) == (norbitals, norbitals)
            @test size(file["basis_centers"]) == (norbitals, 3)
            @test Int.(file["dims_per_shell"]) == expected_dims
            @test Int.(file["orders"]) == expected_orders
            @test Matrix{Float64}(file["H1"]) ≈ expected_h1 atol = 0.0 rtol = 0.0
            @test Matrix{Float64}(file["Vee"]) ≈ expected_vee atol = 0.0 rtol = 0.0
            @test Matrix{Float64}(file["basis_centers"]) == zeros(norbitals, 3)

            @test String(file["bridge/format"]) == "fullida_dense_v1"
            @test Int(file["bridge/version"]) == 1
            @test String(file["bridge/site_type"]) == "Electron"
            @test String(file["bridge/interaction_model"]) == "density_density"
            @test String(file["bridge/model_detail"]) == "two_index_ida"
            @test String(file["bridge/source_branch"]) == "atomic_ida"
            @test String(file["bridge/onebody_key"]) == "H1"
            @test String(file["bridge/interaction_key"]) == "Vee"
            @test String(file["bridge/optional_interaction_key"]) == "Vps"
            @test !Bool(file["bridge/has_optional_interaction"])
            @test Int(file["bridge/norb"]) == norbitals
            @test Int(file["bridge/nelec"]) == 2
            @test Bool(file["bridge/has_nelec"])
            @test String(file["bridge/order/convention"]) == "shell_major_by_radial_index"
            @test String(file["bridge/order/within_shell"]) == "ylm_channel_order"
            @test Int.(file["bridge/order/dims_per_shell"]) == expected_dims
            @test Int.(file["bridge/order/shell_offsets"]) == collect(1:nchannels:(norbitals + 1))
            @test Int.(file["bridge/order/shell_index"]) == vcat((fill(shell, nchannels) for shell in 1:radial_dim)...)
            @test Int.(file["bridge/order/orders_per_shell"]) == expected_orders
            @test String(file["bridge/order/basis_centers_kind"]) == "origin_only_atomic_orbitals"
            @test Float64.(file["bridge/order/shell_centers_r"]) == shell_centers_r
            @test !any(isnan, Float64.(file["bridge/order/shell_centers_r"]))
            @test Float64(file["bridge/order/basis_radius"]) == maximum(shell_centers_r)
            @test Int.(file["bridge/order/permutation_from_in_memory"]) == perm
            @test file["bridge/order/shell_centers_r"] == payload_data.bridge_meta["order/shell_centers_r"]
            @test Matrix{Float64}(file["H1"]) == payload_data.payload["H1"]
            @test Matrix{Float64}(file["Vee"]) == payload_data.payload["Vee"]

            @test String(file["meta/producer"]) == "GaussletBases.write_fullida_dense_jld2"
            @test String(file["meta/producer_type"]) == "AtomicIDAOperators"
            @test Float64(file["meta/Z"]) == 2.0
            @test String(file["meta/source_branch"]) == "atomic_ida"
            @test String(file["meta/source_model"]) == "radial_atomic_ida"
            @test String(file["meta/interaction_model"]) == "density_density_ida"
            @test String(file["meta/interaction_detail"]) == "two_index_ida"
            @test Int(file["meta/nchannels"]) == nchannels
            @test Int(file["meta/nradial"]) == radial_dim
            @test String(file["meta/public_extent_kind"]) == "count"
            @test isnan(Float64(file["meta/public_rmax"]))
            @test Int(file["meta/public_count"]) == 6
            @test !Bool(file["meta/has_public_rmax"])
            @test Bool(file["meta/has_public_count"])
            @test String(file["meta/manifest/producer/package"]) == "GaussletBases"
            @test String(file["meta/manifest/producer/version"]) == package_version
            @test String(file["meta/manifest/producer/entrypoint"]) == "GaussletBases.write_fullida_dense_jld2"
            @test String(file["meta/manifest/producer/object_type"]) == "AtomicIDAOperators"
            @test String(file["meta/manifest/interaction/model"]) == "density_density_ida"
            @test String(file["meta/manifest/interaction/detail"]) == "two_index_ida"
            @test String(file["meta/manifest/source/branch"]) == "atomic_ida"
            @test String(file["meta/manifest/source/model"]) == "radial_atomic_ida"
            @test Float64(file["meta/manifest/source/atomic_charge"]) == 2.0
            @test String(file["meta/manifest/source/basis_spec_type"]) == "RadialBasisSpec"
            @test String(file["meta/manifest/source/basis_family"]) == "G10"
            @test String(file["meta/manifest/source/public_extent_kind"]) == "count"
            @test isnan(Float64(file["meta/manifest/source/public_rmax"]))
            @test Int(file["meta/manifest/source/public_count"]) == 6
            @test !Bool(file["meta/manifest/source/has_public_rmax"])
            @test Bool(file["meta/manifest/source/has_public_count"])
            @test Float64(file["meta/manifest/source/reference_spacing"]) == 1.0
            @test Int(file["meta/manifest/source/tails"]) == 3
            @test Int(file["meta/manifest/source/odd_even_kmax"]) == 2
            @test String(file["meta/manifest/source/supplement_kind"]) == "xgaussian"
            @test Int(file["meta/manifest/source/supplement_count"]) == 1
            @test Float64.(file["meta/manifest/source/supplement/xgaussian_alphas"]) == [0.2]
            @test String(file["meta/manifest/source/mapping/type"]) == "AsinhMapping"
            @test !Bool(file["meta/manifest/source/mapping/is_identity"])
            @test Float64(file["meta/manifest/source/mapping/a"]) == 1.0
            @test Float64(file["meta/manifest/source/mapping/c"]) == 0.15
            @test Float64(file["meta/manifest/source/mapping/s"]) == 0.15
            @test Float64(file["meta/manifest/source/mapping/tail_spacing"]) == 10.0
            @test Int(file["meta/manifest/source/radial_dimension"]) == radial_dim
            @test Int(file["meta/manifest/source/channel_count"]) == nchannels
            @test Int(file["meta/manifest/source/channel_lmax"]) == channels.lmax
            @test Int.(file["meta/manifest/source/channel_l"]) == expected_channel_l
            @test Int.(file["meta/manifest/source/channel_m"]) == expected_channel_m
            @test String(file["meta/manifest/source/channel_convention"]) == "ylm_channels_increasing_l_then_m"
            @test String(file["meta/example"]) == "test_atomic_fullida_dense_export"
        end
    end
end

@testset "Atomic sliced Hamiltonian export" begin
    _, _, _, channels, _, ida = _quick_radial_atomic_fixture()
    nchannels = length(channels)
    radial_dim = size(ida.radial_operators.overlap, 1)
    norbitals = length(orbitals(ida))
    shell_centers_r = Float64[Float64(value) for value in ida.radial_operators.shell_centers_r]
    expected_channel_l = Int[channel.l for channel in channels]
    expected_channel_m = Int[channel.m for channel in channels]
    package_version = string(Base.pkgversion(GaussletBases))
    orbital_perm, channel_perm = GaussletBases._atomic_sliced_permutation(ida)
    ordered_channels = ida.one_body.channels.channel_data[channel_perm]
    expected_dims = fill(nchannels, radial_dim)
    expected_offs = collect(1:nchannels:(norbitals + 1))
    expected_m = Int[channel.m for channel in ordered_channels]
    expected_l = Int[channel.l for channel in ordered_channels]
    expected_labels = String["r=1,l=$(channel.l),m=$(channel.m)" for channel in ordered_channels]
    components = GaussletBases._atomic_onebody_component_matrices(ida)
    expected_h1 = Matrix{Float64}(components.H1[orbital_perm, orbital_perm])
    expected_t = Matrix{Float64}(components.T[orbital_perm, orbital_perm])
    expected_vnuc = Matrix{Float64}(components.Vnuc[orbital_perm, orbital_perm])
    expected_vee = GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[orbital_perm])
    payload_data = sliced_ham_payload(
        ida;
        nelec = 2,
        meta = (example = "test_atomic_sliced_export",),
    )

    mktempdir() do dir
        path = joinpath(dir, "atomic_sliced_ham_test.jld2")
        @test write_sliced_ham_jld2(
            path,
            ida;
            nelec = 2,
            meta = (example = "test_atomic_sliced_export",),
        ) == path

        jldopen(path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "layout" in top_keys
            @test "basis" in top_keys
            @test "ordering" in top_keys
            @test "onebody" in top_keys
            @test "twobody" in top_keys
            @test "meta" in top_keys

            @test Int(file["layout/nslices"]) == radial_dim
            @test Int.(file["layout/dims"]) == expected_dims
            @test Int.(file["layout/offs"]) == expected_offs
            @test Float64.(file["layout/slice_coord"]) == shell_centers_r
            @test !any(isnan, Float64.(file["layout/slice_coord"]))
            @test Int.(file["layout/slice_index"]) == collect(1:radial_dim)
            @test file["layout/slice_coord"] == payload_data.layout_values["slice_coord"]

            m_by_slice = [Int.(collect(v)) for v in file["basis/m_by_slice"]]
            l_by_slice = [Int.(collect(v)) for v in file["basis/l_by_slice"]]
            labels_by_slice = [String.(collect(v)) for v in file["basis/labels_by_slice"]]
            @test length(m_by_slice) == radial_dim
            @test length(l_by_slice) == radial_dim
            @test length(labels_by_slice) == radial_dim
            @test m_by_slice[1] == expected_m
            @test l_by_slice[1] == expected_l
            @test Int.(file["basis/m_flat"]) == vcat(fill(expected_m, radial_dim)...)
            @test Int.(file["basis/l_flat"]) == vcat(fill(expected_l, radial_dim)...)
            @test labels_by_slice[1] == expected_labels

            @test String(file["ordering/major"]) == "slice_major"
            @test String(file["ordering/slice_meaning"]) == "radial_shell"
            @test String(file["ordering/within_slice"]) == "l0_desc_mzigzag"
            @test occursin("slice-major by radial index", String(file["ordering/description"]))
            @test Int.(file["ordering/permutation_from_in_memory"]) == orbital_perm

            @test String(file["onebody/stored"]) == "coo"
            @test Bool(file["onebody/is_hermitian"])
            @test occursin("centrifugal", String(file["onebody/decomposition"]))
            @test String(file["twobody/stored"]) == "coo_all"
            @test String(file["twobody/convention"]) == "density_density_pairdiag_v1"
            @test String(file["twobody/symmetry"]) == "pair_diagonal_density_density"
            @test occursin("density-density interaction model", String(file["twobody/description"]))

            H1blocks = file["onebody/H1blocks"]
            Tblocks = file["onebody/Tblocks"]
            Vnucblocks = file["onebody/Vnucblocks"]
            Vblocks = file["twobody/Vblocks"]
            @test length(H1blocks) == radial_dim
            @test length(H1blocks[1]) == radial_dim
            @test length(Vblocks) == radial_dim
            @test length(Vblocks[1]) == radial_dim

            @test GaussletBases._coo_blocks_to_dense(H1blocks, expected_dims) ≈ expected_h1 atol = 0.0 rtol = 0.0
            @test GaussletBases._coo_blocks_to_dense(Tblocks, expected_dims) ≈ expected_t atol = 0.0 rtol = 0.0
            @test GaussletBases._coo_blocks_to_dense(Vnucblocks, expected_dims) ≈ expected_vnuc atol = 0.0 rtol = 0.0
            @test GaussletBases._pairdiag_blocks_to_density_matrix(Vblocks, expected_dims) ≈ expected_vee atol = 0.0 rtol = 0.0

            @test String(file["meta/format"]) == "atomic_ida_sliced_v1"
            @test String(file["meta/consumer_shape"]) == "slicedmrgutils.HamIO"
            @test Float64(file["meta/Z"]) == 2.0
            @test String(file["meta/producer"]) == "GaussletBases.write_sliced_ham_jld2"
            @test String(file["meta/producer_type"]) == "AtomicIDAOperators"
            @test String(file["meta/source_branch"]) == "atomic_ida"
            @test String(file["meta/source_model"]) == "radial_atomic_ida"
            @test String(file["meta/interaction_model"]) == "density_density_ida"
            @test String(file["meta/interaction_detail"]) == "two_index_ida_pair_diagonal"
            @test String(file["meta/twobody_encoding"]) == "pair_diagonal_density_density"
            @test String(file["meta/slice_kind"]) == "radial_shell"
            @test String(file["meta/slice_coord_kind"]) == "physical_radial_center"
            @test String(file["meta/slice_index_kind"]) == "radial_index"
            @test String(file["meta/orbital_ordering"]) == "slice_major_by_radial_index_then_l0_desc_mzigzag"
            @test Int(file["meta/nchannels"]) == nchannels
            @test Int(file["meta/nradial"]) == radial_dim
            @test String(file["meta/public_extent_kind"]) == "count"
            @test isnan(Float64(file["meta/public_rmax"]))
            @test Int(file["meta/public_count"]) == 6
            @test !Bool(file["meta/has_public_rmax"])
            @test Bool(file["meta/has_public_count"])
            @test Int(file["meta/norb"]) == norbitals
            @test Int(file["meta/nelec"]) == 2
            @test Bool(file["meta/has_nelec"])
            @test Int.(file["meta/permutation_from_in_memory"]) == orbital_perm
            @test String(file["meta/manifest/producer/package"]) == "GaussletBases"
            @test String(file["meta/manifest/producer/version"]) == package_version
            @test String(file["meta/manifest/producer/entrypoint"]) == "GaussletBases.write_sliced_ham_jld2"
            @test String(file["meta/manifest/producer/object_type"]) == "AtomicIDAOperators"
            @test String(file["meta/manifest/interaction/model"]) == "density_density_ida"
            @test String(file["meta/manifest/interaction/detail"]) == "two_index_ida_pair_diagonal"
            @test String(file["meta/manifest/source/branch"]) == "atomic_ida"
            @test String(file["meta/manifest/source/model"]) == "radial_atomic_ida"
            @test Float64(file["meta/manifest/source/atomic_charge"]) == 2.0
            @test String(file["meta/manifest/source/basis_spec_type"]) == "RadialBasisSpec"
            @test String(file["meta/manifest/source/basis_family"]) == "G10"
            @test String(file["meta/manifest/source/public_extent_kind"]) == "count"
            @test isnan(Float64(file["meta/manifest/source/public_rmax"]))
            @test Int(file["meta/manifest/source/public_count"]) == 6
            @test !Bool(file["meta/manifest/source/has_public_rmax"])
            @test Bool(file["meta/manifest/source/has_public_count"])
            @test Float64(file["meta/manifest/source/reference_spacing"]) == 1.0
            @test Int(file["meta/manifest/source/tails"]) == 3
            @test Int(file["meta/manifest/source/odd_even_kmax"]) == 2
            @test String(file["meta/manifest/source/supplement_kind"]) == "xgaussian"
            @test Int(file["meta/manifest/source/supplement_count"]) == 1
            @test Float64.(file["meta/manifest/source/supplement/xgaussian_alphas"]) == [0.2]
            @test String(file["meta/manifest/source/mapping/type"]) == "AsinhMapping"
            @test !Bool(file["meta/manifest/source/mapping/is_identity"])
            @test Float64(file["meta/manifest/source/mapping/a"]) == 1.0
            @test Float64(file["meta/manifest/source/mapping/c"]) == 0.15
            @test Float64(file["meta/manifest/source/mapping/s"]) == 0.15
            @test Float64(file["meta/manifest/source/mapping/tail_spacing"]) == 10.0
            @test Int(file["meta/manifest/source/radial_dimension"]) == radial_dim
            @test Int(file["meta/manifest/source/channel_count"]) == nchannels
            @test Int(file["meta/manifest/source/channel_lmax"]) == channels.lmax
            @test Int.(file["meta/manifest/source/channel_l"]) == expected_channel_l
            @test Int.(file["meta/manifest/source/channel_m"]) == expected_channel_m
            @test String(file["meta/manifest/source/channel_convention"]) == "ylm_channels_increasing_l_then_m"
            @test String(file["meta/example"]) == "test_atomic_sliced_export"
        end
    end
end

@testset "Atomic HamV6 compatibility export and interaction accessor" begin
    _, _, _, channels, _, ida = _quick_radial_atomic_fixture()
    nchannels = length(channels)
    radial_dim = size(ida.radial_operators.overlap, 1)
    norbitals = length(orbitals(ida))
    shell_centers_r = Float64[Float64(value) for value in ida.radial_operators.shell_centers_r]
    expected_dims = fill(nchannels, radial_dim)
    expected_offs = collect(1:nchannels:(norbitals + 1))
    package_version = string(Base.pkgversion(GaussletBases))
    shell_perm = GaussletBases._atomic_shell_major_permutation(ida)
    native_perm, _ = GaussletBases._atomic_sliced_permutation(ida)
    orbital_perm, channel_perm = GaussletBases._atomic_hamv6_permutation(ida)
    ordered_channels = ida.one_body.channels.channel_data[channel_perm]
    expected_m = Int[channel.m for channel in ordered_channels]
    expected_l = Int[channel.l for channel in ordered_channels]
    expected_labels = String["r=1,l=$(channel.l),m=$(channel.m)" for channel in ordered_channels]
    components = GaussletBases._atomic_onebody_component_matrices(ida)
    expected_h1 = Matrix{Float64}(components.H1[orbital_perm, orbital_perm])
    expected_t = Matrix{Float64}(components.T[orbital_perm, orbital_perm])
    expected_vnuc = Matrix{Float64}(components.Vnuc[orbital_perm, orbital_perm])
    expected_vee = GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[orbital_perm])

    @test atomic_ida_density_interaction_matrix(ida) ≈
        GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)) atol = 0.0 rtol = 0.0
    @test atomic_ida_density_interaction_matrix(ida; ordering = :shell_major) ≈
        GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[shell_perm]) atol = 0.0 rtol = 0.0
    @test atomic_ida_density_interaction_matrix(ida; ordering = :sliced_native) ≈
        GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[native_perm]) atol = 0.0 rtol = 0.0
    @test atomic_ida_density_interaction_matrix(ida; ordering = :hamv6) ≈ expected_vee atol = 0.0 rtol = 0.0
    @test atomic_ida_density_interaction_matrix(ida; ordering = orbital_perm) ≈ expected_vee atol = 0.0 rtol = 0.0

    payload_data = atomic_hamv6_payload(
        ida;
        nelec = 2,
        meta = (example = "test_atomic_hamv6_export",),
    )

    mktempdir() do dir
        path = joinpath(dir, "atomic_hamv6_test.jld2")
        @test write_atomic_hamv6_jld2(
            path,
            ida;
            nelec = 2,
            meta = (example = "test_atomic_hamv6_export",),
        ) == path

        jldopen(path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "layout" in top_keys
            @test "basis" in top_keys
            @test "ordering" in top_keys
            @test "onebody" in top_keys
            @test "twobody" in top_keys
            @test "meta" in top_keys

            @test Int(file["layout/nslices"]) == radial_dim
            @test Int.(file["layout/dims"]) == expected_dims
            @test Int.(file["layout/offs"]) == expected_offs
            @test Float64.(file["layout/slice_coord"]) == shell_centers_r
            @test Int.(file["layout/slice_index"]) == collect(1:radial_dim)

            m_by_slice = [Int.(collect(v)) for v in file["basis/m_by_slice"]]
            l_by_slice = [Int.(collect(v)) for v in file["basis/l_by_slice"]]
            labels_by_slice = [String.(collect(v)) for v in file["basis/labels_by_slice"]]
            @test m_by_slice[1] == expected_m
            @test l_by_slice[1] == expected_l
            @test Int.(file["basis/m_flat"]) == vcat(fill(expected_m, radial_dim)...)
            @test Int.(file["basis/l_flat"]) == vcat(fill(expected_l, radial_dim)...)
            @test labels_by_slice[1] == expected_labels

            @test String(file["ordering/major"]) == "slice_major"
            @test String(file["ordering/slice_meaning"]) == "radial_shell"
            @test String(file["ordering/within_slice"]) == "mzigzag_then_l"
            @test occursin("for each m, l increasing", String(file["ordering/description"]))
            @test Int.(file["ordering/permutation_from_in_memory"]) == orbital_perm

            H1blocks = file["onebody/H1blocks"]
            Tblocks = file["onebody/Tblocks"]
            Vnucblocks = file["onebody/Vnucblocks"]
            Vblocks = file["twobody/Vblocks"]
            @test GaussletBases._coo_blocks_to_dense(H1blocks, expected_dims) ≈ expected_h1 atol = 0.0 rtol = 0.0
            @test GaussletBases._coo_blocks_to_dense(Tblocks, expected_dims) ≈ expected_t atol = 0.0 rtol = 0.0
            @test GaussletBases._coo_blocks_to_dense(Vnucblocks, expected_dims) ≈ expected_vnuc atol = 0.0 rtol = 0.0
            @test GaussletBases._pairdiag_blocks_to_density_matrix(Vblocks, expected_dims) ≈ expected_vee atol = 0.0 rtol = 0.0

            @test String(file["meta/format"]) == "atomic_hamv6_v1"
            @test String(file["meta/consumer_shape"]) == "slicedmrgutils.HamIO/HamV6"
            @test Float64(file["meta/Z"]) == 2.0
            @test String(file["meta/source_model"]) == "radial_atomic_ida"
            @test String(file["meta/interaction_model"]) == "density_density_ida"
            @test String(file["meta/interaction_detail"]) == "two_index_ida_pair_diagonal"
            @test String(file["meta/orbital_ordering"]) == "slice_major_by_radial_index_then_mzigzag_then_l"
            @test Int(file["meta/norb"]) == norbitals
            @test Int(file["meta/nelec"]) == 2
            @test String(file["meta/public_extent_kind"]) == "count"
            @test Int.(file["meta/permutation_from_in_memory"]) == orbital_perm
            @test String(file["meta/manifest/producer/package"]) == "GaussletBases"
            @test String(file["meta/manifest/producer/version"]) == package_version
            @test String(file["meta/manifest/producer/entrypoint"]) == "GaussletBases.write_atomic_hamv6_jld2"
            @test String(file["meta/example"]) == "test_atomic_hamv6_export"

            @test file["layout/slice_coord"] == payload_data.layout_values["slice_coord"]
            @test file["basis/m_flat"] == payload_data.basis_values["m_flat"]
            @test file["ordering/permutation_from_in_memory"] == payload_data.ordering_values["permutation_from_in_memory"]
            @test file["meta/Z"] == payload_data.meta_values["Z"]
        end
    end
end

@testset "Angular benchmark exact HamV6 bridge export" begin
    benchmark = _atomic_injected_angular_small_ed_benchmark_fixture()
    hf_payload = angular_benchmark_exact_hamv6_payload(
        benchmark.hf_style;
        nelec = 2,
        meta = (example = "test_angular_bridge_export",),
    )
    payload = angular_benchmark_exact_hamv6_payload(
        benchmark;
        nelec = 2,
        meta = (example = "test_angular_bridge_export",),
    )
    reference_payload = atomic_hamv6_payload(
        benchmark.hf_style.exact_ida_reference;
        nelec = 2,
        meta = (example = "test_angular_bridge_export",),
    )

    @test payload.layout_values == reference_payload.layout_values
    @test payload.basis_values == reference_payload.basis_values
    @test payload.ordering_values == reference_payload.ordering_values
    @test payload.onebody_values == reference_payload.onebody_values
    @test payload.twobody_values == reference_payload.twobody_values
    @test payload.layout_values == hf_payload.layout_values
    @test payload.basis_values == hf_payload.basis_values
    @test payload.ordering_values == hf_payload.ordering_values
    @test payload.onebody_values == hf_payload.onebody_values
    @test payload.twobody_values == hf_payload.twobody_values

    @test payload.meta_values["Z"] == 2.0
    @test payload.meta_values["consumer_shape"] == "slicedmrgutils.HamIO/HamV6"
    @test payload.meta_values["angular_bridge_kind"] == "exact_common_low_l_reference"
    @test payload.meta_values["angular_bridge_scope"] == "exact_common_low_l_reference_only"
    @test payload.meta_values["angular_bridge_consumer_language"] == "slicedmrgutils.HamIO/HamV6"
    @test payload.meta_values["angular_exact_common_lmax"] == benchmark.hf_style.one_body.exact_common_lmax
    @test payload.meta_values["angular_shell_orders"] == benchmark.hf_style.one_body.angular_assembly.shell_orders
    @test payload.meta_values["producer"] == "GaussletBases.write_angular_benchmark_exact_hamv6_jld2"
    @test payload.meta_values["producer_type"] == "AtomicInjectedAngularSmallEDBenchmark"
    @test hf_payload.meta_values["producer_type"] == "AtomicInjectedAngularHFStyleBenchmark"

    mktempdir() do dir
        path = joinpath(dir, "angular_exact_hamv6_bridge_test.jld2")
        @test write_angular_benchmark_exact_hamv6_jld2(
            path,
            benchmark;
            nelec = 2,
            meta = (example = "test_angular_bridge_export",),
        ) == path

        jldopen(path, "r") do file
            @test String(file["ordering/within_slice"]) == "mzigzag_then_l"
            @test String(file["meta/producer"]) == "GaussletBases.write_angular_benchmark_exact_hamv6_jld2"
            @test String(file["meta/producer_type"]) == "AtomicInjectedAngularSmallEDBenchmark"
            @test String(file["meta/angular_bridge_kind"]) == "exact_common_low_l_reference"
            @test String(file["meta/angular_bridge_scope"]) == "exact_common_low_l_reference_only"
            @test Int(file["meta/angular_exact_common_lmax"]) == benchmark.hf_style.one_body.exact_common_lmax
            @test Int.(file["meta/angular_shell_orders"]) == benchmark.hf_style.one_body.angular_assembly.shell_orders
            @test Float64(file["meta/Z"]) == 2.0
            @test file["layout/dims"] == payload.layout_values["dims"]
            @test file["basis/l_flat"] == payload.basis_values["l_flat"]
        end
    end
end

@testset "Atomic export source metadata supports rmax-based recipes" begin
    _, _, radial_ops, atom = _quick_hydrogen_ylm_fixture()
    ida = atomic_ida_operators(radial_ops; lmax = atom.channels.lmax)
    dense_payload = fullida_dense_payload(ida)
    sliced_payload = sliced_ham_payload(ida)

    @test String(dense_payload.meta_values["manifest/source/public_extent_kind"]) == "rmax"
    @test dense_payload.meta_values["Z"] == 1.0
    @test String(dense_payload.meta_values["source_model"]) == "radial_atomic_ida"
    @test String(dense_payload.meta_values["public_extent_kind"]) == "rmax"
    @test dense_payload.meta_values["public_rmax"] == 30.0
    @test dense_payload.meta_values["manifest/source/public_rmax"] == 30.0
    @test !Bool(dense_payload.meta_values["manifest/source/has_public_count"])
    @test Bool(dense_payload.meta_values["manifest/source/has_public_rmax"])
    @test String(dense_payload.meta_values["manifest/source/mapping/type"]) == "AsinhMapping"
    @test sliced_payload.meta_values["manifest/source/public_rmax"] == dense_payload.meta_values["manifest/source/public_rmax"]
    @test sliced_payload.meta_values["manifest/source/atomic_charge"] == dense_payload.meta_values["manifest/source/atomic_charge"]
    @test sliced_payload.meta_values["Z"] == dense_payload.meta_values["Z"]
    @test sliced_payload.meta_values["manifest/source/channel_l"] == dense_payload.meta_values["manifest/source/channel_l"]
    @test sliced_payload.meta_values["manifest/source/channel_m"] == dense_payload.meta_values["manifest/source/channel_m"]
end

@testset "Atomic IDA direct matrix" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    radial_dim = length(rb)
    nchannels = length(channels)
    norbitals = length(orbitals(ida))
    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    density = zeros(Float64, norbitals, norbitals)
    density[orbital_index(1, 1), orbital_index(1, 1)] = 1.2
    density[orbital_index(3, 1), orbital_index(3, 1)] = 0.4
    density[orbital_index(1, 1), orbital_index(3, 1)] = 0.3
    density[orbital_index(3, 1), orbital_index(1, 1)] = 0.3
    density[orbital_index(7, 2), orbital_index(7, 2)] = 0.5
    density[orbital_index(7, 2), orbital_index(3, 2)] = -0.15
    density[orbital_index(3, 2), orbital_index(7, 2)] = -0.15
    density[orbital_index(4, 1), orbital_index(8, 3)] = 0.9
    density[orbital_index(8, 3), orbital_index(4, 1)] = 0.9

    direct = direct_matrix(ida, density)
    reference = _dense_direct_reference(ida, density)

    @test size(direct) == (norbitals, norbitals)
    @test direct ≈ reference atol = 1.0e-12 rtol = 1.0e-12
    @test direct ≈ transpose(direct) atol = 1.0e-12 rtol = 1.0e-12

    radial_projected_density = copy(density)
    for left_channel in 1:nchannels, right_channel in 1:nchannels, left_radial in 1:radial_dim, right_radial in 1:radial_dim
        left_radial == right_radial && continue
        radial_projected_density[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)] = 0.0
    end
    @test direct ≈ direct_matrix(ida, radial_projected_density) atol = 1.0e-12 rtol = 1.0e-12

    for left_channel in 1:nchannels, right_channel in 1:nchannels, left_radial in 1:radial_dim, right_radial in 1:radial_dim
        left_radial == right_radial && continue
        @test direct[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)] == 0.0
    end

    selection_density = zeros(Float64, norbitals, norbitals)
    selection_density[orbital_index(3, 1), orbital_index(1, 1)] = 1.0
    selection_density[orbital_index(1, 1), orbital_index(3, 1)] = 1.0
    selection_direct = direct_matrix(ida, selection_density)
    target_delta_m = channels[1].m - channels[3].m

    for left_channel in 1:nchannels, right_channel in 1:nchannels, radial_index in 1:radial_dim
        delta_m = channels[left_channel].m - channels[right_channel].m
        value = selection_direct[orbital_index(left_channel, radial_index), orbital_index(right_channel, radial_index)]
        if delta_m != target_delta_m && delta_m != -target_delta_m
            @test abs(value) ≤ 1.0e-12
        end
    end
end

@testset "Atomic IDA exchange matrix" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    radial_dim = length(rb)
    nchannels = length(channels)
    norbitals = length(orbitals(ida))
    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    density = zeros(Float64, norbitals, norbitals)
    density[orbital_index(1, 1), orbital_index(1, 1)] = 0.9
    density[orbital_index(3, 1), orbital_index(1, 1)] = 0.25
    density[orbital_index(1, 1), orbital_index(3, 1)] = 0.25
    density[orbital_index(1, 1), orbital_index(3, 2)] = -0.2
    density[orbital_index(3, 2), orbital_index(1, 1)] = -0.2
    density[orbital_index(7, 3), orbital_index(4, 3)] = 0.35
    density[orbital_index(4, 3), orbital_index(7, 3)] = 0.35

    exchange = exchange_matrix(ida, density)
    reference = _dense_exchange_reference(ida, density)

    @test size(exchange) == (norbitals, norbitals)
    @test exchange ≈ reference atol = 1.0e-12 rtol = 1.0e-12
    @test exchange ≈ transpose(exchange) atol = 1.0e-12 rtol = 1.0e-12

    zeroed_density = copy(density)
    zeroed_density[orbital_index(1, 1), orbital_index(3, 2)] = 0.0
    zeroed_density[orbital_index(3, 2), orbital_index(1, 1)] = 0.0
    zeroed_exchange = exchange_matrix(ida, zeroed_density)
    @test exchange != zeroed_exchange

    selection_density = zeros(Float64, norbitals, norbitals)
    selection_density[orbital_index(3, 1), orbital_index(1, 2)] = 1.0
    selection_density[orbital_index(1, 2), orbital_index(3, 1)] = 1.0
    selection_exchange = exchange_matrix(ida, selection_density)
    target_mdiff = channels[3].m - channels[1].m

    for left_channel in 1:nchannels, right_channel in 1:nchannels, left_radial in 1:radial_dim, right_radial in 1:radial_dim
        value = selection_exchange[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)]
        if channels[left_channel].m - channels[right_channel].m != target_mdiff
            @test abs(value) ≤ 1.0e-12
        end
    end

    channel_index_by_lm(l::Int, m::Int) =
        findfirst(channel -> channel.l == l && channel.m == m, channels.channel_data)
    s_channel = channel_index_by_lm(0, 0)
    radial_index = min(3, radial_dim)
    s_orbital = orbital_index(s_channel, radial_index)
    for shell_l in 1:2
        shell_density = zeros(Float64, norbitals, norbitals)
        direct_values = Float64[]
        exchange_values = Float64[]
        for m in -shell_l:shell_l
            shell_channel = channel_index_by_lm(shell_l, m)
            shell_orbital = orbital_index(shell_channel, radial_index)
            density_m = zeros(Float64, norbitals, norbitals)
            density_m[shell_orbital, shell_orbital] = 1.0
            shell_density[shell_orbital, shell_orbital] = 1.0
            push!(direct_values, direct_matrix(ida, density_m)[s_orbital, s_orbital])
            push!(exchange_values, exchange_matrix(ida, density_m)[s_orbital, s_orbital])
        end
        @test direct_values ≈ fill(direct_values[1], length(direct_values)) atol = 1.0e-12 rtol = 1.0e-12
        @test exchange_values ≈ fill(exchange_values[1], length(exchange_values)) atol = 1.0e-12 rtol = 1.0e-12
        @test minimum(abs.(exchange_values)) > 1.0e-12
        @test direct_matrix(ida, shell_density)[s_orbital, s_orbital] ≈ sum(direct_values) atol = 1.0e-12 rtol = 1.0e-12
        @test exchange_matrix(ida, shell_density)[s_orbital, s_orbital] ≈ sum(exchange_values) atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testset "Atomic IDA Fock matrix" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    radial_dim = length(rb)
    norbitals = length(orbitals(ida))
    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    density = zeros(Float64, norbitals, norbitals)
    density[orbital_index(1, 1), orbital_index(1, 1)] = 0.9
    density[orbital_index(3, 1), orbital_index(1, 1)] = 0.25
    density[orbital_index(1, 1), orbital_index(3, 1)] = 0.25
    density[orbital_index(1, 1), orbital_index(3, 2)] = -0.2
    density[orbital_index(3, 2), orbital_index(1, 1)] = -0.2
    density[orbital_index(4, 2), orbital_index(4, 2)] = 0.4

    direct = direct_matrix(ida, density)
    exchange = exchange_matrix(ida, density)
    fock = fock_matrix(ida, density)
    reference = ida.one_body.hamiltonian + _dense_direct_reference(ida, density) - _dense_exchange_reference(ida, density)
    reference = 0.5 .* (reference .+ transpose(reference))

    @test size(fock) == size(ida.one_body.hamiltonian)
    @test fock ≈ ida.one_body.hamiltonian + direct - exchange atol = 1.0e-12 rtol = 1.0e-12
    @test fock ≈ reference atol = 1.0e-12 rtol = 1.0e-12
    @test fock ≈ transpose(fock) atol = 1.0e-12 rtol = 1.0e-12

    block = channel_range(ida, YlmChannel(0, 0))
    @test size(fock[block, block]) == (radial_dim, radial_dim)
end

@testset "Atomic IDA spin-aware Fock matrices" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    radial_dim = length(rb)
    norbitals = length(orbitals(ida))
    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    density_alpha = zeros(Float64, norbitals, norbitals)
    density_beta = zeros(Float64, norbitals, norbitals)

    density_alpha[orbital_index(1, 1), orbital_index(1, 1)] = 0.8
    density_alpha[orbital_index(3, 1), orbital_index(1, 1)] = 0.2
    density_alpha[orbital_index(1, 1), orbital_index(3, 1)] = 0.2
    density_alpha[orbital_index(1, 1), orbital_index(3, 2)] = -0.15
    density_alpha[orbital_index(3, 2), orbital_index(1, 1)] = -0.15

    density_beta[orbital_index(4, 2), orbital_index(4, 2)] = 0.5
    density_beta[orbital_index(2, 1), orbital_index(2, 1)] = 0.3
    density_beta[orbital_index(2, 1), orbital_index(4, 2)] = 0.1
    density_beta[orbital_index(4, 2), orbital_index(2, 1)] = 0.1

    total_density = density_alpha + density_beta
    direct_total = direct_matrix(ida, total_density)
    exchange_alpha = exchange_matrix(ida, density_alpha)
    exchange_beta = exchange_matrix(ida, density_beta)

    fock_alpha = fock_matrix_alpha(ida, density_alpha, density_beta)
    fock_beta = fock_matrix_beta(ida, density_alpha, density_beta)

    @test fock_alpha ≈ ida.one_body.hamiltonian + direct_total - exchange_alpha atol = 1.0e-12 rtol = 1.0e-12
    @test fock_beta ≈ ida.one_body.hamiltonian + direct_total - exchange_beta atol = 1.0e-12 rtol = 1.0e-12
    @test fock_alpha ≈ transpose(fock_alpha) atol = 1.0e-12 rtol = 1.0e-12
    @test fock_beta ≈ transpose(fock_beta) atol = 1.0e-12 rtol = 1.0e-12

    beta_shift = copy(density_beta)
    beta_shift[orbital_index(1, 2), orbital_index(1, 2)] += 0.4
    beta_shift[orbital_index(1, 2), orbital_index(3, 2)] -= 0.1
    beta_shift[orbital_index(3, 2), orbital_index(1, 2)] -= 0.1
    fock_alpha_shift = fock_matrix_alpha(ida, density_alpha, beta_shift)
    @test fock_alpha_shift - fock_alpha ≈ direct_matrix(ida, beta_shift - density_beta) atol = 1.0e-12 rtol = 1.0e-12

    alpha_shift = copy(density_alpha)
    alpha_shift[orbital_index(1, 3), orbital_index(1, 3)] += 0.25
    alpha_shift[orbital_index(1, 3), orbital_index(3, 3)] += 0.05
    alpha_shift[orbital_index(3, 3), orbital_index(1, 3)] += 0.05
    fock_beta_shift = fock_matrix_beta(ida, alpha_shift, density_beta)
    @test fock_beta_shift - fock_beta ≈ direct_matrix(ida, alpha_shift - density_alpha) atol = 1.0e-12 rtol = 1.0e-12

    closed_shell_density = density_alpha
    closed_alpha = fock_matrix_alpha(ida, closed_shell_density, closed_shell_density)
    closed_beta = fock_matrix_beta(ida, closed_shell_density, closed_shell_density)
    @test closed_alpha ≈ closed_beta atol = 1.0e-12 rtol = 1.0e-12
    @test closed_alpha ≈ ida.one_body.hamiltonian + direct_matrix(ida, 2.0 .* closed_shell_density) - exchange_matrix(ida, closed_shell_density) atol = 1.0e-12 rtol = 1.0e-12

    spinless_helper = fock_matrix(ida, density_alpha)
    @test spinless_helper ≈ ida.one_body.hamiltonian + direct_matrix(ida, density_alpha) - exchange_matrix(ida, density_alpha) atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Atomic IDA UHF" begin
    _rb, _grid, _radial_ops, ida, exact_problem, scf = _tiny_atomic_ida_uhf_fixture()
    exact_energy = ground_state_energy(exact_problem)
    norbitals = length(orbitals(ida))

    coeffs = [1.0 0.0; 0.0 1.0; 0.0 0.0]
    @test density_matrix(coeffs) ≈ [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0] atol = 1.0e-12 rtol = 1.0e-12
    @test density_matrix(view(coeffs, :, 1)) ≈ [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0] atol = 1.0e-12 rtol = 1.0e-12

    initial_alpha = zeros(Float64, norbitals, norbitals)
    initial_beta = zeros(Float64, norbitals, norbitals)
    initial_alpha[1, 1] = 1.0
    initial_beta[1, 1] = 1.0
    step = uhf_step(ida, initial_alpha, initial_beta; nalpha = 1, nbeta = 1)

    @test step.fock_alpha ≈ fock_matrix_alpha(ida, initial_alpha, initial_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test step.fock_beta ≈ fock_matrix_beta(ida, initial_alpha, initial_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test step.density_alpha ≈ density_matrix(step.occupied_coefficients_alpha) atol = 1.0e-12 rtol = 1.0e-12
    @test step.density_beta ≈ density_matrix(step.occupied_coefficients_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test step.energy ≈ uhf_energy(ida, step.density_alpha, step.density_beta) atol = 1.0e-12 rtol = 1.0e-12

    @test scf.converged
    @test 1 <= scf.iterations <= 80
    @test !isempty(scf.energies)
    @test length(scf.energies) == scf.iterations
    @test length(scf.residuals) == scf.iterations
    @test scf.fock_alpha ≈ transpose(scf.fock_alpha) atol = 1.0e-12 rtol = 1.0e-12
    @test scf.fock_beta ≈ transpose(scf.fock_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test scf.fock_alpha ≈ fock_matrix_alpha(ida, scf.density_alpha, scf.density_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test scf.fock_beta ≈ fock_matrix_beta(ida, scf.density_alpha, scf.density_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test scf.energy ≈ uhf_energy(ida, scf.density_alpha, scf.density_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test norm(scf.fock_alpha - scf.fock_beta, Inf) ≤ 1.0e-8
    @test norm(scf.density_alpha - scf.density_beta, Inf) ≤ 1.0e-8
    @test scf.energy >= exact_energy - 1.0e-9
    @test scf.energy < -2.7
end

@testset "Atomic IDA two-electron problem" begin
    rb, grid, radial_ops, ida, problem = _tiny_atomic_ida_two_electron_fixture()
    norbitals = length(orbitals(problem))
    states = two_electron_states(problem)

    @test problem isa AtomicIDATwoElectronProblem
    @test size(problem.overlap) == (norbitals^2, norbitals^2)
    @test size(problem.one_body) == size(problem.overlap)
    @test size(problem.two_body) == size(problem.overlap)
    @test size(problem.hamiltonian) == size(problem.overlap)
    @test length(states) == norbitals^2
    @test states[1].index == 1
    @test states[1].up_orbital == orbitals(problem)[1]
    @test states[1].down_orbital == orbitals(problem)[1]
    @test states[2].up_orbital == orbitals(problem)[1]
    @test states[2].down_orbital == orbitals(problem)[2]
    @test states[norbitals + 1].up_orbital == orbitals(problem)[2]
    @test states[norbitals + 1].down_orbital == orbitals(problem)[1]
    @test norm(problem.orbital_overlap - I, Inf) ≤ 1.0e-5
    @test norm(problem.overlap - I, Inf) ≤ 2.0e-5
    @test problem.hamiltonian ≈ transpose(problem.hamiltonian) atol = 1.0e-12 rtol = 1.0e-12
    @test problem.overlap ≈ transpose(problem.overlap) atol = 1.0e-12 rtol = 1.0e-12
    @test problem.two_body ≈ Diagonal(diag(problem.two_body)) atol = 1.0e-12 rtol = 1.0e-12

    coefficients = collect(range(1.0, length = length(states)))
    @test apply_hamiltonian(problem, coefficients) ≈ problem.hamiltonian * coefficients atol = 1.0e-12 rtol = 1.0e-12
    @test apply_overlap(problem, coefficients) ≈ problem.overlap * coefficients atol = 1.0e-12 rtol = 1.0e-12

    E0 = ground_state_energy(problem)
    @test isfinite(E0)
    @test -3.2 < E0 < -2.5
end

@testset "Atomic IDA Lanczos" begin
    _rb, _grid, _radial_ops, _ida, problem = _tiny_atomic_ida_lanczos_fixture()
    dense = ground_state_energy(problem)
    lanczos = lanczos_ground_state(problem; krylovdim = 200, maxiter = 200, tol = 1.0e-7)

    @test lanczos.converged
    @test lanczos.iterations <= 200
    @test lanczos.residual ≤ 1.0e-7
    @test abs(lanczos.value - dense) ≤ 1.0e-10
    @test abs(norm(lanczos.vector) - 1.0) ≤ 1.0e-10
    @test norm(problem.orbital_overlap - I, Inf) ≤ 5.0e-6
    @test norm(problem.overlap - I, Inf) ≤ 1.0e-5
    @test -2.95 < dense < -2.85
end

@testset "Minimal sliced hydrogen chain" begin
    function laplace_quadrature(integrand; step = 0.02, limit = 24.0)
        count = round(Int, 2limit / step)
        total = 0.0
        for index in 0:count
            s = -limit + index * step
            value = exp(s) * integrand(exp(2s))
            total += (index == 0 || index == count ? 0.5 : 1.0) * value
        end
        return 2step * total / sqrt(pi)
    end

    function direct_transverse_density_factor(chain, exponent, a, b)
        transverse = chain.transverse
        value = 0.0
        for p in eachindex(transverse.pair_exponents),
            q in eachindex(transverse.pair_exponents)
            ep = transverse.pair_exponents[p]
            eq = transverse.pair_exponents[q]
            value += transverse.density_coefficients[p, a] *
                transverse.density_coefficients[q, b] * pi^2 /
                (ep * eq + exponent * (ep + eq))
        end
        return value
    end

    function direct_vee(chain, i, j)
        longitudinal = chain.longitudinal
        a = longitudinal.class_indices[i]
        b = longitudinal.class_indices[j]
        delta = longitudinal.centers[i] - longitudinal.centers[j]
        return laplace_quadrature() do exponent
            direct_transverse_density_factor(chain, exponent, a, b) *
                GaussletBases._sliced_longitudinal_pair_factor(
                    longitudinal, delta, exponent) / longitudinal.integral_weight^2
        end
    end

    function direct_transverse_product_factor(chain, exponent, a, b, prefactors)
        transverse = chain.transverse
        value = 0.0
        for p in eachindex(transverse.exponents), q in eachindex(transverse.exponents)
            value += transverse.coefficients[p, a] * transverse.coefficients[q, b] *
                (prefactors[p] * prefactors[q])^2 * pi /
                (transverse.exponents[p] + transverse.exponents[q] + exponent)
        end
        return value
    end

    function direct_h1(chain, i, j)
        longitudinal = chain.longitudinal
        transverse = chain.transverse
        a = longitudinal.class_indices[i]
        b = longitudinal.class_indices[j]
        delta = longitudinal.centers[i] - longitudinal.centers[j]
        distance = abs(longitudinal.lattice_indices[i] - longitudinal.lattice_indices[j])
        midpoint = (longitudinal.centers[i] + longitudinal.centers[j]) / 2
        value = GaussletBases._sliced_longitudinal_pair(longitudinal, delta, true) *
            GaussletBases._sliced_transverse_bilinear(transverse, transverse.overlap, a, b) +
            GaussletBases._sliced_longitudinal_pair(longitudinal, delta, false) *
            GaussletBases._sliced_transverse_bilinear(transverse, transverse.kinetic, a, b)
        prefactors = [GaussletBases.GaussianAnalyticIntegrals.
            polynomial_gaussian_shell_prefactor(exponent, 0)
            for exponent in transverse.exponents]
        for nucleus in chain.nuclei
            value -= laplace_quadrature() do exponent
                direct_transverse_product_factor(chain, exponent, a, b, prefactors) *
                    GaussletBases._sliced_longitudinal_nuclear(
                        longitudinal, distance, midpoint, exponent, nucleus)
            end
        end
        return value
    end

    scalar_allocations(operator, row, column) = @allocated operator[row, column]
    row_allocations(destination, operator, row) = @allocated sliced_row!(destination, operator, row)
    action_allocations(destination, operator, source) = @allocated mul!(destination, operator, source)

    common = (; R = 3.6, phase = 0.0, left_padding = 2.0,
        alignment = :atom_centered, one_body_tolerance = 1e-8,
        interaction_tolerance = 1e-5)
    finite = sliced_hydrogen_chain(2; common..., spacing = 0.2,
        transverse_mode = :finite)
    periodic = sliced_hydrogen_chain(2; common..., spacing = 0.2,
        transverse_mode = :periodic_template)
    periodic_by_count = sliced_hydrogen_chain(2; common..., sites_per_atom = 18,
        transverse_mode = :periodic_template)
    singleton = sliced_hydrogen_chain(1; R = 3.6, phase = 0.0,
        left_padding = 0.0, right_padding = 0.0, alignment = :atom_centered,
        spacing = 0.2, transverse_mode = :finite,
        one_body_tolerance = 1e-8, interaction_tolerance = 1e-5)

    @test finite isa SlicedHydrogenChain
    @test length(finite) == length(periodic) == 39
    @test finite.nuclei == periodic.nuclei == [-1.8, 1.8]
    @test finite.sites_per_atom == periodic.sites_per_atom == 18
    @test finite.longitudinal.centers == periodic.longitudinal.centers
    @test periodic.longitudinal.centers == periodic_by_count.longitudinal.centers
    @test periodic.transverse.coefficients == periodic_by_count.transverse.coefficients
    @test periodic.h1_bands == periodic_by_count.h1_bands
    @test sliced_h1_bandwidth(singleton) == 0
    @test length(singleton) == 1
    @test periodic.longitudinal.prototype_errors[1] < 1e-12
    @test periodic.longitudinal.prototype_errors[2] < 1e-12
    @test periodic.longitudinal.prototype_errors[3] < 1e-12
    @test periodic.diagnostics.longitudinal_overlap_error < 1e-12
    @test periodic.diagnostics.transverse_normalization_error < 1e-12
    @test periodic.diagnostics.periodic_class_error < 1e-12
    @test periodic.diagnostics.omitted_projected_norm < 1e-7
    @test all(periodic.transverse.ranks .>= 1)
    @test all(periodic.transverse.gaps .> periodic.one_body_tolerance)
    @test all(periodic.transverse.coefficients[
        periodic.transverse.phase_pivots[class], class] > 0
        for class in axes(periodic.transverse.coefficients, 2))
    @test finite.coulomb.vee_tail_bound > 0
    @test finite.coulomb.vee_transition_error > 0
    @test finite.coulomb.vee_tail_bound + finite.coulomb.vee_transition_error <
        finite.interaction_tolerance

    finite_h1 = sliced_h1(finite)
    finite_vee = sliced_vee(finite)
    finite_mid = cld(length(finite), 2)
    finite_band = sliced_h1_bandwidth(finite)
    periodic_band = sliced_h1_bandwidth(periodic)
    @test 0 <= finite_band < length(finite)
    @test 0 <= periodic_band < length(periodic)
    @test finite_band < length(finite) - 1
    @test periodic_band < length(periodic) - 1
    @test all(finite_h1[i, j] == 0.0 for j in 1:length(finite),
        i in 1:length(finite) if abs(i - j) > finite_band)
    @test all(sliced_h1(periodic)[i, j] == 0.0 for j in 1:length(periodic),
        i in 1:length(periodic) if abs(i - j) > periodic_band)
    sliced_h1_bandwidth(finite)
    @test @allocated(sliced_h1_bandwidth(finite)) == 0
    @test finite_h1[finite_mid, finite_mid] ≈ direct_h1(finite, finite_mid, finite_mid) atol = 1e-8
    @test finite_vee[finite_mid, finite_mid] ≈ direct_vee(finite, finite_mid, finite_mid) atol = 1e-8
    @test finite_vee[finite_mid, finite_mid + 1] ≈ direct_vee(finite, finite_mid, finite_mid + 1) atol = 1e-8

    h10 = sliced_hydrogen_chain(10; common..., spacing = 0.2,
        transverse_mode = :periodic_template)
    h1 = sliced_h1(h10)
    vee = sliced_vee(h10)
    n = length(h10)
    mid = cld(n, 2)
    dense_h1 = Matrix(h1)
    dense_vee = Matrix(vee)
    @test n == 183
    @test all(isfinite, dense_h1)
    @test all(isfinite, dense_vee)
    @test dense_h1 ≈ dense_h1' atol = 1e-12 rtol = 0
    @test dense_vee ≈ dense_vee' atol = 1e-12 rtol = 0
    @test minimum(diag(dense_vee)) > 0
    @test h1[1, 1] ≈ direct_h1(h10, 1, 1) atol = 1e-8
    @test h1[mid, mid] ≈ direct_h1(h10, mid, mid) atol = 1e-8
    @test vee[mid, mid] ≈ direct_vee(h10, mid, mid) atol = 1e-8
    @test vee[mid, mid + 1] ≈ direct_vee(h10, mid, mid + 1) atol = 1e-8
    @test vee[1, n] ≈ direct_vee(h10, 1, n) atol = 1e-8
    @test vee[mid, mid + 18] == vee[mid + 18, mid + 36]
    @test h10.coulomb.h1_offband_bound + h10.coulomb.h1_tail_bound +
        h10.coulomb.h1_transition_error < h10.one_body_tolerance

    row = zeros(n)
    source = sin.(collect(1:n))
    destination = zeros(n)
    sliced_row!(row, h1, mid)
    @test row == dense_h1[mid, :]
    sliced_row!(row, vee, mid)
    @test row == dense_vee[mid, :]
    mul!(destination, h1, source)
    @test destination ≈ dense_h1 * source atol = 1e-11 rtol = 1e-12
    mul!(destination, vee, source)
    @test destination ≈ dense_vee * source atol = 1e-11 rtol = 1e-12
    @test scalar_allocations(h1, mid, mid) == 0
    @test scalar_allocations(vee, mid, mid) == 0
    @test row_allocations(row, h1, mid) == 0
    @test row_allocations(row, vee, mid) == 0
    @test action_allocations(destination, h1, source) == 0
    @test action_allocations(destination, vee, source) == 0

    h1000 = sliced_hydrogen_chain(1000; common..., spacing = 0.2,
        transverse_mode = :periodic_template)
    @test length(h1000) == 18_003
    @test Base.summarysize(h1000) < 10_000_000
    @test size(h1000.h1_bands, 2) == length(h1000)
    @test size(h1000.coulomb.periodic_blocks) == (18, 18, 1003)
    @test h1000.diagnostics.retained_bytes < Base.summarysize(h1000)
    @test maximum(h1000.longitudinal.centers) - minimum(h1000.longitudinal.centers) > 3500

    @test_throws ArgumentError sliced_hydrogen_chain(0; common..., spacing = 0.2)
    @test_throws ArgumentError sliced_hydrogen_chain(2; common..., spacing = 0.2,
        sites_per_atom = 18)
    @test_throws ArgumentError sliced_hydrogen_chain(2; common..., spacing = 0.19)
    @test_throws ArgumentError sliced_hydrogen_chain(2; common..., spacing = 0.2,
        phase = 0.1)
    @test_throws ArgumentError sliced_hydrogen_chain(2; common..., spacing = 0.2,
        boundary = :periodic)
    @test_throws ArgumentError sliced_hydrogen_chain(2; common..., spacing = 0.2,
        basis_name = "STO-3G")
    @test_throws ArgumentError sliced_hydrogen_chain(2; common..., spacing = 0.2,
        one_body_tolerance = 0.0)
    identity2 = Matrix{Float64}(I, 2, 2)
    @test_throws ArgumentError GaussletBases._sliced_transverse_vector(0.0, Float64[],
        [1.0, 2.0], ones(2), ones(2), identity2, identity2, identity2, 1e-8)
    @test_throws ArgumentError GaussletBases._sliced_longitudinal(0.2,
        [-0.2, 0.0, 0.2], [-1, 0, 1], [1, 2, 3], 1.0)
    @test_throws ArgumentError sliced_hydrogen_chain(2; common..., spacing = 0.2,
        interaction_tolerance = 1e-30)
end
