@testset "Vendored angular sphere-point access" begin
    orders = sphere_point_set_orders()
    curated_orders = curated_sphere_point_set_orders()

    @test first(orders) == 10
    @test last(orders) == 580
    @test length(orders) == 122
    @test issorted(orders)
    @test 15 in orders
    @test 32 in orders
    @test 51 in orders
    @test 100 in orders
    @test 580 in orders
    @test_throws ArgumentError sphere_point_set(9)
    @test_throws ArgumentError sphere_point_set(101)

    set10 = sphere_point_set(10)
    set15_full = sphere_point_set(15)
    set100 = sphere_point_set(100)
    set580 = sphere_point_set(580)

    for set in (set10, set15_full, set100, set580)
        @test set isa SpherePointSet
        @test set.cardinality == set.order
        @test size(set.coordinates) == (set.cardinality, 3)
        @test all(isfinite, set.coordinates)
        @test set.nn_ratio >= 1.0
        @test set.provenance.source_tag == "optimized_sphere_points_full_vendor"
        @test set.provenance.source_project == "GaussletModules/Radial"
        @test occursin("SpherePoints.jld2", set.provenance.source_artifact)
        @test occursin("xyzsets", set.provenance.source_note)

        norms = sqrt.(sum(abs2, set.coordinates; dims = 2))
        @test maximum(abs.(norms .- 1.0)) < 1.0e-12
    end

    @test all(order in orders for order in curated_orders)
    @test set15_full.coordinates ≈ curated_sphere_point_set(15).coordinates atol = 0.0 rtol = 1.0e-14
    @test set15_full.nn_ratio ≈ curated_sphere_point_set(15).nn_ratio atol = 0.0 rtol = 1.0e-14
end

@testset "Curated angular sphere-point access" begin
    orders = curated_sphere_point_set_orders()

    @test orders == [15, 32, 51]
    @test_throws ArgumentError curated_sphere_point_set(14)

    set15 = curated_sphere_point_set(15)
    set32 = curated_sphere_point_set(32)
    set51 = curated_sphere_point_set(51)

    for set in (set15, set32, set51)
        @test set isa CuratedSpherePointSet
        @test set.cardinality == set.order
        @test size(set.coordinates) == (set.cardinality, 3)
        @test all(isfinite, set.coordinates)
        @test set.nn_ratio >= 1.0
        @test set.provenance.source_tag == "optimized_sphere_points_curated_subset"
        @test set.provenance.source_project == "GaussletModules/Radial"
        @test occursin("SpherePoints.jld2", set.provenance.source_artifact)
        @test occursin("xyzsets", set.provenance.source_note)

        norms = sqrt.(sum(abs2, set.coordinates; dims = 2))
        @test maximum(abs.(norms .- 1.0)) < 1.0e-12
    end

    @test set15.nn_ratio ≈ 1.1233641689852316 atol = 0.0 rtol = 1.0e-14
    @test set32.nn_ratio ≈ 1.0000000159730313 atol = 0.0 rtol = 1.0e-14
    @test set51.nn_ratio ≈ 1.0802896662822246 atol = 0.0 rtol = 1.0e-14
    @test set15.coordinates[1, :] ≈ [0.5449408412377406, -0.022962554521322145, 0.838160008972606] atol = 0.0 rtol = 1.0e-14
    @test set32.coordinates[1, :] ≈ [0.40002938831494556, -0.0010525966270966628, 0.9165017078678638] atol = 0.0 rtol = 1.0e-14
end

@testset "Explicit angular Fibonacci and optimization paths" begin
    fib10_a = fibonacci_sphere_point_set(10)
    fib10_b = fibonacci_sphere_point_set(10)
    fib15 = fibonacci_sphere_point_set(15)
    optimized10 = optimize_sphere_point_set(fib10_a; beta = 2.0, iters = 20, gtol = 1.0e-8)

    @test fib10_a isa SpherePointSet
    @test fib15 isa SpherePointSet
    @test fib10_a.order == fib10_a.cardinality == 10
    @test fib15.order == fib15.cardinality == 15
    @test fib10_a.coordinates ≈ fib10_b.coordinates atol = 0.0 rtol = 0.0
    @test fib10_a.nn_ratio == fib10_b.nn_ratio
    @test fib10_a.provenance.source_tag == "deterministic_fibonacci_seed"
    @test fib10_a.provenance.source_project == "GaussletBases"
    @test fib10_a.provenance.source_artifact == "in_memory_generation"
    @test occursin("Fibonacci sphere seed", fib10_a.provenance.source_note)
    @test occursin("no randomization", fib10_a.provenance.source_note)

    fib10_kappa = GaussletBases._sphere_point_kappa_from_beta(fib10_a.order, 2.0)
    fib10_initial_logdet =
        GaussletBases._sphere_point_logdet(
            GaussletBases._sphere_point_gaussian_gram(fib10_a.coordinates, fib10_kappa),
        )
    fib10_optimized_logdet =
        GaussletBases._sphere_point_logdet(
            GaussletBases._sphere_point_gaussian_gram(optimized10.coordinates, fib10_kappa),
        )

    @test optimized10 isa SpherePointSet
    @test optimized10.order == fib10_a.order
    @test optimized10.cardinality == fib10_a.cardinality
    @test optimized10.provenance.source_tag == "optimized_from_input_point_set"
    @test optimized10.provenance.source_project == "GaussletBases"
    @test optimized10.provenance.source_artifact == "in_memory_optimization"
    @test occursin("deterministic_fibonacci_seed", optimized10.provenance.source_note)
    @test occursin("beta=2.0", optimized10.provenance.source_note)
    @test occursin("iters=20", optimized10.provenance.source_note)
    @test occursin("gtol=1.0e-8", optimized10.provenance.source_note)
    @test maximum(abs.(sqrt.(sum(abs2, optimized10.coordinates; dims = 2)) .- 1.0)) < 1.0e-12
    @test optimized10.nn_ratio >= 1.0
    @test fib10_optimized_logdet ≥ fib10_initial_logdet - 1.0e-10
end

@testset "Shell-local injected angular basis" begin
    shell15 = _shell_local_injected_angular_fixture(15)
    shell32 = _shell_local_injected_angular_fixture(32)
    shell51 = _shell_local_injected_angular_fixture(51)

    @test shell15 isa ShellLocalInjectedAngularBasis
    @test shell32 isa ShellLocalInjectedAngularBasis
    @test shell51 isa ShellLocalInjectedAngularBasis

    @test shell15.l_inject == 1
    @test shell32.l_inject == 3
    @test shell51.l_inject == 4

    for shell in (shell15, shell32, shell51)
        diagnostics = shell_local_injected_angular_diagnostics(shell)
        expected_kinetic = Diagonal(diagnostics.expected_injected_kinetic_eigenvalues)

        @test shell.prototype_count == shell.point_set.order
        @test shell.final_count == shell.prototype_count
        @test shell.injected_count == (shell.l_inject + 1)^2
        @test shell.whitened_complement_count ≥ shell.final_count - shell.injected_count
        @test size(shell.prototype_overlap) == (shell.prototype_count, shell.prototype_count)
        @test size(shell.final_overlap) == (shell.final_count, shell.final_count)
        @test size(shell.final_kinetic) == (shell.final_count, shell.final_count)
        @test size(shell.injected_overlap) == (shell.injected_count, shell.final_count)
        @test sum(shell.shell_weights) ≈ 4 * pi atol = 1.0e-12 rtol = 1.0e-12
        @test minimum(shell.theta_nn) > 0.0
        @test minimum(shell.kappa) > 0.0
        @test diagnostics.overlap_error ≤ 1.0e-10
        @test diagnostics.injected_exactness_error ≤ 1.0e-10
        @test diagnostics.injected_kinetic_error ≤ 1.0e-9
        @test norm(shell.injected_kinetic - expected_kinetic, Inf) ≤ 1.0e-9
        @test issorted(diagnostics.expected_injected_kinetic_eigenvalues)
        @test sort(diagnostics.injected_kinetic_eigenvalues) ≈
              diagnostics.expected_injected_kinetic_eigenvalues atol = 1.0e-9 rtol = 1.0e-9
    end
end

@testset "Shell-local angular profiles" begin
    profile15 = _shell_local_angular_profile_fixture(15)
    profile15_again = shell_local_angular_profile(15)
    profile15_uncached =
        GaussletBases._build_shell_local_angular_profile_uncached(sphere_point_set(15))
    profile32 = _shell_local_angular_profile_fixture(32)
    overlap15_32 = _shell_local_angular_profile_overlap_fixture(15, 32)
    overlap15_32_again = adjacent_shell_local_angular_profile_overlap(profile15, profile32)

    @test profile15 isa ShellLocalAngularProfile
    @test profile15.key.order == 15
    @test profile15.key.point_set_source_tag == sphere_point_set(15).provenance.source_tag
    @test profile15.key.gauge_version == :v1_seed_order_dominant_positive
    @test profile15 === profile15_again
    @test build_shell_local_injected_angular_basis(15) === profile15.basis
    @test profile15.profile_id == profile15_uncached.profile_id
    @test profile15.labels == profile15_uncached.labels
    @test profile15.block_kinds == profile15_uncached.block_kinds
    @test profile15.diagnostics.grand_coefficients_checksum ==
          profile15_uncached.diagnostics.grand_coefficients_checksum

    expected_exact_labels = [
        "exact_l$(channel.l)_m$(channel.m)" for channel in profile15.basis.injected_channels.channel_data
    ]
    @test profile15.exact_labels == expected_exact_labels
    @test profile15.labels[1:profile15.basis.injected_count] == expected_exact_labels
    @test profile15.mixed_labels == ["mixed_$(i)" for i in 1:(profile15.basis.final_count - profile15.basis.injected_count)]
    @test profile15.block_kinds ==
          vcat(fill(:exact, profile15.basis.injected_count), fill(:mixed, profile15.basis.final_count - profile15.basis.injected_count))
    @test all(diag(profile15.basis.injected_overlap[:, 1:profile15.basis.injected_count]) .> 0.0)

    mixed_offset = profile15.basis.injected_count
    mixed_metadata = profile15.gauge_metadata
    @test mixed_metadata.mixed_orientation_strategy == :seed_order
    @test length(mixed_metadata.mixed_dominant_grand_indices) == profile15.basis.final_count - profile15.basis.injected_count
    @test length(mixed_metadata.mixed_signs) == profile15.basis.final_count - profile15.basis.injected_count
    for j in 1:(profile15.basis.final_count - profile15.basis.injected_count)
        idx = mixed_metadata.mixed_dominant_grand_indices[j]
        @test profile15.basis.grand_coefficients[idx, mixed_offset + j] > 0.0
    end

    @test overlap15_32 isa ShellLocalAngularProfileOverlap
    @test overlap15_32 === overlap15_32_again
    @test size(overlap15_32.overlap) == (profile15.basis.final_count, profile32.basis.final_count)
    @test overlap15_32.source_labels == profile15.labels
    @test overlap15_32.target_labels == profile32.labels
    @test overlap15_32.source_exact_count == profile15.basis.injected_count
    @test overlap15_32.target_exact_count == profile32.basis.injected_count
    @test overlap15_32.shell_independent
    @test overlap15_32.diagnostics.min_singular_value > 1.0e-8
    @test isfinite(overlap15_32.diagnostics.exact_block_inf_norm)
end

@testset "Atomic shell-local angular assembly" begin
    assembly = _atomic_shell_local_angular_fixture()
    diagnostics = atomic_shell_local_angular_diagnostics(assembly)
    scheduled_orders = assign_atomic_angular_shell_orders(
        [0.05, 0.5, 1.2, 3.0, 8.0];
        ord_min = 15,
        ord_max = 51,
        r_lo = 0.2,
        r_hi = 4.5,
        w_lo = 0.2,
        w_hi = 0.7,
    )
    scheduled_orders_full = assign_atomic_angular_shell_orders(
        [0.05, 0.5, 1.2, 3.0, 8.0];
        ord_min = 24,
        ord_max = 58,
        r_lo = 0.2,
        r_hi = 4.5,
        w_lo = 0.2,
        w_hi = 0.7,
    )

    @test assembly isa AtomicShellLocalInjectedAngularAssembly
    @test assembly.shell_orders == [15, 32, 51, 32]
    @test assembly.shell_dimensions == [15, 32, 51, 32]
    @test assembly.shell_offsets == [1, 16, 48, 99]
    @test assembly.shell_exact_lmax == [1, 3, 4, 3]
    @test assembly.profiles[1] === shell_local_angular_profile(15)
    @test assembly.profiles[2] === shell_local_angular_profile(32)
    @test assembly.profiles[2] === assembly.profiles[4]
    @test assembly.shells[2] === assembly.profiles[2].basis
    @test assembly.shells[2] === assembly.shells[4]
    @test size(assembly.overlap) == (130, 130)
    @test size(assembly.kinetic) == (130, 130)
    @test diagnostics.nshells == 4
    @test diagnostics.total_dim == 130
    @test diagnostics.max_shell_overlap_error ≤ 1.0e-10
    @test diagnostics.max_shell_injected_exactness_error ≤ 1.0e-10
    @test diagnostics.max_shell_injected_kinetic_error ≤ 1.0e-9
    @test diagnostics.max_diagonal_overlap_error ≤ 1.0e-10
    @test diagnostics.max_diagonal_injected_kinetic_error ≤ 1.0e-9
    @test diagnostics.max_pair_overlap_symmetry_error ≤ 1.0e-10
    @test diagnostics.max_pair_kinetic_symmetry_error ≤ 1.0e-10
    @test all(
        diagnostics.shell_interaction_lcap[i] ≥
        diagnostics.shell_interaction_lexpand[i] ≥
        diagnostics.shell_exact_lmax[i] for i in eachindex(diagnostics.shell_exact_lmax)
    )
    @test diagnostics.interaction_pair_plan ≥ diagnostics.interaction_exact_lower_bound
    @test diagnostics.max_shell_interaction_tail ≤ 1.0e-6

    @test scheduled_orders[1] == 15
    @test scheduled_orders[end] == 15
    @test all(order in sphere_point_set_orders() for order in scheduled_orders)
    @test any(order ∉ curated_sphere_point_set_orders() for order in scheduled_orders)
    @test maximum(scheduled_orders) ≥ 32
    @test all(order in sphere_point_set_orders() for order in scheduled_orders_full)
    @test any(order ∉ curated_sphere_point_set_orders() for order in scheduled_orders_full)
end

@testset "Atomic fixed-radial angular sequence" begin
    sequence = _paper_style_fixed_radial_angular_sequence_fixture()
    level10 = sequence.levels[1]
    level15 = sequence.levels[2]
    level32 = sequence.levels[3]
    sidecar10_15 = sequence.adjacent_overlaps[1]
    sidecar15_32 = sequence.adjacent_overlaps[2]
    sidecar10_32 = sequence.direct_overlaps[1]
    nr = length(sequence.shell_ids)

    @test sequence isa AtomicFixedRadialAngularSequence
    @test sequence.N_sph_values == [10, 15, 32]
    @test length(sequence.levels) == 3
    @test length(sequence.adjacent_overlaps) == 2
    @test length(sequence.direct_overlaps) == 1
    @test sequence.shell_ids == collect(1:nr)
    @test issorted(sequence.shell_centers_r)
    @test all(level -> level.radial_basis_id == sequence.radial_basis_id, sequence.levels)
    @test all(level -> level.shell_ids == sequence.shell_ids, sequence.levels)
    @test all(level -> level.shell_centers_r == sequence.shell_centers_r, sequence.levels)
    @test level10.profile === shell_local_angular_profile(10)
    @test level15.profile === shell_local_angular_profile(15)
    @test level32.profile === shell_local_angular_profile(32)
    @test all(profile -> profile === level10.profile, level10.payload.one_body.angular_assembly.profiles)
    @test all(profile -> profile === level15.profile, level15.payload.one_body.angular_assembly.profiles)
    @test all(profile -> profile === level32.profile, level32.payload.one_body.angular_assembly.profiles)
    @test level10.shell_dimensions == fill(level10.profile.basis.final_count, nr)
    @test level15.shell_dimensions == fill(level15.profile.basis.final_count, nr)
    @test level32.shell_dimensions == fill(level32.profile.basis.final_count, nr)

    level_payload = atomic_fixed_radial_angular_level_dense_payload(level10)
    @test level_payload.payload["H1"] ≈ level10.payload.hamiltonian atol = 0.0 rtol = 0.0
    @test level_payload.payload["Vee"] ≈ level10.payload.interaction atol = 0.0 rtol = 0.0
    @test level_payload.payload["shell_ids"] == sequence.shell_ids
    @test level_payload.payload["shell_centers_r"] == sequence.shell_centers_r
    @test level_payload.payload["shell_dimensions"] == level10.shell_dimensions
    @test level_payload.payload["within_shell_labels"] == level10.profile.labels
    @test level_payload.bridge_meta["sequence_id"] == sequence.sequence_id
    @test level_payload.bridge_meta["N_sph"] == 10
    @test level_payload.bridge_meta["angular_profile_id"] == level10.profile.profile_id
    @test level_payload.bridge_meta["gauge_version"] == string(level10.profile.key.gauge_version)

    overlap_payload = atomic_fixed_radial_angular_overlap_sidecar_payload(sidecar10_15)
    @test overlap_payload.payload["overlap"] ≈ sidecar10_15.overlap atol = 0.0 rtol = 0.0
    @test overlap_payload.payload["source_labels"] == level10.profile.labels
    @test overlap_payload.payload["target_labels"] == level15.profile.labels
    @test overlap_payload.bridge_meta["sequence_id"] == sequence.sequence_id
    @test overlap_payload.bridge_meta["source_N_sph"] == 10
    @test overlap_payload.bridge_meta["target_N_sph"] == 15
    @test overlap_payload.bridge_meta["pair_kind"] == "adjacent"
    @test overlap_payload.bridge_meta["shell_independent"]
    @test sidecar10_15.source_profile_id == level10.profile.profile_id
    @test sidecar10_15.target_profile_id == level15.profile.profile_id
    @test sidecar10_15.source_gauge_version == level10.profile.key.gauge_version
    @test sidecar10_15.target_gauge_version == level15.profile.key.gauge_version
    @test sidecar10_15.source_labels == level10.profile.labels
    @test sidecar10_15.target_labels == level15.profile.labels
    @test sidecar10_15.shell_independent
    @test sidecar10_15.pair_kind == :adjacent
    @test sidecar15_32.source_profile_id == level15.profile.profile_id
    @test sidecar15_32.target_profile_id == level32.profile.profile_id
    @test sidecar15_32.pair_kind == :adjacent
    @test sidecar10_32.source_N_sph == 10
    @test sidecar10_32.target_N_sph == 32
    @test sidecar10_32.source_profile_id == level10.profile.profile_id
    @test sidecar10_32.target_profile_id == level32.profile.profile_id
    @test sidecar10_32.source_gauge_version == level10.profile.key.gauge_version
    @test sidecar10_32.target_gauge_version == level32.profile.key.gauge_version
    @test sidecar10_32.source_labels == level10.profile.labels
    @test sidecar10_32.target_labels == level32.profile.labels
    @test sidecar10_32.shell_independent
    @test sidecar10_32.pair_kind == :direct

    legacy_payload = atomic_fixed_radial_legacy_dmrgatom_payload(level10)
    @test legacy_payload.payload["H1"] ≈ level10.payload.hamiltonian atol = 0.0 rtol = 0.0
    @test legacy_payload.payload["Vee"] ≈ level10.payload.interaction atol = 0.0 rtol = 0.0
    @test legacy_payload.payload["dims_per_shell"] == level10.shell_dimensions
    @test legacy_payload.meta_values["Z"] == 2
    @test legacy_payload.bridge_meta["angular_profile_id"] == level10.profile.profile_id
    @test legacy_payload.bridge_meta["gauge_version"] == string(level10.profile.key.gauge_version)
    @test legacy_payload.bridge_meta["basis_centers_are_representative"]
    @test legacy_payload.bridge_meta["basis_centers_center_policy_name"] ==
          "representative_dominant_prototype_direction_v1"
    dominant_indices =
        Int.(legacy_payload.bridge_meta["within_shell_dominant_prototype_indices"])
    @test length(dominant_indices) == level10.profile.basis.final_count
    @test all(index -> 1 <= index <= level10.profile.basis.prototype_count, dominant_indices)
    basis_centers = Matrix{Float64}(legacy_payload.payload["basis_centers"])
    @test size(basis_centers) ==
          (size(level10.payload.hamiltonian, 1), 3)
    point_coordinates = level10.profile.basis.point_set.coordinates
    shell_size = level10.profile.basis.final_count
    first_shell_centers = basis_centers[1:shell_size, :]
    expected_first_shell_centers =
        level10.shell_centers_r[1] .* point_coordinates[dominant_indices, :]
    @test first_shell_centers ≈ expected_first_shell_centers atol = 1.0e-12 rtol = 1.0e-12

    mktempdir() do dir
        level_path = joinpath(dir, "he_fixed_radial_level10.jld2")
        sidecar_path = joinpath(dir, "he_fixed_radial_10_15_overlap.jld2")
        direct_sidecar_path = joinpath(dir, "he_fixed_radial_10_32_overlap.jld2")
        legacy_path = joinpath(dir, "he_fixed_radial_level10.legacy_dmrgatom.jld2")
        @test write_atomic_fixed_radial_angular_level_jld2(level_path, level10) == level_path
        @test write_atomic_fixed_radial_angular_overlap_sidecar_jld2(sidecar_path, sidecar10_15) == sidecar_path
        @test write_atomic_fixed_radial_angular_overlap_sidecar_jld2(direct_sidecar_path, sidecar10_32) == direct_sidecar_path
        @test write_atomic_fixed_radial_legacy_dmrgatom_jld2(legacy_path, level10) == legacy_path
        jldopen(level_path, "r") do file
            @test String(file["bridge/format"]) == "angular_fixed_radial_dense_v1"
            @test Int(file["bridge/level_index"]) == 1
            @test Int(file["bridge/N_sph"]) == 10
            @test String(file["bridge/sequence_id"]) == sequence.sequence_id
            @test String(file["bridge/angular_profile_id"]) == level10.profile.profile_id
            @test Int.(file["shell_ids"]) == sequence.shell_ids
            @test Float64.(file["shell_centers_r"]) == sequence.shell_centers_r
            @test String.(file["within_shell_labels"]) == level10.profile.labels
        end
        jldopen(sidecar_path, "r") do file
            @test String(file["bridge/format"]) == "angular_fixed_radial_profile_overlap_v1"
            @test Int(file["bridge/source_N_sph"]) == 10
            @test Int(file["bridge/target_N_sph"]) == 15
            @test String(file["bridge/source_level_id"]) == level10.level_id
            @test String(file["bridge/target_level_id"]) == level15.level_id
            @test String(file["bridge/source_gauge_version"]) == string(level10.profile.key.gauge_version)
            @test String(file["bridge/target_gauge_version"]) == string(level15.profile.key.gauge_version)
            @test String(file["bridge/pair_kind"]) == "adjacent"
            @test Bool(file["bridge/shell_independent"])
            @test String.(file["source_labels"]) == level10.profile.labels
            @test String.(file["target_labels"]) == level15.profile.labels
            @test Matrix{Float64}(file["overlap"]) ≈ sidecar10_15.overlap atol = 0.0 rtol = 0.0
        end
        jldopen(direct_sidecar_path, "r") do file
            @test String(file["bridge/format"]) == "angular_fixed_radial_profile_overlap_v1"
            @test Int(file["bridge/source_N_sph"]) == 10
            @test Int(file["bridge/target_N_sph"]) == 32
            @test String(file["bridge/source_level_id"]) == level10.level_id
            @test String(file["bridge/target_level_id"]) == level32.level_id
            @test String(file["bridge/source_profile_id"]) == level10.profile.profile_id
            @test String(file["bridge/target_profile_id"]) == level32.profile.profile_id
            @test String(file["bridge/source_gauge_version"]) == string(level10.profile.key.gauge_version)
            @test String(file["bridge/target_gauge_version"]) == string(level32.profile.key.gauge_version)
            @test String(file["bridge/pair_kind"]) == "direct"
            @test Bool(file["bridge/shell_independent"])
            @test String.(file["source_labels"]) == level10.profile.labels
            @test String.(file["target_labels"]) == level32.profile.labels
            @test Matrix{Float64}(file["overlap"]) ≈ sidecar10_32.overlap atol = 0.0 rtol = 0.0
        end
        jldopen(legacy_path, "r") do file
            @test String(file["bridge/format"]) == "legacy_dmrgatom_dense_v1"
            @test Int(file["bridge/N_sph"]) == 10
            @test Int(file["meta/Z"]) == 2
            @test Int.(file["dims_per_shell"]) == level10.shell_dimensions
            @test Matrix{Float64}(file["H1"]) ≈ level10.payload.hamiltonian atol = 0.0 rtol = 0.0
            @test Matrix{Float64}(file["Vee"]) ≈ level10.payload.interaction atol = 0.0 rtol = 0.0
            @test String(file["bridge/basis_centers_center_policy_name"]) ==
                  "representative_dominant_prototype_direction_v1"
            @test Bool(file["bridge/basis_centers_are_representative"])
            @test Int.(file["bridge/within_shell_dominant_prototype_indices"]) ==
                  dominant_indices
            @test Matrix{Float64}(file["basis_centers"]) ≈ basis_centers atol = 0.0 rtol = 0.0
        end
    end
end

@testset "Atomic injected angular one-body benchmark" begin
    benchmark = _atomic_injected_angular_one_body_benchmark_fixture()
    diagnostics = atomic_injected_angular_one_body_diagnostics(benchmark)
    _, _, radial_ops, _, _, _ = _quick_radial_atomic_fixture()
    exact_atom = atomic_one_body_operators(radial_ops; lmax = diagnostics.exact_common_lmax)
    exact_benchmark_spectrum =
        GaussletBases._generalized_spectrum(benchmark.exact_hamiltonian, benchmark.exact_overlap)
    exact_atom_spectrum =
        GaussletBases._generalized_spectrum(exact_atom.hamiltonian, exact_atom.overlap)

    @test benchmark isa AtomicInjectedAngularOneBodyBenchmark
    @test benchmark.angular_assembly isa AtomicShellLocalInjectedAngularAssembly
    @test benchmark.angular_assembly.shell_radii == radial_ops.shell_centers_r
    @test benchmark.exact_common_lmax ≥ 1
    @test length(benchmark.exact_channels) == (benchmark.exact_common_lmax + 1)^2
    @test size(benchmark.overlap) == size(benchmark.hamiltonian)
    @test size(benchmark.overlap, 1) == sum(benchmark.angular_assembly.shell_dimensions)
    @test size(benchmark.exact_overlap) == size(exact_atom.overlap)
    @test size(benchmark.exact_hamiltonian) == size(exact_atom.hamiltonian)
    @test exact_benchmark_spectrum ≈ exact_atom_spectrum atol = 1.0e-10 rtol = 1.0e-10

    @test diagnostics.overlap_symmetry_error ≤ 1.0e-10
    @test diagnostics.hamiltonian_symmetry_error ≤ 1.0e-10
    @test diagnostics.min_overlap_eigenvalue > 1.0e-6
    @test diagnostics.projected_exact_overlap_error ≤ 1.0e-9
    @test diagnostics.projected_exact_hamiltonian_error ≤ 1.0e-8
    @test diagnostics.projected_exact_low_eigenvalue_count == 4
    @test diagnostics.projected_exact_low_eigenvalue_error ≤ 1.0e-8
    @test diagnostics.benchmark_ground_state_error ≤ 1.0e-8
    @test diagnostics.benchmark_ground_state_energy ≈ diagnostics.exact_ground_state_energy atol = 1.0e-8 rtol = 1.0e-8
end

@testset "Atomic injected angular Cartesian moments" begin
    rb, grid, radial_ops, benchmark, bundle, bundle_reconstructed =
        _paper_style_angular_cartesian_moments_fixture(15; Z = 2.0, lmax = 2)

    @test bundle isa AtomicInjectedAngularCartesianMomentBundle
    @test bundle.radial_operators === radial_ops
    @test bundle.angular_assembly === benchmark.angular_assembly
    @test bundle.shell_ranges == benchmark.shell_ranges
    @test bundle.S ≈ benchmark.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test bundle.radial_moment_source == :explicit_basis_and_grid
    @test bundle_reconstructed.radial_moment_source == :reconstructed_from_source_manifest
    @test bundle_reconstructed.S ≈ bundle.S atol = 0.0 rtol = 0.0
    @test bundle_reconstructed.X ≈ bundle.X atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.Y ≈ bundle.Y atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.Z ≈ bundle.Z atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.X2 ≈ bundle.X2 atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.Y2 ≈ bundle.Y2 atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.Z2 ≈ bundle.Z2 atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.R2 ≈ bundle.R2 atol = 1.0e-10 rtol = 1.0e-10

    for matrix in (bundle.S, bundle.X, bundle.Y, bundle.Z, bundle.X2, bundle.Y2, bundle.Z2, bundle.R2)
        @test all(isfinite, matrix)
        @test opnorm(matrix - transpose(matrix), Inf) ≤ 1.0e-8
    end

    @test opnorm(bundle.X, Inf) > 1.0e-8
    @test opnorm(bundle.Y, Inf) > 1.0e-8
    @test opnorm(bundle.Z, Inf) > 1.0e-8
    @test opnorm(bundle.X2, Inf) > 1.0e-8
    @test opnorm(bundle.Y2, Inf) > 1.0e-8
    @test opnorm(bundle.Z2, Inf) > 1.0e-8
    @test bundle.R2 ≈ bundle.X2 + bundle.Y2 + bundle.Z2 atol = 1.0e-8 rtol = 1.0e-8

    s_exact = build_atomic_injected_angular_one_body_benchmark(
        radial_ops;
        shell_orders = fill(10, length(radial_ops.shell_centers_r)),
        l_inject = 0,
    )
    s_bundle = build_atomic_injected_angular_cartesian_moments(
        s_exact;
        radial_basis = rb,
        radial_grid = grid,
    )
    exact_transform = s_exact.exact_transform
    projected_X = exact_transform * s_bundle.X * transpose(exact_transform)
    projected_Y = exact_transform * s_bundle.Y * transpose(exact_transform)
    projected_Z = exact_transform * s_bundle.Z * transpose(exact_transform)
    projected_X2 = exact_transform * s_bundle.X2 * transpose(exact_transform)
    projected_Y2 = exact_transform * s_bundle.Y2 * transpose(exact_transform)
    projected_Z2 = exact_transform * s_bundle.Z2 * transpose(exact_transform)
    projected_R2 = exact_transform * s_bundle.R2 * transpose(exact_transform)

    @test s_exact.exact_common_lmax == 0
    @test projected_X ≈ zeros(size(projected_X)) atol = 1.0e-10 rtol = 0.0
    @test projected_Y ≈ zeros(size(projected_Y)) atol = 1.0e-10 rtol = 0.0
    @test projected_Z ≈ zeros(size(projected_Z)) atol = 1.0e-10 rtol = 0.0
    @test projected_X2 ≈ projected_Y2 atol = 1.0e-9 rtol = 1.0e-9
    @test projected_Y2 ≈ projected_Z2 atol = 1.0e-9 rtol = 1.0e-9
    @test projected_R2 ≈ projected_X2 + projected_Y2 + projected_Z2 atol = 1.0e-9 rtol = 1.0e-9
end

@testset "Atomic injected angular HF-style benchmark" begin
    benchmark = _atomic_injected_angular_hf_style_benchmark_fixture()
    diagnostics = atomic_injected_angular_hf_style_diagnostics(benchmark)
    exact_ida = atomic_ida_operators(benchmark.one_body.radial_operators; lmax = benchmark.one_body.exact_common_lmax)
    exact_shell_major_interaction = atomic_ida_density_interaction_matrix(exact_ida; ordering = :shell_major)

    @test benchmark isa AtomicInjectedAngularHFStyleBenchmark
    @test benchmark.one_body isa AtomicInjectedAngularOneBodyBenchmark
    @test benchmark.exact_ida_reference isa AtomicIDAOperators
    @test size(benchmark.interaction) == size(benchmark.one_body.overlap)
    @test size(benchmark.exact_interaction) == size(benchmark.one_body.exact_overlap)
    @test benchmark.exact_interaction ≈ exact_shell_major_interaction atol = 1.0e-12 rtol = 1.0e-12

    @test diagnostics.interaction_symmetry_error ≤ 1.0e-10
    @test diagnostics.full_converged
    @test diagnostics.exact_converged
    @test diagnostics.full_residual ≤ 1.0e-8
    @test diagnostics.exact_residual ≤ 1.0e-8 || diagnostics.exact_iterations == 13
    @test diagnostics.full_electron_count_error ≤ 1.0e-10
    @test diagnostics.exact_electron_count_error ≤ 1.0e-10
    @test isfinite(diagnostics.full_energy)
    @test isfinite(diagnostics.exact_energy)
    @test abs(diagnostics.energy_difference_to_exact_reference) ≤ 1.0e-8
    @test diagnostics.ground_orbital_energy_error ≤ 1.0e-8

    hfdmrg = _local_hfdmrg_module()
    if hfdmrg !== nothing
        payload = build_atomic_injected_angular_hfdmrg_payload(benchmark)
        @test payload isa AtomicInjectedAngularHFDMRGHFAdapter
        @test payload.route == :dense_density_density
        @test payload.solver_mode == :restricted_closed_shell
        @test payload.nup == 1
        @test payload.ndn == 1
        hfdmrg_result = run_atomic_injected_angular_hfdmrg_hf(
            payload;
            hfmod = hfdmrg,
            nblockcenter = 1,
            maxiter = 40,
            cutoff = 1.0e-10,
            scf_cutoff = 1.0e-11,
            verbose = false,
        )
        @test hfdmrg_result.route == payload.route
        @test hfdmrg_result.solver_mode == :restricted_closed_shell
        @test hfdmrg_result.nblockcenter == 1
        @test hfdmrg_result.blocksize == min(size(payload.hamiltonian, 1), 64)
        @test abs(diagnostics.full_energy - hfdmrg_result.energy) ≤ 1.0e-8
    end
end

@testset "Atomic injected angular HFDMRG-facing HF adapter" begin
    adapter = _atomic_injected_angular_hfdmrg_hf_adapter_fixture()
    diagnostics = atomic_injected_angular_hfdmrg_hf_adapter_diagnostics(adapter)
    benchmark = _atomic_injected_angular_hf_style_benchmark_fixture()
    from_benchmark = build_atomic_injected_angular_hfdmrg_hf_adapter(benchmark)
    payload_from_benchmark = build_atomic_injected_angular_hfdmrg_payload(benchmark)
    from_small_ed =
        build_atomic_injected_angular_hfdmrg_hf_adapter(_atomic_injected_angular_small_ed_benchmark_fixture())
    benchmark_diagnostics = atomic_injected_angular_hfdmrg_hf_adapter_diagnostics(from_benchmark)
    open_shell_seeds =
        build_atomic_injected_angular_hfdmrg_hf_seeds(benchmark; nup = 2, ndn = 1)
    explicit_psiup0 = open_shell_seeds.psiup0[:, [2, 1]]
    explicit_psidn0 = open_shell_seeds.psidn0
    open_shell_adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = 2,
        ndn = 1,
    )
    explicit_open_shell_adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = 2,
        ndn = 1,
        psiup0 = explicit_psiup0,
        psidn0 = explicit_psidn0,
    )
    mixed_seed_adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = 2,
        ndn = 1,
        psiup0 = explicit_psiup0,
    )

    @test adapter isa AtomicInjectedAngularHFDMRGHFAdapter
    @test adapter.one_body isa AtomicInjectedAngularOneBodyBenchmark
    @test isnothing(adapter.hf_style)
    @test adapter.route == :dense_density_density
    @test size(adapter.hamiltonian) == size(adapter.interaction)
    @test size(adapter.hamiltonian, 1) == diagnostics.basis_dim
    @test size(adapter.psiup0) == (diagnostics.basis_dim, diagnostics.nup)
    @test size(adapter.psidn0) == (diagnostics.basis_dim, diagnostics.ndn)
    @test diagnostics.nup == 1
    @test diagnostics.ndn == 1
    @test diagnostics.solver_mode == :restricted_closed_shell
    @test diagnostics.overlap_identity_error ≤ 2.0e-6
    @test diagnostics.hamiltonian_symmetry_error ≤ 1.0e-8
    @test diagnostics.interaction_symmetry_error ≤ 1.0e-10
    @test diagnostics.psiup0_orthogonality_error ≤ 1.0e-12
    @test diagnostics.psidn0_orthogonality_error ≤ 1.0e-12
    @test diagnostics.psiup0_source == :default_one_body_orbitals
    @test diagnostics.psidn0_source == :default_one_body_orbitals
    @test !diagnostics.has_benchmark_reference
    @test ismissing(diagnostics.benchmark_full_energy)
    @test ismissing(diagnostics.benchmark_exact_energy)
    @test from_benchmark.hf_style === benchmark
    @test payload_from_benchmark.hf_style === benchmark
    @test benchmark_diagnostics.has_benchmark_reference
    @test benchmark_diagnostics.solver_mode == :restricted_closed_shell
    @test isfinite(benchmark_diagnostics.benchmark_full_energy)
    @test isfinite(benchmark_diagnostics.benchmark_exact_energy)
    @test from_benchmark.route == adapter.route
    @test from_benchmark.nup == adapter.nup
    @test from_benchmark.ndn == adapter.ndn
    @test from_benchmark.hamiltonian ≈ adapter.hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test from_benchmark.interaction ≈ adapter.interaction atol = 1.0e-12 rtol = 1.0e-12
    @test from_benchmark.psiup0 ≈ adapter.psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test from_benchmark.psidn0 ≈ adapter.psidn0 atol = 1.0e-12 rtol = 1.0e-12
    @test payload_from_benchmark.route == from_benchmark.route
    @test payload_from_benchmark.solver_mode == from_benchmark.solver_mode
    @test payload_from_benchmark.hamiltonian ≈ from_benchmark.hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test payload_from_benchmark.interaction ≈ from_benchmark.interaction atol = 1.0e-12 rtol = 1.0e-12
    @test payload_from_benchmark.psiup0 ≈ from_benchmark.psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test payload_from_benchmark.psidn0 ≈ from_benchmark.psidn0 atol = 1.0e-12 rtol = 1.0e-12
    @test from_small_ed.route == adapter.route
    @test from_small_ed.nup == adapter.nup
    @test from_small_ed.ndn == adapter.ndn
    @test from_small_ed.hamiltonian ≈ adapter.hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test from_small_ed.interaction ≈ adapter.interaction atol = 1.0e-12 rtol = 1.0e-12
    @test from_small_ed.psiup0 ≈ adapter.psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test from_small_ed.psidn0 ≈ adapter.psidn0 atol = 1.0e-12 rtol = 1.0e-12
    @test size(open_shell_seeds.psiup0) == (diagnostics.basis_dim, 2)
    @test size(open_shell_seeds.psidn0) == (diagnostics.basis_dim, 1)
    @test open_shell_adapter.nup == 2
    @test open_shell_adapter.ndn == 1
    @test open_shell_adapter.solver_mode == :unrestricted
    @test open_shell_adapter.psiup0_source == :default_one_body_orbitals
    @test open_shell_adapter.psidn0_source == :default_one_body_orbitals
    @test opnorm(transpose(open_shell_adapter.psiup0) * open_shell_adapter.psiup0 - Matrix{Float64}(I, 2, 2), Inf) ≤ 1.0e-12
    @test opnorm(transpose(open_shell_adapter.psidn0) * open_shell_adapter.psidn0 - Matrix{Float64}(I, 1, 1), Inf) ≤ 1.0e-12
    @test explicit_open_shell_adapter.solver_mode == :unrestricted
    @test explicit_open_shell_adapter.psiup0_source == :explicit_seed
    @test explicit_open_shell_adapter.psidn0_source == :explicit_seed
    @test explicit_open_shell_adapter.psiup0 ≈ explicit_psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test explicit_open_shell_adapter.psidn0 ≈ explicit_psidn0 atol = 1.0e-12 rtol = 1.0e-12
    @test mixed_seed_adapter.solver_mode == :unrestricted
    @test mixed_seed_adapter.psiup0_source == :explicit_seed
    @test mixed_seed_adapter.psidn0_source == :default_one_body_orbitals
    @test mixed_seed_adapter.psiup0 ≈ explicit_psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test mixed_seed_adapter.psidn0 ≈ open_shell_seeds.psidn0 atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Angular He legacy-trim HF payload anchors" begin
    hfdmrg = _local_hfdmrg_module()
    hfdmrg === nothing && return

    he_anchor_reference = -2.861679990485
    anchor_settings = (
        nblockcenter = 2,
        blocksize = 100,
        maxiter = 100,
        cutoff = 1.0e-8,
        scf_cutoff = 1.0e-9,
        verbose = false,
    )

    payload10 = _paper_style_angular_hfdmrg_payload_fixture(10)
    payload15 = _paper_style_angular_hfdmrg_payload_fixture(15)
    result10 = _solve_hfdmrg_from_payload_direct(payload10, hfdmrg; anchor_settings...)
    result15 = _solve_hfdmrg_from_payload_direct(payload15, hfdmrg; anchor_settings...)

    @test payload10.solver_mode == :restricted_closed_shell
    @test payload15.solver_mode == :restricted_closed_shell
    @test payload10.nup == 1
    @test payload10.ndn == 1
    @test payload15.nup == 1
    @test payload15.ndn == 1
    @test first(payload10.one_body.angular_assembly.shell_orders) == 10
    @test first(payload15.one_body.angular_assembly.shell_orders) == 15
    @test all(==(10), payload10.one_body.angular_assembly.shell_orders)
    @test all(==(15), payload15.one_body.angular_assembly.shell_orders)
    @test result10.energy ≈ he_anchor_reference atol = 5.0e-6 rtol = 1.0e-6
    @test result15.energy ≈ he_anchor_reference atol = 5.0e-6 rtol = 1.0e-6
    @test abs(result15.energy - he_anchor_reference) ≤ abs(result10.energy - he_anchor_reference)
end

@testset "Angular Be legacy-trim one-body anchors" begin
    for order in (10, 15)
        benchmark = _paper_style_angular_one_body_benchmark_fixture(order)
        full_spectrum = sort(real(eigvals(Hermitian(benchmark.hamiltonian), Hermitian(benchmark.overlap))))
        exact_spectrum =
            sort(real(eigvals(Hermitian(benchmark.exact_hamiltonian), Hermitian(benchmark.exact_overlap))))
        low_count = 8
        closed_shell_noninteracting = 2.0 * (full_spectrum[1] + full_spectrum[2])
        one_s_like_count = count(<(-4.0), full_spectrum[1:10])

        @test benchmark.exact_common_lmax == 1
        @test benchmark.angular_assembly.shell_orders == fill(order, length(benchmark.angular_assembly.shell_orders))
        @test all(lcap ≥ benchmark.exact_common_lmax + 4 for lcap in benchmark.angular_assembly.shell_kinetic_lcap)
        @test one_s_like_count == 1
        @test full_spectrum[1] < -7.0
        @test full_spectrum[2] > -3.0
        @test closed_shell_noninteracting > -22.0
        @test closed_shell_noninteracting < -18.0
        @test full_spectrum[1:low_count] ≈ exact_spectrum[1:low_count] atol = 1.0e-8 rtol = 1.0e-8
    end
end

@testset "Angular Ne legacy-trim one-body discriminator" begin
    for order in (10, 15, 32)
        benchmark = _paper_style_angular_one_body_benchmark_fixture(order; Z = 10.0, lmax = 6)
        full_spectrum = sort(real(eigvals(Hermitian(benchmark.hamiltonian), Hermitian(benchmark.overlap))))
        exact_spectrum =
            sort(real(eigvals(Hermitian(benchmark.exact_hamiltonian), Hermitian(benchmark.exact_overlap))))
        low_count = 8
        nocc = 5
        closed_shell_noninteracting = 2.0 * sum(full_spectrum[1:nocc])
        exact_noninteracting = 2.0 * sum(exact_spectrum[1:nocc])
        one_s_like_count = count(<(-20.0), full_spectrum[1:nocc])

        @test benchmark.angular_assembly.shell_orders == fill(order, length(benchmark.angular_assembly.shell_orders))
        if order in (10, 15)
            @test benchmark.exact_common_lmax == 1
        else
            @test benchmark.exact_common_lmax ≥ 3
        end
        @test one_s_like_count == 1
        @test full_spectrum[1] < -40.0
        @test full_spectrum[2] > -20.0
        @test closed_shell_noninteracting > -201.0
        @test closed_shell_noninteracting < -199.0
        @test closed_shell_noninteracting ≈ exact_noninteracting atol = 1.0e-5 rtol = 1.0e-7
        @test full_spectrum[1:low_count] ≈ exact_spectrum[1:low_count] atol = 1.0e-6 rtol = 1.0e-8
    end
end

@testset "Angular Ne legacy-trim interaction moment-span discriminator" begin
    for order in (10, 15)
        benchmark = _paper_style_angular_one_body_benchmark_fixture(order; Z = 10.0, lmax = 6)
        assembly = benchmark.angular_assembly
        bare_moment_lmax =
            maximum(maximum(keys(blocks)) for blocks in assembly.shell_moment_blocks)
        interaction_moment_lmax =
            maximum(maximum(keys(blocks)) for blocks in assembly.shell_interaction_moment_blocks)
        required_product_lmax = 2 * benchmark.exact_common_lmax
        current_interaction =
            GaussletBases._assemble_atomic_injected_angular_interaction(
                benchmark.radial_operators,
                assembly,
            )
        interaction_pair_plan = 2 * maximum(assembly.shell_interaction_lexpand)

        @test benchmark.exact_common_lmax == 1
        @test bare_moment_lmax == benchmark.exact_common_lmax
        @test bare_moment_lmax < required_product_lmax
        @test interaction_moment_lmax > bare_moment_lmax
        @test interaction_moment_lmax ≥ required_product_lmax
        @test interaction_pair_plan ≥ required_product_lmax
        @test opnorm(current_interaction - transpose(current_interaction), Inf) ≤ 1.0e-10
        @test minimum(diag(current_interaction)) > 0.0
    end
end

@testset "Angular Ne legacy-trim HF branch repair" begin
    benchmark10 = _paper_style_angular_hf_style_benchmark_fixture(10)
    benchmark15 = _paper_style_angular_hf_style_benchmark_fixture(15)
    benchmark32 = _paper_style_angular_hf_style_benchmark_fixture(32)
    energy10 = benchmark10.scf_result.energy
    energy15 = benchmark15.scf_result.energy
    energy32 = benchmark32.scf_result.energy

    @test energy10 < -128.4
    @test energy15 < -128.5
    @test abs(energy10 - energy32) < 0.1
    @test abs(energy15 - energy32) < 0.05
end

@testset "Atomic injected angular small-ED benchmark" begin
    benchmark = _atomic_injected_angular_small_ed_benchmark_fixture()
    diagnostics = atomic_injected_angular_small_ed_diagnostics(benchmark)
    exact_problem = benchmark.exact_reference_problem

    @test benchmark isa AtomicInjectedAngularSmallEDBenchmark
    @test benchmark.hf_style isa AtomicInjectedAngularHFStyleBenchmark
    @test benchmark.orbital_count == size(benchmark.hf_style.one_body.overlap, 1)
    @test benchmark.state_count == benchmark.orbital_count^2
    @test size(benchmark.orbital_overlap) == (benchmark.orbital_count, benchmark.orbital_count)
    @test size(benchmark.orbital_one_body) == size(benchmark.orbital_overlap)
    @test size(benchmark.orbital_interaction) == size(benchmark.orbital_overlap)
    @test exact_problem isa AtomicIDATwoElectronProblem

    @test diagnostics.orbital_overlap_symmetry_error ≤ 1.0e-10
    @test diagnostics.orbital_one_body_symmetry_error ≤ 1.0e-8
    @test diagnostics.orbital_interaction_symmetry_error ≤ 1.0e-10
    @test diagnostics.orbital_overlap_identity_error ≤ 2.0e-6
    @test diagnostics.state_overlap_identity_error_estimate ≤ 3.0e-6
    @test diagnostics.min_orbital_overlap_eigenvalue > 1.0e-6
    @test diagnostics.min_state_overlap_eigenvalue_estimate > 1.0e-6
    @test diagnostics.state_interaction_diagonal_min > 0.0
    @test diagnostics.state_interaction_diagonal_max > diagnostics.state_interaction_diagonal_min
    @test diagnostics.full_converged
    @test diagnostics.full_residual ≤ 1.0e-7
    @test diagnostics.exact_reference_energy ≈ 13.020668426715936 atol = 1.0e-8 rtol = 1.0e-8
    @test diagnostics.full_energy ≈ 12.97749161121589 atol = 1.0e-5 rtol = 1.0e-8
    @test diagnostics.energy_difference_to_exact_reference < -1.0e-3
end

@testset "Atomic Ylm one-body layer" begin
    rb, grid, radial_ops, channels, atom = _quick_radial_atomic_fixture()[1:5]

    @test length(channels) == 9
    @test channels[1] == YlmChannel(0, 0)
    @test channels[2] == YlmChannel(1, -1)
    @test channels[4] == YlmChannel(1, 1)
    @test channels[end] == YlmChannel(2, 2)
    @test size(atom.overlap) == (9 * length(rb), 9 * length(rb))
    @test size(atom.hamiltonian) == size(atom.overlap)

    for (i, channel) in enumerate(channels)
        block = channel_range(atom, i)
        direct_block = radial_ops.kinetic + radial_ops.nuclear + centrifugal(radial_ops, channel.l)
        @test block == ((i - 1) * length(rb) + 1):(i * length(rb))
        @test channel_range(atom, channel) == block
        @test channel_overlap(atom, i) ≈ radial_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
        @test channel_overlap(atom, channel) ≈ radial_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
        @test channel_hamiltonian(atom, i) ≈ direct_block atol = 1.0e-12 rtol = 1.0e-12
        @test channel_hamiltonian(atom, channel) ≈ direct_block atol = 1.0e-12 rtol = 1.0e-12

        for j in 1:length(channels)
            i == j && continue
            other = channel_range(atom, j)
            @test atom.overlap[block, other] == zeros(Float64, length(rb), length(rb))
            @test atom.hamiltonian[block, other] == zeros(Float64, length(rb), length(rb))
        end
    end
end

@testset "Gaunt table backend" begin
    table = GaussletBases.build_gaunt_table(2; Lmax = 4, atol = 1.0e-14, basis = :complex)
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()

    @test table isa GaussletBases.GauntTable{Float64}
    @test GaussletBases.gaunt_lmax(table) == 2
    @test GaussletBases.gaunt_Lmax(table) == 4
    @test GaussletBases.gaunt_nnz(table) > 0
    @test GaussletBases.gaunt_hasblock(table, 0, 0, 0)
    @test !GaussletBases.gaunt_hasblock(table, 1, 0, 0)
    @test GaussletBases.gaunt_value(table, 0, 0, 0, 0, 0, 0) ≈ inv(sqrt(4 * pi)) atol = 1.0e-12 rtol = 1.0e-12

    for L in 0:GaussletBases.gaunt_Lmax(table)
        counted = 0
        for (l1, l2, entries) in GaussletBases.gaunt_each_block(table, L)
            @test GaussletBases.gaunt_legal_triple(l1, l2, L)
            previous = nothing
            for entry in entries
                counted += 1
                @test GaussletBases.gaunt_legal_ms(l1, entry.m1, l2, entry.m2, L, entry.M)
                @test entry.M == entry.m1 - entry.m2
                current = (entry.M, entry.m1, entry.m2)
                previous === nothing || @test previous <= current
                previous = current
            end
        end
        @test counted == GaussletBases.gaunt_nnz(table, L)

        expected_tensor = [
            GaussletBases.gaunt_value(
                table,
                L,
                channels[alpha].l,
                channels[alpha].m,
                channels[alphap].l,
                channels[alphap].m,
                M,
            )
            for alpha in 1:length(channels), alphap in 1:length(channels), M in -L:L
        ]
        @test gaunt_tensor(ida, L) ≈ expected_tensor atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testset "Angular kernel sectorization" begin
    _, _, _, channels, _, ida = _quick_radial_atomic_fixture()
    sectors = ida.angular_sectors
    nchannels = length(channels)

    @test sectors.nchannels == nchannels
    @test length(sectors.pair_to_sector) == nchannels^2
    @test length(sectors.pair_to_local) == nchannels^2
    @test sum(length(sector.pair_indices) for sector in sectors.sectors) == nchannels^2

    for (sector_index, sector) in enumerate(sectors.sectors)
        @test issorted(sector.pair_indices)
        matrix_by_L = sectors.sector_matrices

        for local_index in eachindex(sector.pair_indices)
            pair_index = sector.pair_indices[local_index]
            alpha = sector.left_channel_indices[local_index]
            beta = sector.right_channel_indices[local_index]
            @test channels[alpha].m + channels[beta].m == sector.msum
            @test sectors.pair_to_sector[pair_index] == sector_index
            @test sectors.pair_to_local[pair_index] == local_index
        end

        @test length(matrix_by_L) == GaussletBases.gaunt_Lmax(ida.gaunt_table) + 1
        for level in matrix_by_L
            @test size(level[sector_index]) == (length(sector.pair_indices), length(sector.pair_indices))
        end
    end

    for L in 0:GaussletBases.gaunt_Lmax(ida.gaunt_table)
        expected_kernel = _direct_dense_angular_kernel(ida.gaunt_table, channels, L)
        @test angular_kernel(ida, L) ≈ expected_kernel atol = 1.0e-12 rtol = 1.0e-12

        for alpha in 1:nchannels, beta in 1:nchannels, alphap in 1:nchannels, betap in 1:nchannels
            value = GaussletBases._angular_kernel_sector_value(sectors, L, alpha, beta, alphap, betap)
            @test value ≈ expected_kernel[alpha, alphap, beta, betap] atol = 1.0e-12 rtol = 1.0e-12

            if channels[alpha].m + channels[beta].m != channels[alphap].m + channels[betap].m
                @test value == 0.0
            end
        end
    end
end

@testset "Hydrogen Ylm spectrum" begin
    rb, grid, radial_ops, atom = _quick_hydrogen_ylm_fixture()

    @test norm(atom.overlap - I, Inf) ≤ 1.0e-5
    spectrum = sort(real(eigen(Hermitian(atom.hamiltonian)).values))
    E0 = spectrum[1]

    l1_channels = [YlmChannel(1, m) for m in -1:1]
    l1_energies = [
        begin
            @test norm(channel_overlap(atom, channel) - I, Inf) ≤ 1.0e-5
            minimum(real(eigen(Hermitian(channel_hamiltonian(atom, channel))).values))
        end
        for channel in l1_channels
    ]
    l2_channels = [YlmChannel(2, m) for m in -2:2]
    l2_energies = [
        begin
            @test norm(channel_overlap(atom, channel) - I, Inf) ≤ 1.0e-5
            minimum(real(eigen(Hermitian(channel_hamiltonian(atom, channel))).values))
        end
        for channel in l2_channels
    ]

    @test abs(E0 + 0.5) ≤ 5.0e-5
    @test maximum(abs.(l1_energies .- l1_energies[1])) ≤ 1.0e-10
    @test maximum(abs.(l2_energies .- l2_energies[1])) ≤ 1.0e-10
    @test abs(l1_energies[1] + 0.125) ≤ 5.0e-3
    @test abs(l2_energies[1] + 1.0 / 18.0) ≤ 5.0e-3
end
