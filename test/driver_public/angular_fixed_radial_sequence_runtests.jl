@testset "Atomic fixed-radial angular public sequence" begin
    sequence = _paper_style_fixed_radial_angular_sequence_fixture()
    levels = sequence.levels
    overlaps = (
        (sequence.adjacent_overlaps[1], 1, 2, :adjacent),
        (sequence.adjacent_overlaps[2], 2, 3, :adjacent),
        (sequence.direct_overlaps[1], 1, 3, :direct),
    )
    nr = length(sequence.shell_ids)
    @test sequence isa AtomicFixedRadialAngularSequence
    @test (sequence.N_sph_values, length(levels), length(sequence.adjacent_overlaps),
        length(sequence.direct_overlaps)) == ([10, 15, 32], 3, 2, 1)
    @test sequence.shell_ids == collect(1:nr)
    @test issorted(sequence.shell_centers_r)
    for (level, order) in zip(levels, (10, 15, 32))
        profile = shell_local_angular_profile(order)
        @test level isa AtomicFixedRadialAngularSequenceLevel
        @test level.N_sph == order
        @test level.radial_basis_id == sequence.radial_basis_id
        @test level.shell_ids == sequence.shell_ids
        @test level.shell_centers_r == sequence.shell_centers_r
        @test level.profile === profile
        @test level.shell_dimensions == fill(profile.basis.final_count, nr)
        @test all(p -> p === profile, level.payload.one_body.angular_assembly.profiles)
    end
    for (sidecar, source_index, target_index, pair_kind) in overlaps
        source, target = levels[source_index], levels[target_index]
        @test sidecar isa AtomicFixedRadialAngularSequenceOverlapSidecar
        @test (sidecar.source_N_sph, sidecar.target_N_sph) ==
              (source.N_sph, target.N_sph)
        @test (sidecar.source_level_id, sidecar.target_level_id) ==
              (source.level_id, target.level_id)
        @test (sidecar.source_profile_id, sidecar.target_profile_id) ==
              (source.profile.profile_id, target.profile.profile_id)
        @test (sidecar.source_gauge_version, sidecar.target_gauge_version) ==
              (source.profile.key.gauge_version, target.profile.key.gauge_version)
        @test sidecar.source_labels == source.profile.labels
        @test sidecar.target_labels == target.profile.labels
        @test size(sidecar.overlap) ==
              (source.profile.basis.final_count, target.profile.basis.final_count)
        @test sidecar.shell_independent
        @test sidecar.pair_kind == pair_kind
    end
    level_payload = atomic_fixed_radial_angular_level_dense_payload(levels[1])
    @test level_payload.payload["H1"] ≈ levels[1].payload.hamiltonian atol = 0.0 rtol = 0.0
    @test level_payload.payload["Vee"] ≈ levels[1].payload.interaction atol = 0.0 rtol = 0.0
    @test level_payload.payload["shell_ids"] == sequence.shell_ids
    @test level_payload.payload["shell_centers_r"] == sequence.shell_centers_r
    @test level_payload.payload["shell_dimensions"] == levels[1].shell_dimensions
    @test level_payload.payload["within_shell_labels"] == levels[1].profile.labels
    @test level_payload.bridge_meta["sequence_id"] == sequence.sequence_id
    @test level_payload.bridge_meta["level_id"] == levels[1].level_id
    @test level_payload.bridge_meta["angular_profile_id"] == levels[1].profile.profile_id
    @test level_payload.bridge_meta["gauge_version"] ==
          string(levels[1].profile.key.gauge_version)
    direct = sequence.direct_overlaps[1]
    overlap_payload = atomic_fixed_radial_angular_overlap_sidecar_payload(direct)
    @test overlap_payload.payload["overlap"] ≈ direct.overlap atol = 0.0 rtol = 0.0
    @test overlap_payload.payload["source_labels"] == direct.source_labels
    @test overlap_payload.payload["target_labels"] == direct.target_labels
    mktempdir() do dir
        path = joinpath(dir, "he_fixed_radial_10_32_overlap.jld2")
        @test write_atomic_fixed_radial_angular_overlap_sidecar_jld2(path, direct) == path
        jldopen(path, "r") do file
            @test String(file["bridge/sequence_id"]) == sequence.sequence_id
            @test String(file["bridge/source_level_id"]) == levels[1].level_id
            @test String(file["bridge/target_level_id"]) == levels[3].level_id
            @test String(file["bridge/source_profile_id"]) == levels[1].profile.profile_id
            @test String(file["bridge/target_profile_id"]) == levels[3].profile.profile_id
            @test String(file["bridge/source_gauge_version"]) ==
                  string(levels[1].profile.key.gauge_version)
            @test String(file["bridge/target_gauge_version"]) ==
                  string(levels[3].profile.key.gauge_version)
            @test String(file["bridge/pair_kind"]) == "direct"
            @test String.(file["source_labels"]) == direct.source_labels
            @test String.(file["target_labels"]) == direct.target_labels
            @test Matrix{Float64}(file["overlap"]) ≈ direct.overlap atol = 0.0 rtol = 0.0
        end
    end
end
