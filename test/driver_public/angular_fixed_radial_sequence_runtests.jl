@testset "Atomic contact geometry" begin
    rb, _ = _quick_radial_operator_fixture()
    # Actual frozen McLaren24 coordinates/order, not the vendored 24-point rule.
    a, b, c = 0.86624681810782056, 0.42251865376111153, 0.26663540151670473
    points = [-a -b c; -a b -c; a -b -c; a b c; -a -c -b; -a c b;
        a -c b; a c -b; -b -a -c; -b a c; b -a c; b a -c;
        -b -c a; -b c -a; b -c -a; b c a; -c -a b; -c a -b;
        c -a -b; c a b; -c -b -a; -c b a; c -b a; c b -a]
    points ./= sqrt.(sum(abs2, points; dims=2))
    custom = SpherePointSet(24, 24, points, GaussletBases._sphere_point_set_nn_ratio(points),
        SpherePointSetProvenance("mcl24", "frozen contact input", "", "ordered coordinates"))
    p24, p32 = shell_local_angular_profile(custom), shell_local_angular_profile(32)
    shells = [p32.basis; fill(p24.basis, length(rb)-1)]
    dims = [s.final_count for s in shells]
    offsets = cumsum([1; dims[1:end-1]])
    ranges = [o:o+d-1 for (o, d) in zip(offsets, dims)]
    labels = ["$g:$k" for (g, d) in enumerate(dims) for k in 1:d]
    e = eigen(SymTridiagonal(zeros(160), [k/sqrt(4k*k-1) for k in 1:159]))
    r, w = 20 .* (e.values .+ 1), 40 .* e.vectors[1, :].^2
    grid = RadialQuadratureGrid(r, w; mapping=mapping(rb))
    ez = eigen(SymTridiagonal(zeros(96), [k/sqrt(4k*k-1) for k in 1:95]))
    xyz = reduce(vcat, [[sqrt(1-z*z)*cos(2pi*k/192) sqrt(1-z*z)*sin(2pi*k/192) z]
        for z in ez.values for k in 0:191])
    wa = [v*2pi/192 for v in 2 .* ez.vectors[1, :].^2 for k in 0:191]
    aq = (directions=xyz, weights=wa)
    source = (producer_commit="frozen-test", basis_manifest_sha256="contact-test")
    kw = (radial_grid=grid, angular_quadrature=aq, labels=labels, provenance=source)
    result = atomic_contact_integrals(rb, shells, ranges; kw...)
    dense = materialize_atomic_contact_integrals(result)
    @test size(dense.U) == (sum(dims), sum(dims))
    @test minimum(dense.U) >= 0
    @test dense.U ≈ dense.U' rtol=1e-13
    @test diag(dense.T) ≈ diag(dense.U) rtol=1e-13
    @test norm(dense.T-dense.T') > 1e-6
    @test result.diagnostics.row_sums ≈ dense.w .* vec(sum(dense.U; dims=2)) rtol=1e-13
    @test evaluate_angular_basis(p24, xyz) == evaluate_angular_basis(p24.basis, xyz)
    @test result.site_map[last(first(ranges))].label == labels[last(first(ranges))]
    @test_throws ArgumentError evaluate_angular_basis(p24, ones(2, 3))
    @test_throws ArgumentError evaluate_angular_basis(p24, fill(NaN, 2, 3))
    @test_throws ArgumentError evaluate_angular_basis(p24, ones(3, 2))
    @test_throws ArgumentError atomic_contact_integrals(rb, shells, reverse(ranges); kw...)
    @test_throws ArgumentError atomic_contact_integrals(rb, shells, ranges; merge(kw,
        (angular_quadrature=(directions=xyz, weights=-wa),))...)
    assembly = build_atomic_shell_local_angular_assembly(centers(rb); shell_orders=fill(32, length(rb)))
    from_assembly = atomic_contact_integrals(rb, assembly; radial_grid=grid, angular_quadrature=aq)
    @test length(from_assembly.site_map) == 32length(rb)
    @test from_assembly.site_map[1].label == assembly.profiles[1].labels[1]
    @test from_assembly.factors.radial.R == result.factors.radial.R
    Z = 2.0
    hydrogenic = 2Z^1.5 .* r .* exp.(-Z .* r)
    channels = [YlmChannel(0, 0), YlmChannel(1, -1)]
    samples = hcat(hydrogenic, hydrogenic)
    supplement = (grid=grid, reduced_values=samples, channels=channels)
    augmented = atomic_contact_integrals(rb, shells, ranges; kw...,
        supplements=supplement, supplement_labels=["1s", "p_y"])
    full = materialize_atomic_contact_integrals(augmented)
    n = length(dense.w)
    @test full.U[1:n, 1:n] == dense.U && full.T[1:n, 1:n] == dense.T
    @test full.w[n+1] ≈ 64pi/Z^3 rtol=2e-13
    @test full.U[n+1, n+1] ≈ Z^3/(8pi) rtol=2e-13
    @test full.T[n+1, n+1] ≈ full.U[n+1, n+1] rtol=1e-14
    @test full.w[n+2] < 1e-27
    @test norm(full.U[1:n, n+1:n+2]) > 0
    @test norm(full.T[1:n, n+1:n+2]-full.T[n+1:n+2, 1:n]') > 1e-6
    flipped = atomic_contact_integrals(rb, shells, ranges; kw...,
        supplements=(grid=grid, reduced_values=-samples, channels=channels), supplement_labels=["1s", "p_y"])
    dflip = materialize_atomic_contact_integrals(flipped)
    @test dflip.U == full.U && dflip.w == full.w
    @test dflip.T[1:n, n+1:end] == -full.T[1:n, n+1:end]
    @test dflip.T[n+1:end, 1:n] == -full.T[n+1:end, 1:n]
    @test_throws ArgumentError atomic_contact_integrals(rb, shells, ranges; kw..., supplements=supplement)
    @test_throws ArgumentError atomic_contact_integrals(rb, shells, ranges; kw...,
        supplements=(grid=grid, reduced_values=complex.(samples), channels=channels), supplement_labels=["s", "p"])
    @test_throws ArgumentError atomic_contact_integrals(rb, shells, ranges; kw...,
        supplements=(grid=RadialQuadratureGrid(r, 2w), reduced_values=samples, channels=channels), supplement_labels=["s", "p"])
    # Stream a few actual product columns, checking all base/residual directions.
    select = [1, dims[1]+1, n+1, n+2]
    V = GaussletBases._basis_values_matrix(rb, r)
    F = hcat(evaluate_angular_basis(p32, xyz)[:, 1], evaluate_angular_basis(p24, xyz)[:, 1],
        fill(1/sqrt(4pi), length(wa)), [GaussletBases._real_spherical_harmonic(channels[2], x) for x in eachrow(xyz)])
    reduced = hcat(V[:, 1:2], samples)
    oracleU, oracleT = zeros(4, 4), zeros(4, 4)
    for j in eachindex(r)
        values = F .* transpose(reduced[j, :] ./ r[j])
        measure = w[j]*r[j]^2 .* wa
        oracleU .+= transpose(values.^2) * (measure .* values.^2)
        oracleT .+= transpose(values.^3) * (measure .* values)
    end
    @test full.U[select, select] ≈ oracleU rtol=2e-12
    @test full.T[select, select] ≈ oracleT rtol=2e-12
    # Six original B24 s/p samples (receipt 00c68015...), a discrete-input check,
    # not a new continuum quadrature or a reconstruction of the residuals.
    raw = [0.007094191204987262 5.242491753709607e-6 0.039312807473682654 0.02270981658381326 0.003969817681926234;
        0.04571752661928699 6.02827369999803e-5 -0.010138340484549957 0.011099507089308479 0.012589912106637978;
        0.2950504139145793 0.00011398384007470904 0.1582651604492071 0.051940157573943324 0.005870534644926449;
        1.918050634131147 0.001987820169373544 -0.3428949832362452 -0.48889702279902286 -0.5115507298335632;
        12.341714270269183 0.014992198375384985 -0.044687087391415255 -0.13310901461255187 -0.020696428637050485;
        79.9976322026874 0.006074747064896995 1.7647846385034586e-8 5.543979941427574e-8 9.613182214490913e-8]
    sg = RadialQuadratureGrid(raw[:, 1], raw[:, 2]; mapping=mapping(rb))
    sv = raw[:, [3, 4, 5, 5, 5]]
    sc = [YlmChannel(0, 0), YlmChannel(0, 0), YlmChannel(1, -1), YlmChannel(1, 0), YlmChannel(1, 1)]
    sr = atomic_contact_integrals(rb, shells, ranges; merge(kw, (radial_grid=sg,))...,
        supplements=(grid=sg, reduced_values=sv, channels=sc), supplement_labels=string.(sc))
    sd = materialize_atomic_contact_integrals(sr)
    @test sd.U[n+1, n+1] ≈ sum(raw[:, 2].*sv[:, 1].^4 ./ raw[:, 1].^2)/(4pi) rtol=1e-13
    @test sd.U[n+3, n+3] ≈ sum(raw[:, 2].*sv[:, 3].^4 ./ raw[:, 1].^2)*9/(20pi) rtol=1e-13
    @test sd.U[n+1, n+3] ≈ sum(raw[:, 2].*sv[:, 1].^2 .* sv[:, 3].^2 ./ raw[:, 1].^2)/(4pi) rtol=1e-13
    @test sr.site_map[end].label == string(sc[end])
    mktempdir() do dir
        path = joinpath(dir, "contacts.jld2")
        @test write_atomic_contact_integrals_jld2(path, augmented) == path
        jldopen(path, "r") do file
            @test file["U"] == full.U && file["T"] == full.T && file["w"] == full.w
            @test file["site_map"] == augmented.site_map && file["provenance"] == source
            @test file["factors"].angular.profile_index == augmented.factors.angular.profile_index
        end
        @test_throws ArgumentError write_atomic_contact_integrals_jld2(path, result)
    end
end

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

@testset "Fixed-radial sequence identity records the interaction cap" begin
    _, _, radial_ops = _paper_style_angular_anchor_radial_fixture(; Z = 2.0, lmax = 2)
    auto_seq = build_atomic_fixed_radial_angular_sequence(radial_ops, [10])
    stored_seq = build_atomic_fixed_radial_angular_sequence(radial_ops, [10]; interaction_lmax = :stored)
    @test auto_seq.radial_basis_id == stored_seq.radial_basis_id
    @test auto_seq.levels[1].payload.interaction != stored_seq.levels[1].payload.interaction
    @test auto_seq.sequence_id != stored_seq.sequence_id
    @test auto_seq.levels[1].level_id != stored_seq.levels[1].level_id
    level = auto_seq.levels[1]
    legacy = typeof(level.payload)(ntuple(i -> getfield(level.payload, i), 13)...)
    unknown = AtomicFixedRadialAngularSequenceLevel(ntuple(i -> getfield(level, i), 10)..., legacy)
    cap = level.payload.interaction_lmax_plan.interaction_lmax
    mktempdir() do dir
        for (entry, expected) in ((level, (cap, cap, "auto", false)),
                                 (unknown, ("unknown", "unknown", "unknown", "unknown")))
            path = joinpath(dir, string(expected[3], ".jld2"))
            write_atomic_fixed_radial_angular_level_jld2(path, entry)
            jldopen(path, "r") do file
                @test tuple((file["meta/manifest/interaction/" * key] for key in
                    ("multipole_Lmax_used", "multipole_Lmax_required", "multipole_Lmax_mode", "truncated"))...) == expected
            end
        end
    end
end
