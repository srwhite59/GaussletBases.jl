@testset "Positive spherical cubature catalog" begin
    names = cubature_sphere_point_set_names()
    @test length(names) == 15
    copied_names = cubature_sphere_point_set_names()
    pop!(copied_names)
    @test cubature_sphere_point_set_names() == names
    @test_throws ArgumentError cubature_sphere_point_set(:unknown_cubature)
    for name in names
        rule = cubature_sphere_point_set(name)
        @test rule isa CubatureSpherePointSet
        @test rule === cubature_sphere_point_set(Symbol(name))
        @test rule.name == name && rule.family in (:ih, :oh)
        points, weights = rule.point_set.coordinates, rule.cubature_weights
        @test size(points) == (rule.point_set.order, 3)
        @test rule.point_set.cardinality == length(weights) == size(points, 1)
        @test maximum(abs, vec(sum(abs2, points; dims=2)) .- 1) < 1e-13
        @test all(isfinite, weights) && minimum(weights) > 0
        @test abs(sum(weights) - 4pi) < 1e-13
        @test 2 * rule.default_l_flat <= rule.cubature_degree
        @test rule.default_l_inject >= rule.default_l_flat >= 0
        @test isfinite(rule.default_beta) && rule.default_beta > 0
        channels = ylm_channels(rule.cubature_degree)
        Y = [GaussletBases._real_spherical_harmonic(ch, view(points, a, :))
            for ch in channels.channel_data, a in eachindex(weights)]
        moments = Y * weights
        @test abs(moments[1] - sqrt(4pi)) < 1e-13
        @test maximum(abs, moments[2:end]) < 1e-13
        nlow = (rule.cubature_degree ÷ 2 + 1)^2
        low = Y[1:nlow, :]
        @test opnorm(low * Diagonal(weights) * low' - I, Inf) < 1e-13
    end
end

@testset "Polar coiflet moments and locality" begin
    flat_moments(b, l) = reduce(vcat, [GaussletBases._shell_local_ylm_moment_block(b, k) for k in 0:l])
    target(x, w, l) = [sqrt(w[a]) * GaussletBases._real_spherical_harmonic(ch, view(x, a, :))
        for ch in ylm_channels(l).channel_data, a in eachindex(w)]
    gl = eigen(SymTridiagonal(zeros(48), [k/sqrt(4k*k-1) for k in 1:47]))
    directions = [sqrt(1-z*z)*f(2pi*(j-0.5)/96) for z in gl.values, j in 1:96, f in (cos, sin)]
    xyz = hcat(vec(directions[:, :, 1]), vec(directions[:, :, 2]), repeat(gl.values, 96))
    qw = repeat(2 .* gl.vectors[1, :].^2 .* (2pi/96), 96)
    function locality(b, points)
        rho = evaluate_angular_basis(b, xyz).^2
        d2 = 2 .* (1 .- clamp.(xyz * points', -1.0, 1.0))
        return sqrt(sum(qw .* rho .* d2)/size(points, 1)), sum(qw .* rho .* d2.^3)/size(points, 1)
    end
    for name in ("Ih_12", "McLaren_24", "Ih_32", "Oh_50_Unif")
        cs = cubature_sphere_point_set(name)
        parent = shell_local_angular_profile(cs.point_set; beta=cs.default_beta, l_inject=cs.default_l_inject)
        p = coiflet_angular_profile(cs)
        actual = flat_moments(p.basis, cs.default_l_flat)
        reference = target(cs.point_set.coordinates, cs.cubature_weights, cs.default_l_flat)
        @test maximum(abs, actual-reference) < 1e-13
        @test opnorm(p.basis.final_overlap-I, Inf) < 1e-11
        @test minimum(actual[1, :]) > 0
        @test p.gauge_metadata.min_site_alignment >= 0.98
        @test p.basis.shell_weights == parent.basis.shell_weights
        @test p.gauge_metadata.mixed_signs == parent.gauge_metadata.mixed_signs
        @test p === coiflet_angular_profile(Symbol(name))
        @test maximum(abs, eigvals(Symmetric(p.basis.final_kinetic)) -
            eigvals(Symmetric(parent.basis.final_kinetic))) < 1e-11
        lp, lc = locality(parent.basis, cs.point_set.coordinates), locality(p.basis, cs.point_set.coordinates)
        @test lc[1] <= 1.02lp[1]
        name == "Ih_32" && @test lc[2] <= 1.02lp[2]
        @info "Polar coiflet qualification" rule=name moment_error=maximum(abs, actual-reference) alignment=p.gauge_metadata.min_site_alignment sigma_parent=lp[1] sigma_coiflet=lc[1] tau6_parent=lp[2] tau6_coiflet=lc[2]
    end
    cs = cubature_sphere_point_set("Oh_14")
    x, w = cs.point_set.coordinates, cs.cubature_weights
    axes = [a for a in eachindex(w) if count(v -> abs(v) > 0.99, view(x, a, :)) == 1]
    corners = setdiff(collect(eachindex(w)), axes)
    alt = copy(w)
    alt[axes] .= 1.1 .* w[axes]
    alt[corners] .= (4pi-sum(alt[axes]))/length(corners)
    opts = (l_flat=1, beta=cs.default_beta, l_inject=cs.default_l_inject)
    original = coiflet_angular_profile(cs; l_flat=1)
    alternative = coiflet_angular_profile(cs.point_set, alt; opts...)
    whitened = coiflet_angular_profile(cs.point_set, w; opts..., whiten_target=true)
    @test original.key != alternative.key && original.key != whitened.key
    @test original.profile_id != alternative.profile_id
    @test original === coiflet_angular_profile(cs.point_set, BigFloat.(w); opts...)
    @test alternative === coiflet_angular_profile(cs.point_set, alt; opts...)
    @test whitened === coiflet_angular_profile(cs.point_set, w; opts..., whiten_target=true)
    @test maximum(abs, flat_moments(original.basis, 1)-target(x, w, 1)) < 1e-13
    @test maximum(abs, flat_moments(alternative.basis, 1)-target(x, alt, 1)) < 1e-13
    @test maximum(abs, flat_moments(whitened.basis, 1)-target(x, w, 1)) < 1e-13
    ps = sphere_point_set(50)
    A = target(ps.coordinates, ones(50), 6)
    w0 = fill(4pi/50, 50)
    wr = w0 + pinv(A; rtol=1e-12) * (vcat(sqrt(4pi), zeros(size(A, 1)-1))-A*w0)
    @test minimum(wr) > 0
    @test_throws ArgumentError coiflet_angular_profile(ps, wr; l_flat=3)
    rounded = coiflet_angular_profile(ps, wr; l_flat=3, whiten_target=true)
    raw = target(ps.coordinates, wr, 3)
    effective = inv(sqrt(Symmetric(raw*raw'))) * raw
    @test maximum(abs, flat_moments(rounded.basis, 3)-effective) < 1e-11
    @test maximum(abs, effective-raw) > 1e-11
    @test opnorm(rounded.basis.final_overlap-I, Inf) < 1e-11
    @test minimum(flat_moments(rounded.basis, 0)) > 0
    @test rounded.gauge_metadata.min_site_alignment >= 0.95
    @test rounded === coiflet_angular_profile(ps, wr; l_flat=3, whiten_target=true)
    @test_throws ArgumentError coiflet_angular_profile("Ih_32"; l_inject=-1)
    @test_throws ArgumentError coiflet_angular_profile("Ih_32"; l_inject=2)
    @test_throws ArgumentError coiflet_angular_profile("Ih_12"; l_flat=3)
    @test_throws ArgumentError coiflet_angular_profile("Ih_12"; l_flat=-1)
    @test_throws ArgumentError coiflet_angular_profile(cs; whiten=:invalid)
    for bad in (w[1:end-1], -w, fill(0.0, 14), fill(Inf, 14), fill(NaN, 14), 2w)
        @test_throws ArgumentError coiflet_angular_profile(cs.point_set, bad; opts...)
    end
    tiny = BigFloat.(w)
    tiny[2] += tiny[1] - big"1e-1000"
    tiny[1] = big"1e-1000"
    @test_throws ArgumentError GaussletBases._coiflet_weights_digest(tiny, 14, false)
    @test_throws ArgumentError coiflet_angular_profile(cs.point_set, w .* (1+1e-7); opts...)
    admitted = coiflet_angular_profile(cs.point_set, w .* (1+1e-7); opts..., whiten_target=true)
    @test opnorm(admitted.basis.final_overlap-I, Inf) < 1e-11
    @test_throws ArgumentError coiflet_angular_profile(cs.point_set, w .* (1+1e-4); opts..., whiten_target=true)
    bad_gram = copy(w)
    bad_gram[1] += 1e-3
    bad_gram[2] -= 1e-3
    @test_throws ArgumentError coiflet_angular_profile(cs.point_set, bad_gram; opts..., whiten_target=true)
    @test_throws ArgumentError GaussletBases._polar_coiflet_isometry(0.99Matrix{Float64}(I, 2, 2), Matrix{Float64}(I, 2, 2))
    @test_throws DimensionMismatch GaussletBases._polar_coiflet_isometry(ones(2, 3), ones(2, 2))
end

@testset "Explicit angular profile assembly" begin
    radii = collect(range(0.08, 4.5; length=8))
    for orders in (fill(32, 8), [15, 24, 32, 15, 24, 32, 15, 32])
        profiles = shell_local_angular_profile.(orders)
        legacy = build_atomic_shell_local_angular_assembly(radii; shell_orders=orders)
        actual = build_atomic_shell_local_angular_assembly(radii, profiles)
        @test actual.overlap == legacy.overlap
        @test actual.kinetic == legacy.kinetic
        @test actual.shell_orders == actual.shell_dimensions == orders
        @test actual.shell_offsets == legacy.shell_offsets
        @test actual.shell_interaction_lcap == legacy.shell_interaction_lcap
        @test actual.shell_interaction_lexpand == legacy.shell_interaction_lexpand
        @test actual.shell_radii == radii && actual.shell_radii !== radii
        @test actual.profiles !== profiles && all(actual.profiles[i] === profiles[i] for i in 1:8)
    end
    cs = cubature_sphere_point_set("Oh_14")
    parent = shell_local_angular_profile(cs.point_set; beta=cs.default_beta, l_inject=cs.default_l_inject)
    polar, p32 = coiflet_angular_profile(cs), coiflet_angular_profile("Ih_32")
    @test parent.key != polar.key && parent.basis.final_count == polar.basis.final_count
    for profiles in (fill(polar, 8), [parent, polar, p32, parent, p32, polar, parent, polar])
        actual = build_atomic_shell_local_angular_assembly(radii, profiles)
        keys = [p.key for p in profiles]
        nk = length(unique(keys))
        @test length(unique(objectid.(actual.shell_kinetic_moment_blocks))) == nk
        @test length(unique(objectid.(actual.shell_interaction_moment_blocks))) == nk
        @test length(unique(objectid.(actual.shell_moment_blocks))) == 8
        @test length(unique(objectid.(actual.shells))) == nk
        @test length(unique(objectid.(vec(actual.overlap_blocks)))) == nk^2
        @test length(unique(objectid.(vec(actual.kinetic_blocks)))) == nk^2
        @test all(actual.shells[i] === profiles[i].basis for i in 1:8)
        @test actual.shell_orders == [p.basis.final_count for p in profiles]
        @test maximum(abs, actual.overlap-actual.overlap') < 1e-13
        @test maximum(abs, actual.kinetic-actual.kinetic') < 1e-13
        representatives = [findfirst(==(key), keys) for key in unique(keys)]
        for i in representatives
            kinetic = GaussletBases._build_shell_local_kinetic_moment_blocks(actual.shells[i])
            interaction = GaussletBases._build_shell_local_interaction_moment_blocks(actual.shells[i])
            @test all(actual.shell_kinetic_moment_blocks[i][l] == block for (l, block) in kinetic.blocks)
            @test all(actual.shell_interaction_moment_blocks[i][l] == block for (l, block) in interaction.blocks)
            for j in representatives
                @test actual.overlap_blocks[i, j] == GaussletBases._shell_local_overlap(actual.shells[i], actual.shells[j])
                @test actual.kinetic_blocks[i, j] == GaussletBases._shell_local_pair_kinetic(
                    actual.shells[i], actual.shells[j], actual.shell_kinetic_moment_blocks[i],
                    actual.shell_kinetic_moment_blocks[j], actual.shell_kinetic_lcap[i], actual.shell_kinetic_lcap[j])
            end
        end
        if nk > 1
            @test norm(actual.overlap_blocks[1, 1]-actual.overlap_blocks[1, 2]) > 1e-6
        end
    end
    @test_throws ArgumentError build_atomic_shell_local_angular_assembly(Float64[], ShellLocalAngularProfile[])
    @test_throws ArgumentError build_atomic_shell_local_angular_assembly(radii, [polar])
    @test_throws ArgumentError build_atomic_shell_local_angular_assembly([1.0], ShellLocalAngularProfile[])
    # Independent Gaunt quadrature versus projected IDA using ALL returned moments.
    e = eigen(SymTridiagonal(zeros(32), [k/sqrt(4k*k-1) for k in 1:31]))
    xyz = reduce(vcat, [[sqrt(1-z*z)*cos(2pi*(j-0.5)/64) sqrt(1-z*z)*sin(2pi*(j-0.5)/64) z]
        for z in e.values for j in 1:64])
    weights = [w*2pi/64 for w in 2 .* e.vectors[1, :].^2 for j in 1:64]
    channels = ylm_channels(4).channel_data
    Y = [GaussletBases._real_spherical_harmonic(ch, view(xyz, q, :)) for ch in channels, q in eachindex(weights)]
    products(M) = reduce(hcat, [M[i, :] .* M[j, :] for i in axes(M, 1) for j in axes(M, 1)])
    for (name_a, la, name_b, lb) in (("Oh_6", 0, "Oh_38", 2), ("Ih_12", 1, "Ih_32", 1), ("McLaren_24", 1, "Oh_38", 2))
        assembly = build_atomic_shell_local_angular_assembly([0.25, 1.25], coiflet_angular_profile.([name_a, name_b]))
        da, db = assembly.shell_interaction_moment_blocks
        na, nb = assembly.shell_dimensions
        interaction = zeros(na, nb)
        for l in 0:maximum(assembly.shell_interaction_lcap)
            ma, mb = get(da, l, zeros(2l+1, na)), get(db, l, zeros(2l+1, nb))
            interaction .+= (4pi/(2l+1)*0.25^l/1.25^(l+1)) .* (ma' * mb)
        end
        wa = sqrt(4pi) .* vec(assembly.shell_moment_blocks[1][0])
        wb = sqrt(4pi) .* vec(assembly.shell_moment_blocks[2][0])
        interaction ./= wa * wb'
        ma = reduce(vcat, [assembly.shell_moment_blocks[1][l] for l in 0:la])
        mb = reduce(vcat, [assembly.shell_moment_blocks[2][l] for l in 0:lb])
        projected = products(ma)' * interaction * products(mb)
        lm = (2min(la, lb)+1)^2
        scale = [sqrt(4pi/(2ch.l+1)*0.25^ch.l/1.25^(ch.l+1)) for ch in channels[1:lm]]
        weighted = Y[1:lm, :] .* weights'
        ga = scale .* (weighted * products(Y[1:(la+1)^2, :]))
        gb = scale .* (weighted * products(Y[1:(lb+1)^2, :]))
        error = maximum(abs, projected-ga'*gb)
        @info "Explicit-profile four-index Coulomb oracle" name_a name_b carried_lmax=maximum(assembly.shell_interaction_lcap) max_error_Ha=error
        @test error < 1e-13
    end
end

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
    for name in (:U, :T)
        reference = getproperty(dense, name)
        delta = abs.(getproperty(full, name)[1:n, 1:n] - reference)
        limits = ifelse.(abs.(reference) .>= 1e-6 .* diag(dense.U),
            1e-10 .* abs.(reference), 1e-14 .* diag(dense.U))
        @info "Contact base-block preservation" matrix=name max_absolute_delta=maximum(delta) max_error_gate_ratio=maximum(delta ./ limits)
        @test all(delta .<= limits)
    end
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
        link = joinpath(dir, "dangling.jld2")
        symlink(joinpath(dir, "missing.jld2"), link)
        @test_throws ArgumentError write_atomic_contact_integrals_jld2(link, result)
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
