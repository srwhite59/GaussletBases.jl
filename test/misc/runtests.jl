using SHA

@testset "occupied-first capture contract" begin
    CRG = GaussletBases.CartesianResidualGaussians
    S = Matrix{Float64}(I, 3, 3)
    X = [0.5 0.0 0.0; 0.0 0.8 0.0; 0.0 0.0 0.2; 0.0 0.0 0.0]
    Y = reshape([1.0, 0.0, 0.0], 3, 1)
    result = CRG.occupied_first_injection_geometry(
        4, X, S, Y; optional_capture_cutoff = 0.5)
    diagnostics = result.diagnostics
    @test diagnostics.occupied_base_capture_singular_values == [0.5]
    @test diagnostics.occupied_base_capture_min == 0.25
    @test diagnostics.occupied_recovery_after_mandatory_inclusion_singular_values ≈ [1.0]
    @test diagnostics.occupied_recovery_after_mandatory_inclusion_loss <=
        diagnostics.capture_tol
    @test diagnostics.complement_metric_minimum_eigenvalue ≈ 0.36
    @test diagnostics.raw_full_capture_range[1] >= -diagnostics.capture_tol
    @test diagnostics.raw_full_capture_range[2] <= 1 + diagnostics.capture_tol
    @test collect(diagnostics.raw_complement_capture_range) ≈ [0.04, 0.64]
    @test result.capture_eigenvalues ≈ [0.64, 0.04]

    malformed_X = copy(X)
    malformed_X[1, 1] = 1.1
    @test_throws ArgumentError CRG.occupied_first_injection_geometry(
        4, malformed_X, S, Y)
    @test_throws ArgumentError CRG._rg_validate_capture_eigenvalues(
        [-0.01, 0.5], diagnostics.capture_tol, "synthetic capture")
    @test_throws ArgumentError CRG._rg_validate_capture_eigenvalues(
        [0.5, 1.01], diagnostics.capture_tol, "synthetic capture")

    Y1 = reshape([1.0, 0.0, 0.0], 3, 1)
    Y2 = reshape([0.6, 0.8, 0.0], 3, 1)
    union_result = CRG.protected_occupied_union(S, [Y1, Y2])
    @test union_result.retained_rank == 2
    @test union_result.gram_eigenvalues ≈ [0.4, 1.6]
    @test maximum(union_result.recovery_losses) < 1.0e-12
    duplicate_result = CRG.protected_occupied_union(S, [Y1, Y1])
    @test duplicate_result.retained_rank == 1
    @test maximum(duplicate_result.recovery_losses) < 1.0e-12
    @test_throws ArgumentError CRG.protected_occupied_union(S, [1.1 .* Y1, Y2])
end

@testset "numerical-complete residual metric contract" begin
    CRG = GaussletBases.CartesianResidualGaussians
    labels = ["a_s1", "a_s2"]
    centers = NTuple{3,Float64}[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
    owners = [1, 1]
    S_AA = Matrix{Float64}(I, 2, 2)
    X = Diagonal([sqrt(1.0 - 2.0e-10), 1.0]) |> Matrix
    residual = CRG.build_residual_gaussian_basis(
        2, X, S_AA, labels, centers, owners;
        residual_occupation_cutoff = 1.0e-10,
        residual_injection_cutoff = 0.0,
        residual_compactness = nothing)
    @test residual.candidate_count == 2
    @test residual.residual_dimension == 1
    @test residual.residual_occupations[1] > 1.0e-10
    @test isnothing(residual.injected_G)
    @test norm(residual.T_G + X * residual.T_A, Inf) <= 1.0e-10
    @test norm(CRG.residual_gaussian_overlap(
        residual.T_G, residual.T_A, X, S_AA) - I, Inf) <= 1.0e-10

    malformed_X = copy(X)
    malformed_X[1, 1] = sqrt(1.0 + 1.0e-6)
    @test_throws ArgumentError CRG.build_residual_gaussian_basis(
        2, malformed_X, S_AA, labels, centers, owners;
        residual_occupation_cutoff = 1.0e-10,
        residual_injection_cutoff = 0.0,
        residual_compactness = nothing)
end

@testset "vendored legacy BasisSets provenance" begin
    path = joinpath(_PROJECT_ROOT, "data", "legacy", "BasisSets")
    text = read(path, String)
    body_start = findfirst("#BASIS SET:", text)
    @test body_start !== nothing
    body = text[first(body_start):end]
    normalized_body = join((replace(line, r"[ \t]+$" => "")
        for line in split(body, '\n'; keepempty = true)), '\n')
    @test bytes2hex(sha256(codeunits(normalized_body))) ==
        "b83883f4589234dd512eb760c95280854a2f42d007ab6e3533abda39a2829051"
    @test count(line -> startswith(line, "#BASIS SET:"), eachline(path)) == 60

    expected = [
        ("H", "cc-pVTZ", 6, 8),
        ("H", "cc-pVQZ", 10, 12),
        ("Be", "cc-pV5Z", 21, 42),
        ("Ne", "cc-pV5Z", 21, 54),
        ("Cr", "cc-pV5Z", 32, 434),
    ]
    for (atom, basis, shell_count, primitive_count) in expected
        shells, parsed_path = GaussletBases._legacy_basis_shells(
            atom, basis; basisfile = path)
        @test parsed_path == path
        @test length(shells) == shell_count
        @test sum(length(shell[2]) for shell in shells) == primitive_count
    end
end

@testset "REPL displays" begin
    (
        family,
        map,
        ub_spec,
        rb_spec,
        ub,
        rb,
        grid,
        ops,
        rep,
        channels,
        atom,
        ida,
        _tiny_rb,
        _tiny_grid,
        _tiny_radial_ops,
        tiny_ida,
        tiny_problem,
    ) = _quick_display_fixture()

    @test sprint(show, family) == "GaussletFamily(:G10)"
    @test occursin("AsinhMapping(", sprint(show, map))
    @test occursin("UniformBasisSpec(", sprint(show, ub_spec))
    @test occursin("RadialBasisSpec(", sprint(show, rb_spec))
    @test occursin("UniformBasis(length=3", sprint(show, ub))
    @test occursin("RadialBasis(length=6", sprint(show, rb))
    @test occursin("RadialQuadratureGrid(length=", sprint(show, grid))
    @test occursin("RadialAtomicOperators(size=(6, 6)", sprint(show, ops))
    @test occursin("BasisRepresentation1D(kind=:uniform", sprint(show, rep))
    @test occursin("YlmChannel(l=1, m=0)", sprint(show, YlmChannel(1, 0)))
    @test occursin("YlmChannelSet(lmax=2, nchannels=9)", sprint(show, channels))
    @test occursin("AtomicOneBodyOperators(nchannels=9", sprint(show, atom))
    @test occursin("AtomicOrbital(index=1, channel=YlmChannel(l=0, m=0), radial_index=1)", sprint(show, orbitals(ida)[1]))
    @test occursin("AtomicIDAOperators(nchannels=9", sprint(show, ida))
    @test occursin("AtomicIDATwoElectronState(index=1", sprint(show, two_electron_states(tiny_problem)[1]))
    @test occursin("AtomicIDATwoElectronProblem(norbitals=", sprint(show, tiny_problem))
end

if _RUN_SLOW_TESTS
    @testset "REPL displays (slow advanced objects)" begin
        (
            _family,
            _map,
            _ub_spec,
            hb_spec,
            _rb_spec,
            _ub,
            hb,
            _rb,
            _grid,
            _ops,
            _rep,
            partition,
            hierarchy,
            pgdg,
            augmented_pgdg,
            spec,
            global_layer,
            contracted_layer,
            _channels,
            _atom,
            _ida,
        ) = _slow_display_fixture()

        @test occursin("HalfLineBasisSpec(", sprint(show, hb_spec))
        @test occursin("HalfLineBasis(length=", sprint(show, hb))
        @test occursin("BasisPartition1D(nbasis=3, nboxes=3)", sprint(show, partition))
        @test occursin("BasisBox1D(index=1", sprint(show, boxes(partition)[1]))
        @test occursin("HierarchicalBasisPartition1D(nbasis=3, nboxes=5, nleaves=4)", sprint(show, hierarchy))
        @test occursin("HierarchicalBasisBox1D(index=4", sprint(show, boxes(hierarchy)[4]))
        @test occursin("LeafLocalPGDG1D(nleaves=4, nbasis=8, primitives_per_leaf=2, naugmented=0)", sprint(show, pgdg))
        @test occursin("LeafLocalPGDG1D(nleaves=4, nbasis=9, primitives_per_leaf=2, naugmented=1)", sprint(show, augmented_pgdg))
        @test occursin("LeafGaussianSpec1D(relative_position=0.5, width_scale=0.2)", sprint(show, spec))
        @test occursin("GlobalMappedPrimitiveLayer1D(nbasis=", sprint(show, global_layer))
        @test occursin("LeafBoxContractionLayer1D(nleaves=4", sprint(show, contracted_layer))
    end
end
