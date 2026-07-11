using LinearAlgebra
using Test
using JLD2

using GaussletBases

function _synthetic_external_gto_probes()
    orbitals = CartesianGaussianShellOrbitalRepresentation3D[
        CartesianGaussianShellOrbitalRepresentation3D(
            "s_core",
            (0, 0, 0),
            (0.0, 0.0, 0.0),
            [0.7],
            [1.0],
            :axiswise_normalized_cartesian_gaussian,
        ),
        CartesianGaussianShellOrbitalRepresentation3D(
            "pz_valence",
            (0, 0, 1),
            (0.0, 0.0, 0.0),
            [0.45],
            [1.0],
            :axiswise_normalized_cartesian_gaussian,
        ),
        CartesianGaussianShellOrbitalRepresentation3D(
            "s_shifted",
            (0, 0, 0),
            (0.35, 0.0, 0.0),
            [0.9],
            [1.0],
            :axiswise_normalized_cartesian_gaussian,
        ),
    ]
    return CartesianGaussianShellSupplementRepresentation3D(
        :external_gto_import_test,
        orbitals,
        (; source = :cartesian_external_gto_import_test),
    )
end

function _synthetic_working_basis()
    return build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 4.0),
        reference_spacing = 1.0,
    ))
end

function _orthonormal_source_coefficients(S_GG)
    return inv(cholesky(Symmetric(S_GG)).U)
end

function _source_packet(; alpha_spin = :restricted, beta = nothing, fingerprint = nothing)
    probes = _synthetic_external_gto_probes()
    S_GG = Matrix{Float64}(GaussletBases._cartesian_supplement_cross_overlap(probes, probes))
    C = _orthonormal_source_coefficients(S_GG)
    alpha = ExternalGTOOrbitalSpinBlock(alpha_spin, C[:, 1:2], [1.0, 1.0])
    return ExternalGTOOrbitalPacket(
        probes,
        S_GG,
        alpha;
        beta = beta,
        S_GG_fingerprint = something(fingerprint, external_gto_overlap_fingerprint(S_GG)),
        provenance = (; source = :synthetic_test_packet),
    )
end

function _synthetic_protected_member(working, supplement; final_dimension = 8)
    nG = length(working)^3
    G_L = zeros(Float64, nG, final_dimension)
    G_L[1:final_dimension, :] .= I(final_dimension)
    A_L = zeros(Float64, length(supplement.orbitals), final_dimension)
    H = Matrix(Diagonal(collect(1.0:final_dimension)))
    V = Matrix(Diagonal(fill(0.5, final_dimension)))
    recipe = (; nesting = :pqs, ns = 5, core_spacing = 0.3, s_factor = 1.0,
        basisname = "synthetic", lmax = 1, atom_symbols = ["H"],
        nuclear_charges = [1.0], atom_locations = [(0.0, 0.0, 0.0)],
        nup = 1, ndn = 1, source_artifact = "synthetic_source.jld2",
        source_commit = "synthetic-source-commit")
    expansion = GaussletBases._cartesian_coulomb_expansion_summary(
        :compact, GaussletBases.coulomb_gaussian_expansion(doacc = false))
    return (; recipe, inputs = (; base = working, supplement),
        raw = (; G_L, A_L), H, V, coulomb_expansion = expansion)
end

function _write_synthetic_protected_artifact(path, member; H1_L = member.H)
    n = size(H1_L, 1)
    diagnostics = (; B_min = 1.0, B_median = 1.0, B_max = 1.0,
        F_S_F_identity_error = 0.0, Z_S_M_Qperp_error = 0.0, Qperp_identity_error = 0.0,
        protected_span_min_sv = 1.0, L_identity_error = 0.0, M_L_diag_delta_max = 0.0,
        M_L_offdiag_max = 0.0, M_L_fro_delta = 0.0, H1_L_symmetry_error = 0.0,
        Vee_L_symmetry_error = 0.0, B_lt_0p999 = 0, B_lt_0p99 = 0,
        B_lt_0p98 = 0, B_lt_0p95 = 0, B_lt_0p9 = 0)
    recipe = member.recipe
    GaussletBases.write_protected_localized_ida_hamiltonian(path;
        H1_L, Vee_L = member.V, nup = recipe.nup, ndn = recipe.ndn,
        nuclear_charges = recipe.nuclear_charges,
        nuclear_positions = reduce(vcat,
            [reshape(collect(location), 1, :) for location in recipe.atom_locations]),
        sector_counts = (; base = n, compact_R = 0, protected_Z = 0, broad_Z = 0),
        diagnostics,
        provenance = (; source_artifact = recipe.source_artifact,
            source_commit = recipe.source_commit,
            current_commit = GaussletBases._plb_current_commit()),
        basis_controls = (; nesting = recipe.nesting, ns = recipe.ns,
            core_spacing = recipe.core_spacing, s_factor = recipe.s_factor,
            basisname = recipe.basisname, lmax = recipe.lmax),
        geometry_inputs = (; atom_symbols = recipe.atom_symbols,
            nuclear_charges = recipe.nuclear_charges,
            atom_locations = recipe.atom_locations),
        coulomb_expansion = member.coulomb_expansion)
    return path
end

function _mutated_sidecar(path, mutate!)
    target = tempname() * ".jld2"
    cp(path, target; force = true)
    jldopen(target, "r+") do file
        mutate!(file)
    end
    return target
end

_replace_jld2!(file, key, value) = (delete!(file, key); file[key] = value)

@testset "External GTO orbital import" begin
    working = _synthetic_working_basis()
    packet = _source_packet()
    result = import_external_gto_orbitals(working, packet)
    S_FG = gto_overlap_matrix(working, packet.probes)

    @test result.cross_overlap_size == (length(working)^3, length(packet.probes.orbitals))
    @test result.cross_overlap_finite
    @test result.ordering_fingerprint_valid
    @test result.S_GG_fingerprint_valid
    @test result.S_GG_symmetry_error < 1.0e-12
    @test result.S_GG_expected_error < 1.0e-12
    @test size(result.alpha.imported_coefficients) == (length(working)^3, 2)
    @test result.alpha.imported_coefficients ≈ S_FG * packet.alpha.coefficients
    @test result.alpha.source_orthogonality_error < 1.0e-12
    @test result.alpha.capture_matrix ≈
        transpose(result.alpha.imported_coefficients) * result.alpha.imported_coefficients
    @test result.alpha.density_trace_source == 2.0
    @test result.alpha.density_trace_capture ≈
        dot(packet.alpha.occupations, result.alpha.orbital_captures)
    @test result.alpha.density_trace_loss ≈
        result.alpha.density_trace_source - result.alpha.density_trace_capture
    @test result.alpha.worst_orbital_capture == minimum(result.alpha.orbital_captures)

    theta = 0.37
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    rotated_alpha = ExternalGTOOrbitalSpinBlock(
        :restricted,
        packet.alpha.coefficients * rotation,
        [1.0, 1.0],
    )
    rotated_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        packet.S_GG,
        rotated_alpha;
        provenance = (; source = :rotated_synthetic_test_packet),
    )
    rotated_result = import_external_gto_orbitals(working, rotated_packet)
    @test rotated_result.alpha.source_orthogonality_error < 1.0e-12
    @test rotated_result.alpha.density_trace_capture ≈
        result.alpha.density_trace_capture atol = 1.0e-12 rtol = 1.0e-12

    C = _orthonormal_source_coefficients(packet.S_GG)
    beta = ExternalGTOOrbitalSpinBlock(:beta, C[:, 2:2], [1.0])
    spin_packet = _source_packet(alpha_spin = :alpha, beta = beta)
    spin_result = import_external_gto_orbitals(working, spin_packet)
    @test spin_result.alpha.spin == :alpha
    @test spin_result.beta !== nothing
    @test spin_result.beta.spin == :beta
    @test size(spin_result.beta.imported_coefficients, 2) == 1
    @test spin_result.beta.source_orthogonality_error < 1.0e-12

    bad_fingerprint_packet = _source_packet(fingerprint = "not-the-source-overlap")
    @test_throws ArgumentError import_external_gto_orbitals(working, bad_fingerprint_packet)

    wrong_metric = copy(packet.S_GG)
    wrong_metric[1, 1] += 0.05
    wrong_metric_alpha = ExternalGTOOrbitalSpinBlock(
        :restricted,
        _orthonormal_source_coefficients(wrong_metric)[:, 1:2],
        [1.0, 1.0],
    )
    wrong_metric_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        wrong_metric,
        wrong_metric_alpha;
        S_GG_fingerprint = external_gto_overlap_fingerprint(wrong_metric),
        provenance = (; source = :wrong_but_refingerprinted_metric_test_packet),
    )
    @test_throws ArgumentError import_external_gto_orbitals(working, wrong_metric_packet)

    nonsymmetric_metric = copy(packet.S_GG)
    nonsymmetric_metric[1, 2] += 1.0e-3
    nonsymmetric_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        nonsymmetric_metric,
        packet.alpha;
        S_GG_fingerprint = external_gto_overlap_fingerprint(nonsymmetric_metric),
        provenance = (; source = :nonsymmetric_metric_test_packet),
    )
    @test_throws ArgumentError import_external_gto_orbitals(working, nonsymmetric_packet)

    bad_order_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        packet.S_GG,
        packet.alpha;
        ordering_fingerprint = "not-the-ao-order",
        provenance = (; source = :bad_ordering_test_packet),
    )
    @test_throws ArgumentError import_external_gto_orbitals(working, bad_order_packet)
    @test_throws ArgumentError import_external_gto_orbitals(nothing, bad_order_packet)

    @test_throws ArgumentError _source_packet(alpha_spin = :beta)
    @test_throws ArgumentError _source_packet(alpha_spin = :restricted, beta = beta)
    bad_beta = ExternalGTOOrbitalSpinBlock(:alpha, C[:, 2:2], [1.0])
    @test_throws ArgumentError _source_packet(alpha_spin = :alpha, beta = bad_beta)

    bad_alpha = ExternalGTOOrbitalSpinBlock(
        :restricted,
        packet.alpha.coefficients .* 2.0,
        [1.0, 1.0],
    )
    bad_orthogonality_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        packet.S_GG,
        bad_alpha;
        provenance = (; source = :bad_orthogonality_test_packet),
    )
    @test_throws ArgumentError import_external_gto_orbitals(working, bad_orthogonality_packet)
end

@testset "Protected external GTO representation sidecar" begin
    working = _synthetic_working_basis()
    base_packet = _source_packet()
    C = _orthonormal_source_coefficients(base_packet.S_GG)
    alpha = ExternalGTOOrbitalSpinBlock(:alpha, C[:, 1:2], [1.0, 1.0])
    beta = ExternalGTOOrbitalSpinBlock(:beta, C[:, 3:3], [1.0])
    packet = ExternalGTOOrbitalPacket(base_packet.probes, base_packet.S_GG,
        alpha; beta, provenance = (; source = :protected_sidecar_test))
    member = _synthetic_protected_member(working, packet.probes)
    raw_to_L = [member.raw.G_L; member.raw.A_L]
    direct = GaussletBases._cartesian_final_gto_cross_overlap_handoff(
        member.inputs.base, member.inputs.supplement, raw_to_L, packet.probes)
    value = GaussletBases.protected_localized_external_gto_import(member, packet)

    @test direct.diagnostics.orientation == :final_by_gto
    @test value.identity.orientation == :final_by_external
    @test value.identity.site_order_kind == :native
    @test value.cross_overlap == direct.cross_overlap
    @test size(value.cross_overlap) == (size(member.H, 1), length(packet.ao_labels))
    @test value.imported.alpha.imported_coefficients ≈
        value.cross_overlap * packet.alpha.coefficients
    @test value.imported.beta.imported_coefficients ≈
        value.cross_overlap * packet.beta.coefficients
    wide_packet = ExternalGTOOrbitalPacket(packet.probes, packet.S_GG, ExternalGTOOrbitalSpinBlock(:restricted, C, ones(3)))
    small_member = _synthetic_protected_member(working, packet.probes; final_dimension = 2)
    wide = GaussletBases.protected_localized_external_gto_import(small_member, wide_packet)
    @test (length(wide.source_metric.principal_singular_values), length(wide.source_metric.projected_gram_eigenvalues)) == (2, 3)
    @test (length(wide.occupied_capture.alpha.principal_singular_values), length(wide.occupied_capture.alpha.projected_gram_eigenvalues)) == (2, 3)
    wide_path = tempname() * ".jld2"
    GaussletBases.write_protected_localized_external_gto_representation(wide_path, wide)
    @test size(GaussletBases.read_protected_localized_external_gto_representation(wide_path; packet = wide_packet, member = small_member).cross_overlap) == (2, 3)

    theta = 0.29
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    rotated_packet = ExternalGTOOrbitalPacket(packet.probes, packet.S_GG,
        ExternalGTOOrbitalSpinBlock(:alpha,
            packet.alpha.coefficients * rotation, packet.alpha.occupations);
        beta = packet.beta, provenance = (; source = :rotated_sidecar_test))
    rotated = GaussletBases.protected_localized_external_gto_import(
        member, rotated_packet)
    @test rotated.imported.alpha.density_trace_capture ≈
        value.imported.alpha.density_trace_capture atol = 1.0e-12

    member_path = tempname() * ".jld2"
    _write_synthetic_protected_artifact(member_path, member)
    value = GaussletBases.protected_localized_external_gto_import(
        member, packet; member_artifact = member_path)
    path = tempname() * ".jld2"
    GaussletBases.write_protected_localized_external_gto_representation(path, value)
    loaded = GaussletBases.read_protected_localized_external_gto_representation(
        path; packet, member, protected_artifact = member_path)
    @test loaded.cross_overlap == value.cross_overlap
    @test loaded.imported.alpha.imported_coefficients == value.imported.alpha.imported_coefficients
    fresh = GaussletBases.external_gto_import_from_saved_representation(
        loaded, rotated_packet)
    @test fresh.alpha.imported_coefficients ≈ loaded.cross_overlap * rotated_packet.alpha.coefficients
    @test fresh.alpha.density_trace_capture ≈
        value.imported.alpha.density_trace_capture atol = 1.0e-12

    jldopen(path, "r") do file
        required = String["artifact_kind", "format_version", "convention_id",
            "convention_version", "site_order_kind", "orientation", "external/provenance/repr"]
        groups = (
            "cross_overlap" => (:S_LG, :fingerprint_sha256, :final_dimension, :external_dimension),
            "external" => (:ao_count, :ao_labels, :ordering_fingerprint_sha256,
                :S_GG_fingerprint_sha256, :alpha_coefficients_fingerprint_sha256,
                :beta_coefficients_fingerprint_sha256),
            "protected" => (:final_dimension, :artifact_kind, :convention_id,
                :recipe_fingerprint_sha256, :H1_L_fingerprint_sha256,
                :Vee_L_fingerprint_sha256, :source_artifact, :member_artifact,
                :source_commit, :current_commit),
            "protected/basis_controls" => (:nesting, :ns, :core_spacing, :s_factor, :basisname, :lmax),
            "protected/geometry_inputs" => (:atom_symbols, :nuclear_charges, :atom_locations, :nup, :ndn),
            "diagnostics/cross_overlap" => (:final_dimension, :external_dimension, :finite, :max_abs, :fingerprint_sha256),
            "diagnostics/source_metric" => (:retained_rank, :discarded_count, :tau, :metric_eigenvalue_min, :metric_eigenvalue_max, :projected_gram_eigenvalues,
                :principal_singular_values, :S_GG_symmetry_error, :S_GG_expected_error),
            "coulomb_expansion" => GaussletBases._CARTESIAN_COULOMB_EXPANSION_SUMMARY_KEYS)
        foreach(group -> append!(required, ["$(group.first)/$(name)" for name in group.second]), groups)
        for spin in (:alpha, :beta)
            append!(required, ["imported/$(spin)/coefficients", "imported/$(spin)/occupations"])
            append!(required, ["diagnostics/$(spin)/$(name)" for name in
                (:spin, :source_orthogonality_error, :capture_matrix, :orbital_captures,
                    :density_trace_source, :density_trace_capture, :density_trace_loss,
                    :worst_orbital_capture, :projected_gram_eigenvalues, :principal_singular_values)])
        end
        @test all(key -> haskey(file, key), required)
        for key in ("G_L", "A_L", "H1_L", "Vee_L", "protected/G_L",
                "protected/A_L", "protected/H1_L", "protected/Vee_L")
            @test !haskey(file, key)
        end
    end

    mutations = [
        file -> _replace_jld2!(file, "artifact_kind", :wrong),
        file -> _replace_jld2!(file, "format_version", 2),
        file -> _replace_jld2!(file, "convention_id", :wrong),
        file -> _replace_jld2!(file, "convention_version", 2),
        file -> _replace_jld2!(file, "site_order_kind", :z_order),
        file -> _replace_jld2!(file, "orientation", :external_by_final),
        file -> _replace_jld2!(file, "cross_overlap/final_dimension", size(member.H, 1) + 1),
        file -> _replace_jld2!(file, "cross_overlap/fingerprint_sha256", "wrong"),
        file -> _replace_jld2!(file, "diagnostics/alpha/density_trace_capture",
            file["diagnostics/alpha/density_trace_capture"] + 0.1),
        file -> _replace_jld2!(file, "diagnostics/alpha/principal_singular_values",
            file["diagnostics/alpha/principal_singular_values"] .+ 0.1),
        file -> _replace_jld2!(file, "diagnostics/alpha/source_orthogonality_error", 0.1),
        file -> _replace_jld2!(file,
            "diagnostics/source_metric/principal_singular_values",
            reverse(file["diagnostics/source_metric/principal_singular_values"])),
        file -> _replace_jld2!(file, "diagnostics/beta/spin", :alpha),
        file -> _replace_jld2!(file, "protected/artifact_kind", :wrong),
    ]
    for mutate! in mutations
        bad = _mutated_sidecar(path, mutate!)
        @test_throws Exception GaussletBases.read_protected_localized_external_gto_representation(bad)
    end
    bad_beta = _mutated_sidecar(path, file -> begin
        delete!(file, "external/beta_coefficients_fingerprint_sha256"); delete!(file, "imported/beta/coefficients"); delete!(file, "imported/beta/occupations")
        foreach(name -> String(name) == "spin" || delete!(file, "diagnostics/beta/$(name)"), collect(keys(file["diagnostics/beta"])))
    end)
    @test_throws ArgumentError GaussletBases.read_protected_localized_external_gto_representation(bad_beta)
    missing_coulomb = _mutated_sidecar(path,
        file -> delete!(file, "coulomb_expansion/policy"))
    @test_throws ArgumentError GaussletBases.read_protected_localized_external_gto_representation(missing_coulomb)

    wrong_packet = ExternalGTOOrbitalPacket(packet.probes, packet.S_GG,
        packet.alpha; beta = packet.beta, ordering_fingerprint = "wrong")
    @test_throws ArgumentError GaussletBases.read_protected_localized_external_gto_representation(
        path; packet = wrong_packet)
    @test_throws ArgumentError GaussletBases.read_protected_localized_external_gto_representation(
        path; packet = rotated_packet)
    wrong_occupations = ExternalGTOOrbitalPacket(packet.probes, packet.S_GG,
        ExternalGTOOrbitalSpinBlock(:alpha, packet.alpha.coefficients, [0.5, 0.5]);
        beta = packet.beta)
    @test_throws ArgumentError GaussletBases.read_protected_localized_external_gto_representation(
        path; packet = wrong_occupations)
    wrong_member = merge(member, (; H = member.H + I))
    @test_throws ArgumentError GaussletBases.read_protected_localized_external_gto_representation(
        path; member = wrong_member)

    unusable_member = merge(member, (; raw = nothing))
    @test_throws ArgumentError GaussletBases.protected_localized_external_gto_import(
        unusable_member, wrong_packet)

    tampered = GaussletBases.read_protected_localized_external_gto_representation(path)
    tampered.cross_overlap[1, 1] += 1.0e-3
    @test_throws ArgumentError GaussletBases.external_gto_import_from_saved_representation(
        tampered, packet)
    @test_throws ArgumentError GaussletBases.write_protected_localized_external_gto_representation(
        tempname() * ".jld2", tampered)

    relocated_member = tempname() * ".jld2"
    cp(member_path, relocated_member)
    @test GaussletBases.read_protected_localized_external_gto_representation(path; protected_artifact = relocated_member).identity.protected.member_artifact == member_path
    jldopen(member_path, "r+") do file
        _replace_jld2!(file, "H1_L", file["H1_L"] + I)
    end
    @test_throws ArgumentError GaussletBases.protected_localized_external_gto_import(member, packet; member_artifact = member_path)
end
