using LinearAlgebra
using Test

using GaussletBases

const CRD = GaussletBases.CartesianReferenceDensity

function _unconverged_packet(packet)
    diagnostics = merge(packet.rhf_diagnostics, (; converged = false))
    return typeof(packet)(packet.spec, packet.supplement, packet.overlap,
        packet.overlap_fingerprint, packet.occupied_coefficients,
        packet.occupations, packet.orbital_energies, packet.density_matrix,
        diagnostics, packet.density_fit, packet.potential_fit,
        packet.validation, packet.provenance)
end

function _packet_with_occupations(packet, occupations)
    return typeof(packet)(packet.spec, packet.supplement, packet.overlap,
        packet.overlap_fingerprint, packet.occupied_coefficients,
        occupations, packet.orbital_energies, packet.density_matrix,
        packet.rhf_diagnostics, packet.density_fit, packet.potential_fit,
        packet.validation, packet.provenance)
end

function _packet_with_overlap_fingerprint(packet, fingerprint)
    return typeof(packet)(packet.spec, packet.supplement, packet.overlap,
        String(fingerprint), packet.occupied_coefficients,
        packet.occupations, packet.orbital_energies, packet.density_matrix,
        packet.rhf_diagnostics, packet.density_fit, packet.potential_fit,
        packet.validation, packet.provenance)
end

function _be_core_fixture()
    spec = CRD.be_core_reference_packet_spec(ns = 3, core_spacing = 0.15)
    packet = CRD.build_atomic_hf_reference_packet(spec)
    system = (;
        atom_symbols = ["Be"],
        nuclear_charges = [4.0],
        atom_locations = [(0.0, 0.0, 0.0)],
        nup = 2,
        ndn = 2,
    )
    basis = (;
        q = spec.ns,
        core_spacing = spec.core_spacing,
        radius = spec.radius,
        reference_spacing = spec.reference_spacing,
        d = spec.core_spacing,
        s_factor = spec.s_factor,
    )
    base = GaussletBases.cartesian_base_working_basis(
        system; basis, supplemented = true)
    ham = GaussletBases.cartesian_base_hamiltonian_assembly(base)
    J0 = CRD.atomic_reference_packet_terminal_hartree_gg(
        base, packet; source = :potential_fit)
    exact = CRD.atomic_reference_packet_terminal_hartree_gg(
        base, packet; source = :density_fit)
    return (; spec, packet, base, ham, J0, exact)
end

function _orthonormal_packet_reference(packet)
    S = Symmetric(packet.overlap)
    values, vectors = eigen(S)
    Ssqrt = vectors * Diagonal(sqrt.(values)) * transpose(vectors)
    Sinvsqrt = vectors * Diagonal(1 ./ sqrt.(values)) * transpose(vectors)
    C_final = Ssqrt * packet.occupied_coefficients
    pair = gaussian_coulomb_pair_matrix(
        packet.supplement; max_orbitals = length(packet.supplement.orbitals))
    J_AA, _ = CRD._coulomb_exchange_matrices(pair, packet.density_matrix)
    J_final = Sinvsqrt' * J_AA * Sinvsqrt
    E0 = sum(packet.density_matrix .* J_AA)
    V = Matrix{Float64}(I, size(J_final, 1), size(J_final, 2)) .+
        0.01 .* abs.(0.5 .* (J_final .+ transpose(J_final)))
    V = 0.5 .* (V .+ transpose(V))
    return (; C_final, occupations = packet.occupations, J_final, V, E0)
end

function _terminal_control_reference(J)
    values, vectors = eigen(Symmetric(0.5 .* (J .+ transpose(J))))
    return reshape(vectors[:, argmax(values)], :, 1), [2.0]
end

@testset "Screened Hartree correction assembly" begin
    fixture = _be_core_fixture()
    packet_reference = _orthonormal_packet_reference(fixture.packet)

    correction = CRD.build_screened_hartree_correction(
        packet_reference.V,
        packet_reference.J_final,
        packet_reference.E0,
        packet_reference.C_final,
        packet_reference.occupations;
        packet = fixture.packet,
        source = :density_fit,
    )

    @test size(correction.delta_one_body) == size(packet_reference.V)
    @test all(isfinite, correction.delta_one_body)
    @test correction.diagnostics.V_IDA_input_symmetry_error < 1.0e-12
    @test correction.diagnostics.J0_G_input_symmetry_error < 1.0e-12
    @test correction.diagnostics.Delta_J0_symmetry_error < 1.0e-12
    @test correction.diagnostics.J0_G_finite
    @test correction.diagnostics.Delta_J0_finite
    @test correction.diagnostics.q0_charge ≈
        sum(packet_reference.occupations) atol = 1.0e-10
    @test correction.diagnostics.P0_trace ≈
        sum(packet_reference.occupations) atol = 1.0e-10
    @test correction.diagnostics.represented_orthogonality_error < 1.0e-10
    @test correction.diagnostics.direct_hartree_energy_anchor_error ≈
        0.0 atol = 1.0e-10
    @test correction.diagnostics.derivative_anchor_error ≈ 0.0 atol = 1.0e-12
    @test correction.energy_accounting ==
        "Delta_J0 + C is accounted as screened direct electron-electron/Hartree interaction, not as a physical kinetic/nuclear one-body change."
    @test correction.packet_summary.atom == "Be"
    @test correction.packet_summary.density_fit_self_energy_relative_error < 1.0e-8
    @test packet_reference.E0 ≈
        fixture.packet.density_fit.row.exact_self_energy atol = 1.0e-8

    packet_count = length(fixture.packet.supplement.orbitals)
    packet_embedding = CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, fixture.packet.supplement, fixture.packet.overlap,
        ones(Int, packet_count); owner_index = 1, center = fixture.spec.center,
        supplement_indices = 1:packet_count)
    @test sum(packet_embedding.occupations) == fixture.spec.electron_count
    exact_mapping = packet_embedding.overlap_mapping
    @test exact_mapping.stored_packet_fingerprint == fixture.packet.overlap_fingerprint
    @test exact_mapping.recomputed_packet_fingerprint == fixture.packet.overlap_fingerprint
    @test exact_mapping.mapped_block_fingerprint == fixture.packet.overlap_fingerprint
    @test exact_mapping.mapped_fingerprint_exact_match
    @test exact_mapping.max_abs_error == 0.0
    @test exact_mapping.inf_error == 0.0
    @test exact_mapping.overlap_atol == 1.0e-10
    corrupt_fingerprint = _packet_with_overlap_fingerprint(fixture.packet, "corrupt")
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        corrupt_fingerprint, fixture.packet.supplement, fixture.packet.overlap,
        ones(Int, packet_count); owner_index = 1, center = fixture.spec.center,
        supplement_indices = 1:packet_count)
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, fixture.packet.supplement, fixture.packet.overlap,
        ones(Int, packet_count); owner_index = 1, center = fixture.spec.center,
        supplement_indices = 1:packet_count, overlap_atol = NaN)
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, fixture.packet.supplement, fixture.packet.overlap,
        ones(Int, packet_count); owner_index = 1, center = fixture.spec.center,
        supplement_indices = 1:packet_count, overlap_atol = Inf)
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, fixture.packet.supplement, fixture.packet.overlap,
        ones(Int, packet_count); owner_index = 1, center = fixture.spec.center,
        supplement_indices = 1:packet_count, overlap_atol = -1.0)

    orbital_type = GaussletBases.CartesianGaussianShellOrbitalRepresentation3D
    centers = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))
    molecular_orbitals = orbital_type[]
    for (prefix, center) in zip(("a_", "b_"), centers), orbital in fixture.packet.supplement.orbitals
        push!(molecular_orbitals, orbital_type(prefix * orbital.label,
            orbital.angular_powers, center, copy(orbital.exponents),
            copy(orbital.coefficients), orbital.primitive_normalization))
    end
    molecular_supplement = typeof(fixture.packet.supplement)(
        fixture.packet.supplement.supplement_kind, molecular_orbitals,
        fixture.packet.supplement.metadata)
    molecular_overlap = Matrix{Float64}(GaussletBases._cartesian_supplement_cross_overlap(
        molecular_supplement, molecular_supplement))
    molecular_owners = vcat(fill(1, packet_count), fill(2, packet_count))
    owner2 = (packet_count + 1):(2packet_count)
    mapped_fingerprint = CRD._matrix_fingerprint(molecular_overlap[owner2, owner2])
    equivalent_overlap = copy(molecular_overlap)
    if mapped_fingerprint == fixture.packet.overlap_fingerprint
        # Exact translation can reproduce the same bytes; one symmetric ULP makes
        # the numerical-equivalence branch deterministic without changing physics.
        equivalent_overlap[first(owner2), first(owner2)] =
            nextfloat(equivalent_overlap[first(owner2), first(owner2)])
    end
    translated_embedding = CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, molecular_supplement, equivalent_overlap, molecular_owners;
        owner_index = 2, center = centers[2], supplement_indices = owner2)
    translated_mapping = translated_embedding.overlap_mapping
    @test !translated_mapping.mapped_fingerprint_exact_match
    @test translated_mapping.mapped_block_fingerprint !=
        translated_mapping.stored_packet_fingerprint
    @test translated_mapping.inf_error < 1.0e-10
    @test translated_mapping.max_abs_error <= translated_mapping.inf_error
    @test norm(transpose(translated_embedding.Y) * equivalent_overlap *
        translated_embedding.Y - I, Inf) < 1.0e-10
    @test sum(translated_embedding.occupations) == fixture.spec.electron_count
    too_different = copy(equivalent_overlap)
    too_different[first(owner2), first(owner2)] += 2.0e-10
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, molecular_supplement, too_different, molecular_owners;
        owner_index = 2, center = centers[2], supplement_indices = owner2)
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, molecular_supplement, equivalent_overlap, molecular_owners;
        owner_index = 1, center = centers[2], supplement_indices = owner2)
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, molecular_supplement, equivalent_overlap, molecular_owners;
        owner_index = 2, center = (0.0, 0.0, 2.1), supplement_indices = owner2)
    wrong_basis_metadata = merge(molecular_supplement.metadata,
        (; basis_name = "wrong-basis"))
    wrong_basis_supplement = typeof(molecular_supplement)(
        molecular_supplement.supplement_kind, molecular_orbitals,
        wrong_basis_metadata)
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, wrong_basis_supplement, equivalent_overlap, molecular_owners;
        owner_index = 2, center = centers[2], supplement_indices = owner2)
    packet_labels = string.(getproperty.(fixture.packet.supplement.orbitals, :label))
    px = findfirst(label -> startswith(label, "px"), packet_labels)
    py = findfirst(==("py" * packet_labels[px][3:end]), packet_labels)
    reordered_orbitals = copy(fixture.packet.supplement.orbitals)
    reordered_orbitals[px], reordered_orbitals[py] =
        reordered_orbitals[py], reordered_orbitals[px]
    reordered_supplement = typeof(fixture.packet.supplement)(
        fixture.packet.supplement.supplement_kind, reordered_orbitals,
        fixture.packet.supplement.metadata)
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        fixture.packet, reordered_supplement, fixture.packet.overlap,
        ones(Int, packet_count); owner_index = 1, center = fixture.spec.center,
        supplement_indices = 1:packet_count)
    bad_occupations = _packet_with_occupations(fixture.packet, [1.5])
    @test_throws ArgumentError CRD.atomic_reference_packet_occupied_embedding(
        bad_occupations, fixture.packet.supplement, fixture.packet.overlap,
        ones(Int, packet_count); owner_index = 1, center = fixture.spec.center,
        supplement_indices = 1:packet_count)

    C_A = reshape([1.0, 0.0, 0.0], 3, 1)
    C_B = reshape([0.6, 0.8, 0.0], 3, 1)
    occupations = [[2.0], [2.0]]
    additive_reference = CRD.represented_additive_reference_p0_q0(
        [C_A, C_B], occupations)
    additive_J = [1.2 0.1 0.0; 0.1 0.9 0.2; 0.0 0.2 0.7]
    additive_V = [1.0 0.2 0.1; 0.2 1.1 0.3; 0.1 0.3 0.8]
    additive_E0 = sum(additive_reference.P0 .* additive_J)
    V_before, J_before = copy(additive_V), copy(additive_J)
    additive = CRD.build_additive_screened_hartree_correction(
        additive_V, additive_J, additive_E0, [C_A, C_B], occupations)
    @test additive.diagnostics.P0_trace ≈ 4.0 atol = 1.0e-12
    @test additive.diagnostics.q0_charge ≈ 4.0 atol = 1.0e-12
    @test additive.diagnostics.additive_reference.block_traces ≈ [2.0, 2.0]
    @test additive.diagnostics.additive_reference.interpacket_occupied_overlap_max ≈ 0.6
    @test additive.diagnostics.direct_hartree_energy_anchor_error ≈ 0.0 atol = 1.0e-12
    @test additive.diagnostics.derivative_anchor_error ≈ 0.0 atol = 1.0e-12
    @test additive_V == V_before
    @test additive_J == J_before
    @test_throws ArgumentError CRD.build_additive_screened_hartree_correction(
        additive_V, additive_J, additive_E0, [C_A, 1.1 .* C_B], occupations)

    placements = [(; packet = fixture.packet, center = (-2.0, 0.0, 0.0)),
        (; packet = fixture.packet, center = (2.0, 0.0, 0.0))]
    cloud_energy = CRD.atomic_reference_packet_additive_density_energy(placements)
    @test cloud_energy.pair_energies[1, 2] ≈ cloud_energy.pair_energies[2, 1] atol = 1.0e-12
    @test cloud_energy.pair_energies[1, 1] ≈
        fixture.packet.density_fit.row.fit_self_energy atol = 1.0e-10
    @test cloud_energy.total ≈ sum(cloud_energy.pair_energies) atol = 1.0e-12
    @test cloud_energy.total ≈ cloud_energy.combined_cloud_oracle atol = 1.0e-12

    packet_raw = CRD.atomic_reference_packet_fitted_potential_raw_blocks(
        fixture.base, fixture.packet.supplement, fixture.packet)
    @test all(isfinite, packet_raw.GG) && all(isfinite, packet_raw.GA) &&
        all(isfinite, packet_raw.AA)
    @test norm(packet_raw.GG - transpose(packet_raw.GG), Inf) < 1.0e-10
    @test norm(packet_raw.AA - transpose(packet_raw.AA), Inf) < 1.0e-10
    @test abs(sum(fixture.packet.density_matrix .* packet_raw.AA) -
        fixture.packet.density_fit.row.fit_self_energy) <= 1.0e-9

    terminal_C, terminal_occ = _terminal_control_reference(fixture.J0)
    terminal_diagnostic = CRD.build_atomic_packet_screened_hartree_correction(
        fixture.base,
        fixture.ham,
        fixture.packet;
        reference_coefficients = terminal_C,
        occupations = terminal_occ,
        check_density_fit_matrix = true,
        diagnostic_only = true,
    )
    @test terminal_diagnostic.E0_G ==
        fixture.packet.density_fit.row.fit_self_energy
    @test all(isfinite, terminal_diagnostic.J0_G)
    @test all(isfinite, terminal_diagnostic.delta_one_body)
    @test terminal_diagnostic.diagnostics.derivative_anchor_error ≈
        0.0 atol = 1.0e-12
    @test terminal_diagnostic.diagnostics.potential_fit_exact_relative_fro < 1.0e-3
    @test abs(terminal_diagnostic.diagnostics.direct_hartree_energy_anchor_error) > 1.0e-6

    @test_throws ArgumentError CRD.build_atomic_packet_screened_hartree_correction(
        fixture.base,
        fixture.ham,
        fixture.packet;
        reference_coefficients = terminal_C,
        occupations = terminal_occ,
        check_density_fit_matrix = true,
    )

    shifted_center = (0.2, 0.0, 0.0)
    shifted_potential = CRD.atomic_reference_packet_terminal_hartree_gg(
        fixture.base, fixture.packet; source = :potential_fit, center = shifted_center)
    shifted_density = CRD.atomic_reference_packet_terminal_hartree_gg(
        fixture.base, fixture.packet; source = :density_fit, center = shifted_center)
    @test norm(shifted_potential - shifted_density) /
        max(norm(shifted_density), eps(Float64)) < 1.0e-3
    @test norm(shifted_potential - fixture.J0, Inf) > 1.0e-8

    cloud, cloud_density = CRD._packet_cloud_from_readback(fixture.packet)
    compact_expansion = CRD._atomic_reference_coulomb_expansion(:compact)
    explicit_density = GaussletBases.CartesianGaussianRawBlocks.atomic_reference_hartree_gg_block(
        fixture.base.terminal_basis,
        fixture.base.parent.parent_axis_bundle_object,
        cloud,
        cloud_density;
        expansion = compact_expansion)
    @test fixture.exact == explicit_density.GG
    @test explicit_density.diagnostics.coulomb_expansion_term_count == 45

    unconverged = _unconverged_packet(fixture.packet)
    @test_throws ArgumentError CRD.atomic_reference_packet_additive_density_energy(
        [(; packet = unconverged, center = fixture.spec.center)])
    @test_throws ArgumentError CRD.build_screened_hartree_correction(
        packet_reference.V,
        packet_reference.J_final,
        packet_reference.E0,
        packet_reference.C_final,
        packet_reference.occupations;
        packet = unconverged,
        diagnostic_only = true)
    @test_throws ArgumentError CRD.build_atomic_packet_screened_hartree_correction(
        fixture.base,
        fixture.ham,
        unconverged;
        reference_coefficients = terminal_C,
        occupations = terminal_occ,
        diagnostic_only = true)

    packet_path = joinpath(mktempdir(), "screened_hartree_packet.jld2")
    CRD.write_atomic_hf_reference_packet(packet_path, fixture.packet)
    readback = CRD.read_atomic_hf_reference_packet(packet_path)
    readback_embedding = CRD.atomic_reference_packet_occupied_embedding(
        readback, fixture.packet.supplement, fixture.packet.overlap,
        ones(Int, packet_count); owner_index = 1, center = fixture.spec.center,
        supplement_indices = 1:packet_count)
    @test readback_embedding.Y == packet_embedding.Y
    @test readback_embedding.overlap_mapping == packet_embedding.overlap_mapping
    unconverged_readback = merge(readback, (; rhf_converged = false))
    @test_throws ArgumentError CRD.build_atomic_packet_screened_hartree_correction(
        fixture.base,
        fixture.ham,
        unconverged_readback;
        reference_coefficients = terminal_C,
        occupations = terminal_occ,
        diagnostic_only = true)

    bad_reference = copy(packet_reference.C_final)
    bad_reference[:, 1] .*= 1.1
    @test_throws ArgumentError CRD.build_screened_hartree_correction(
        packet_reference.V,
        packet_reference.J_final,
        packet_reference.E0,
        bad_reference,
        packet_reference.occupations,
    )

    @test_throws ArgumentError CRD.build_screened_hartree_correction(
        packet_reference.V,
        packet_reference.J_final,
        1.1 * packet_reference.E0,
        packet_reference.C_final,
        packet_reference.occupations,
    )

    asymmetric_J = copy(packet_reference.J_final)
    asymmetric_J[1, 2] += 1.0e-4
    @test_throws ArgumentError CRD.build_screened_hartree_correction(
        packet_reference.V,
        asymmetric_J,
        packet_reference.E0,
        packet_reference.C_final,
        packet_reference.occupations,
    )
    asym_diagnostic = CRD.build_screened_hartree_correction(
        packet_reference.V,
        asymmetric_J,
        packet_reference.E0,
        packet_reference.C_final,
        packet_reference.occupations;
        diagnostic_only = true,
    )
    @test asym_diagnostic.diagnostics.J0_G_input_symmetry_error > 0.0
    @test asym_diagnostic.diagnostics.Delta_J0_input_symmetry_error > 0.0

    @test_throws DimensionMismatch CRD.build_screened_hartree_correction(
        packet_reference.V[1:(end - 1), 1:(end - 1)],
        packet_reference.J_final,
        packet_reference.E0,
        packet_reference.C_final,
        packet_reference.occupations,
    )
    @test_throws DimensionMismatch CRD.build_atomic_packet_screened_hartree_correction(
        fixture.base,
        fixture.ham,
        fixture.packet;
        reference_coefficients = terminal_C[1:(end - 1), :],
        occupations = terminal_occ,
    )
end
