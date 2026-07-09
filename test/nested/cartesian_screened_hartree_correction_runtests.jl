using LinearAlgebra
using Test

using GaussletBases

const CRD = GaussletBases.CartesianReferenceDensity

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
