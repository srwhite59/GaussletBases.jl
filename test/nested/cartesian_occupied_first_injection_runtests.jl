using LinearAlgebra
using Test

using GaussletBases

const CRD = GaussletBases.CartesianReferenceDensity
const CRG = GaussletBases.CartesianResidualGaussians
const CFR = GaussletBases.CartesianFinalBasisRealization

function _be_full_reference_packet()
    spec = CRD.atomic_hf_reference_packet_spec(;
        atom = "Be",
        nuclear_charge = 4.0,
        electron_count = 4,
        basis_name = "cc-pV5Z",
        lmax = 1,
        ns = 5,
        core_spacing = 1.2 / (4.0 * (5 - 1)),
        fill_shell_convention = "Be closed-shell all-electron occupied-first injection test",
    )
    return CRD.build_atomic_hf_reference_packet(spec)
end

function _channel(orbital)
    l = sum(orbital.angular_powers)
    l == 0 && return "s"
    l == 1 && return "p"
    return "l$(l)"
end

function _packet_metadata(packet)
    labels = String[orbital.label for orbital in packet.supplement.orbitals]
    channels = String[_channel(orbital) for orbital in packet.supplement.orbitals]
    owners = ones(Int, length(labels))
    return labels, owners, channels
end

function _atom_system(packet)
    electrons = round(Int, packet.spec.nuclear_charge)
    return (; atom_symbols = [packet.spec.atom],
        nuclear_charges = [packet.spec.nuclear_charge],
        atom_locations = [packet.spec.center],
        nup = div(electrons, 2), ndn = div(electrons, 2))
end

function _atom_basis(packet)
    spec = packet.spec
    return (; ns = spec.ns, core_spacing = spec.core_spacing,
        radius = spec.radius, reference_spacing = spec.reference_spacing,
        s_factor = spec.s_factor)
end

function _due_diligence_summary(base)
    due = base.terminal_due_diligence
    rows = isempty(due.terminal_rows) ? base.terminal_inventory.rows :
        due.terminal_rows
    retained(row) = hasproperty(row, :retained_count) ? row.retained_count :
        row.final_cols
    row_summary = [(;
        region_kind = row.region_kind,
        source_mode_shape = hasproperty(row, :source_mode_shape) ?
            row.source_mode_shape : nothing,
        retained_count = retained(row),
        warning_flags = hasproperty(row, :warning_flags) ?
            row.warning_flags : (:none,)) for row in rows]
    return (; parent_bounds = due.geometry.parent_physical_bounds,
        parent_axis_counts = due.geometry.parent_axis_counts,
        radius = due.geometry.radius,
        padding = due.geometry.xmax_transverse,
        final_dimension = due.dimensions.base_final_dimension,
        retained_counts = [row.retained_count for row in row_summary],
        shell_slab_topology = [row.region_kind for row in row_summary],
        row_summary,
        warnings = due.warnings)
end

function _run_packet_case(packet, case_label, expected_optional)
    base = GaussletBases.cartesian_base_working_basis(
        _atom_system(packet); basis = _atom_basis(packet), supplemented = true)
    supplement = packet.supplement
    X = CFR._terminal_residual_mixed_overlap(
        base.terminal_basis, base.parent.parent_axis_bundle_object, supplement)
    S = Matrix{Float64}(GaussletBases._cartesian_supplement_cross_overlap(
        supplement, supplement))
    Y = packet.occupied_coefficients
    labels, owners, channels = _packet_metadata(packet)
    result = CRG.occupied_first_injection_geometry(
        base.terminal_basis.final_dimension,
        X,
        S,
        Y;
        optional_capture_cutoff = 0.99,
        labels,
        owners,
        channels,
        provenance = (; source = :atomic_hf_reference_packet, case_label,
            overlap_fingerprint = packet.overlap_fingerprint),
    )
    k = size(Y, 2)
    diagnostics = result.diagnostics
    due = _due_diligence_summary(base)
    @test result.diagnostics.provenance.case_label == case_label
    @test S ≈ packet.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test diagnostics.occupied_orthogonality_error < 1.0e-10
    @test length(diagnostics.occupied_base_capture_singular_values) == k
    @test 0.0 < diagnostics.occupied_base_capture_min <= 1.0 + diagnostics.capture_tol
    @test diagnostics.occupied_recovery_after_mandatory_inclusion_loss < 1.0e-10
    @test minimum(diagnostics.occupied_recovery_after_mandatory_inclusion_singular_values) ≈
        1.0 atol = 1.0e-10
    @test diagnostics.complement_metric_minimum_eigenvalue >= -diagnostics.capture_tol
    @test diagnostics.raw_full_capture_range[1] >= -diagnostics.capture_tol
    @test diagnostics.raw_full_capture_range[2] <= 1.0 + diagnostics.capture_tol
    @test diagnostics.raw_complement_capture_range[1] >= -diagnostics.capture_tol
    @test diagnostics.raw_complement_capture_range[2] <= 1.0 + diagnostics.capture_tol
    @test size(result.mandatory_A) == (size(S, 1), k)
    @test size(result.optional_A, 2) == expected_optional
    @test size(result.injected_A, 2) == k + expected_optional
    @test result.injected_G ≈ X * result.injected_A atol = 1.0e-12
    @test length(result.full_capture_eigenvalues) == size(S, 1)
    @test length(result.capture_eigenvalues) == size(S, 1) - k
    @test diagnostics.optional_kept_count == expected_optional
    @test diagnostics.optional_rejected_count == size(S, 1) - k - expected_optional
    @test minimum(result.capture_eigenvalues[result.optional_kept_indices]) >= 0.99
    @test maximum(result.capture_eigenvalues[result.optional_rejected_indices]) < 0.99
    @test diagnostics.weakest_kept_optional_capture >= 0.99
    @test diagnostics.strongest_rejected_weak_capture < 0.99
    @test diagnostics.final_rank == k + expected_optional
    @test diagnostics.rejected_weak_directions_not_mwg_residual_channels
    @test !isempty(diagnostics.kept_direction_summary)
    @test !isempty(diagnostics.rejected_direction_summary)
    @test due.final_dimension == size(X, 1)

    @test_throws ArgumentError CRG.occupied_first_injection_geometry(
        size(X, 1), X, S, 1.1 .* Y)
    @test_throws DimensionMismatch CRG.occupied_first_injection_geometry(
        size(X, 1), X[:, 1:(end - 1)], S, Y)
    println("occupied_first_case=", case_label,
        " base_capture_singular_values=", diagnostics.occupied_base_capture_singular_values,
        " base_capture_min=", diagnostics.occupied_base_capture_min,
        " post_recovery_singular_values=",
        diagnostics.occupied_recovery_after_mandatory_inclusion_singular_values,
        " post_recovery_loss=",
        diagnostics.occupied_recovery_after_mandatory_inclusion_loss,
        " capture_tol=", diagnostics.capture_tol,
        " R_min=", diagnostics.complement_metric_minimum_eigenvalue,
        " full_range=", diagnostics.raw_full_capture_range,
        " complement_range=", diagnostics.raw_complement_capture_range,
        " optional_kept=", diagnostics.optional_kept_count,
        " optional_rejected=", diagnostics.optional_rejected_count)
    println("occupied_first_due_diligence_", case_label, "=", due)
    return (; result, due)
end

@testset "Occupied-first injection geometry" begin
    be = _be_full_reference_packet()
    ne = CRD.build_atomic_hf_reference_packet(
        CRD.ne_all_electron_reference_packet_spec(ns = 5))

    be_result = _run_packet_case(be, "Be", 15)
    ne_result = _run_packet_case(ne, "Ne", 12)

    @test size(be_result.result.mandatory_A, 2) == 2
    @test size(ne_result.result.mandatory_A, 2) == 5
end
