using LinearAlgebra
using Test

using GaussletBases

const CRD = GaussletBases.CartesianReferenceDensity
const CRG = GaussletBases.CartesianResidualGaussians

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

function _metric_complement(S, Y)
    values, vectors = eigen(Symmetric(0.5 .* (S .+ transpose(S))))
    U = vectors * Diagonal(1.0 ./ sqrt.(values)) * transpose(vectors)
    Q0 = U - Y * (transpose(Y) * S * U)
    qvalues, qvectors = eigen(Symmetric(0.5 .* (transpose(Q0) * S * Q0 .+
        transpose(transpose(Q0) * S * Q0))))
    keep = findall(>(1.0e-10), qvalues)
    return Q0 * qvectors[:, keep] * Diagonal(1.0 ./ sqrt.(qvalues[keep]))
end

function _controlled_mixed_overlap(S, Y; strong_optional = 2)
    Z = _metric_complement(S, Y)
    W = hcat(Y, Z)
    nA = size(S, 1)
    nG = nA + 4
    k = size(Y, 2)
    lambda = fill(0.35, nA)
    lambda[1:k] .= 1.0
    for offset in 1:min(strong_optional, nA - k)
        lambda[k + offset] = offset == 1 ? 0.9995 : 0.995
    end
    B = zeros(Float64, nG, nA)
    for column in 1:nA
        B[column, column] = sqrt(lambda[column])
    end
    return Matrix{Float64}(B * inv(W)), lambda
end

function _run_packet_case(packet, case_label)
    S = packet.overlap
    Y = packet.occupied_coefficients
    X, lambda = _controlled_mixed_overlap(S, Y)
    labels, owners, channels = _packet_metadata(packet)
    result = CRG.occupied_first_injection_geometry(
        size(X, 1),
        X,
        S,
        Y;
        optional_capture_cutoff = 0.99,
        labels,
        owners,
        channels,
        provenance = (; source = :atomic_hf_reference_packet, case_label),
    )
    k = size(Y, 2)
    @test result.diagnostics.provenance.case_label == case_label
    @test result.diagnostics.occupied_orthogonality_error < 1.0e-10
    @test result.occupied_recovery_error < 1.0e-10
    @test result.diagnostics.weakest_occupied_capture ≈ 1.0 atol = 1.0e-10
    @test size(result.mandatory_A) == (size(S, 1), k)
    @test size(result.optional_A, 2) == 2
    @test size(result.injected_A, 2) == k + 2
    @test result.injected_G ≈ X * result.injected_A atol = 1.0e-12
    @test length(result.full_capture_eigenvalues) == size(S, 1)
    @test length(result.capture_eigenvalues) == size(S, 1) - k
    @test result.diagnostics.optional_kept_count == 2
    @test result.diagnostics.optional_rejected_count == size(S, 1) - k - 2
    @test minimum(result.capture_eigenvalues[result.optional_kept_indices]) >= 0.99
    @test maximum(result.capture_eigenvalues[result.optional_rejected_indices]) < 0.99
    @test result.diagnostics.weakest_kept_optional_capture >= 0.99
    @test result.diagnostics.strongest_rejected_weak_capture < 0.99
    @test result.diagnostics.final_rank == k + 2
    @test result.diagnostics.rejected_weak_directions_not_mwg_residual_channels
    @test !isempty(result.diagnostics.kept_direction_summary)
    @test !isempty(result.diagnostics.rejected_direction_summary)
    @test lambda[k + 1] > result.diagnostics.strongest_rejected_weak_capture

    @test_throws ArgumentError CRG.occupied_first_injection_geometry(
        size(X, 1), X, S, 1.1 .* Y)
    @test_throws DimensionMismatch CRG.occupied_first_injection_geometry(
        size(X, 1), X[:, 1:(end - 1)], S, Y)
    return result
end

function _indefinite_metric_case()
    S = Diagonal([1.0, 1.0, -1.0])
    Y = reshape([1.0, 0.0, 0.0], 3, 1)
    X = [1.0 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.25; 0.0 0.0 0.0]
    return @test_throws ArgumentError CRG.occupied_first_injection_geometry(4, X, S, Y)
end

@testset "Occupied-first injection geometry" begin
    be = _be_full_reference_packet()
    ne = CRD.build_atomic_hf_reference_packet(
        CRD.ne_all_electron_reference_packet_spec(ns = 5))

    be_result = _run_packet_case(be, "Be")
    ne_result = _run_packet_case(ne, "Ne")

    @test size(be_result.mandatory_A, 2) == 2
    @test size(ne_result.mandatory_A, 2) == 5
    _indefinite_metric_case()
end
