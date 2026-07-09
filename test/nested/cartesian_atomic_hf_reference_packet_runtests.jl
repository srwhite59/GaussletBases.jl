using LinearAlgebra
using Test

using GaussletBases

const CRD = GaussletBases.CartesianReferenceDensity

function _packet_roundtrip_smoke(spec, label)
    packet = CRD.build_atomic_hf_reference_packet(spec)
    path = joinpath(mktempdir(), "$(label)_atomic_hf_reference_packet.jld2")
    CRD.write_atomic_hf_reference_packet(path, packet)
    readback = CRD.read_atomic_hf_reference_packet(path)
    validation = CRD.validate_atomic_hf_reference_packet(path)
    p0 = CRD.atomic_reference_packet_p0_q0(readback)

    @test readback.artifact_kind == :atomic_hf_reference_density_fit
    @test readback.convention_id == :atomic_hf_reference_density_fit_v1
    @test readback.electron_count == spec.electron_count
    @test size(readback.C_occ, 2) == div(spec.electron_count, 2)
    @test validation.occupied_orthogonality_error < 1.0e-10
    @test abs(validation.density_trace_error) < 1.0e-10
    @test validation.density_matrix_error < 1.0e-12
    @test abs(validation.density_fit_charge_error) < 1.0e-10
    @test validation.density_fit_self_energy_relative_error < 1.0e-8
    @test validation.potential_fit_radial_relmax < 1.0e-4
    @test packet.rhf_diagnostics.coulomb_expansion_doacc === true
    @test packet.rhf_diagnostics.coulomb_expansion_terms >= 100
    @test packet.rhf_diagnostics.coulomb_expansion_maxu >= 100.0
    @test readback.rhf_coulomb_expansion_doacc === true
    @test readback.rhf_coulomb_expansion_terms == packet.rhf_diagnostics.coulomb_expansion_terms
    @test readback.rhf_coulomb_expansion_maxu == packet.rhf_diagnostics.coulomb_expansion_maxu
    @test abs(p0.trace - spec.electron_count) < 1.0e-10
    @test sum(p0.q_AA) > 0.0
    return readback
end

@testset "Atomic HF reference packet roundtrip" begin
    _packet_roundtrip_smoke(CRD.be_core_reference_packet_spec(), "be_core")
    _packet_roundtrip_smoke(
        CRD.ne_all_electron_reference_packet_spec(), "ne_all_electron")
end
