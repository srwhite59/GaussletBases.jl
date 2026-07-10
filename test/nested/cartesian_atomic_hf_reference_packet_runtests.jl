using LinearAlgebra
using Test
using JLD2

using GaussletBases

const CRD = GaussletBases.CartesianReferenceDensity

function _packet_with_convergence(packet, converged::Bool)
    diagnostics = merge(packet.rhf_diagnostics, (; converged))
    return typeof(packet)(packet.spec, packet.supplement, packet.overlap,
        packet.overlap_fingerprint, packet.occupied_coefficients,
        packet.occupations, packet.orbital_energies, packet.density_matrix,
        diagnostics, packet.density_fit, packet.potential_fit,
        packet.validation, packet.provenance)
end

function _packet_roundtrip_smoke(spec, label)
    packet = CRD.build_atomic_hf_reference_packet(spec)
    accepted_potential = CRD.fit_atomic_reference_potential(packet.density_fit)
    polish = packet.potential_fit.row.moment_polish
    path = joinpath(mktempdir(), "$(label)_atomic_hf_reference_packet.jld2")
    CRD.write_atomic_hf_reference_packet(path, packet)
    readback = CRD.read_atomic_hf_reference_packet(path)
    validation = CRD.validate_atomic_hf_reference_packet(path)
    in_memory_validation = CRD.validate_atomic_hf_reference_packet(packet)
    p0 = CRD.atomic_reference_packet_p0_q0(readback)

    @test readback.artifact_kind == :atomic_hf_reference_density_fit
    @test readback.convention_id == :atomic_hf_reference_density_fit_v1
    @test readback.electron_count == spec.electron_count
    @test packet.rhf_diagnostics.converged === true
    @test readback.rhf_converged === true
    @test validation.rhf_converged === true
    @test in_memory_validation.rhf_converged === true
    @test size(readback.C_occ, 2) == div(spec.electron_count, 2)
    @test validation.occupied_orthogonality_error < 1.0e-10
    @test abs(validation.density_trace_error) < 1.0e-10
    @test validation.density_matrix_error < 1.0e-12
    @test abs(validation.density_fit_charge_error) < 1.0e-10
    @test validation.density_fit_self_energy_relative_error < 1.0e-8
    @test validation.potential_fit_radial_relmax < 1.0e-4
    @test length(packet.potential_fit.coefficients) == 33
    @test packet.potential_fit.exponents == accepted_potential.exponents
    @test packet.potential_fit.coefficients[1:5] == accepted_potential.coefficients[1:5]
    @test polish.policy_id == :determinant_densityfit_coulomb_moment_v1
    @test polish.retained_rank == 28
    @test polish.moment_max_abs_error <= 1.0e-9
    @test packet.potential_fit.row.absmax <= accepted_potential.row.absmax + 1.0e-9
    @test packet.potential_fit.row.tail_charge_error <=
        accepted_potential.row.tail_charge_error + 1.0e-9
    @test readback.potential_fit.row.moment_polish == polish
    @test packet.rhf_diagnostics.coulomb_expansion_doacc === true
    @test packet.rhf_diagnostics.coulomb_expansion_terms >= 100
    @test packet.rhf_diagnostics.coulomb_expansion_maxu >= 100.0
    @test readback.rhf_coulomb_expansion_doacc === true
    @test readback.rhf_coulomb_expansion_terms == packet.rhf_diagnostics.coulomb_expansion_terms
    @test readback.rhf_coulomb_expansion_maxu == packet.rhf_diagnostics.coulomb_expansion_maxu
    @test readback.coulomb_expansions.rhf.policy === :high
    @test readback.coulomb_expansions.rhf.term_count == 135
    @test readback.coulomb_expansions.density_self_energy.policy === :compact
    @test readback.coulomb_expansions.density_self_energy.term_count == 45
    @test readback.coulomb_expansions.potential_tail_scaffold.policy === :compact
    @test readback.coulomb_expansions.potential_tail_scaffold.term_count == 45
    legacy_path = joinpath(dirname(path), "$(label)_legacy_packet.jld2")
    cp(path, legacy_path; force = true)
    JLD2.jldopen(legacy_path, "r+") do file
        for role in (:rhf, :density_self_energy, :potential_tail_scaffold),
                name in GaussletBases._CARTESIAN_COULOMB_EXPANSION_SUMMARY_KEYS
            delete!(file, "coulomb_expansion/$(role)/$(name)")
        end
        for name in (:policy_id, :retained_rank, :coefficient_delta_max,
                :moment_max_abs_error)
            delete!(file, "potential_fit/moment_polish/$(name)")
        end
    end
    legacy = CRD.read_atomic_hf_reference_packet(legacy_path)
    @test isnothing(legacy.coulomb_expansions.rhf)
    @test isnothing(legacy.coulomb_expansions.density_self_energy)
    @test isnothing(legacy.coulomb_expansions.potential_tail_scaffold)
    @test isnothing(legacy.potential_fit.row.moment_polish)
    @test abs(p0.trace - spec.electron_count) < 1.0e-10
    @test sum(p0.q_AA) > 0.0

    unconverged = _packet_with_convergence(packet, false)
    @test CRD.validate_atomic_hf_reference_packet(unconverged).rhf_converged === false
    @test_throws ArgumentError CRD.write_atomic_hf_reference_packet(
        joinpath(dirname(path), "$(label)_unconverged_packet.jld2"), unconverged)
    @test_throws ArgumentError CRD._require_atomic_reference_converged(
        (; converged = false), "atomic reference packet construction")

    altered_path = joinpath(dirname(path), "$(label)_altered_unconverged_packet.jld2")
    cp(path, altered_path; force = true)
    JLD2.jldopen(altered_path, "r+") do file
        delete!(file, "hf/converged")
        file["hf/converged"] = false
    end
    @test CRD.validate_atomic_hf_reference_packet(altered_path).rhf_converged === false
    return readback
end

@testset "Atomic HF reference packet roundtrip" begin
    _packet_roundtrip_smoke(CRD.be_core_reference_packet_spec(), "be_core")
    _packet_roundtrip_smoke(
        CRD.ne_all_electron_reference_packet_spec(), "ne_all_electron")
end
