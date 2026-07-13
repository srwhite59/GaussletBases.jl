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

function _packet_with_potential_row(packet, row)
    fit0 = packet.potential_fit
    fit = typeof(fit0)(fit0.coefficients, fit0.exponents, fit0.radial_grid,
        fit0.radial_exact, fit0.radial_fit, fit0.radial_error, row)
    return typeof(packet)(packet.spec, packet.supplement, packet.overlap,
        packet.overlap_fingerprint, packet.occupied_coefficients,
        packet.occupations, packet.orbital_energies, packet.density_matrix,
        packet.rhf_diagnostics, packet.density_fit, fit,
        packet.validation, packet.provenance)
end

function _packet_without_potential_row_field(packet, field)
    names = Tuple(filter(!=(field), propertynames(packet.potential_fit.row)))
    row = NamedTuple{names}(Tuple(getproperty(packet.potential_fit.row, name)
        for name in names))
    return _packet_with_potential_row(packet, row)
end

function _packet_with_potential_arrays(packet, coefficients, exponents)
    fit0 = packet.potential_fit
    fit = typeof(fit0)(coefficients, exponents, fit0.radial_grid,
        fit0.radial_exact, fit0.radial_fit, fit0.radial_error, fit0.row)
    return typeof(packet)(packet.spec, packet.supplement, packet.overlap,
        packet.overlap_fingerprint, packet.occupied_coefficients,
        packet.occupations, packet.orbital_energies, packet.density_matrix,
        packet.rhf_diagnostics, packet.density_fit, fit,
        packet.validation, packet.provenance)
end

function _packet_roundtrip_smoke(spec, label)
    packet = CRD.build_atomic_hf_reference_packet(spec)
    accepted_potential = CRD.fit_atomic_reference_potential(packet.density_fit)
    path = joinpath(mktempdir(), "$(label)_atomic_hf_reference_packet.jld2")
    CRD.write_atomic_hf_reference_packet(path, packet)
    readback = CRD.read_atomic_hf_reference_packet(path)
    validation = CRD.validate_atomic_hf_reference_packet(path)
    in_memory_validation = CRD.validate_atomic_hf_reference_packet(packet)
    p0 = CRD.atomic_reference_packet_p0_q0(readback)
    p0_in_memory = CRD.atomic_reference_packet_p0_q0(packet)

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
    @test packet.potential_fit.coefficients == accepted_potential.coefficients
    @test packet.potential_fit.row.source_coulomb_terms == 45
    @test packet.potential_fit.row.total_terms == 33
    @test packet.potential_fit.row.absmax == accepted_potential.row.absmax
    @test packet.potential_fit.row.tail_charge_error ==
        accepted_potential.row.tail_charge_error
    @test isfinite(packet.potential_fit.row.consistency_error)
    @test validation.potential_fit_consistency_error ==
        packet.potential_fit.row.consistency_error
    @test readback.potential_fit.row.consistency_error ==
        packet.potential_fit.row.consistency_error
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
    end
    legacy = CRD.read_atomic_hf_reference_packet(legacy_path)
    @test isnothing(legacy.coulomb_expansions.rhf)
    @test isnothing(legacy.coulomb_expansions.density_self_energy)
    @test isnothing(legacy.coulomb_expansions.potential_tail_scaffold)
    @test abs(p0.trace - spec.electron_count) < 1.0e-10
    @test p0_in_memory.P_AA == p0.P_AA
    @test sum(p0.q_AA) > 0.0

    retired = _packet_with_potential_row(packet,
        merge(packet.potential_fit.row, (;
            moment_polish = (; policy_id = :retired_test,))))
    @test_throws ArgumentError CRD.write_atomic_hf_reference_packet(
        joinpath(dirname(path), "$(label)_retired_polish_packet.jld2"), retired)
    polished_path = joinpath(dirname(path), "$(label)_file_polished_packet.jld2")
    cp(path, polished_path; force = true)
    JLD2.jldopen(polished_path, "r+") do file
        file["potential_fit/moment_polish/policy_id"] = :retired_test
    end
    @test_throws ArgumentError CRD.read_atomic_hf_reference_packet(polished_path)
    incomplete_path = joinpath(dirname(path), "$(label)_incomplete_ordinary_packet.jld2")
    cp(path, incomplete_path; force = true)
    JLD2.jldopen(incomplete_path, "r+") do file
        delete!(file, "potential_fit/consistency_error")
    end
    incomplete_error = try
        CRD.read_atomic_hf_reference_packet(incomplete_path)
        nothing
    catch error
        error
    end
    @test incomplete_error isa ArgumentError
    @test occursin("regenerate", sprint(showerror, incomplete_error))

    unconverged = _packet_with_convergence(packet, false)
    @test CRD.validate_atomic_hf_reference_packet(unconverged).rhf_converged === false
    unconverged_path = joinpath(dirname(path), "$(label)_unconverged_packet.jld2")
    @test_throws ArgumentError CRD.write_atomic_hf_reference_packet(
        unconverged_path, unconverged)
    @test !isfile(unconverged_path)
    @test_throws ArgumentError CRD.atomic_reference_packet_p0_q0(unconverged)
    @test_throws ArgumentError CRD._require_atomic_reference_converged(
        (; converged = false), "atomic reference packet construction")

    if label == "be_core"
        missing = _packet_without_potential_row_field(packet, :retained_rank)
        missing_path = joinpath(dirname(path), "missing_ordinary_field.jld2")
        @test_throws ArgumentError CRD.write_atomic_hf_reference_packet(
            missing_path, missing)
        @test !isfile(missing_path)

        nonfinite = _packet_with_potential_row(packet,
            merge(packet.potential_fit.row, (; core_relmax = NaN)))
        nonfinite_path = joinpath(dirname(path), "nonfinite_ordinary_field.jld2")
        @test_throws ArgumentError CRD.write_atomic_hf_reference_packet(
            nonfinite_path, nonfinite)
        @test !isfile(nonfinite_path)

        mismatched = _packet_with_potential_arrays(packet,
            packet.potential_fit.coefficients[1:(end - 1)],
            packet.potential_fit.exponents)
        mismatch_path = joinpath(dirname(path), "mismatched_potential_arrays.jld2")
        @test_throws ArgumentError CRD.write_atomic_hf_reference_packet(
            mismatch_path, mismatched)
        @test !isfile(mismatch_path)

        inconsistent = _packet_with_potential_row(packet,
            merge(packet.potential_fit.row, (;
                determinant_field_expectation =
                    packet.potential_fit.row.determinant_field_expectation + 1.0e-4)))
        inconsistent_path = joinpath(dirname(path), "inconsistent_ordinary_energy.jld2")
        @test_throws ArgumentError CRD.write_atomic_hf_reference_packet(
            inconsistent_path, inconsistent)
        @test !isfile(inconsistent_path)

        invalid_rank = _packet_with_potential_row(packet,
            merge(packet.potential_fit.row, (;
                retained_rank = packet.potential_fit.row.source_coulomb_terms)))
        invalid_rank_path = joinpath(dirname(path), "invalid_potential_rank.jld2")
        @test_throws ArgumentError CRD.validate_atomic_hf_reference_packet(invalid_rank)
        @test_throws ArgumentError CRD.write_atomic_hf_reference_packet(
            invalid_rank_path, invalid_rank)
        @test !isfile(invalid_rank_path)
    end

    altered_path = joinpath(dirname(path), "$(label)_altered_unconverged_packet.jld2")
    cp(path, altered_path; force = true)
    JLD2.jldopen(altered_path, "r+") do file
        delete!(file, "hf/converged")
        file["hf/converged"] = false
    end
    @test CRD.validate_atomic_hf_reference_packet(altered_path).rhf_converged === false
    unconverged_readback = CRD.read_atomic_hf_reference_packet(altered_path)
    @test_throws ArgumentError CRD.atomic_reference_packet_p0_q0(unconverged_readback)
    return readback
end

@testset "Atomic HF reference packet roundtrip" begin
    _packet_roundtrip_smoke(CRD.be_core_reference_packet_spec(), "be_core")
    _packet_roundtrip_smoke(
        CRD.ne_all_electron_reference_packet_spec(), "ne_all_electron")
end
