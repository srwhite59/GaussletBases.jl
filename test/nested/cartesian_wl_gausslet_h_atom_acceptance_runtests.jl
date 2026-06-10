# Runtime role: WL/Cartesian gausslet-only H atom acceptance path.
#
# This is a bounded scientific acceptance check for the existing gausslet-only
# one-electron Cartesian route. It assembles H1 = kinetic + nuclear attraction
# for one proton at the origin, solves against the carried final overlap, and
# records a compact report. It does not add GTO supplements, PQS retained
# transforms, CPB/GTO bundle consumption, route/global refactors, or broad
# metadata assertions.

using Test
using LinearAlgebra
using GaussletBases

function _wl_h_atom_acceptance_result(;
    fixture_label::Symbol,
    axis_count::Int,
    mapping,
    reference_spacing::Real,
    xmax::Real,
    spacing_parameter_name::Symbol,
    spacing_parameter,
    acceptance_upper_bound::Real,
    acceptance_distance_upper_bound::Real,
)
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = axis_count,
        mapping,
        reference_spacing,
    ))
    operator_build_elapsed_s = @elapsed operators = ordinary_cartesian_ida_operators(
        basis;
        expansion,
        Z = 1.0,
        backend = :pgdg_localized_experimental,
    )
    overlap = Hermitian(operators.overlap_3d)
    hamiltonian = Hermitian(operators.one_body_hamiltonian)
    solve_elapsed_s = @elapsed solve = eigen(hamiltonian, overlap)
    lowest_energy = Float64(first(solve.values))
    helper_energy = mapped_cartesian_hydrogen_energy(
        basis;
        expansion,
        Z = 1.0,
        backend = :pgdg_localized_experimental,
    )
    axis_centers = centers(basis)
    central_axis_spacing =
        minimum(diff(sort(Float64[Float64(center) for center in axis_centers])))
    return (;
        basis,
        operators,
        report = (;
            fixture_label,
            route = :wl_cartesian_gausslet_only,
            proton_count = 1,
            proton_location = (0.0, 0.0, 0.0),
            nuclear_charge = 1.0,
            backend = :pgdg_localized_experimental,
            gaussian_expansion_term_count = length(expansion),
            axis_count,
            parent_support_size = axis_count^3,
            basis_dimension = size(operators.one_body_hamiltonian, 1),
            spacing_parameter_name,
            spacing_parameter,
            reference_spacing = Float64(reference_spacing),
            requested_xmax = Float64(xmax),
            mapping_a = mapping.a,
            mapping_s = mapping.s,
            mapping_c = mapping.a * mapping.s,
            tail_spacing = mapping.tail_spacing,
            distortion_strength = mapping.s,
            central_axis_spacing,
            outer_half_width = Float64(xmax),
            axis_extent = (minimum(axis_centers), maximum(axis_centers)),
            solve_kind = :generalized,
            lowest_one_electron_energy = lowest_energy,
            distance_from_minus_half_hartree = lowest_energy + 0.5,
            acceptance_upper_bound = Float64(acceptance_upper_bound),
            acceptance_distance_upper_bound = Float64(acceptance_distance_upper_bound),
            public_helper_energy = helper_energy,
            operator_build_elapsed_s,
            solve_elapsed_s,
        ),
    )
end

@testset "WL gausslet-only H atom scientific acceptance" begin
    fixtures = (
        _wl_h_atom_acceptance_result(
            fixture_label = :small_fitted_strength,
            axis_count = 7,
            mapping = fit_asinh_mapping_for_strength(
                s = 0.5,
                npoints = 7,
                xmax = 6.0,
            ),
            reference_spacing = 1.0,
            xmax = 6.0,
            spacing_parameter_name = :reference_spacing,
            spacing_parameter = 1.0,
            acceptance_upper_bound = -0.45,
            acceptance_distance_upper_bound = 0.05,
        ),
        _wl_h_atom_acceptance_result(
            fixture_label = :coarse_distorted_wl_style,
            axis_count = 11,
            mapping = AsinhMapping(c = 0.25, s = 1.0, tail_spacing = 10.0),
            reference_spacing = 1.0,
            xmax = 8.0,
            spacing_parameter_name = :asinh_c_core_spacing,
            spacing_parameter = 0.25,
            acceptance_upper_bound = -0.49,
            acceptance_distance_upper_bound = 0.01,
        ),
    )

    for result in fixtures
        report = result.report
        println("WL gausslet-only H atom acceptance report: ", report)

        @test report.route == :wl_cartesian_gausslet_only
        @test report.proton_count == 1
        @test report.proton_location == (0.0, 0.0, 0.0)
        @test report.nuclear_charge == 1.0
        @test report.parent_support_size == report.axis_count^3
        @test report.basis_dimension == report.parent_support_size
        @test report.solve_kind == :generalized
        @test size(result.operators.overlap_3d) ==
            (report.basis_dimension, report.basis_dimension)
        @test size(result.operators.one_body_hamiltonian) ==
            (report.basis_dimension, report.basis_dimension)
        @test all(isfinite, result.operators.overlap_3d)
        @test all(isfinite, result.operators.one_body_hamiltonian)
        @test isapprox(
            result.operators.one_body_hamiltonian,
            transpose(result.operators.one_body_hamiltonian);
            atol = 1.0e-10,
            rtol = 1.0e-10,
        )
        @test minimum(eigvals(Hermitian(result.operators.overlap_3d))) > 0.0
        @test isfinite(report.lowest_one_electron_energy)
        @test report.lowest_one_electron_energy > -0.5
        @test report.lowest_one_electron_energy < report.acceptance_upper_bound
        @test report.distance_from_minus_half_hartree > 0.0
        @test report.distance_from_minus_half_hartree <
            report.acceptance_distance_upper_bound
        @test abs(report.lowest_one_electron_energy - report.public_helper_energy) < 1.0e-10
    end

    coarse_report = fixtures[2].report
    @test coarse_report.fixture_label == :coarse_distorted_wl_style
    @test coarse_report.spacing_parameter_name == :asinh_c_core_spacing
    @test coarse_report.spacing_parameter == 0.25
    @test coarse_report.distortion_strength == 1.0
    @test coarse_report.requested_xmax == 8.0
    @test coarse_report.axis_count == 11
    @test coarse_report.basis_dimension == 1331
end
