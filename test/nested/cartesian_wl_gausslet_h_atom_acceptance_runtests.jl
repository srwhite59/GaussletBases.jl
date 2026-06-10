# Runtime role: first WL/Cartesian gausslet-only H atom acceptance path.
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

function _wl_h_atom_acceptance_result()
    axis_count = 7
    mapping_strength = 0.5
    outer_half_width = 6.0
    reference_spacing = 1.0
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = axis_count,
        mapping = fit_asinh_mapping_for_strength(
            s = mapping_strength,
            npoints = axis_count,
            xmax = outer_half_width,
        ),
        reference_spacing,
    ))
    operators = ordinary_cartesian_ida_operators(
        basis;
        expansion,
        Z = 1.0,
        backend = :pgdg_localized_experimental,
    )
    overlap = Hermitian(operators.overlap_3d)
    hamiltonian = Hermitian(operators.one_body_hamiltonian)
    solve = eigen(hamiltonian, overlap)
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
            route = :wl_cartesian_gausslet_only,
            proton_count = 1,
            proton_location = (0.0, 0.0, 0.0),
            nuclear_charge = 1.0,
            backend = :pgdg_localized_experimental,
            gaussian_expansion_term_count = length(expansion),
            axis_count,
            parent_support_size = axis_count^3,
            basis_dimension = size(operators.one_body_hamiltonian, 1),
            reference_spacing,
            mapping_strength,
            central_axis_spacing,
            outer_half_width,
            axis_extent = (minimum(axis_centers), maximum(axis_centers)),
            solve_kind = :generalized,
            lowest_one_electron_energy = lowest_energy,
            distance_from_minus_half_hartree = lowest_energy + 0.5,
            public_helper_energy = helper_energy,
        ),
    )
end

@testset "WL gausslet-only H atom scientific acceptance" begin
    result = _wl_h_atom_acceptance_result()
    report = result.report
    println("WL gausslet-only H atom acceptance report: ", report)

    @test report.route == :wl_cartesian_gausslet_only
    @test report.proton_count == 1
    @test report.proton_location == (0.0, 0.0, 0.0)
    @test report.nuclear_charge == 1.0
    @test report.axis_count == 7
    @test report.parent_support_size == 343
    @test report.basis_dimension == 343
    @test report.solve_kind == :generalized
    @test size(result.operators.overlap_3d) == (343, 343)
    @test size(result.operators.one_body_hamiltonian) == (343, 343)
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
    @test report.lowest_one_electron_energy < -0.45
    @test report.distance_from_minus_half_hartree > 0.0
    @test report.distance_from_minus_half_hartree < 0.05
    @test abs(report.lowest_one_electron_energy - report.public_helper_energy) < 1.0e-10
end
