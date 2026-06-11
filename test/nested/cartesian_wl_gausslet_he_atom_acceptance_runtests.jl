# Runtime role: decomposed White-Lindsey He acceptance readiness audit.
#
# This file exercises the gausslet-only one-center Z = 2 decomposed WL path
# through the first restricted Hartree-Fock acceptance checkpoint. It must not
# fall back to full-parent CPBs, direct Cartesian product operators, or ordinary
# Cartesian IDA operators.

using Test
using LinearAlgebra
using GaussletBases

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const WLHeAcceptanceCPBM = GaussletBases.CartesianPairBlockMaterialization
const WLHeAcceptanceCPB = GaussletBases.CartesianCPB

function _wl_he_fixture_diagnostics(seed_report, axis_inputs)
    basis = seed_report.fixture.basis
    axis_centers = basis.center_data
    reference_centers = basis.reference_center_data
    coordinate_steps = diff(axis_centers)
    retained_shell_range_value = seed_report.inventory.retained_ranges.shell
    retained_shell_ranges =
        retained_shell_range_value isa Tuple ?
        retained_shell_range_value :
        (retained_shell_range_value,)
    return (;
        axis_count = length(axis_centers),
        source_side_count = seed_report.inventory.source_side_count,
        parent_side_count = seed_report.fixture.parent_side_count,
        explicit_mapping_supplied =
            hasproperty(seed_report.fixture, :explicit_mapping_supplied) ?
            seed_report.fixture.explicit_mapping_supplied :
            false,
        shell_layer_count = length(seed_report.fixture.sequence.shell_layers),
        core_side_count = seed_report.fixture.structure.core_side_count,
        reference_coordinate_endpoints =
            (first(reference_centers), last(reference_centers)),
        coordinate_endpoints = (first(axis_centers), last(axis_centers)),
        coordinate_minimum = minimum(axis_centers),
        coordinate_maximum = maximum(axis_centers),
        coordinate_tail_extent = maximum(abs, axis_centers),
        coordinate_spacing_minimum = minimum(coordinate_steps),
        coordinate_spacing_maximum = maximum(coordinate_steps),
        fixture_box_status =
            maximum(abs, axis_centers) < 1.0 ?
            :compact_z2_seed_box_under_one_bohr :
            :z2_seed_box_at_least_one_bohr,
        core_spacing = axis_inputs.core_spacing,
        mapping_a = axis_inputs.mapping_a,
        mapping_s = axis_inputs.mapping_s,
        tail_spacing = axis_inputs.tail_spacing,
        retained_core_range = seed_report.inventory.retained_ranges.core,
        retained_shell_ranges,
        retained_shell_range =
            length(retained_shell_ranges) == 1 ?
            only(retained_shell_ranges) :
            retained_shell_ranges,
    )
end

function _wl_he_parent_axis_inputs(seed_report, expansion)
    return _wl_he_parent_axis_inputs_from_basis(seed_report.fixture.basis, expansion)
end

function _wl_he_parent_axis_inputs_from_basis(basis, expansion)
    doside_source_1d = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    parent_axis_bundle_object = (;
        x = doside_source_1d,
        y = doside_source_1d,
        z = doside_source_1d,
    )
    pgdg = doside_source_1d.pgdg_intermediate
    mapping = basis.spec.mapping_value
    core_spacing = mapping.a * mapping.s
    return (;
        parent_axis_counts = ntuple(_ -> length(basis.center_data), 3),
        parent_axis_bundle_object,
        overlap_1d = (; x = pgdg.overlap, y = pgdg.overlap, z = pgdg.overlap),
        kinetic_1d = (; x = pgdg.kinetic, y = pgdg.kinetic, z = pgdg.kinetic),
        core_spacing,
        mapping_a = mapping.a,
        mapping_s = mapping.s,
        tail_spacing = mapping.tail_spacing,
        spacing_rule = :white_lindsey_atomic_mapping,
        spacing_rule_formula = :asinh_mapping_c_equals_d_s_equals_sqrt_dZ,
        spacing_rule_status = :standard_z_dependent_spacing,
        previous_shared_fixture_spacing_rule_status =
            :violated_he_z_dependent_spacing_rule,
    )
end

function _wl_he_center_records()
    return (
        (;
            center_key = :helium_nucleus,
            center_index = 1,
            nuclear_charge = 2.0,
            location = (0.0, 0.0, 0.0),
        ),
    )
end

function _wl_he_overlap_metric_diagnostics(overlap_matrix, decomposed_inventory)
    sym_s = Symmetric((overlap_matrix + transpose(overlap_matrix)) ./ 2)
    eigenvalues = eigvals(sym_s)
    minimum_eigenvalue = minimum(eigenvalues)
    maximum_eigenvalue = maximum(eigenvalues)
    return (;
        blocker = minimum_eigenvalue <= 0.0 ?
            :decomposed_wl_overlap_metric_not_positive_definite :
            nothing,
        minimum_eigenvalue,
        maximum_eigenvalue,
        condition_estimate = maximum_eigenvalue / minimum_eigenvalue,
        symmetry_error = maximum(abs.(overlap_matrix - transpose(overlap_matrix))),
        rank_estimate = count(eigenvalues .> 1.0e-10),
        retained_dimension = decomposed_inventory.retained_dimension,
    )
end

function _wl_he_lowest_one_electron_energy(
    hamiltonian_matrix,
    overlap_matrix,
    overlap_diagnostics,
)
    if !isnothing(overlap_diagnostics.blocker)
        return (;
            status = :blocked_decomposed_wl_he_one_electron_solve,
            blocker = overlap_diagnostics.blocker,
            energy = NaN,
            solve_kind = :blocked_overlap_metric,
        )
    end
    overlap_identity =
        isapprox(overlap_matrix, Matrix{Float64}(I, size(overlap_matrix)...);
            atol = 1.0e-9,
            rtol = 1.0e-9,
        )
    sym_h = Symmetric((hamiltonian_matrix + transpose(hamiltonian_matrix)) ./ 2)
    eig = overlap_identity ?
        eigen(sym_h) :
        eigen(
            sym_h,
            Symmetric((overlap_matrix + transpose(overlap_matrix)) ./ 2),
        )
    lowest_index = argmin(eig.values)
    return (;
        status = :materialized_decomposed_wl_he_one_electron_solve,
        blocker = nothing,
        energy = eig.values[lowest_index],
        lowest_orbital_coefficients = eig.vectors[:, lowest_index],
        solve_kind =
            overlap_identity ? :ordinary_symmetric : :generalized_symmetric,
        overlap_identity,
    )
end

function _wl_he_restricted_density_density_fock(one_body, interaction, density)
    h = 0.5 .* (Matrix{Float64}(one_body) .+ transpose(Matrix{Float64}(one_body)))
    v =
        0.5 .* (Matrix{Float64}(interaction) .+ transpose(Matrix{Float64}(interaction)))
    rho = 0.5 .* (Matrix{Float64}(density) .+ transpose(Matrix{Float64}(density)))
    occupations = vec(diag(rho))
    fock = h + 2.0 .* Diagonal(v * occupations) - rho .* v
    return 0.5 .* (fock .+ transpose(fock))
end

function _wl_he_restricted_density_density_energy(one_body, interaction, density)
    h = 0.5 .* (Matrix{Float64}(one_body) .+ transpose(Matrix{Float64}(one_body)))
    v =
        0.5 .* (Matrix{Float64}(interaction) .+ transpose(Matrix{Float64}(interaction)))
    rho = 0.5 .* (Matrix{Float64}(density) .+ transpose(Matrix{Float64}(density)))
    occupations = vec(diag(rho))
    direct = 2.0 * dot(occupations, v * occupations)
    exchange = dot(vec(rho), vec(v .* rho))
    return 2.0 * tr(rho * h) + direct - exchange
end

function _wl_he_scalar_density_density_convention_check()
    h = -1.25
    v = 0.375
    density = reshape([1.0], 1, 1)
    energy = _wl_he_restricted_density_density_energy(
        reshape([h], 1, 1),
        reshape([v], 1, 1),
        density,
    )
    return (;
        status =
            isapprox(energy, 2h + v; atol = 1.0e-14, rtol = 0.0) ?
            :passed_closed_shell_density_density_scalar_convention :
            :failed_closed_shell_density_density_scalar_convention,
        energy,
        expected = 2h + v,
        convention = :one_spatial_orbital_E_equals_2h_plus_V,
    )
end

function _wl_he_restricted_density_density_interaction_energy(
    interaction,
    density,
)
    v =
        0.5 .* (Matrix{Float64}(interaction) .+ transpose(Matrix{Float64}(interaction)))
    rho = 0.5 .* (Matrix{Float64}(density) .+ transpose(Matrix{Float64}(density)))
    occupations = vec(diag(rho))
    direct = 2.0 * dot(occupations, v * occupations)
    exchange = dot(vec(rho), vec(v .* rho))
    return direct - exchange
end

function _wl_he_core_unit(inventory)
    matches = Tuple(
        unit for unit in inventory.retained_units
        if unit.metadata.stratum_kind === :direct_core
    )
    length(matches) == 1 || return nothing
    return only(matches)
end

function _wl_he_density_diagnostics(density, interaction, inventory, axis_centers)
    rho = 0.5 .* (Matrix{Float64}(density) .+ transpose(Matrix{Float64}(density)))
    diagonal = vec(diag(rho))
    density_trace = sum(diagonal)
    peak_column = argmax(diagonal)
    core_unit = _wl_he_core_unit(inventory)
    direct_core_fraction = :unavailable
    boundary_fraction = :unavailable
    direct_core_rms_radius = :unavailable
    direct_core_peak_coordinate = :unavailable
    if !isnothing(core_unit)
        core_range = core_unit.column_range
        direct_core_weight = sum(diagonal[core_range])
        direct_core_fraction = direct_core_weight / density_trace
        boundary_fraction = 1.0 - direct_core_fraction
        source_cpb = only(core_unit.source_cpbs)
        intervals = WLHeAcceptanceCPB.intervals(source_cpb)
        coordinates = NTuple{3,Float64}[]
        radii2 = Float64[]
        for ix in intervals[1], iy in intervals[2], iz in intervals[3]
            coordinate = (
                axis_centers[ix],
                axis_centers[iy],
                axis_centers[iz],
            )
            push!(coordinates, coordinate)
            push!(
                radii2,
                coordinate[1]^2 + coordinate[2]^2 + coordinate[3]^2,
            )
        end
        if direct_core_weight > 0.0
            weights = diagonal[core_range]
            direct_core_rms_radius = sqrt(dot(weights, radii2) / direct_core_weight)
        end
        if peak_column in core_range
            direct_core_peak_coordinate =
                coordinates[peak_column - first(core_range) + 1]
        else
            direct_core_peak_coordinate = :outside_direct_core
        end
    end
    return (;
        density_trace,
        electron_count = 2.0 * density_trace,
        peak_density_diagonal = maximum(diagonal),
        peak_density_column = peak_column,
        direct_core_fraction,
        boundary_fraction,
        direct_core_rms_radius,
        direct_core_peak_coordinate,
        representative_positive_coulomb_energy =
            _wl_he_restricted_density_density_interaction_energy(
                interaction,
                rho,
            ),
    )
end

function _wl_he_restricted_hartree_fock(
    one_body,
    interaction;
    maxiter = 80,
    damping = 0.25,
    tol = 1.0e-10,
)
    h = 0.5 .* (Matrix{Float64}(one_body) .+ transpose(Matrix{Float64}(one_body)))
    v =
        0.5 .* (Matrix{Float64}(interaction) .+ transpose(Matrix{Float64}(interaction)))
    size(h) == size(v) ||
        throw(DimensionMismatch("He RHF one-body and interaction dimensions must match"))

    initial = eigen(Symmetric(h))
    occupied = initial.vectors[:, 1:1]
    density = occupied * transpose(occupied)
    density = 0.5 .* (density .+ transpose(density))
    energies = Float64[]
    residuals = Float64[]
    converged = false
    iterations = 0
    energy_change = Inf
    fock = h
    eig = initial

    for iteration in 1:maxiter
        iterations = iteration
        fock = _wl_he_restricted_density_density_fock(h, v, density)
        eig = eigen(Symmetric(fock))
        occupied_new = eig.vectors[:, 1:1]
        density_new = occupied_new * transpose(occupied_new)
        density_new = 0.5 .* (density_new .+ transpose(density_new))
        mixed_density = (1 - damping) .* density_new .+ damping .* density
        mixed_density = 0.5 .* (mixed_density .+ transpose(mixed_density))
        energy = _wl_he_restricted_density_density_energy(h, v, mixed_density)
        residual = norm(mixed_density - density, Inf)
        energy_change = isempty(energies) ? Inf : abs(energy - energies[end])
        push!(energies, energy)
        push!(residuals, residual)
        density = mixed_density
        if residual <= tol && energy_change <= tol
            converged = true
            break
        end
    end

    fock = _wl_he_restricted_density_density_fock(h, v, density)
    eig = eigen(Symmetric(fock))
    occupied = eig.vectors[:, 1:1]
    projector = occupied * transpose(occupied)
    projector = 0.5 .* (projector .+ transpose(projector))
    energy = _wl_he_restricted_density_density_energy(h, v, density)
    interaction_energy =
        _wl_he_restricted_density_density_interaction_energy(v, density)
    one_electron_energy = 2.0 * tr(density * h)
    return (;
        status =
            converged ?
            :converged_decomposed_wl_he_restricted_hartree_fock :
            :blocked_decomposed_wl_he_restricted_hartree_fock,
        blocker = converged ? nothing : :decomposed_wl_hartree_fock_not_converged,
        energy,
        one_electron_energy,
        electron_electron_energy = interaction_energy,
        fock,
        density,
        occupations = 2.0 .* vec(diag(density)),
        occupied_coefficients = occupied,
        orbital_energies = Vector{Float64}(eig.values),
        iterations,
        converged,
        residual = isempty(residuals) ? Inf : residuals[end],
        energy_change,
        solve_kind = :restricted_closed_shell_density_density_hartree_fock,
        spin_convention = :restricted_closed_shell_one_alpha_one_beta,
        spatial_occupation_count = 1,
        electron_count = sum(2.0 .* vec(diag(density))),
        projector_residual = norm(projector - density, Inf),
    )
end

function _wl_decomposed_he_atom_acceptance_audit(;
    parent_side_count::Int = 7,
    d::Real = 0.2,
    tail_spacing::Real = 10.0,
)
    expansion = coulomb_gaussian_expansion(doacc = false)
    seed_report =
        GaussletBases._white_lindsey_low_order_materialized_seed_report(
            parent_side_count = parent_side_count,
            Z = 2.0,
            d = d,
            tail_spacing = tail_spacing,
        )
    axis_inputs = _wl_he_parent_axis_inputs(seed_report, expansion)
    center_records = _wl_he_center_records()
    parent_axes = ntuple(_ -> seed_report.fixture.basis.center_data, 3)
    shellification_inventory_source =
        WLHeAcceptanceCPBM.white_lindsey_shellification_decomposed_unit_pair_inventory(
            parent_axes,
            ((0.0, 0.0, 0.0),);
            metadata = (; q = 5),
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
        )
    decomposed_inventory = shellification_inventory_source.inventory
    fixture_diagnostics = _wl_he_fixture_diagnostics(seed_report, axis_inputs)
    scalar_convention = _wl_he_scalar_density_density_convention_check()
    overlap_global =
        WLHeAcceptanceCPBM.route_global_decomposed_wl_overlap_matrix(
            decomposed_inventory;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object =
                axis_inputs.parent_axis_bundle_object,
            overlap_1d = axis_inputs.overlap_1d,
        )
    kinetic_global =
        WLHeAcceptanceCPBM.route_global_decomposed_wl_kinetic_matrix(
            decomposed_inventory;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object =
                axis_inputs.parent_axis_bundle_object,
            overlap_1d = axis_inputs.overlap_1d,
            kinetic_1d = axis_inputs.kinetic_1d,
        )
    by_center_global =
        WLHeAcceptanceCPBM.route_global_electron_nuclear_by_center_matrices(
            decomposed_inventory;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object =
                axis_inputs.parent_axis_bundle_object,
            coulomb_expansion = expansion,
            center_records,
        )
    hamiltonian =
        WLHeAcceptanceCPBM.white_lindsey_decomposed_one_electron_hamiltonian(
            kinetic_global,
            by_center_global,
        )
    overlap_diagnostics =
        _wl_he_overlap_metric_diagnostics(
            overlap_global.matrix,
            decomposed_inventory,
        )
    solve = _wl_he_lowest_one_electron_energy(
        hamiltonian.matrix,
        overlap_global.matrix,
        overlap_diagnostics,
    )
    density_density =
        WLHeAcceptanceCPBM.route_global_decomposed_wl_density_density_matrix(
            decomposed_inventory;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object =
                axis_inputs.parent_axis_bundle_object,
            coulomb_expansion = expansion,
        )
    if density_density.status ===
       :materialized_route_global_density_density_interaction_matrix
        hartree_fock = _wl_he_restricted_hartree_fock(
            hamiltonian.matrix,
            density_density.matrix,
        )
    else
        hartree_fock = (;
            status = :blocked_decomposed_wl_he_restricted_hartree_fock,
            blocker = :missing_decomposed_wl_density_density_interaction_route,
            energy = :unavailable,
            one_electron_energy = :unavailable,
            electron_electron_energy = :unavailable,
            iterations = 0,
            converged = false,
            residual = :unavailable,
            energy_change = :unavailable,
            solve_kind = :blocked_missing_density_density_interaction,
            spin_convention = :restricted_closed_shell_one_alpha_one_beta,
            spatial_occupation_count = 1,
            electron_count = :unavailable,
            projector_residual = :unavailable,
        )
    end
    if hartree_fock.status ===
       :converged_decomposed_wl_he_restricted_hartree_fock
        density_diagnostics = _wl_he_density_diagnostics(
            hartree_fock.density,
            density_density.matrix,
            decomposed_inventory,
            seed_report.fixture.basis.center_data,
        )
    else
        density_diagnostics = (;
            density_trace = :unavailable,
            electron_count = :unavailable,
            peak_density_diagonal = :unavailable,
            peak_density_column = :unavailable,
            direct_core_fraction = :unavailable,
            boundary_fraction = :unavailable,
            direct_core_rms_radius = :unavailable,
            direct_core_peak_coordinate = :unavailable,
            representative_positive_coulomb_energy = :unavailable,
        )
    end
    two_electron_route_available =
        density_density.status ===
        :materialized_route_global_density_density_interaction_matrix
    hf_available =
        hartree_fock.status ===
        :converged_decomposed_wl_he_restricted_hartree_fock
    one_electron_contribution = 2.0 * solve.energy
    lowest_density =
        solve.status === :materialized_decomposed_wl_he_one_electron_solve ?
        solve.lowest_orbital_coefficients *
        transpose(solve.lowest_orbital_coefficients) :
        zeros(Float64, 0, 0)
    hydrogenic_self_coulomb =
        two_electron_route_available ?
        _wl_he_restricted_density_density_interaction_energy(
            density_density.matrix,
            lowest_density,
        ) :
        :unavailable
    hydrogenic_self_coulomb_reference = 5.0 * 2.0 / 8.0
    return (;
        object_kind = :decomposed_wl_he_atom_acceptance_audit,
        acceptance_suite = :decomposed_wl_gausslet_two_electron_acceptance,
        acceptance_fixture = :he_atom,
        status =
            hf_available ?
            :accepted_decomposed_wl_he_atom_restricted_hartree_fock :
            :blocked_decomposed_wl_he_atom_two_electron_acceptance,
        blocker =
            two_electron_route_available ?
            hartree_fock.blocker :
            :missing_decomposed_wl_density_density_interaction_route,
        q = 5,
        ns = 5,
        center_count = length(center_records),
        center_keys = Tuple(record.center_key for record in center_records),
        nuclear_charges = Tuple(record.nuclear_charge for record in center_records),
        core_spacing = axis_inputs.core_spacing,
        parent_side_count = fixture_diagnostics.parent_side_count,
        coordinate_endpoints = fixture_diagnostics.coordinate_endpoints,
        coordinate_tail_extent = fixture_diagnostics.coordinate_tail_extent,
        parent_axis_counts = axis_inputs.parent_axis_counts,
        decomposed_unit_count = decomposed_inventory.unit_count,
        decomposed_unit_pair_count = decomposed_inventory.pair_count,
        shellification_backed_decomposed_wl_inventory =
            decomposed_inventory.source_kind ===
            :cartesian_shellification_retained_unit_pair_plan,
        low_order_materialized_seed_inventory_used =
            decomposed_inventory.source_kind ===
            :white_lindsey_low_order_materialized_seed_ranges,
        retained_dimension = decomposed_inventory.retained_dimension,
        route_global_overlap_status = overlap_global.status,
        route_global_kinetic_status = kinetic_global.status,
        route_global_by_center_status = by_center_global.status,
        route_global_by_center_center_count = by_center_global.center_count,
        route_global_by_center_nuclear_charge_applied =
            by_center_global.nuclear_charge_applied,
        route_global_by_center_centers_summed = by_center_global.centers_summed,
        hamiltonian_status = hamiltonian.status,
        nuclear_charge_application_stage =
            hamiltonian.charge_application_stage,
        centers_summed_at_hamiltonian_assembly =
            hamiltonian.centers_summed_at_hamiltonian_assembly,
        overlap_minimum_eigenvalue = overlap_diagnostics.minimum_eigenvalue,
        overlap_condition_estimate = overlap_diagnostics.condition_estimate,
        overlap_rank_estimate = overlap_diagnostics.rank_estimate,
        one_electron_solve_kind = solve.solve_kind,
        lowest_one_electron_orbital_energy = solve.energy,
        closed_shell_one_electron_energy_contribution = one_electron_contribution,
        electron_electron_status = density_density.status,
        electron_electron_blocker = density_density.blocker,
        electron_electron_matrix_shape = density_density.matrix_shape,
        electron_electron_matrix_dimension =
            density_density.retained_dimension,
        electron_electron_matrix_symmetry_error =
            two_electron_route_available ?
            maximum(
                abs.(
                    density_density.matrix -
                    transpose(density_density.matrix),
                ),
            ) :
            :unavailable,
        electron_electron_matrix_values_finite =
            two_electron_route_available ?
            all(isfinite, density_density.matrix) :
            false,
        source_weight_division_stage =
            density_density.source_weight_division_stage,
        raw_pair_numerator_contracted =
            density_density.raw_pair_numerator_contracted,
        final_retained_weight_division_applied =
            density_density.final_retained_weight_division_applied,
        axis_integral_weights_deferred =
            density_density.axis_integral_weights_deferred,
        ida_density_density_convention =
            density_density.ida_density_density_convention,
        hydrogenic_1s_self_coulomb_reference =
            hydrogenic_self_coulomb_reference,
        hydrogenic_1s_h1_eigenvalue_measured = solve.energy,
        hydrogenic_1s_h1_eigenvalue_reference = -2.0,
        hydrogenic_1s_h1_eigenvalue_error = abs(solve.energy - -2.0),
        hydrogenic_1s_self_coulomb_measured =
            hydrogenic_self_coulomb,
        hydrogenic_1s_self_coulomb_error =
            two_electron_route_available ?
            abs(hydrogenic_self_coulomb - hydrogenic_self_coulomb_reference) :
            :unavailable,
        hartree_fock_status = hartree_fock.status,
        hartree_fock_solve_kind = hartree_fock.solve_kind,
        hartree_fock_iterations = hartree_fock.iterations,
        hartree_fock_converged = hartree_fock.converged,
        hartree_fock_residual = hartree_fock.residual,
        hartree_fock_energy_change = hartree_fock.energy_change,
        occupancy_spin_convention = hartree_fock.spin_convention,
        spatial_occupation_count = hartree_fock.spatial_occupation_count,
        hartree_fock_electron_count = hartree_fock.electron_count,
        bare_closed_shell_one_electron_energy = one_electron_contribution,
        rhf_one_electron_energy = hartree_fock.one_electron_energy,
        rhf_electron_electron_energy = hartree_fock.electron_electron_energy,
        rhf_total_energy = hartree_fock.energy,
        electron_electron_contribution = hartree_fock.electron_electron_energy,
        total_he_energy = hartree_fock.energy,
        rhf_total_decomposition_error =
            hf_available ?
            abs(
                hartree_fock.energy -
                (
                    hartree_fock.one_electron_energy +
                    hartree_fock.electron_electron_energy
                ),
            ) :
            :unavailable,
        rhf_density_trace = density_diagnostics.density_trace,
        rhf_density_electron_count = density_diagnostics.electron_count,
        converged_density_positive_coulomb_energy =
            density_diagnostics.representative_positive_coulomb_energy,
        scalar_density_density_convention_status = scalar_convention.status,
        scalar_density_density_convention_energy = scalar_convention.energy,
        scalar_density_density_convention_expected = scalar_convention.expected,
        exact_he_reference_energy = -2.9037243770341196,
        decomposed_wl_units_consumed = true,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
    )
end

@testset "decomposed WL gausslet-only He atom acceptance readiness" begin
    report = _wl_decomposed_he_atom_acceptance_audit()

    @test report.status ==
          :accepted_decomposed_wl_he_atom_restricted_hartree_fock
    @test isnothing(report.blocker)
    @test report.acceptance_fixture == :he_atom
    @test report.q == 5
    @test report.ns == 5
    @test report.center_count == 1
    @test report.center_keys == (:helium_nucleus,)
    @test report.nuclear_charges == (2.0,)
    @test report.core_spacing ≈ 0.2 atol = 1.0e-14 rtol = 0.0
    @test report.parent_side_count == 7
    @test report.parent_axis_counts == (7, 7, 7)
    @test report.coordinate_endpoints[1] ≈ -report.coordinate_endpoints[2] atol =
          1.0e-12 rtol = 0.0
    @test report.coordinate_tail_extent ≈ 0.9666200087560217 atol =
          1.0e-12 rtol = 0.0
    @test report.retained_dimension == 223
    @test report.decomposed_unit_count == 27
    @test report.decomposed_unit_pair_count == 378
    @test report.shellification_backed_decomposed_wl_inventory
    @test !report.low_order_materialized_seed_inventory_used
    @test report.route_global_overlap_status ==
          :materialized_route_global_overlap_matrix
    @test report.route_global_kinetic_status ==
          :materialized_route_global_kinetic_matrix
    @test report.route_global_by_center_status ==
          :materialized_route_global_electron_nuclear_by_center_matrix_set
    @test !report.route_global_by_center_nuclear_charge_applied
    @test !report.route_global_by_center_centers_summed
    @test report.hamiltonian_status ==
          :materialized_white_lindsey_decomposed_one_electron_hamiltonian
    @test report.nuclear_charge_application_stage ==
          :white_lindsey_hamiltonian_assembly
    @test report.centers_summed_at_hamiltonian_assembly
    @test report.overlap_minimum_eigenvalue > 0.0
    @test report.overlap_condition_estimate < 1.0000001
    @test report.overlap_rank_estimate == report.retained_dimension
    @test report.one_electron_solve_kind == :ordinary_symmetric
    @test isfinite(report.lowest_one_electron_orbital_energy)
    @test report.lowest_one_electron_orbital_energy < -1.0
    @test isfinite(report.closed_shell_one_electron_energy_contribution)
    @test report.electron_electron_status ==
          :materialized_route_global_density_density_interaction_matrix
    @test isnothing(report.electron_electron_blocker)
    @test report.electron_electron_matrix_shape == (223, 223)
    @test report.electron_electron_matrix_dimension == 223
    @test report.electron_electron_matrix_symmetry_error < 1.0e-10
    @test report.electron_electron_matrix_values_finite
    @test report.raw_pair_numerator_contracted
    @test report.source_weight_division_stage ==
          :after_retained_raw_numerator_assembly
    @test report.final_retained_weight_division_applied
    @test !report.axis_integral_weights_deferred
    @test report.ida_density_density_convention ==
          :full_retained_two_index_interaction_matrix
    @test report.hydrogenic_1s_self_coulomb_reference ≈ 1.25 atol =
          1.0e-14 rtol = 0.0
    @test report.hydrogenic_1s_h1_eigenvalue_reference ≈ -2.0 atol =
          1.0e-14 rtol = 0.0
    @test isfinite(report.hydrogenic_1s_h1_eigenvalue_measured)
    @test report.hydrogenic_1s_h1_eigenvalue_measured ==
          report.lowest_one_electron_orbital_energy
    @test report.hydrogenic_1s_h1_eigenvalue_error < 0.2
    @test isfinite(report.hydrogenic_1s_self_coulomb_measured)
    @test report.hydrogenic_1s_self_coulomb_measured > 0.0
    @test report.hydrogenic_1s_self_coulomb_error < 0.5
    @test report.hartree_fock_status ==
          :converged_decomposed_wl_he_restricted_hartree_fock
    @test report.hartree_fock_solve_kind ==
          :restricted_closed_shell_density_density_hartree_fock
    @test report.occupancy_spin_convention ==
          :restricted_closed_shell_one_alpha_one_beta
    @test report.spatial_occupation_count == 1
    @test report.hartree_fock_converged
    @test report.hartree_fock_iterations <= 80
    @test report.hartree_fock_residual < 1.0e-8
    @test report.hartree_fock_electron_count ≈ 2.0 atol = 1.0e-8 rtol = 0.0
    @test isfinite(report.electron_electron_contribution)
    @test report.electron_electron_contribution > 0.0
    @test isfinite(report.total_he_energy)
    @test report.rhf_total_energy ≈
          report.rhf_one_electron_energy + report.rhf_electron_electron_energy atol =
          1.0e-10 rtol = 0.0
    @test report.rhf_total_energy ≈ -2.3944175346639884 atol =
          1.0e-10 rtol = 0.0
    @test report.rhf_electron_electron_energy ≈ 1.3106054775285387 atol =
          1.0e-10 rtol = 0.0
    @test report.rhf_total_decomposition_error < 1.0e-10
    @test report.rhf_density_trace ≈ 1.0 atol = 1.0e-8 rtol = 0.0
    @test report.rhf_density_electron_count ≈ 2.0 atol = 1.0e-8 rtol = 0.0
    @test report.converged_density_positive_coulomb_energy ≈
          report.rhf_electron_electron_energy atol = 1.0e-10 rtol = 0.0
    @test report.converged_density_positive_coulomb_energy > 0.0
    @test report.scalar_density_density_convention_status ==
          :passed_closed_shell_density_density_scalar_convention
    @test report.scalar_density_density_convention_energy ≈
          report.scalar_density_density_convention_expected atol = 1.0e-14 rtol = 0.0
    @test report.total_he_energy > report.exact_he_reference_energy
    @test report.total_he_energy < -2.0
    @test report.decomposed_wl_units_consumed
    @test !report.full_parent_window_cpb_used
    @test !report.direct_cartesian_product_assembly_used
    @test !report.ordinary_cartesian_ida_operators_used
end
