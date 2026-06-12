using Test
using LinearAlgebra
using GaussletBases

const PQSH1CPBM = GaussletBases.CartesianPairBlockMaterialization
const PQSH1CPOP = GaussletBases.CartesianPairOperatorPlans
const PQSH1CUP = GaussletBases.CartesianUnitPairs
const PQSH1CRTC = GaussletBases.CartesianRetainedUnitTransformContracts
const PQSH1CRU = GaussletBases.CartesianRetainedUnits
const PQSH1CTL = GaussletBases.CartesianTerminalLowering
const PQSH1CPB = GaussletBases.CartesianCPB
const PQSH1CRPS = GaussletBases.CartesianRawProductSources
const PQSH1CCPM = GaussletBases.CartesianContractedParentMetrics

function _pqs_h1_test_bundle(count::Int)
    xmax = 8.0
    tail = 10.0
    endpoint = (count - 1) / 2
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count,
            mapping = AsinhMapping(
                a = 0.25,
                s = asinh(xmax / 0.25) / (endpoint - xmax / tail),
                tail_spacing = tail,
            ),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
end

function _pqs_h1_axis_metrics(bundles)
    pgdg_x = GaussletBases._nested_axis_pgdg(bundles, :x)
    pgdg_y = GaussletBases._nested_axis_pgdg(bundles, :y)
    pgdg_z = GaussletBases._nested_axis_pgdg(bundles, :z)
    return (;
        x = (overlap = pgdg_x.overlap, kinetic = pgdg_x.kinetic),
        y = (overlap = pgdg_y.overlap, kinetic = pgdg_y.kinetic),
        z = (overlap = pgdg_z.overlap, kinetic = pgdg_z.kinetic),
    )
end

function _pqs_h1_minimal_lowering_plan()
    return PQSH1CTL.TerminalLoweringPlan(
        PQSH1CTL.PQSLowering(q = 5),
        (),
        (),
        (;
            object_kind = :pqs_direct_retained_final_h1_lowering_summary,
            status = :available_terminal_lowering_plan,
            policy_kind = :pqs_direct_retained_final_h1,
            terminal_region_count = 0,
            available_contract_count = 0,
            selected_contract_count = 0,
            selected_contract_kinds = (),
            all_terminal_regions_have_selected_contract = true,
            materialized = false,
            retained_spaces_materialized = false,
            final_retained_units_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :pqs_direct_retained_final_h1),
    )
end

function _pqs_h1_retained_unit(unit_key::Symbol, source_cpb, raw_plan)
    return PQSH1CRU.RetainedUnitRecord(
        unit_key,
        1,
        :pqs_shell_retained_unit,
        Symbol(unit_key, "_contract"),
        Symbol(unit_key, "_terminal_region"),
        :synthetic_terminal_region,
        :synthetic_terminal_region,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
        nothing,
        (source_cpb,),
        nothing,
        :not_materialized,
        nothing,
        :not_materialized,
        nothing,
        nothing,
        false,
        (;
            q = 5,
            source_mode_shape = raw_plan.source_mode_dims,
            raw_product_source_axis_transform_facts = raw_plan.axis_transform_facts,
            route_core_sidecar_status = :not_materialized,
            transform_source = :descriptor_axis_local_coefficients,
        ),
    )
end

function _pqs_h1_self_pair_materialization_record(source_cpb, raw_plan)
    unit = _pqs_h1_retained_unit(:pqs_direct_retained_final_h1_unit, source_cpb, raw_plan)
    retained_plan = PQSH1CRU.RetainedUnitPlan(
        PQSH1CRU.MetadataOnlyRetainedUnits(),
        _pqs_h1_minimal_lowering_plan(),
        (unit,),
        (;
            object_kind = :pqs_direct_retained_final_h1_retained_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = 1,
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :pqs_direct_retained_final_h1),
    )
    unit_pair_plan = PQSH1CUP.unit_pair_plan(retained_plan)
    transform_plan = PQSH1CRTC.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan = PQSH1CPOP.pair_operator_plan(
        unit_pair_plan,
        transform_plan;
        route_core_sidecars = false,
    )
    materialization_plan = PQSH1CPBM.pair_block_materialization_plan(pair_operator_plan)
    return only(PQSH1CPBM.pair_block_materialization_records(materialization_plan))
end

function _pqs_h1_axis_projected_operator(axis_operator, fact, interval)
    support_operator = axis_operator[interval, interval]
    coefficients = fact.coefficient_matrix
    return transpose(coefficients) * support_operator * coefficients
end

function _pqs_h1_support_product_matrix(states, mx, my, mz)
    result = Matrix{Float64}(undef, length(states), length(states))
    for column in eachindex(states), row in eachindex(states)
        ix, iy, iz = states[row]
        jx, jy, jz = states[column]
        result[row, column] = mx[ix, jx] * my[iy, jy] * mz[iz, jz]
    end
    return result
end

function _pqs_h1_shell_support_kinetic(states, metrics)
    return (
        _pqs_h1_support_product_matrix(
            states,
            metrics.x.kinetic,
            metrics.y.overlap,
            metrics.z.overlap,
        ) +
        _pqs_h1_support_product_matrix(
            states,
            metrics.x.overlap,
            metrics.y.kinetic,
            metrics.z.overlap,
        ) +
        _pqs_h1_support_product_matrix(
            states,
            metrics.x.overlap,
            metrics.y.overlap,
            metrics.z.kinetic,
        )
    )
end

function _pqs_h1_shell_support_nuclear(states, axis_layers, expansion, center)
    x_matrices = gaussian_factor_matrices(
        axis_layers.x;
        exponents = expansion.exponents,
        center = center.location[1],
    )
    y_matrices = gaussian_factor_matrices(
        axis_layers.y;
        exponents = expansion.exponents,
        center = center.location[2],
    )
    z_matrices = gaussian_factor_matrices(
        axis_layers.z;
        exponents = expansion.exponents,
        center = center.location[3],
    )
    result = zeros(Float64, length(states), length(states))
    for term_index in eachindex(expansion.coefficients)
        result .+=
            -Float64(expansion.coefficients[term_index]) *
            _pqs_h1_support_product_matrix(
                states,
                x_matrices[term_index],
                y_matrices[term_index],
                z_matrices[term_index],
            )
    end
    return result
end

function _pqs_h1_boundary_matrix_result(term::Symbol, block)
    return (;
        object_kind = :pqs_retained_source_one_body_matrix,
        status = :materialized_pqs_retained_source_one_body_matrix,
        blocker = nothing,
        term,
        matrix = block,
        matrix_space = :retained_pqs_source_modes,
        retained_dimension = size(block, 1),
        matrix_materialized = true,
        retained_source_matrix_materialized = true,
    )
end

@testset "PQS direct-retained final H1 gate" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundle5 = _pqs_h1_test_bundle(5)
    bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle5, bundle5, bundle5)
    layer = GaussletBases._nested_projected_q_shell_layer(
        bundles,
        (1:5, 1:5, 1:5),
        (2:4, 2:4, 2:4);
        bond_axis = :z,
        q = 5,
        L = 5,
        term_coefficients = Float64.(expansion.coefficients),
    )
    descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(layer)
    metrics = _pqs_h1_axis_metrics(bundles)
    shell_plan = PQSH1CCPM._pqs_shell_realization_plan(descriptor, metrics)
    shell_overlap = PQSH1CCPM._pqs_product_box_support_overlap_matrix(
        descriptor.support_states,
        metrics,
    )
    source_cpb = PQSH1CPB.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :pqs_direct_retained_final_h1_source,
    )
    raw_plan = PQSH1CRPS.raw_product_box_plan(
        source_cpb;
        source_key = :pqs_direct_retained_final_h1_source,
        source_mode_dims = (5, 5, 5),
        axis_transform_matrices = descriptor.axis_local_coefficients,
    )
    retained_rule = PQSH1CRPS.pqs_boundary_product_mode_retained_rule(raw_plan)
    final_basis = PQSH1CPBM.pqs_source_shell_realization_final_basis(
        raw_plan,
        retained_rule;
        shell_support_indices = descriptor.support_indices,
        shell_overlap,
        shell_projection = shell_plan.shell_projection_matrix,
        lowdin_cleanup = shell_plan.lowdin_cleanup,
    )
    record = _pqs_h1_self_pair_materialization_record(source_cpb, raw_plan)
    pgdg = bundle5.pgdg_intermediate
    overlap_1d = _pqs_h1_axis_projected_operator(
        pgdg.overlap,
        raw_plan.axis_transform_facts[1],
        1:5,
    )
    kinetic_1d = _pqs_h1_axis_projected_operator(
        pgdg.kinetic,
        raw_plan.axis_transform_facts[1],
        1:5,
    )
    retained_overlap_block = PQSH1CPBM.pqs_source_pair_retained_overlap_block(
        record;
        overlap_1d = (; x = overlap_1d, y = overlap_1d, z = overlap_1d),
    )
    retained_kinetic_block = PQSH1CPBM.pqs_source_pair_retained_kinetic_block(
        record;
        overlap_1d = (; x = overlap_1d, y = overlap_1d, z = overlap_1d),
        kinetic_1d = (; x = kinetic_1d, y = kinetic_1d, z = kinetic_1d),
    )
    center = (;
        center_key = :origin,
        center_index = 1,
        location = (0.0, 0.0, 0.0),
        charge = 1.0,
    )
    axis_layers = (;
        x = bundle5.pgdg_intermediate.base_layer,
        y = bundle5.pgdg_intermediate.base_layer,
        z = bundle5.pgdg_intermediate.base_layer,
    )
    retained_nuclear =
        PQSH1CPBM.pqs_source_pair_retained_centered_electron_nuclear_by_center_block(
            record;
            axis_layers,
            coulomb_expansion = expansion,
            center_record = center,
        )

    final_overlap =
        PQSH1CPBM.pqs_source_shell_final_one_body_from_boundary_matrix(
            final_basis,
            _pqs_h1_boundary_matrix_result(
                :retained_source_overlap,
                retained_overlap_block.block,
            );
            term = :overlap,
        )
    final_kinetic =
        PQSH1CPBM.pqs_source_shell_final_one_body_from_boundary_matrix(
            final_basis,
            _pqs_h1_boundary_matrix_result(
                :retained_source_kinetic,
                retained_kinetic_block.block,
            );
            term = :kinetic,
        )
    final_nuclear =
        PQSH1CPBM.pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(
            final_basis,
            retained_nuclear,
        )
    final_hamiltonian = PQSH1CPBM.pqs_source_shell_final_one_electron_hamiltonian(
        final_kinetic,
        (final_nuclear,),
    )

    shell_hamiltonian =
        _pqs_h1_shell_support_kinetic(descriptor.support_states, metrics) .+
        center.charge .* _pqs_h1_shell_support_nuclear(
            descriptor.support_states,
            axis_layers,
            expansion,
            center,
        )
    R = final_basis.final_shell_coefficients
    oracle_hamiltonian = transpose(R) * shell_hamiltonian * R
    hamiltonian = final_hamiltonian.hamiltonian_matrix
    h1 = first(eigvals(Symmetric((hamiltonian + transpose(hamiltonian)) ./ 2)))
    oracle_h1 = first(
        eigvals(Symmetric((oracle_hamiltonian + transpose(oracle_hamiltonian)) ./ 2)),
    )

    @test raw_plan.source_mode_dims == (5, 5, 5)
    @test raw_plan.source_mode_count == 125
    @test retained_rule.retained_count == 98
    @test final_basis.final_retained_count == 98
    @test final_overlap.final_operator ≈ final_basis.final_overlap atol = 1.0e-12
    @test final_basis.final_overlap_identity_error < 1.0e-10
    @test retained_overlap_block.metadata.retained_direct_boundary_product_used
    @test retained_kinetic_block.metadata.retained_direct_boundary_product_used
    @test retained_nuclear.metadata.retained_direct_boundary_product_used
    @test !retained_overlap_block.metadata.raw_source_operator_block_materialized
    @test !retained_kinetic_block.metadata.raw_source_operator_block_materialized
    @test !retained_nuclear.metadata.raw_source_operator_block_materialized
    @test !retained_overlap_block.metadata.source_space_input_used
    @test !retained_kinetic_block.metadata.source_space_input_used
    @test !retained_nuclear.metadata.source_space_input_used
    @test all(isfinite, hamiltonian)
    @test hamiltonian ≈ transpose(hamiltonian) atol = 1.0e-12 rtol = 0.0
    @test hamiltonian ≈ oracle_hamiltonian atol = 1.0e-12 rtol = 0.0
    @test h1 ≈ oracle_h1 atol = 1.0e-12 rtol = 0.0
    @test final_hamiltonian.hamiltonian_data_materialized
    @test final_hamiltonian.h1_solve_materialized == false
    @test final_hamiltonian.eigensolve_materialized == false
    @test final_hamiltonian.generalized_overlap_solve_materialized == false
    @test final_hamiltonian.ida_data_materialized == false
    @test final_hamiltonian.density_density_materialized == false
    @test final_hamiltonian.rhf_materialized == false
    @test final_hamiltonian.driver_route_materialized == false
    @test final_hamiltonian.exports_materialized == false
    @test final_hamiltonian.artifacts_materialized == false
end
