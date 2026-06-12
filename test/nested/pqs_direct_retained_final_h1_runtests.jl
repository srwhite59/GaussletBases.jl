using Test
using LinearAlgebra
using GaussletBases

const PQSH1CFBR = GaussletBases.CartesianFinalBasisRealization
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

function _pqs_h1_support_product_matrix(left_states, right_states, mx, my, mz)
    result = Matrix{Float64}(undef, length(left_states), length(right_states))
    for (left_index, (ix, iy, iz)) in pairs(left_states)
        for (right_index, (jx, jy, jz)) in pairs(right_states)
            result[left_index, right_index] =
                mx[ix, jx] * my[iy, jy] * mz[iz, jz]
        end
    end
    return result
end

function _pqs_h1_support_kinetic_matrix(states, metrics)
    return (
        _pqs_h1_support_product_matrix(
            states,
            states,
            metrics.x.kinetic,
            metrics.y.overlap,
            metrics.z.overlap,
        ) +
        _pqs_h1_support_product_matrix(
            states,
            states,
            metrics.x.overlap,
            metrics.y.kinetic,
            metrics.z.overlap,
        ) +
        _pqs_h1_support_product_matrix(
            states,
            states,
            metrics.x.overlap,
            metrics.y.overlap,
            metrics.z.kinetic,
        )
    )
end

function _pqs_h1_support_nuclear_matrix(states, gaussian_factor_terms, expansion)
    result = zeros(Float64, length(states), length(states))
    for term_index in eachindex(expansion.coefficients)
        result .+=
            -Float64(expansion.coefficients[term_index]) *
            _pqs_h1_support_product_matrix(
                states,
                states,
                @view(gaussian_factor_terms[term_index, :, :]),
                @view(gaussian_factor_terms[term_index, :, :]),
                @view(gaussian_factor_terms[term_index, :, :]),
            )
    end
    return result
end

function _pqs_h1_core_indices_and_states(inner_box, dims)
    states = NTuple{3,Int}[]
    indices = Int[]
    for ix in inner_box[1], iy in inner_box[2], iz in inner_box[3]
        push!(states, (ix, iy, iz))
        push!(indices, GaussletBases._cartesian_flat_index(ix, iy, iz, dims))
    end
    return indices, states
end

function _pqs_h1_complete_fixture()
    expansion = coulomb_gaussian_expansion(doacc = false)
    current_box = (1:7, 1:7, 1:7)
    inner_box = (2:6, 2:6, 2:6)
    dims = (7, 7, 7)
    bundle7 = _pqs_h1_test_bundle(7)
    bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle7, bundle7, bundle7)
    metrics = _pqs_h1_axis_metrics(bundles)
    layer = GaussletBases._nested_projected_q_shell_layer(
        bundles,
        current_box,
        inner_box;
        bond_axis = :z,
        q = 5,
        L = 5,
        raw_source_dims = (5, 5, 5),
        selected_q = 5,
        term_coefficients = Float64.(expansion.coefficients),
    )
    descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(layer)
    shell_plan = PQSH1CCPM._pqs_shell_realization_plan(descriptor, metrics)
    core_indices, core_states = _pqs_h1_core_indices_and_states(inner_box, dims)
    shell_coefficients = shell_plan.shell_projection_matrix * shell_plan.lowdin_cleanup
    shell_states = descriptor.support_states
    core_overlap = _pqs_h1_support_product_matrix(
        core_states,
        core_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    core_shell_overlap = _pqs_h1_support_product_matrix(
        core_states,
        shell_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    shell_overlap = _pqs_h1_support_product_matrix(
        shell_states,
        shell_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    final_basis = PQSH1CFBR.pqs_complete_core_shell_final_basis(
        core_support_indices = core_indices,
        shell_support_indices = descriptor.support_indices,
        core_overlap = core_overlap,
        core_shell_overlap = core_shell_overlap,
        shell_overlap = shell_overlap,
        shell_final_coefficients = shell_coefficients,
        metadata = (;
            fixture = :pqs_direct_retained_final_h1_complete_core_shell,
            current_box,
            inner_box,
            raw_source_dims = (5, 5, 5),
        ),
    )
    return (;
        expansion,
        bundle7,
        metrics,
        current_box,
        inner_box,
        raw_source_dims = (5, 5, 5),
        core_states,
        shell_states,
        all_states = vcat(core_states, shell_states),
        final_basis,
    )
end

function _pqs_h1_fixed_block_oracle_energy(bundle, expansion)
    fixed_block = one_center_atomic_full_parent_fixed_block(
        bundle;
        expansion,
        nside = 5,
    )
    hamiltonian = Matrix{Float64}(fixed_block.kinetic) - Matrix{Float64}(fixed_block.gaussian_sum)
    return first(eigvals(Symmetric((hamiltonian + transpose(hamiltonian)) ./ 2)))
end

@testset "PQS complete core-shell final H1 gate" begin
    fixture = _pqs_h1_complete_fixture()
    center = (;
        center_key = :origin,
        center_index = 1,
        location = (0.0, 0.0, 0.0),
        charge = 1.0,
    )
    states = fixture.all_states
    metrics = fixture.metrics
    support_overlap = _pqs_h1_support_product_matrix(
        states,
        states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    support_kinetic = _pqs_h1_support_kinetic_matrix(states, metrics)
    support_nuclear = _pqs_h1_support_nuclear_matrix(
        states,
        fixture.bundle7.pgdg_intermediate.gaussian_factor_terms,
        fixture.expansion,
    )

    final_overlap = PQSH1CFBR.pqs_complete_core_shell_final_one_body_matrix(
        fixture.final_basis,
        support_overlap;
        term = :overlap,
    )
    final_kinetic = PQSH1CFBR.pqs_complete_core_shell_final_one_body_matrix(
        fixture.final_basis,
        support_kinetic;
        term = :kinetic,
    )
    final_nuclear = PQSH1CFBR.pqs_complete_core_shell_final_one_body_matrix(
        fixture.final_basis,
        support_nuclear;
        term = :electron_nuclear_by_center,
        center_record = center,
        metadata = (;
            nuclear_factor_source = :pgdg_intermediate_gaussian_factor_terms,
            raw_base_layer_gaussian_factor_matrices_used = false,
        ),
    )
    final_hamiltonian = PQSH1CFBR.pqs_complete_core_shell_final_one_electron_hamiltonian(
        final_kinetic,
        (final_nuclear,),
    )
    h1 = PQSH1CFBR.pqs_complete_core_shell_final_h1_solve(final_hamiltonian)
    fixed_h1 =
        _pqs_h1_fixed_block_oracle_energy(fixture.bundle7, fixture.expansion)

    @test fixture.current_box == (1:7, 1:7, 1:7)
    @test fixture.inner_box == (2:6, 2:6, 2:6)
    @test fixture.raw_source_dims == (5, 5, 5)
    @test fixture.final_basis.core_support_count == 125
    @test fixture.final_basis.shell_support_count == 218
    @test fixture.final_basis.shell_final_retained_count == 98
    @test fixture.final_basis.final_retained_count == 223
    @test fixture.final_basis.final_overlap_identity_error < 1.0e-10
    @test final_overlap.final_operator ≈ fixture.final_basis.final_overlap atol = 1.0e-12 rtol = 0.0
    @test final_overlap.final_operator ≈
          Matrix{Float64}(I, h1.final_dimension, h1.final_dimension) atol = 1.0e-10 rtol = 0.0
    @test final_nuclear.metadata.nuclear_factor_source ===
          :pgdg_intermediate_gaussian_factor_terms
    @test final_nuclear.metadata.raw_base_layer_gaussian_factor_matrices_used == false
    @test final_hamiltonian.hamiltonian_matrix_finite
    @test final_hamiltonian.hamiltonian_matrix ≈
          transpose(final_hamiltonian.hamiltonian_matrix) atol = 1.0e-12 rtol = 0.0
    @test h1.solve_kind === :ordinary_symmetric
    @test h1.lowest_energy ≈ -0.48047934800387226 atol = 1.0e-12 rtol = 0.0
    @test abs(h1.lowest_energy - fixed_h1) < 1.0e-6
    @test final_kinetic.old_fixed_block_matrix_authority_used == false
    @test final_nuclear.old_fixed_block_matrix_authority_used == false
    @test final_hamiltonian.metadata.old_fixed_block_matrix_authority_used == false
    @test final_hamiltonian.metadata.current_route_safe_term_matrices_used == false
    @test h1.generalized_overlap_solve_materialized == false
    @test h1.ida_data_materialized == false
    @test h1.density_density_materialized == false
    @test h1.rhf_materialized == false
    @test h1.driver_route_materialized == false
    @test h1.exports_materialized == false
    @test h1.artifacts_materialized == false
end
