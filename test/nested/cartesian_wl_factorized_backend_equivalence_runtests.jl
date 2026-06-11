using Test
using LinearAlgebra
using GaussletBases

const WLFactorizedEquivalenceCPBM =
    GaussletBases.CartesianPairBlockMaterialization

function _factorized_equivalence_basis()
    return build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 7,
            mapping = AsinhMapping(c = 0.2, s = sqrt(0.2 * 2.0), tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
end

function _factorized_equivalence_axis_inputs(basis, expansion)
    source_1d = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = source_1d.pgdg_intermediate
    return (;
        parent_axis_counts = ntuple(_ -> length(basis.center_data), 3),
        parent_axis_bundle_object = (; x = source_1d, y = source_1d, z = source_1d),
        overlap_1d = (; x = pgdg.overlap, y = pgdg.overlap, z = pgdg.overlap),
        kinetic_1d = (; x = pgdg.kinetic, y = pgdg.kinetic, z = pgdg.kinetic),
    )
end

function _factorized_equivalence_inventory(basis, axis_inputs)
    source =
        WLFactorizedEquivalenceCPBM.white_lindsey_shellification_decomposed_unit_pair_inventory(
            ntuple(_ -> basis.center_data, 3),
            ((0.0, 0.0, 0.0),);
            shellification_policy =
                GaussletBases.CartesianShellification.OneCenterShellification(
                    core_side = 5,
                    q = 5,
                ),
            metadata = (; q = 5, ns = 5),
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
        )
    return source.inventory
end

function _factorized_equivalence_lowest_orbital(hamiltonian, overlap)
    sym_h = Symmetric((hamiltonian + transpose(hamiltonian)) ./ 2)
    overlap_identity = isapprox(
        overlap,
        Matrix{Float64}(I, size(overlap)...);
        atol = 1.0e-9,
        rtol = 1.0e-9,
    )
    eig = overlap_identity ?
        eigen(sym_h) :
        eigen(sym_h, Symmetric((overlap + transpose(overlap)) ./ 2))
    index = argmin(eig.values)
    return (; energy = eig.values[index], coefficients = eig.vectors[:, index])
end

function _factorized_equivalence_interaction_energy(interaction, density)
    v =
        0.5 .* (Matrix{Float64}(interaction) .+ transpose(Matrix{Float64}(interaction)))
    rho = 0.5 .* (Matrix{Float64}(density) .+ transpose(Matrix{Float64}(density)))
    occupations = vec(diag(rho))
    direct = 2.0 * dot(occupations, v * occupations)
    exchange = dot(vec(rho), vec(v .* rho))
    return direct - exchange
end

function _factorized_equivalence_streaming_density_density(
    inventory,
    axis_inputs,
    expansion,
)
    axis_counts =
        WLFactorizedEquivalenceCPBM._axis_counts_tuple(axis_inputs.parent_axis_counts)
    raw_pair_factor_terms =
        WLFactorizedEquivalenceCPBM._white_lindsey_density_density_pair_factor_terms(
            axis_inputs.parent_axis_bundle_object,
        )
    coefficients =
        WLFactorizedEquivalenceCPBM._white_lindsey_density_density_expansion_coefficients(
            expansion,
        )
    coefficient_cache =
        WLFactorizedEquivalenceCPBM._white_lindsey_decomposed_operator_unit_coefficient_cache(
            inventory,
        )
    retained_weights =
        WLFactorizedEquivalenceCPBM._white_lindsey_density_density_retained_weights(
            inventory,
            axis_inputs.parent_axis_bundle_object,
            axis_counts,
            coefficient_cache,
        )
    @test retained_weights.status === :materialized_retained_density_weights
    matrix = zeros(Float64, inventory.retained_dimension, inventory.retained_dimension)
    stream_state =
        WLFactorizedEquivalenceCPBM._route_global_density_density_streaming_fill!(
            matrix,
            inventory.unit_pairs,
            axis_counts,
            raw_pair_factor_terms,
            coefficients;
            unit_coefficient_cache = coefficient_cache,
            prepared_unit_cache =
                WLFactorizedEquivalenceCPBM._white_lindsey_decomposed_operator_prepared_unit_cache(
                    inventory,
                ),
        )
    @test isnothing(stream_state.blocker)
    WLFactorizedEquivalenceCPBM._white_lindsey_apply_retained_density_weight_division!(
        matrix,
        retained_weights.weights,
    )
    return matrix
end

@testset "decomposed WL factorized backend matches pair-streaming reference" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = _factorized_equivalence_basis()
    axis_inputs = _factorized_equivalence_axis_inputs(basis, expansion)
    inventory = _factorized_equivalence_inventory(basis, axis_inputs)
    center_record = (;
        center_key = :helium_nucleus,
        center_index = 1,
        nuclear_charge = 2.0,
        location = (0.0, 0.0, 0.0),
    )

    factorized_overlap =
        WLFactorizedEquivalenceCPBM._route_global_decomposed_wl_factorized_one_body_matrix(
            inventory,
            :overlap;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            overlap_1d = axis_inputs.overlap_1d,
        ).global_matrix_result.matrix
    streaming_overlap =
        WLFactorizedEquivalenceCPBM._route_global_decomposed_wl_streaming_one_body_matrix(
            inventory,
            :overlap;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            overlap_1d = axis_inputs.overlap_1d,
        ).global_matrix_result.matrix

    factorized_kinetic =
        WLFactorizedEquivalenceCPBM._route_global_decomposed_wl_factorized_one_body_matrix(
            inventory,
            :kinetic;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            overlap_1d = axis_inputs.overlap_1d,
            kinetic_1d = axis_inputs.kinetic_1d,
        ).global_matrix_result.matrix
    streaming_kinetic =
        WLFactorizedEquivalenceCPBM._route_global_decomposed_wl_streaming_one_body_matrix(
            inventory,
            :kinetic;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            overlap_1d = axis_inputs.overlap_1d,
            kinetic_1d = axis_inputs.kinetic_1d,
        ).global_matrix_result.matrix

    factorized_nuclear =
        WLFactorizedEquivalenceCPBM._route_global_decomposed_wl_factorized_electron_nuclear_by_center_matrix(
            inventory;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
            coulomb_expansion = expansion,
            center_record,
        ).global_matrix_result.matrix
    streaming_nuclear =
        WLFactorizedEquivalenceCPBM._route_global_decomposed_wl_streaming_electron_nuclear_by_center_matrix(
            inventory;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
            coulomb_expansion = expansion,
            center_record,
        ).global_matrix_result.matrix

    factorized_density =
        WLFactorizedEquivalenceCPBM.route_global_decomposed_wl_density_density_matrix(
            inventory;
            parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
            coulomb_expansion = expansion,
        ).matrix
    streaming_density = _factorized_equivalence_streaming_density_density(
        inventory,
        axis_inputs,
        expansion,
    )

    overlap_difference = maximum(abs.(factorized_overlap .- streaming_overlap))
    kinetic_difference = maximum(abs.(factorized_kinetic .- streaming_kinetic))
    nuclear_difference = maximum(abs.(factorized_nuclear .- streaming_nuclear))
    density_difference = maximum(abs.(factorized_density .- streaming_density))

    @test overlap_difference <= 1.0e-12
    @test kinetic_difference <= 1.0e-12
    @test nuclear_difference <= 1.0e-12
    @test density_difference <= 1.0e-12

    factorized_h1 = factorized_kinetic .- 2.0 .* factorized_nuclear
    streaming_h1 = streaming_kinetic .- 2.0 .* streaming_nuclear
    factorized_solve =
        _factorized_equivalence_lowest_orbital(factorized_h1, factorized_overlap)
    streaming_solve =
        _factorized_equivalence_lowest_orbital(streaming_h1, streaming_overlap)
    factorized_density_orbital =
        factorized_solve.coefficients * transpose(factorized_solve.coefficients)
    streaming_density_orbital =
        streaming_solve.coefficients * transpose(streaming_solve.coefficients)

    factorized_j = _factorized_equivalence_interaction_energy(
        factorized_density,
        factorized_density_orbital,
    )
    streaming_j = _factorized_equivalence_interaction_energy(
        streaming_density,
        streaming_density_orbital,
    )

    @test abs(factorized_solve.energy - streaming_solve.energy) <= 1.0e-12
    @test abs(factorized_j - streaming_j) <= 1.0e-12
end
