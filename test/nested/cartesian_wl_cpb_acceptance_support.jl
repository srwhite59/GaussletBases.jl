using LinearAlgebra
using GaussletBases

const WLCBPAcceptanceCPB = GaussletBases.CartesianCPB
const WLCBPAcceptanceCPGB = GaussletBases.CartesianParentGaussletBases
const WLCBPAcceptanceProvider = GaussletBases.CartesianCPBBlockProviders

const _WL_CPB_ACCEPTANCE_Q = 5
const _WL_CPB_ACCEPTANCE_NS = _WL_CPB_ACCEPTANCE_Q
const _WL_CPB_ACCEPTANCE_D = 0.15
const _WL_CPB_ACCEPTANCE_REFERENCE_SPACING = 1.0
const _WL_CPB_ACCEPTANCE_TAIL_SPACING = 10.0
const _WL_CPB_ACCEPTANCE_BACKEND = :pgdg_localized_experimental

function _wl_cpb_acceptance_axis(count::Integer)
    return build_basis(MappedUniformBasisSpec(
        :G10;
        count = Int(count),
        mapping = white_lindsey_atomic_mapping(
            Z = 1.0,
            d = _WL_CPB_ACCEPTANCE_D,
            tail_spacing = _WL_CPB_ACCEPTANCE_TAIL_SPACING,
        ),
        reference_spacing = _WL_CPB_ACCEPTANCE_REFERENCE_SPACING,
    ))
end

function _wl_cpb_acceptance_bundle(axis, expansion)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        axis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = _WL_CPB_ACCEPTANCE_BACKEND,
    )
end

function _wl_cpb_acceptance_parent(axis_counts)
    axis_x = _wl_cpb_acceptance_axis(axis_counts.x)
    axis_y = axis_counts.y == axis_counts.x ?
             axis_x :
             _wl_cpb_acceptance_axis(axis_counts.y)
    axis_z =
        axis_counts.z == axis_counts.x ? axis_x :
        axis_counts.z == axis_counts.y ? axis_y :
        _wl_cpb_acceptance_axis(axis_counts.z)
    axes = (
        x = axis_x,
        y = axis_y,
        z = axis_z,
    )
    parent = WLCBPAcceptanceCPGB.CartesianParentGaussletBasis3D(
        axes;
        metadata = (;
            basis_family = :post_cpb_white_lindsey_acceptance_fixture,
            q = _WL_CPB_ACCEPTANCE_Q,
            n_s = _WL_CPB_ACCEPTANCE_NS,
            core_spacing = _WL_CPB_ACCEPTANCE_D,
            q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    axis_bundle = (;
        x = _wl_cpb_acceptance_bundle(axes.x, expansion),
        y = _wl_cpb_acceptance_bundle(axes.y, expansion),
        z = _wl_cpb_acceptance_bundle(axes.z, expansion),
    )
    return (; parent, axis_bundle, expansion)
end

function _wl_cpb_acceptance_checked_block(block, expected_status)
    block_summary = WLCBPAcceptanceProvider.summary(block)
    block_summary.status === expected_status || error(
        "post-CPB WL acceptance block failed: status=$(block_summary.status), blocker=$(block_summary.blocker)",
    )
    return block_summary
end

function _wl_cpb_acceptance_result(;
    system::Symbol,
    centers,
    charges,
    axis_counts = (x = 7, y = 7, z = 7),
    old_direct_electronic_baseline,
    physical_electronic_reference,
    physical_total_reference = nothing,
    proton_repulsion = 0.0,
)
    setup = _wl_cpb_acceptance_parent(axis_counts)
    parent = setup.parent
    axis_bundle = setup.axis_bundle
    expansion = setup.expansion
    packet = WLCBPAcceptanceCPGB.parent_overlap_axis_factor_packet(parent, axis_bundle)
    cpb = WLCBPAcceptanceCPB.cpb(
        1:axis_counts.x,
        1:axis_counts.y,
        1:axis_counts.z;
        role = :wl_cpb_acceptance_full_parent_window,
    )
    interval_pair = WLCBPAcceptanceProvider.cpb_interval_pair(parent, cpb, cpb)

    assembly_elapsed_s = @elapsed begin
        overlap_block =
            WLCBPAcceptanceProvider.cpb_overlap_operator_block(packet, interval_pair)
        kinetic_block =
            WLCBPAcceptanceProvider.cpb_kinetic_operator_block(packet, interval_pair)
        overlap_summary = _wl_cpb_acceptance_checked_block(
            overlap_block,
            :materialized_cpb_overlap_operator_block,
        )
        kinetic_summary = _wl_cpb_acceptance_checked_block(
            kinetic_block,
            :materialized_cpb_kinetic_operator_block,
        )
        nuclear_blocks = Tuple(
            WLCBPAcceptanceProvider.cpb_electron_nuclear_by_center_local_block(
                axis_bundle,
                expansion,
                (;
                    center_key = Symbol(:center_, center_index),
                    center_index,
                    nuclear_charge = 1.0,
                    location = centers[center_index],
                ),
                interval_pair,
            ) for center_index in eachindex(centers)
        )
        nuclear_summaries = Tuple(
            _wl_cpb_acceptance_checked_block(
                block,
                :materialized_cpb_electron_nuclear_by_center_local_block,
            ) for block in nuclear_blocks
        )
        overlap_matrix = overlap_block.axis_product_block.dense_block
        hamiltonian =
            Matrix{Float64}(kinetic_block.sum_axis_product_block.dense_block)
        for block in nuclear_blocks
            hamiltonian .+= charges[block.center_record.center_index] .* block.dense_block
        end
    end

    solve_elapsed_s = @elapsed solve = eigen(
        Hermitian(hamiltonian),
        Hermitian(overlap_matrix),
    )
    electronic_energy = Float64(first(solve.values))
    total_energy = electronic_energy + proton_repulsion
    axis_centers = centers == ((0.0, 0.0, 0.0),) ?
                   GaussletBases.centers(WLCBPAcceptanceCPGB.axis_basis(parent, :x)) :
                   GaussletBases.centers(WLCBPAcceptanceCPGB.axis_basis(parent, :z))
    central_axis_spacing =
        minimum(diff(sort(Float64[Float64(center) for center in axis_centers])))

    report = (;
        route = :post_cpb_wl_gausslet_only,
        system,
        q = _WL_CPB_ACCEPTANCE_Q,
        n_s = _WL_CPB_ACCEPTANCE_NS,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = _WL_CPB_ACCEPTANCE_D,
        reference_spacing = _WL_CPB_ACCEPTANCE_REFERENCE_SPACING,
        tail_spacing = _WL_CPB_ACCEPTANCE_TAIL_SPACING,
        backend = _WL_CPB_ACCEPTANCE_BACKEND,
        basis_family = :G10,
        axis_counts,
        parent_support_size = axis_counts.x * axis_counts.y * axis_counts.z,
        basis_dimension = size(hamiltonian, 1),
        center_count = length(centers),
        centers,
        nuclear_charges = charges,
        gaussian_expansion_term_count = length(expansion),
        cpb_local_operator_blocks_materialized = true,
        cpb_overlap_block_status = overlap_summary.status,
        cpb_kinetic_block_status = kinetic_summary.status,
        cpb_electron_nuclear_block_statuses =
            Tuple(summary.status for summary in nuclear_summaries),
        cpb_local_electron_nuclear_by_center_blocks_materialized =
            all(summary -> summary.provider_level_electron_nuclear_block_materialized,
                nuclear_summaries),
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        nuclear_charge_application_stage = :wl_acceptance_hamiltonian_assembly,
        solve_kind = :generalized,
        lowest_electronic_energy = electronic_energy,
        proton_proton_repulsion = proton_repulsion,
        born_oppenheimer_total_energy = total_energy,
        distance_from_physical_electronic_reference =
            electronic_energy - physical_electronic_reference,
        distance_from_physical_total_reference =
            isnothing(physical_total_reference) ?
            nothing :
            total_energy - physical_total_reference,
        distance_from_old_direct_electronic_baseline =
            electronic_energy - old_direct_electronic_baseline,
        central_axis_spacing,
        axis_extent = (minimum(axis_centers), maximum(axis_centers)),
        assembly_elapsed_s,
        solve_elapsed_s,
    )
    return (; overlap_matrix, hamiltonian, report)
end
