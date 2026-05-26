using Printf
using LinearAlgebra

using GaussletBases

const CCP = GaussletBases.CartesianContractedParents
const CCPM = GaussletBases.CartesianContractedParentMetrics

function _axis_metric(pgdg)
    return (
        overlap = pgdg.overlap,
        position = pgdg.position,
        weights = pgdg.weights,
        centers = pgdg.centers,
        source = :nested_pgdg_axis,
    )
end

basis = bond_aligned_homonuclear_qw_basis(
    family = :G10,
    bond_length = 1.4,
    core_spacing = 0.7,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
    bond_axis = :z,
)
expansion = coulomb_gaussian_expansion(doacc = false)
bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
source = GaussletBases._nested_bond_aligned_diatomic_source(
    basis,
    bundles;
    bond_axis = :z,
    nside = 5,
    term_coefficients = Float64.(expansion.coefficients),
    packet_kernel = :factorized_direct,
    shared_shell_layer_policy = :endcap_panel_owned,
    shared_shell_endcap_panel_q = 4,
    shared_shell_endcap_panel_L = 4,
)
fixed_block = GaussletBases._nested_fixed_block(source)
contracted = CCP.cartesian_contracted_parent(fixed_block)
sidecar = fixed_block.staged_by_center_sidecar[]
axis_metrics = (
    x = _axis_metric(GaussletBases._nested_axis_pgdg(bundles, :x)),
    y = _axis_metric(GaussletBases._nested_axis_pgdg(bundles, :y)),
    z = _axis_metric(GaussletBases._nested_axis_pgdg(bundles, :z)),
)

# Compile both paths first, then report warm packet-construction timings.
support_warm = CCPM.cartesian_contracted_parent_metric_packet(
    contracted;
    axis_metrics,
    construction_path = :support_local_product,
)
product_warm = CCPM.cartesian_contracted_parent_metric_packet(
    contracted;
    axis_metrics,
    construction_path = :product_staged_metric_contraction,
)
support_timed = @timed CCPM.cartesian_contracted_parent_metric_packet(
    contracted;
    axis_metrics,
    construction_path = :support_local_product,
)
product_timed = @timed CCPM.cartesian_contracted_parent_metric_packet(
    contracted;
    axis_metrics,
    construction_path = :product_staged_metric_contraction,
)
support_packet = support_timed.value
product_packet = product_timed.value

max_overlap_diff = norm(product_packet.overlap - support_packet.overlap, Inf)
max_weight_diff = norm(product_packet.weights - support_packet.weights, Inf)
max_center_diff = norm(product_packet.centers - support_packet.centers, Inf)
max_fixed_overlap_diff = norm(product_packet.overlap - fixed_block.overlap, Inf)

println("Cartesian contracted parent metric packet staged-sidecar benchmark")
@printf("parent_dimension = %d\n", size(fixed_block.coefficient_matrix, 1))
@printf("contracted_dimension = %d\n", size(fixed_block.coefficient_matrix, 2))
@printf("staged_unit_count = %d\n", sidecar.diagnostics.unit_count)
@printf("product_unit_count = %d\n", sidecar.diagnostics.product_unit_count)
@printf("generic_unit_count = %d\n", sidecar.diagnostics.generic_unit_count)
@printf("max_support_count = %d\n", sidecar.diagnostics.max_support_count)
@printf("support_path = %s\n", String(support_packet.diagnostics.construction_path))
@printf("support_dense_parent_matrix_used = %s\n", string(support_packet.diagnostics.dense_parent_matrix_used))
@printf("support_warm_time_s = %.6f\n", support_timed.time)
@printf("support_warm_alloc_mib = %.3f\n", support_timed.bytes / 1024^2)
@printf("product_path = %s\n", String(product_packet.diagnostics.construction_path))
@printf("product_dense_parent_matrix_used = %s\n", string(product_packet.diagnostics.dense_parent_matrix_used))
@printf("product_block_count = %d\n", product_packet.diagnostics.product_block_count)
@printf("fallback_block_count = %d\n", product_packet.diagnostics.fallback_block_count)
@printf("product_warm_time_s = %.6f\n", product_timed.time)
@printf("product_warm_alloc_mib = %.3f\n", product_timed.bytes / 1024^2)
@printf("max_overlap_diff = %.6e\n", max_overlap_diff)
@printf("max_weight_diff = %.6e\n", max_weight_diff)
@printf("max_center_diff = %.6e\n", max_center_diff)
@printf("max_fixed_overlap_diff = %.6e\n", max_fixed_overlap_diff)
@printf("warmup_trace_guard = %.12e\n", tr(support_warm.overlap) + tr(product_warm.overlap))
