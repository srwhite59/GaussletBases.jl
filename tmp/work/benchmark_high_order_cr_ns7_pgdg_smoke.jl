using Printf
using GaussletBases

cold = @timed GaussletBases._experimental_high_order_cr_map_count7_pgdg_smoke_diagnostic()
warm = @timed GaussletBases._experimental_high_order_cr_map_count7_pgdg_smoke_diagnostic()
result = warm.value
diagnostics = result.diagnostics

function _mib(bytes)
    return bytes / 1024^2
end

println("Cr-map count=7 high-order PGDG smoke diagnostic")
println("route                       ", diagnostics.route)
println("classification              ", diagnostics.classification)
println("backend                     ", diagnostics.backend)
println("mapping family              ", diagnostics.mapping_family)
println("high-order 7 meaning        ", diagnostics.high_order_7_meaning)
println("parent count meaning        ", diagnostics.parent_count_meaning)
println("parent side                 ", diagnostics.parent_side)
println("parent dimension            ", diagnostics.parent_dimension)
println("route comparability         ", diagnostics.route_comparability)
println("ordinary ns7 comparable     ", diagnostics.ordinary_ns7_comparable)
println("ordinary ns7 ref side       ", diagnostics.ordinary_ns7_reference_parent_side)
println("ordinary ns7 ref dimension  ", diagnostics.ordinary_ns7_reference_parent_dimension)
println("doside                      ", diagnostics.doside)
println("sides                       ", diagnostics.sides)
println("retained dimension          ", diagnostics.retained_dimension)
println("block labels                ", diagnostics.block_labels)
println("block column ranges         ", diagnostics.block_column_ranges)
println("shell dimensions            ", diagnostics.shell_dimensions)
println("shell cleanup kept ranks    ", [spectrum.kept_rank for spectrum in diagnostics.shell_cleanup_spectra])
println("contracted weight zeroish   ", count(abs.(result.stack.contracted_weights) .<= 1.0e-14))
println("contracted weight negative  ", count(result.stack.contracted_weights .< -1.0e-14))
@printf("contracted weight min/max   %.12e / %.12e\n",
    minimum(result.stack.contracted_weights),
    maximum(result.stack.contracted_weights),
)
@printf("overlap error               %.6e\n", diagnostics.overlap_error)
@printf("overlap eig min/max         %.12e / %.12e\n",
    diagnostics.overlap_spectrum.minimum_eigenvalue,
    diagnostics.overlap_spectrum.maximum_eigenvalue,
)
println("contracted weights finite   ", diagnostics.contracted_weights_finite)
println("axis overlap finite         ", diagnostics.axis_overlap_finite)
println("axis weights finite         ", diagnostics.axis_weight_finite)
println("reaches QW operators        ", diagnostics.reaches_atomic_qw_operators)
println("same-density evaluation     ", diagnostics.same_density_operator_evaluation)
println("same-density route compare  ", diagnostics.same_density_route_comparison)
println("occupied capture status     ", diagnostics.occupied_capture_status)
println("smallest missing interface  ", diagnostics.smallest_missing_interface)
@printf("axis data time/allocation   %.5f s / %.2f MiB\n",
    diagnostics.timing_seconds.axis_data,
    _mib(diagnostics.allocation_bytes.axis_data),
)
@printf("stack time/allocation       %.5f s / %.2f MiB\n",
    diagnostics.timing_seconds.stack,
    _mib(diagnostics.allocation_bytes.stack),
)
@printf("cold entry time/allocation  %.5f s / %.2f MiB\n",
    cold.time,
    _mib(cold.bytes),
)
@printf("warm entry time/allocation  %.5f s / %.2f MiB\n",
    warm.time,
    _mib(warm.bytes),
)
