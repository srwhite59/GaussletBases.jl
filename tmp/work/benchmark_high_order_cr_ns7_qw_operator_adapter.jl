using GaussletBases
using LinearAlgebra
using Printf

function _timed(f, label)
    GC.gc()
    timing = @timed f()
    @printf("%-32s %9.5f s %9.2f MiB\n", label, timing.time, timing.bytes / 2.0^20)
    return timing.value, timing
end

smoke, smoke_timing = _timed("Cr-map count7 smoke") do
    GaussletBases._experimental_high_order_cr_map_count7_pgdg_smoke_diagnostic()
end
diagnostic, adapter_timing = _timed("Cr-map count7 QW adapter") do
    GaussletBases._experimental_high_order_stack_to_atomic_qw_operator_diagnostic(smoke)
end

fields = diagnostic.diagnostics
operators = diagnostic.operators
overlap, overlap_timing = _timed("Cr-map count7 GTO overlap") do
    gto_overlap_matrix(operators, diagnostic.supplement)
end
warm_overlap, warm_overlap_timing = _timed("Cr-map count7 GTO overlap warm") do
    gto_overlap_matrix(operators, diagnostic.supplement)
end
reference_handoff, reference_timing = _timed("Cr-map count7 handoff overlap ref") do
    GaussletBases._cartesian_final_gto_cross_overlap_handoff(
        diagnostic.fixed_block,
        diagnostic.supplement,
        operators.raw_to_final,
        diagnostic.supplement;
        provenance = :benchmark_high_order_cr_ns7_qw_operator_adapter,
    )
end
reference_overlap = reference_handoff.cross_overlap
overlap_error = maximum(abs.(overlap .- reference_overlap))
warm_overlap_error = maximum(abs.(warm_overlap .- reference_overlap))

println()
println("Cr-map count=7 high-order PGDG QW adapter diagnostic")
println("route                           ", fields.route)
println("classification                  ", fields.classification)
println("operator construction scope     ", fields.operator_construction_scope)
println("interaction treatment           ", fields.interaction_treatment)
println("stack backend                   ", fields.stack_backend)
println("fixed backend                   ", fields.fixed_backend)
println("adapter backend                 ", fields.adapter_backend)
println("operator backend                ", fields.operator_backend)
println("parent side                     ", fields.parent_side)
println("parent dimension                ", fields.parent_dimension)
println("route comparability             ", fields.route_comparability)
println("ordinary ns7 comparable         ", fields.ordinary_ns7_comparable)
println("ordinary ns7 ref side           ", fields.ordinary_ns7_reference_parent_side)
println("ordinary ns7 ref dimension      ", fields.ordinary_ns7_reference_parent_dimension)
println("high-order retained dimension   ", fields.high_order_retained_dimension)
println("support count                   ", fields.support_count)
println("supplement dimension            ", fields.supplement_dimension)
println("final operator dimension        ", fields.final_operator_dimension)
println("residual count                  ", fields.residual_count)
println("contracted weight zeroish       ", fields.contracted_weight_zeroish_count)
println("contracted weight negative      ", fields.contracted_weight_negative_count)
@printf("contracted weight min/max       %.12e / %.12e\n",
    fields.contracted_weight_minimum,
    fields.contracted_weight_maximum,
)
@printf("fixed overlap error             %.6e\n", fields.fixed_overlap_error)
@printf("operator overlap error          %.6e\n", fields.operator_overlap_error)
@printf("H symmetry error                %.6e\n", fields.h_symmetry_error)
@printf("V symmetry error                %.6e\n", fields.v_symmetry_error)
println("MWG widths finite positive      ", fields.residual_widths_finite_positive)
println("residual centers finite         ", all(isfinite, operators.residual_centers))
println("same-density status             ", fields.same_density_two_electron_evaluation)
println("same-density route compare      ", fields.same_density_route_comparison)
println("occupied capture status         ", fields.occupied_capture_status)
println("smallest remaining interface    ", fields.smallest_missing_interface)
println("GTO overlap shape               ", size(overlap))
println("GTO handoff contract            ", reference_handoff.diagnostics.cross_overlap_contract)
@printf("GTO dense fallback max error     %.6e\n", overlap_error)
@printf("warm GTO fallback max error      %.6e\n", warm_overlap_error)
@printf("internal smoke time/allocation  %.5f s / %.2f MiB\n",
    fields.timing_seconds.smoke,
    fields.allocation_bytes.smoke / 2.0^20,
)
@printf("axis data time/allocation       %.5f s / %.2f MiB\n",
    fields.timing_seconds.axis_data,
    fields.allocation_bytes.axis_data / 2.0^20,
)
@printf("stack time/allocation           %.5f s / %.2f MiB\n",
    fields.timing_seconds.stack,
    fields.allocation_bytes.stack / 2.0^20,
)
@printf("fixed packet time/allocation    %.5f s / %.2f MiB\n",
    fields.timing_seconds.fixed_block,
    fields.allocation_bytes.fixed_block / 2.0^20,
)
@printf("operator time/allocation        %.5f s / %.2f MiB\n",
    fields.timing_seconds.operators,
    fields.allocation_bytes.operators / 2.0^20,
)
@printf("adapter total time/allocation   %.5f s / %.2f MiB\n",
    fields.timing_seconds.total,
    fields.allocation_bytes.total / 2.0^20,
)
@printf("script smoke wall/allocation    %.5f s / %.2f MiB\n",
    smoke_timing.time,
    smoke_timing.bytes / 2.0^20,
)
@printf("script adapter wall/allocation  %.5f s / %.2f MiB\n",
    adapter_timing.time,
    adapter_timing.bytes / 2.0^20,
)
@printf("script GTO overlap allocation   %.5f s / %.2f MiB\n",
    overlap_timing.time,
    overlap_timing.bytes / 2.0^20,
)
@printf("script warm GTO overlap alloc    %.5f s / %.2f MiB\n",
    warm_overlap_timing.time,
    warm_overlap_timing.bytes / 2.0^20,
)
@printf("script dense ref allocation     %.5f s / %.2f MiB\n",
    reference_timing.time,
    reference_timing.bytes / 2.0^20,
)
