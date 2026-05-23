using GaussletBases
using LinearAlgebra
using Printf

rb = build_basis(RadialBasisSpec(:G10;
    count = 12,
    mapping = AsinhMapping(c = 0.12, s = 0.16),
    reference_spacing = 1.0,
    tails = 3,
    odd_even_kmax = 2,
    xgaussians = [XGaussian(alpha = 0.12)],
))
grid = radial_quadrature(rb; accuracy = :medium, quadrature_rmax = 16.0)
points = quadrature_points(grid)

exponents = collect(exp.(range(log(0.04), log(18.0), length = 24)))
l = 1
design = GaussletBases._radial_ylm_solid_harmonic_gto_design_matrix(points, l, exponents)
target_coefficients = zeros(Float64, length(exponents), 6)
for orbital in axes(target_coefficients, 2)
    for exponent_index in axes(target_coefficients, 1)
        target_coefficients[exponent_index, orbital] =
            sin(0.37 * exponent_index * orbital) / (1.0 + exponent_index)
    end
end
targets = design * target_coefficients

fit_radial_ylm_to_solid_harmonic_gto(
    grid,
    targets,
    exponents;
    l = l,
    m = 0,
    metric_svd_cutoff = 1.0e-12,
)

GC.gc()
timed = @timed fit_radial_ylm_to_solid_harmonic_gto(
    grid,
    targets,
    exponents;
    l = l,
    m = 0,
    metric_svd_cutoff = 1.0e-12,
)
fit = timed.value
radial_ylm_fit_cartesian_gto_adapter(fit)

GC.gc()
adapter_timed = @timed radial_ylm_fit_cartesian_gto_adapter(fit)
adapter = adapter_timed.value
basis = build_basis(MappedUniformBasisSpec(:G10;
    count = 5,
    mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 4.0),
    reference_spacing = 1.0,
))

project_radial_ylm_gto_adapter_to_cartesian(basis, adapter)

GC.gc()
projection_timed = @timed project_radial_ylm_gto_adapter_to_cartesian(basis, adapter)
projection = projection_timed.value

d_fit = fit_radial_ylm_to_solid_harmonic_gto(
    grid,
    targets[:, 1:1],
    exponents;
    l = 2,
    m = 0,
    metric_svd_cutoff = 1.0e-12,
)
d_adapter = radial_ylm_fit_cartesian_gto_adapter(d_fit)
project_cartesian_gto_to_supplement_subspace(d_adapter, d_adapter.supplement)

GC.gc()
subspace_timed = @timed project_cartesian_gto_to_supplement_subspace(
    d_adapter,
    d_adapter.supplement,
)
subspace_projection = subspace_timed.value

@printf("radial_basis_size=%d\n", length(rb))
@printf("radial_grid_size=%d\n", length(points))
@printf("l=%d orbital_count=%d exponent_count=%d\n", l, size(targets, 2), length(exponents))
@printf("cartesian_working_dimension=%d\n", length(basis)^3)
@printf("cartesian_probe_count=%d coefficient_map=%dx%d\n",
    length(adapter.supplement.orbitals),
    size(adapter.coefficient_map, 1),
    size(adapter.coefficient_map, 2),
)
@printf("effective_rank=%d metric_condition=%.6e\n",
    fit.diagnostics.effective_rank,
    fit.diagnostics.metric_condition,
)
@printf("max_relative_residual=%.6e\n", maximum(fit.diagnostics.relative_residual_norms))
@printf("time_seconds=%.6f\n", timed.time)
@printf("alloc_mb=%.3f\n", timed.bytes / 1024^2)
@printf("adapter_time_seconds=%.6f\n", adapter_timed.time)
@printf("adapter_alloc_mb=%.3f\n", adapter_timed.bytes / 1024^2)
@printf("projection_columns=%d max_relative_norm_loss=%.6e\n",
    size(projection.cartesian_coefficients, 2),
    maximum(projection.diagnostics.relative_norm_losses),
)
@printf("projection_raw_max_singular_value=%.6e\n",
    maximum(projection.diagnostics.raw_projected_overlap_singular_values),
)
@printf("projection_source_orthonormal_max_singular_value=%.6e\n",
    maximum(projection.diagnostics.source_orthonormal_projected_singular_values),
)
@printf("projection_time_seconds=%.6f\n", projection_timed.time)
@printf("projection_alloc_mb=%.3f\n", projection_timed.bytes / 1024^2)
@printf("subspace_projection_source_l=%d source_columns=%d subspace_orbitals=%d\n",
    d_fit.l,
    subspace_projection.diagnostics.source_column_count,
    subspace_projection.diagnostics.subspace_orbital_count,
)
@printf("subspace_projection_max_relative_residual=%.6e\n",
    subspace_projection.diagnostics.max_relative_residual_norm,
)
@printf("subspace_projection_max_capture_singular_value=%.6e\n",
    maximum(subspace_projection.diagnostics.normalized_capture_singular_values),
)
@printf("subspace_projection_time_seconds=%.6f\n", subspace_timed.time)
@printf("subspace_projection_alloc_mb=%.3f\n", subspace_timed.bytes / 1024^2)
