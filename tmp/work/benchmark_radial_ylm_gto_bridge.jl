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

@printf("radial_basis_size=%d\n", length(rb))
@printf("radial_grid_size=%d\n", length(points))
@printf("l=%d orbital_count=%d exponent_count=%d\n", l, size(targets, 2), length(exponents))
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
