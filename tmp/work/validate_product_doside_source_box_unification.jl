using LinearAlgebra
using GaussletBases

let CCPM = GaussletBases.CartesianContractedParentMetrics,
    TERMS = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )

function product_doside_probe_metrics()
    return (
        x = (
            overlap = Matrix{Float64}(I, 2, 2),
            position = Matrix(Diagonal([0.25, 1.75])),
            x2 = Matrix(Diagonal([0.25^2, 1.75^2])),
            kinetic = [
                2.0 -0.25
                -0.25 2.5
            ],
            weights = [1.0, 1.0],
            centers = [0.25, 1.75],
            source = :product_doside_source_box_probe,
        ),
        y = (
            overlap = Matrix{Float64}(I, 2, 2),
            position = Matrix(Diagonal([-0.5, 0.5])),
            x2 = Matrix(Diagonal([0.25, 0.25])),
            kinetic = [
                1.5 0.125
                0.125 1.25
            ],
            weights = [1.0, 1.0],
            centers = [-0.5, 0.5],
            source = :product_doside_source_box_probe,
        ),
        z = (
            overlap = Matrix{Float64}(I, 1, 1),
            position = Matrix(Diagonal([3.25])),
            x2 = Matrix(Diagonal([3.25^2])),
            kinetic = Matrix(Diagonal([0.75])),
            weights = [1.0],
            centers = [3.25],
            source = :product_doside_source_box_probe,
        ),
    )
end

function product_coefficients(support_states, axes, axis_function_indices)
    coefficients = zeros(Float64, length(support_states), length(axis_function_indices))
    for col in eachindex(axis_function_indices)
        fx, fy, fz = axis_function_indices[col]
        for row in eachindex(support_states)
            sx, sy, sz = support_states[row]
            coefficients[row, col] =
                axes[1].coefficient_matrix[sx, fx] *
                axes[2].coefficient_matrix[sy, fy] *
                axes[3].coefficient_matrix[sz, fz]
        end
    end
    return coefficients
end

function product_unit(name::Symbol, axes, axis_function_indices)
    support_states =
        NTuple{3,Int}[(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)]
    coefficients = product_coefficients(support_states, axes, axis_function_indices)
    retained_count = length(axis_function_indices)
    return GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        name,
        :product_doside,
        1:retained_count,
        collect(1:4),
        support_states,
        coefficients,
        axes,
        axis_function_indices,
        (source = :validate_product_doside_source_box_unification,),
        (support_count = 4, retained_count = retained_count),
    )
end

function authoritative_block(left_unit, right_unit, metrics, term::Symbol)
    term == :kinetic && return CCPM._product_doside_retained_kinetic_block(
        left_unit,
        right_unit,
        metrics,
    )
    return CCPM._product_doside_retained_low_order_block(
        left_unit,
        right_unit,
        metrics;
        term,
    )
end

function axis_metric_matrix(metrics, axis::Int, kind::Symbol)
    return getproperty(getproperty(metrics, (:x, :y, :z)[axis]), kind)
end

function term_factor_kinds(term::Symbol)
    term == :overlap && return ((:overlap, :overlap, :overlap),)
    term == :position_x && return ((:position, :overlap, :overlap),)
    term == :position_y && return ((:overlap, :position, :overlap),)
    term == :position_z && return ((:overlap, :overlap, :position),)
    term == :x2_x && return ((:x2, :overlap, :overlap),)
    term == :x2_y && return ((:overlap, :x2, :overlap),)
    term == :x2_z && return ((:overlap, :overlap, :x2),)
    term == :kinetic && return (
        (:kinetic, :overlap, :overlap),
        (:overlap, :kinetic, :overlap),
        (:overlap, :overlap, :kinetic),
    )
    throw(ArgumentError("unsupported explicit support reference term $(term)"))
end

function explicit_support_operator(left_unit, right_unit, metrics, term::Symbol)
    operator = zeros(
        Float64,
        length(left_unit.support_states),
        length(right_unit.support_states),
    )
    for factor_kinds in term_factor_kinds(term)
        factors = ntuple(
            axis -> axis_metric_matrix(metrics, axis, factor_kinds[axis]),
            3,
        )
        for col in eachindex(right_unit.support_states)
            xj, yj, zj = right_unit.support_states[col]
            for row in eachindex(left_unit.support_states)
                xi, yi, zi = left_unit.support_states[row]
                operator[row, col] +=
                    factors[1][xi, xj] * factors[2][yi, yj] * factors[3][zi, zj]
            end
        end
    end
    return operator
end

function explicit_support_reference_block(left_unit, right_unit, metrics, term::Symbol)
    support_operator = explicit_support_operator(left_unit, right_unit, metrics, term)
    return transpose(left_unit.coefficient_matrix) *
           support_operator *
           right_unit.coefficient_matrix
end

function local_gaussian_fixture_terms()
    axis_terms_x = Array{Float64,3}(undef, 2, 2, 2)
    axis_terms_x[1, :, :] .= [
        1.0 0.125
        0.125 0.75
    ]
    axis_terms_x[2, :, :] .= [
        0.35 -0.05
        -0.05 0.45
    ]
    axis_terms_y = Array{Float64,3}(undef, 2, 2, 2)
    axis_terms_y[1, :, :] .= [
        0.9 -0.1
        -0.1 0.8
    ]
    axis_terms_y[2, :, :] .= [
        0.4 0.075
        0.075 0.5
    ]
    axis_terms_z = Array{Float64,3}(undef, 2, 1, 1)
    axis_terms_z[1, 1, 1] = 1.25
    axis_terms_z[2, 1, 1] = 0.6
    return (
        term_coefficients = [0.7, 0.2],
        axis_gaussian_terms = (
            x = axis_terms_x,
            y = axis_terms_y,
            z = axis_terms_z,
        ),
    )
end

function centered_local_gaussian_probe_expansion()
    return CoulombGaussianExpansion(
        [0.7, 0.2],
        [0.35, 0.9];
        del = 0.0,
        s = 1.0,
        c = 1.0,
        maxu = 1.0,
    )
end

function centered_local_gaussian_axis_layers()
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 2,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return (x = basis, y = basis, z = basis)
end

function explicit_centered_axis_terms(axis_layers, expansion, center)
    return (
        x = _term_first_axis_array(
            gaussian_factor_matrices(axis_layers.x;
                exponents = expansion.exponents,
                center = center[1],
            ),
        ),
        y = _term_first_axis_array(
            gaussian_factor_matrices(axis_layers.y;
                exponents = expansion.exponents,
                center = center[2],
            ),
        ),
        z = _term_first_axis_array(
            gaussian_factor_matrices(axis_layers.z;
                exponents = expansion.exponents,
                center = center[3],
            ),
        ),
    )
end

function _term_first_axis_array(matrices)
    nterms = length(matrices)
    first_matrix = Matrix{Float64}(matrices[1])
    terms = Array{Float64,3}(undef, nterms, size(first_matrix, 1), size(first_matrix, 2))
    terms[1, :, :] .= first_matrix
    for term in 2:nterms
        terms[term, :, :] .= Matrix{Float64}(matrices[term])
    end
    return terms
end

function explicit_local_gaussian_support_operator(
    left_unit,
    right_unit,
    term_coefficients,
    axis_gaussian_terms,
)
    operator = zeros(
        Float64,
        length(left_unit.support_states),
        length(right_unit.support_states),
    )
    for col in eachindex(right_unit.support_states)
        xj, yj, zj = right_unit.support_states[col]
        for row in eachindex(left_unit.support_states)
            xi, yi, zi = left_unit.support_states[row]
            value = 0.0
            for term in eachindex(term_coefficients)
                value +=
                    term_coefficients[term] *
                    axis_gaussian_terms.x[term, xi, xj] *
                    axis_gaussian_terms.y[term, yi, yj] *
                    axis_gaussian_terms.z[term, zi, zj]
            end
            operator[row, col] = value
        end
    end
    return operator
end

function explicit_local_gaussian_reference_block(
    left_unit,
    right_unit,
    term_coefficients,
    axis_gaussian_terms,
)
    support_operator = explicit_local_gaussian_support_operator(
        left_unit,
        right_unit,
        term_coefficients,
        axis_gaussian_terms,
    )
    return transpose(left_unit.coefficient_matrix) *
           support_operator *
           right_unit.coefficient_matrix
end

identity_axis = Matrix{Float64}(I, 2, 2)
rotation_axis = [
    inv(sqrt(2.0)) inv(sqrt(2.0))
    inv(sqrt(2.0)) -inv(sqrt(2.0))
]
identity_axes = (
    GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
    GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
    GaussletBases._nested_product_staged_fixed_axis(1),
)
rotated_axes = (
    GaussletBases._nested_product_staged_active_axis(1:2, rotation_axis),
    GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
    GaussletBases._nested_product_staged_fixed_axis(1),
)
axis_function_indices =
    GaussletBases._nested_product_axis_function_indices(3, 1, 2, 2, 2)
identity_unit =
    product_unit(:identity_product_doside_source_box_probe, identity_axes, axis_function_indices)
rotated_unit =
    product_unit(:rotated_product_doside_source_box_probe, rotated_axes, axis_function_indices)
metrics = product_doside_probe_metrics()

pair_plan = CCPM._product_doside_source_box_pair_plan(identity_unit, rotated_unit, metrics)
@assert pair_plan.pair_kind == :product_doside_source_box_pair
@assert pair_plan.left_source_dimensions == (2, 2, 1)
@assert pair_plan.right_source_dimensions == (2, 2, 1)
@assert pair_plan.left_retained_count == 4
@assert pair_plan.right_retained_count == 4
@assert pair_plan.diagnostics.operator_factor_source == :explicit_metric_operator_data
@assert pair_plan.diagnostics.operator_metric_sources ==
        (
            :product_doside_source_box_probe,
            :product_doside_source_box_probe,
            :product_doside_source_box_probe,
        )
@assert pair_plan.diagnostics.input_metric_operator_data ==
        :caller_supplied_explicit_data
@assert !pair_plan.diagnostics.input_metric_operator_data_pgdg_checked
@assert !pair_plan.diagnostics.pgdg_analytic_operator_provenance_claimed
@assert !pair_plan.diagnostics.numerical_reference_fallback
@assert pair_plan.diagnostics.existing_product_staged_retained_helpers_authoritative
@assert !pair_plan.diagnostics.packet_adoption
@assert !pair_plan.diagnostics.fixed_block_routing
@assert !pair_plan.diagnostics.qwhamiltonian_consumes

max_errors = Dict{Tuple{Symbol,Symbol,Symbol},Float64}()
max_explicit_support_errors = Dict{Tuple{Symbol,Symbol,Symbol},Float64}()
local_gaussian_fixture = local_gaussian_fixture_terms()
centered_expansion = centered_local_gaussian_probe_expansion()
centered_axis_layers = centered_local_gaussian_axis_layers()
centered_center = (0.15, -0.2, 0.5)
centered_axis_terms = explicit_centered_axis_terms(
    centered_axis_layers,
    centered_expansion,
    centered_center,
)
max_local_gaussian_errors = Dict{Tuple{Symbol,Symbol},Float64}()
max_centered_local_gaussian_errors = Dict{Tuple{Symbol,Symbol},Float64}()
for (label, left_unit, right_unit) in (
    (:self_identity, identity_unit, identity_unit),
    (:self_rotated, rotated_unit, rotated_unit),
    (:cross_identity_rotated, identity_unit, rotated_unit),
)
    for term in TERMS
        reference = CCPM._product_doside_source_box_reference_block(
            left_unit,
            right_unit,
            metrics;
            term,
        )
        expected = authoritative_block(left_unit, right_unit, metrics, term)
        error = norm(reference.block - expected, Inf)
        max_errors[(label, term, :authoritative)] = error
        @assert error <= 1.0e-12
        explicit_support = explicit_support_reference_block(
            left_unit,
            right_unit,
            metrics,
            term,
        )
        explicit_support_error = norm(reference.block - explicit_support, Inf)
        max_explicit_support_errors[(label, term, :explicit_support)] =
            explicit_support_error
        @assert explicit_support_error <= 1.0e-12
        @assert reference.block_error <= 1.0e-12
        @assert reference.diagnostics.authoritative_block_compared
        @assert reference.diagnostics.operator_factor_source ==
                :explicit_metric_operator_data
        @assert reference.diagnostics.input_metric_operator_data ==
                :caller_supplied_explicit_data
        @assert !reference.diagnostics.input_metric_operator_data_pgdg_checked
        @assert !reference.diagnostics.pgdg_analytic_operator_provenance_claimed
        @assert !reference.diagnostics.numerical_reference_fallback
        @assert !reference.diagnostics.packet_adoption
        @assert !reference.diagnostics.fixed_block_routing
        @assert !reference.diagnostics.qwhamiltonian_consumes
        @assert !reference.diagnostics.public_default_consumes
    end
    local_gaussian = CCPM._product_doside_source_box_local_gaussian_sum_block(
        left_unit,
        right_unit;
        term_coefficients = local_gaussian_fixture.term_coefficients,
        axis_gaussian_terms = local_gaussian_fixture.axis_gaussian_terms,
    )
    explicit_local_gaussian = explicit_local_gaussian_reference_block(
        left_unit,
        right_unit,
        local_gaussian_fixture.term_coefficients,
        local_gaussian_fixture.axis_gaussian_terms,
    )
    local_gaussian_error = norm(local_gaussian.block - explicit_local_gaussian, Inf)
    max_local_gaussian_errors[(label, :gaussian_sum)] = local_gaussian_error
    @assert local_gaussian.path == :product_doside_source_box_local_gaussian_sum
    @assert local_gaussian.pair_plan.pair_kind ==
            :product_doside_source_box_local_gaussian_pair
    @assert local_gaussian.pair_plan.supported_terms == (:gaussian_sum,)
    @assert local_gaussian_error <= 1.0e-12
    @assert local_gaussian.diagnostics.source_box_first
    @assert local_gaussian.diagnostics.local_gaussian_source_box_terms
    @assert local_gaussian.diagnostics.positive_gaussian_sum_convention
    @assert !local_gaussian.diagnostics.nuclear_charge_applied
    @assert !local_gaussian.diagnostics.nuclear_attraction_sign_applied
    @assert local_gaussian.diagnostics.axis_gaussian_terms_source ==
            :caller_supplied_explicit_data
    @assert !local_gaussian.diagnostics.input_local_gaussian_data_pgdg_checked
    @assert !local_gaussian.diagnostics.pgdg_analytic_operator_provenance_claimed
    @assert !local_gaussian.diagnostics.numerical_reference_fallback
    @assert !local_gaussian.diagnostics.shell_row_algorithm
    @assert !local_gaussian.diagnostics.support_local_oracle_used
    @assert !local_gaussian.diagnostics.retained_pqs_weight_division_allowed
    @assert !local_gaussian.diagnostics.ida_mwg_semantics_changed
    @assert !local_gaussian.diagnostics.ecp_terms_implemented
    @assert !local_gaussian.diagnostics.electron_electron_terms_implemented
    @assert !local_gaussian.diagnostics.packet_adoption
    @assert !local_gaussian.diagnostics.fixed_block_routing
    @assert !local_gaussian.diagnostics.qwhamiltonian_consumes
    @assert local_gaussian.diagnostics.output_finite
    centered_local_gaussian =
        CCPM._product_doside_source_box_centered_local_gaussian_sum_block(
            left_unit,
            right_unit,
            centered_axis_layers,
            centered_expansion;
            center = centered_center,
        )
    explicit_centered_local_gaussian =
        CCPM._product_doside_source_box_local_gaussian_sum_block(
            left_unit,
            right_unit;
            term_coefficients = centered_expansion.coefficients,
            axis_gaussian_terms = centered_axis_terms,
        )
    centered_local_gaussian_error =
        norm(centered_local_gaussian.block - explicit_centered_local_gaussian.block, Inf)
    max_centered_local_gaussian_errors[(label, :gaussian_sum)] =
        centered_local_gaussian_error
    @assert centered_local_gaussian.path ==
            :product_doside_source_box_centered_local_gaussian_sum
    @assert centered_local_gaussian_error <= 1.0e-12
    @assert centered_local_gaussian.diagnostics.source_box_first
    @assert centered_local_gaussian.diagnostics.local_gaussian_source_box_terms
    @assert centered_local_gaussian.diagnostics.centered_local_gaussian_terms_generated
    @assert centered_local_gaussian.diagnostics.axis_gaussian_terms_source ==
            :analytic_gaussian_factor_matrices
    @assert centered_local_gaussian.diagnostics.analytic_primitive_backend_required
    @assert centered_local_gaussian.diagnostics.analytic_primitive_backend_checked
    @assert centered_local_gaussian.diagnostics.positive_gaussian_sum_convention
    @assert !centered_local_gaussian.diagnostics.nuclear_charge_applied
    @assert !centered_local_gaussian.diagnostics.nuclear_attraction_sign_applied
    @assert !centered_local_gaussian.diagnostics.numerical_reference_fallback
    @assert !centered_local_gaussian.diagnostics.shell_row_algorithm
    @assert !centered_local_gaussian.diagnostics.support_local_oracle_used
    @assert !centered_local_gaussian.diagnostics.retained_pqs_weight_division_allowed
    @assert !centered_local_gaussian.diagnostics.ida_mwg_semantics_changed
    @assert !centered_local_gaussian.diagnostics.ecp_terms_implemented
    @assert !centered_local_gaussian.diagnostics.electron_electron_terms_implemented
    @assert !centered_local_gaussian.diagnostics.packet_adoption
    @assert !centered_local_gaussian.diagnostics.fixed_block_routing
    @assert !centered_local_gaussian.diagnostics.qwhamiltonian_consumes
    @assert centered_local_gaussian.diagnostics.output_finite
end

shadow = CCPM._product_doside_source_box_shadow_blocks(
    identity_unit,
    rotated_unit,
    metrics;
    terms = TERMS,
)
@assert shadow.retained_dimension == 8
@assert shadow.ranges.left == 1:4
@assert shadow.ranges.right == 5:8
@assert shadow.diagnostics.source_box_shadow_only
@assert shadow.diagnostics.private_shadow_only
@assert shadow.diagnostics.existing_product_staged_retained_helpers_authoritative
@assert shadow.diagnostics.operator_factor_source == :explicit_metric_operator_data
@assert shadow.diagnostics.input_metric_operator_data == :caller_supplied_explicit_data
@assert !shadow.diagnostics.input_metric_operator_data_pgdg_checked
@assert !shadow.diagnostics.pgdg_analytic_operator_provenance_claimed
@assert !shadow.diagnostics.numerical_reference_fallback
@assert !shadow.diagnostics.packet_adoption
@assert !shadow.diagnostics.fixed_block_routing
@assert !shadow.diagnostics.qwhamiltonian_consumes
@assert shadow.diagnostics.max_right_left_transpose_error <= 1.0e-12
for term in TERMS
    block = shadow.blocks[term]
    components = shadow.component_blocks[term]
    @assert all(isfinite, block)
    @assert block[1:4, 1:4] ≈ components.left_left
    @assert block[1:4, 5:8] ≈ components.left_right
    @assert block[5:8, 1:4] ≈ components.right_left
    @assert block[5:8, 5:8] ≈ components.right_right
    @assert components.right_left_transpose_error <= 1.0e-12
end

try
    CCPM._product_doside_source_box_reference_block(
        identity_unit,
        identity_unit,
        metrics;
        term = :weights,
    )
    error("unsupported product/doside source-box term was accepted")
catch err
    err isa ArgumentError || rethrow()
end

try
    CCPM._product_doside_source_box_local_gaussian_sum_block(
        identity_unit,
        identity_unit;
        term_coefficients = [1.0],
        axis_gaussian_terms = local_gaussian_fixture.axis_gaussian_terms,
    )
    error("mismatched local-Gaussian term counts were accepted")
catch err
    err isa ArgumentError || rethrow()
end

println(
    "product/doside source-box probe ",
    (
        terms = TERMS,
        retained_dimensions = (identity = 4, rotated = 4, shadow = shadow.retained_dimension),
        max_error = isempty(max_errors) ? 0.0 : maximum(values(max_errors)),
        max_explicit_support_error = isempty(max_explicit_support_errors) ?
                                     0.0 :
                                     maximum(values(max_explicit_support_errors)),
        max_local_gaussian_error = isempty(max_local_gaussian_errors) ?
                                   0.0 :
                                   maximum(values(max_local_gaussian_errors)),
        max_centered_local_gaussian_error =
            isempty(max_centered_local_gaussian_errors) ?
            0.0 :
            maximum(values(max_centered_local_gaussian_errors)),
        max_transpose_error = shadow.diagnostics.max_right_left_transpose_error,
        numerical_reference_fallback = shadow.diagnostics.numerical_reference_fallback,
    ),
)
end
