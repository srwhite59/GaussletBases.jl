function _one_center_atomic_shell_increment(nside::Int)
    nside >= 3 || throw(ArgumentError("one-center atomic shell contract requires nside >= 3"))
    return nside^3 - (nside - 2)^3
end

function _one_center_atomic_complete_shell_retention(nside::Int)
    return _nested_complete_shell_retention_from_nside(nside)
end

function _one_center_atomic_shell_layer_count(working_box_side_count::Int, nside::Int)
    working_box_side_count >= nside || throw(
        ArgumentError("one-center atomic structure diagnostics require working_box_side_count >= nside"),
    )
    current_side = working_box_side_count
    nlayers = 0
    while current_side > nside
        current_side -= 2
        nlayers += 1
    end
    return nlayers, current_side
end

function _one_center_atomic_legacy_profile_working_box(
    parent_side_count::Int,
    working_box::UnitRange{Int},
)
    return _one_center_atomic_legacy_profile_working_box(
        parent_side_count,
        (working_box, working_box, working_box),
    )
end

function _one_center_atomic_legacy_profile_working_box(
    parent_side_count::Int,
    working_box::NTuple{3,UnitRange{Int}},
)
    expected_parent = 1:parent_side_count
    for interval in working_box
        interval == intersect(interval, expected_parent) || throw(
            ArgumentError("one-center atomic legacy profile working box must lie inside 1:$parent_side_count"),
        )
    end
    working_box_sides = Tuple(length.(working_box))
    (working_box_sides[1] == working_box_sides[2] && working_box_sides[2] == working_box_sides[3]) || throw(
        ArgumentError("one-center atomic legacy profile requires a cubic working box"),
    )
    return working_box
end

function _nested_same_term_exponents(
    left::AbstractVector{<:Real},
    right::AbstractVector{<:Real};
    atol::Float64 = 1.0e-12,
)
    length(left) == length(right) || return false
    isempty(left) && return true
    return maximum(abs.(Float64.(left) .- Float64.(right))) <= atol
end

function _one_center_atomic_term_coefficients(
    bundle::_MappedOrdinaryGausslet1DBundle,
    expansion::CoulombGaussianExpansion,
)
    _nested_same_term_exponents(bundle.exponents, expansion.exponents) || throw(
        ArgumentError("compact one-center atomic fixed-block storage requires bundle Gaussian exponents to match the supplied CoulombGaussianExpansion"),
    )
    return Float64[Float64(value) for value in expansion.coefficients]
end

struct OneCenterAtomicNestedLayerStructure
    layer_index::Int
    face_retained_count::Int
    edge_retained_count::Int
    corner_retained_count::Int
    retained_dimension::Int
    provenance::_CartesianNestedShellLayerProvenance3D
end

struct OneCenterAtomicNestedStructureDiagnostics
    parent_side_count::Int
    working_box_side_count::Int
    nside::Int
    core_side_count::Int
    shell_layer_count::Int
    expected_shell_increment::Int
    expected_face_retained_count::Int
    expected_edge_retained_count::Int
    expected_corner_retained_count::Int
    layer_structures::Vector{OneCenterAtomicNestedLayerStructure}
    total_face_retained_count::Int
    total_edge_retained_count::Int
    total_corner_retained_count::Int
    total_expected_gausslet_count::Int
    total_actual_gausslet_count::Int
    layers_match_expected::Bool
end

function _build_one_center_atomic_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle,
    working_box::NTuple{3,UnitRange{Int}};
    nside::Int,
    packet_kernel::Symbol = :factorized_direct,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    bundles = _CartesianNestedAxisBundles3D(bundle, bundle, bundle)
    retention = _one_center_atomic_complete_shell_retention(nside)
    return _nested_complete_shell_sequence_for_box(
        bundles,
        working_box;
        nside = nside,
        retain_xy = retention.retain_xy,
        retain_xz = retention.retain_xz,
        retain_yz = retention.retain_yz,
        retain_x_edge = retention.retain_x_edge,
        retain_y_edge = retention.retain_y_edge,
        retain_z_edge = retention.retain_z_edge,
        packet_kernel = packet_kernel,
        term_coefficients = term_coefficients,
    )
end

"""
    build_one_center_atomic_full_parent_shell_sequence(
        bundle::_MappedOrdinaryGausslet1DBundle;
        nside,
    )

Build the canonical one-center atomic nested shell sequence on the full parent
cube of the supplied mapped 1D gausslet bundle.

This helper is the supported atomic one-center backbone:

- it always uses full parent coverage
- it always uses `working_box = (1:n, 1:n, 1:n)`
- it peels complete shells until the direct inner cube reaches `nside`
- it uses the legacy/W&L complete-shell contract
  - shell increment `= nside^3 - (nside - 2)^3`
  - faces retain `(nside - 2) × (nside - 2)`
  - edges retain `nside - 2`
  - corners are carried directly

It should be used in place of any older central-box atomic diagnostic fixture.
"""
function build_one_center_atomic_full_parent_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    nside::Int,
)
    n = length(bundle.basis)
    term_coefficients = _one_center_atomic_term_coefficients(bundle, expansion)
    return _build_one_center_atomic_shell_sequence(
        bundle,
        (1:n, 1:n, 1:n);
        nside = nside,
        term_coefficients = term_coefficients,
    )
end

function build_one_center_atomic_full_parent_shell_sequence(
    basis::MappedUniformBasis;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    exponents::AbstractVector{<:Real} = expansion.exponents,
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    kwargs...,
)
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = exponents,
        center = center,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    return build_one_center_atomic_full_parent_shell_sequence(
        bundle;
        expansion = expansion,
        kwargs...,
    )
end

"""
    build_one_center_atomic_legacy_profile_shell_sequence(
        bundle::_MappedOrdinaryGausslet1DBundle;
        working_box,
        nside,
    )

Build the explicit legacy-profile one-center atomic nested shell sequence on a
chosen inner working box of the supplied parent lattice.

This helper is intentionally separate from the modern canonical full-parent
path. It keeps the same exact complete-shell retention contract:

- shell increment `= nside^3 - (nside - 2)^3`
- faces retain `(nside - 2) × (nside - 2)`
- edges retain `nside - 2`
- corners are carried directly

but applies it on the explicit inner working box supplied by the caller, for
example `(2:28, 2:28, 2:28)` on a `29^3` parent lattice.
"""
function build_one_center_atomic_legacy_profile_shell_sequence(
    bundle::_MappedOrdinaryGausslet1DBundle;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    working_box::Union{UnitRange{Int},NTuple{3,UnitRange{Int}}},
    nside::Int,
)
    n = length(bundle.basis)
    normalized_working_box = _one_center_atomic_legacy_profile_working_box(n, working_box)
    term_coefficients = _one_center_atomic_term_coefficients(bundle, expansion)
    return _build_one_center_atomic_shell_sequence(
        bundle,
        normalized_working_box;
        nside = nside,
        term_coefficients = term_coefficients,
    )
end

function build_one_center_atomic_legacy_profile_shell_sequence(
    basis::MappedUniformBasis;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    exponents::AbstractVector{<:Real} = expansion.exponents,
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    kwargs...,
)
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = exponents,
        center = center,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    return build_one_center_atomic_legacy_profile_shell_sequence(
        bundle;
        expansion = expansion,
        kwargs...,
    )
end

"""
    one_center_atomic_full_parent_fixed_block(
        bundle::_MappedOrdinaryGausslet1DBundle;
        nside,
        kwargs...,
    )

Build the canonical one-center atomic nested fixed block on the full parent
cube of the supplied mapped 1D gausslet bundle.

This is a thin convenience wrapper around
[`build_one_center_atomic_full_parent_shell_sequence`](@ref) followed by
`_nested_fixed_block(...)`.
"""
function one_center_atomic_full_parent_fixed_block(
    bundle::_MappedOrdinaryGausslet1DBundle;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    timing::Union{Bool,Symbol} = false,
    timing_io::IO = stdout,
    kwargs...,
)
    return _nested_capture_timeg_report(timing, timing_io) do
        term_coefficients = _one_center_atomic_term_coefficients(bundle, expansion)
        sequence = @timeg "fixed_block.sequence_build" begin
            n = length(bundle.basis)
            _build_one_center_atomic_shell_sequence(
                bundle,
                (1:n, 1:n, 1:n);
                term_coefficients = term_coefficients,
                kwargs...,
            )
        end
        return @timeg "fixed_block.adapter" begin
            _nested_fixed_block(sequence, bundle)
        end
    end
end

function one_center_atomic_full_parent_fixed_block(
    basis::MappedUniformBasis;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    exponents::AbstractVector{<:Real} = expansion.exponents,
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    timing::Union{Bool,Symbol} = false,
    timing_io::IO = stdout,
    kwargs...,
)
    return _nested_capture_timeg_report(timing, timing_io) do
        bundle = @timeg "fixed_block.parent_bundle" begin
            _mapped_ordinary_gausslet_1d_bundle(
                basis;
                exponents = isempty(exponents) ? expansion.exponents : exponents,
                center = center,
                backend = gausslet_backend,
                refinement_levels = refinement_levels,
            )
        end
        term_coefficients = _one_center_atomic_term_coefficients(bundle, expansion)
        sequence = @timeg "fixed_block.sequence_build" begin
            n = length(bundle.basis)
            _build_one_center_atomic_shell_sequence(
                bundle,
                (1:n, 1:n, 1:n);
                term_coefficients = term_coefficients,
                kwargs...,
            )
        end
        return @timeg "fixed_block.adapter" begin
            _nested_fixed_block(sequence, bundle)
        end
    end
end

"""
    one_center_atomic_legacy_profile_fixed_block(
        bundle::_MappedOrdinaryGausslet1DBundle;
        working_box,
        nside,
        kwargs...,
    )

Build the legacy-profile one-center atomic nested fixed block on an explicit
inner working box of the supplied parent lattice.

This is a thin convenience wrapper around
[`build_one_center_atomic_legacy_profile_shell_sequence`](@ref) followed by
`_nested_fixed_block(...)`.
"""
function one_center_atomic_legacy_profile_fixed_block(
    bundle::_MappedOrdinaryGausslet1DBundle;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    timing::Union{Bool,Symbol} = false,
    timing_io::IO = stdout,
    working_box::Union{UnitRange{Int},NTuple{3,UnitRange{Int}}},
    nside::Int,
    kwargs...,
)
    return _nested_capture_timeg_report(timing, timing_io) do
        term_coefficients = _one_center_atomic_term_coefficients(bundle, expansion)
        normalized_working_box =
            _one_center_atomic_legacy_profile_working_box(length(bundle.basis), working_box)
        sequence = @timeg "fixed_block.sequence_build" begin
            _build_one_center_atomic_shell_sequence(
                bundle,
                normalized_working_box;
                nside = nside,
                term_coefficients = term_coefficients,
                kwargs...,
            )
        end
        return @timeg "fixed_block.adapter" begin
            _nested_fixed_block(sequence, bundle)
        end
    end
end

function one_center_atomic_legacy_profile_fixed_block(
    basis::MappedUniformBasis;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    exponents::AbstractVector{<:Real} = expansion.exponents,
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    timing::Union{Bool,Symbol} = false,
    timing_io::IO = stdout,
    working_box::Union{UnitRange{Int},NTuple{3,UnitRange{Int}}},
    nside::Int,
    kwargs...,
)
    return _nested_capture_timeg_report(timing, timing_io) do
        bundle = @timeg "fixed_block.parent_bundle" begin
            _mapped_ordinary_gausslet_1d_bundle(
                basis;
                exponents = isempty(exponents) ? expansion.exponents : exponents,
                center = center,
                backend = gausslet_backend,
                refinement_levels = refinement_levels,
            )
        end
        term_coefficients = _one_center_atomic_term_coefficients(bundle, expansion)
        normalized_working_box =
            _one_center_atomic_legacy_profile_working_box(length(bundle.basis), working_box)
        sequence = @timeg "fixed_block.sequence_build" begin
            _build_one_center_atomic_shell_sequence(
                bundle,
                normalized_working_box;
                nside = nside,
                term_coefficients = term_coefficients,
                kwargs...,
            )
        end
        return @timeg "fixed_block.adapter" begin
            _nested_fixed_block(sequence, bundle)
        end
    end
end

function _one_center_atomic_nested_layer_structure(
    shell::_CartesianNestedCompleteShell3D,
    layer_index::Int,
)
    face_retained_count = sum(length, shell.face_column_ranges)
    edge_retained_count = sum(length, shell.edge_column_ranges)
    corner_retained_count = sum(length, shell.corner_column_ranges)
    return OneCenterAtomicNestedLayerStructure(
        layer_index,
        face_retained_count,
        edge_retained_count,
        corner_retained_count,
        size(shell.coefficient_matrix, 2),
        shell.provenance,
    )
end

function _one_center_atomic_layer_provenance(
    working_box::NTuple{3,UnitRange{Int}},
    layer_index::Int,
    retained_dimension::Int,
)
    offset = layer_index - 1
    source_box = ntuple(3) do axis
        (first(working_box[axis]) + offset):(last(working_box[axis]) - offset)
    end
    next_inner_box = ntuple(3) do axis
        (first(source_box[axis]) + 1):(last(source_box[axis]) - 1)
    end
    return _CartesianNestedShellLayerProvenance3D(
        source_box,
        next_inner_box,
        prod(length.(source_box)) - prod(length.(next_inner_box)),
        retained_dimension,
    )
end

function one_center_atomic_nested_structure_diagnostics(
    parent_side_count::Int;
    nside::Int,
    working_box_side_count::Int = parent_side_count,
)
    nlayers, core_side_count = _one_center_atomic_shell_layer_count(working_box_side_count, nside)
    retention = _one_center_atomic_complete_shell_retention(nside)
    canonical_working_box = ntuple(_ -> 1:working_box_side_count, 3)
    layer_structures = OneCenterAtomicNestedLayerStructure[
        OneCenterAtomicNestedLayerStructure(
            layer,
            retention.face_retained_count,
            retention.edge_retained_count,
            retention.corner_retained_count,
            retention.shell_increment,
            _one_center_atomic_layer_provenance(
                canonical_working_box,
                layer,
                retention.shell_increment,
            ),
        ) for layer in 1:nlayers
    ]
    total_face_retained_count = nlayers * retention.face_retained_count
    total_edge_retained_count = nlayers * retention.edge_retained_count
    total_corner_retained_count = nlayers * retention.corner_retained_count
    total_expected = core_side_count^3 + nlayers * retention.shell_increment
    return OneCenterAtomicNestedStructureDiagnostics(
        parent_side_count,
        working_box_side_count,
        nside,
        core_side_count,
        nlayers,
        retention.shell_increment,
        retention.face_retained_count,
        retention.edge_retained_count,
        retention.corner_retained_count,
        layer_structures,
        total_face_retained_count,
        total_edge_retained_count,
        total_corner_retained_count,
        total_expected,
        total_expected,
        true,
    )
end

function one_center_atomic_nested_structure_diagnostics(
    sequence::_CartesianNestedShellSequence3D;
    parent_side_count::Int,
    nside::Int,
)
    working_box_sides = Tuple(length.(sequence.working_box))
    working_box_sides[1] == working_box_sides[2] == working_box_sides[3] || throw(
        ArgumentError("one-center atomic structure diagnostics require a cubic working box"),
    )
    core_side_count = round(Int, cbrt(length(sequence.core_indices)))
    layer_structures = OneCenterAtomicNestedLayerStructure[]
    for (layer_index, shell) in enumerate(sequence.shell_layers)
        shell isa _CartesianNestedCompleteShell3D || throw(
            ArgumentError("one-center atomic structure diagnostics require complete shell layers"),
        )
        push!(layer_structures, _one_center_atomic_nested_layer_structure(shell, layer_index))
    end
    total_face_retained_count = sum(layer.face_retained_count for layer in layer_structures)
    total_edge_retained_count = sum(layer.edge_retained_count for layer in layer_structures)
    total_corner_retained_count = sum(layer.corner_retained_count for layer in layer_structures)
    retention = _one_center_atomic_complete_shell_retention(nside)
    total_expected = core_side_count^3 + length(layer_structures) * retention.shell_increment
    return OneCenterAtomicNestedStructureDiagnostics(
        parent_side_count,
        working_box_sides[1],
        nside,
        core_side_count,
        length(layer_structures),
        retention.shell_increment,
        retention.face_retained_count,
        retention.edge_retained_count,
        retention.corner_retained_count,
        layer_structures,
        total_face_retained_count,
        total_edge_retained_count,
        total_corner_retained_count,
        total_expected,
        size(sequence.coefficient_matrix, 2),
        all(
            layer.face_retained_count == retention.face_retained_count &&
            layer.edge_retained_count == retention.edge_retained_count &&
            layer.corner_retained_count == retention.corner_retained_count &&
            layer.retained_dimension == retention.shell_increment
            for layer in layer_structures
        ),
    )
end

function one_center_atomic_nested_structure_diagnostics(
    fixed_block::_NestedFixedBlock3D;
    nside::Int,
)
    parent_side_count = length(fixed_block.parent_basis)
    return one_center_atomic_nested_structure_diagnostics(
        fixed_block.shell;
        parent_side_count = parent_side_count,
        nside = nside,
    )
end

function one_center_atomic_nested_structure_diagnostics(
    bundle::_MappedOrdinaryGausslet1DBundle;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    nside::Int,
)
    sequence = build_one_center_atomic_full_parent_shell_sequence(
        bundle;
        expansion = expansion,
        nside = nside,
    )
    return one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = length(bundle.basis),
        nside = nside,
    )
end

function one_center_atomic_nested_structure_diagnostics(
    basis::MappedUniformBasis;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    exponents::AbstractVector{<:Real} = expansion.exponents,
    center::Real = 0.0,
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    nside::Int,
)
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = exponents,
        center = center,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    return one_center_atomic_nested_structure_diagnostics(
        bundle;
        expansion = expansion,
        nside = nside,
    )
end

function one_center_atomic_nested_structure_report(
    diagnostics::OneCenterAtomicNestedStructureDiagnostics;
    supplement_orbital_count::Union{Nothing,Int} = nothing,
    total_expected_basis_count::Union{Nothing,Int} = nothing,
    total_actual_basis_count::Union{Nothing,Int} = nothing,
    low_one_body_eigenvalues::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    common_contract = _nested_glass_box_contract(diagnostics)
    lines = String[
        "parent_side_count = $(diagnostics.parent_side_count)",
        "working_box_side_count = $(diagnostics.working_box_side_count)",
        "nside = $(diagnostics.nside)",
        "core_side_count = $(diagnostics.core_side_count)",
        "shell_layers = $(diagnostics.shell_layer_count)",
        "expected_shell_increment = $(diagnostics.expected_shell_increment)",
        "expected_face_retained_count = $(diagnostics.expected_face_retained_count)",
        "expected_edge_retained_count = $(diagnostics.expected_edge_retained_count)",
        "expected_corner_retained_count = $(diagnostics.expected_corner_retained_count)",
        "total_face_retained_count = $(diagnostics.total_face_retained_count)",
        "total_edge_retained_count = $(diagnostics.total_edge_retained_count)",
        "total_corner_retained_count = $(diagnostics.total_corner_retained_count)",
        "total_expected_gausslet_count = $(diagnostics.total_expected_gausslet_count)",
        "total_actual_gausslet_count = $(common_contract.fixed_dimension)",
        "layers_match_expected = $(diagnostics.layers_match_expected)",
    ]
    if !isnothing(supplement_orbital_count)
        push!(lines, "supplement_orbital_count = $(supplement_orbital_count)")
    end
    if !isnothing(total_expected_basis_count)
        push!(lines, "total_expected_basis_count = $(total_expected_basis_count)")
    end
    if !isnothing(total_actual_basis_count)
        push!(lines, "total_actual_basis_count = $(total_actual_basis_count)")
    end
    if !isnothing(low_one_body_eigenvalues)
        push!(lines, "low_one_body_eigenvalues = $(repr(Float64[low_one_body_eigenvalues...]))")
    end
    for (index, layer) in enumerate(diagnostics.layer_structures)
        push!(lines, "layer_$(layer.layer_index)_faces = $(layer.face_retained_count)")
        push!(lines, "layer_$(layer.layer_index)_edges = $(layer.edge_retained_count)")
        push!(lines, "layer_$(layer.layer_index)_corners = $(layer.corner_retained_count)")
        push!(
            lines,
            "layer_$(layer.layer_index)_retained_dimension = $(common_contract.layer_dimensions[index])",
        )
        push!(lines, "layer_$(layer.layer_index)_source_box = $(layer.provenance.source_box)")
        push!(
            lines,
            "layer_$(layer.layer_index)_next_inner_box = $(layer.provenance.next_inner_box)",
        )
        push!(
            lines,
            "layer_$(layer.layer_index)_source_point_count = $(layer.provenance.source_point_count)",
        )
    end
    return join(lines, "\n")
end
