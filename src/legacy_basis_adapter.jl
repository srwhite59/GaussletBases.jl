"""
    LegacyAtomicGaussianShell

One filtered shell from the legacy named-basis loader.

Each shell records:

- angular momentum `l`
- primitive exponents / widths
- shell contraction coefficients

This stays atomic-only and centered on one atom. It is not a molecule-wide
placement object.
"""
struct LegacyAtomicGaussianShell
    l::Int
    exponents::Vector{Float64}
    widths::Vector{Float64}
    coefficients::Vector{Float64}
end

"""
    LegacyAtomicGaussianSupplement

Shared atomic named-basis-plus-`lmax` supplement object for the current
ordinary-QW and nested-QW atomic comparison line.

The object keeps two views at once:

- the full legacy shell list up to the requested `lmax`
- the active centered `l = 0` projection used by the present analytic 1D
  supplement route where that still applies

That mirrors the legacy `basisaddname` / `lmaxadd` intent while staying honest
about the current consumer split:

- the one-dimensional hybrid builder still only has a true active centered `s`
  channel and rejects non-`s` shells
- the ordinary-QW and nested-QW atomic paths can now consume the full shell
  metadata through an explicit atomic-centered 3D Cartesian shell route for
  `lmax <= 2`
- the narrower two-center molecular QW routes remain intentionally smaller for
  now and still stop at `lmax <= 1`

The stored active fields:

- `primitive_exponents`
- `primitive_widths`
- `primitive_gaussians`
- `contraction_matrix`
- `widths`
- `gaussians`

are exactly the active centered `s` supplement seen by the present hybrid/QW
code. The broader shell metadata is kept in `shells`.
"""
struct LegacyAtomicGaussianSupplement
    atom::String
    basis_name::String
    basisfile::String
    lmax::Int
    shells::Vector{LegacyAtomicGaussianShell}
    primitive_exponents::Vector{Float64}
    primitive_widths::Vector{Float64}
    primitive_gaussians::Vector{Gaussian}
    contraction_matrix::Matrix{Float64}
    widths::Vector{Float64}
    gaussians::Vector{Gaussian}
    uncontracted::Bool
    max_width::Union{Nothing, Float64}
end

const LegacySGaussianData = LegacyAtomicGaussianSupplement

"""
    LegacyBondAlignedDiatomicGaussianSupplement

Narrow bond-aligned diatomic molecular supplement assembled from one named
atomic basis and two explicit nuclear centers.

This stays deliberately small in scope:

- one homonuclear atomic shell definition loaded by the legacy named-basis path
- two physical centers on a bond-aligned diatomic
- explicit 3D Cartesian `s/p` shell content for the active QW routes

It is not yet a general arbitrary-molecule placement object.
"""
struct LegacyBondAlignedDiatomicGaussianSupplement
    atomic_source::LegacyAtomicGaussianSupplement
    nuclei::Vector{NTuple{3,Float64}}
    max_width::Union{Nothing, Float64}
end

"""
    LegacyBondAlignedHeteronuclearGaussianSupplement

Narrow bond-aligned diatomic molecular supplement assembled from two explicit
named atomic sources and two nuclear centers.

This is still deliberately small in scope:

- one named atomic shell definition for nucleus A
- one named atomic shell definition for nucleus B
- two physical centers on a bond-aligned diatomic
- explicit 3D Cartesian `s/p` shell content for the active ordinary QW route

It is not yet a general arbitrary-molecule placement object.
"""
struct LegacyBondAlignedHeteronuclearGaussianSupplement
    atomic_sources::NTuple{2,LegacyAtomicGaussianSupplement}
    nuclei::Vector{NTuple{3,Float64}}
    max_width::Union{Nothing, Float64}
end

function _legacy_atomic_has_nonseparable_shells(data::LegacyAtomicGaussianSupplement)
    return any(shell -> shell.l > 0, data.shells)
end

struct _AtomicCartesianShellOrbital3D
    label::String
    lx::Int
    ly::Int
    lz::Int
    exponents::Vector{Float64}
    coefficients::Vector{Float64}
    center::NTuple{3,Float64}
    owner_nucleus_index::Int
end

function _AtomicCartesianShellOrbital3D(
    label::AbstractString,
    lx::Integer,
    ly::Integer,
    lz::Integer,
    exponents::AbstractVector{<:Real},
    coefficients::AbstractVector{<:Real},
    center::NTuple{3,<:Real},
)
    return _AtomicCartesianShellOrbital3D(
        String(label),
        Int(lx),
        Int(ly),
        Int(lz),
        Float64.(collect(exponents)),
        Float64.(collect(coefficients)),
        (Float64(center[1]), Float64(center[2]), Float64(center[3])),
        0,
    )
end

struct _AtomicCartesianShellSupplement3D
    source::LegacyAtomicGaussianSupplement
    orbitals::Vector{_AtomicCartesianShellOrbital3D}
end

struct _BondAlignedDiatomicCartesianShellSupplement3D{S}
    source::S
    orbitals::Vector{_AtomicCartesianShellOrbital3D}
end

function _legacy_atomic_shell_contraction_columns(
    shell::LegacyAtomicGaussianShell,
    uncontracted::Bool,
)
    nprimitive = length(shell.exponents)
    if uncontracted
        return [Vector{Float64}(view(Matrix{Float64}(I, nprimitive, nprimitive), :, column)) for column in 1:nprimitive]
    end
    coefficients = normalize(shell.coefficients)
    return [Float64[coefficients...]]
end

function _atomic_cartesian_shell_labels(l::Int)
    l == 0 && return [("s", (0, 0, 0))]
    l == 1 && return [("px", (1, 0, 0)), ("py", (0, 1, 0)), ("pz", (0, 0, 1))]
    l == 2 && return [
        ("dxx", (2, 0, 0)),
        ("dyy", (0, 2, 0)),
        ("dzz", (0, 0, 2)),
        ("dxy", (1, 1, 0)),
        ("dxz", (1, 0, 1)),
        ("dyz", (0, 1, 1)),
    ]
    throw(ArgumentError("atomic Cartesian shell support currently stops at l = 2"))
end

function _atomic_cartesian_shell_supplement_3d(
    data::LegacyAtomicGaussianSupplement,
)
    any(shell -> shell.l > 2, data.shells) && throw(
        ArgumentError("explicit atomic Cartesian shell supplement currently supports only lmax <= 2"),
    )
    center_value =
        isempty(data.primitive_gaussians) ? 0.0 : Float64(data.primitive_gaussians[1].center_value)
    orbitals = _AtomicCartesianShellOrbital3D[]
    shell_counts = Dict{String,Int}()
    for shell in data.shells
        shell_entries = _atomic_cartesian_shell_labels(shell.l)
        contraction_columns = _legacy_atomic_shell_contraction_columns(shell, data.uncontracted)
        for coefficients in contraction_columns
            for (prefix, (lx, ly, lz)) in shell_entries
                shell_counts[prefix] = get(shell_counts, prefix, 0) + 1
                label = string(prefix, shell_counts[prefix])
                push!(
                    orbitals,
                    _AtomicCartesianShellOrbital3D(
                        label,
                        lx,
                        ly,
                        lz,
                        Float64[shell.exponents...],
                        Float64[coefficients...],
                        (center_value, center_value, center_value),
                    ),
                )
            end
        end
    end
    return _AtomicCartesianShellSupplement3D(data, orbitals)
end

function _bond_aligned_diatomic_cartesian_shell_supplement_3d(
    data::LegacyBondAlignedDiatomicGaussianSupplement,
)
    return _BondAlignedDiatomicCartesianShellSupplement3D(
        data,
        _bond_aligned_two_center_cartesian_orbitals(
            (data.atomic_source, data.atomic_source),
            data.nuclei,
            data.max_width,
        ),
    )
end

function _bond_aligned_diatomic_cartesian_shell_supplement_3d(
    data::LegacyBondAlignedHeteronuclearGaussianSupplement,
)
    return _BondAlignedDiatomicCartesianShellSupplement3D(
        data,
        _bond_aligned_two_center_cartesian_orbitals(data.atomic_sources, data.nuclei, data.max_width),
    )
end

function _legacy_width_filtered_shell(
    shell::LegacyAtomicGaussianShell,
    max_width::Union{Nothing, Real},
)
    max_width === nothing && return shell
    width_limit = Float64(max_width)
    keep = findall(width -> width <= width_limit, shell.widths)
    isempty(keep) && return nothing
    return LegacyAtomicGaussianShell(
        shell.l,
        shell.exponents[keep],
        shell.widths[keep],
        shell.coefficients[keep],
    )
end

function _bond_aligned_two_center_cartesian_orbitals(
    atomic_sources::NTuple{2,LegacyAtomicGaussianSupplement},
    nuclei::AbstractVector{<:NTuple{3,<:Real}},
    max_width::Union{Nothing, Real} = nothing,
)
    any(source -> any(shell -> shell.l > 1, source.shells), atomic_sources) && throw(
        ArgumentError("bond-aligned diatomic Cartesian shell supplement currently supports only lmax <= 1"),
    )
    length(nuclei) == 2 || throw(
        ArgumentError("bond-aligned diatomic Cartesian shell supplement currently expects exactly two nuclear centers"),
    )

    orbitals = _AtomicCartesianShellOrbital3D[]
    shell_counts = Dict{Tuple{Int,String},Int}()
    for (nucleus_index, nucleus_raw) in pairs(nuclei)
        source = atomic_sources[nucleus_index]
        nucleus = (
            Float64(nucleus_raw[1]),
            Float64(nucleus_raw[2]),
            Float64(nucleus_raw[3]),
        )
        for source_shell in source.shells
            # Molecular supplements use `max_width` as a core/locality policy:
            # drop diffuse primitives inside a contraction, and drop the shell
            # only if no primitive survives.
            shell = _legacy_width_filtered_shell(source_shell, max_width)
            shell === nothing && continue
            shell_entries = _atomic_cartesian_shell_labels(shell.l)
            contraction_columns = _legacy_atomic_shell_contraction_columns(
                shell,
                source.uncontracted,
            )
            for coefficients in contraction_columns
                for (prefix, (lx, ly, lz)) in shell_entries
                    key = (nucleus_index, prefix)
                    shell_counts[key] = get(shell_counts, key, 0) + 1
                    label = string(nucleus_index == 1 ? "a_" : "b_", prefix, shell_counts[key])
                    push!(
                        orbitals,
                        _AtomicCartesianShellOrbital3D(
                            label,
                            lx,
                            ly,
                            lz,
                            Float64[shell.exponents...],
                            Float64[coefficients...],
                            nucleus,
                            nucleus_index,
                        ),
                    )
                end
            end
        end
    end
    return orbitals
end

function Base.show(io::IO, data::LegacyAtomicGaussianSupplement)
    shell_ls = unique(sort([shell.l for shell in data.shells]))
    print(
        io,
        "LegacyAtomicGaussianSupplement(atom=\"",
        data.atom,
        "\", basis=\"",
        data.basis_name,
        "\", lmax=",
        data.lmax,
        ", shell_ls=",
        shell_ls,
        ", nshells=",
        length(data.shells),
        ", nactive_primitive=",
        length(data.primitive_gaussians),
        ", nactive_contracted=",
        length(data.gaussians),
        ", uncontracted=",
        data.uncontracted,
    )
    if data.max_width !== nothing
        print(io, ", max_width=", data.max_width)
    end
    print(io, ")")
end

function Base.show(io::IO, data::LegacyBondAlignedDiatomicGaussianSupplement)
    print(
        io,
        "LegacyBondAlignedDiatomicGaussianSupplement(atom=\"",
        data.atomic_source.atom,
        "\", basis=\"",
        data.atomic_source.basis_name,
        "\", lmax=",
        data.atomic_source.lmax,
        ", nuclei=",
        data.nuclei,
        ", nshells=",
        length(data.atomic_source.shells),
        ", uncontracted=",
        data.atomic_source.uncontracted,
    )
    if data.max_width !== nothing
        print(io, ", max_width=", data.max_width)
    end
    print(io, ")")
end

function Base.show(io::IO, data::LegacyBondAlignedHeteronuclearGaussianSupplement)
    print(
        io,
        "LegacyBondAlignedHeteronuclearGaussianSupplement(atoms=(",
        "\"",
        data.atomic_sources[1].atom,
        "\", \"",
        data.atomic_sources[2].atom,
        "\"), bases=(",
        "\"",
        data.atomic_sources[1].basis_name,
        "\", \"",
        data.atomic_sources[2].basis_name,
        "\"), lmax=",
        max(data.atomic_sources[1].lmax, data.atomic_sources[2].lmax),
        ", nuclei=",
        data.nuclei,
        ", nshells=(",
        length(data.atomic_sources[1].shells),
        ", ",
        length(data.atomic_sources[2].shells),
        "), uncontracted=(",
        data.atomic_sources[1].uncontracted,
        ", ",
        data.atomic_sources[2].uncontracted,
        ")",
    )
    if data.max_width !== nothing
        print(io, ", max_width=", data.max_width)
    end
    print(io, ")")
end

function _vendored_legacy_basisfile_path()
    return normpath(joinpath(@__DIR__, "..", "data", "legacy", "BasisSets"))
end

function _legacy_basisfile_path(; basisfile::Union{Nothing, AbstractString} = nothing)
    basisfile !== nothing && return String(basisfile)

    env_path = get(ENV, "GAUSSLETBASES_BASISSETS_PATH", "")
    !isempty(env_path) && return env_path

    vendored_path = _vendored_legacy_basisfile_path()
    isfile(vendored_path) && return vendored_path

    return joinpath(homedir(), "BasisSets")
end

function _legacy_basis_angular_momentum(label::AbstractString)
    lowered = lowercase(String(label))
    lowered == "s" && return 0
    lowered == "p" && return 1
    lowered == "d" && return 2
    lowered == "f" && return 3
    lowered == "g" && return 4
    lowered == "h" && return 5
    lowered == "i" && return 6
    lowered == "j" && return 7
    lowered == "k" && return 8
    throw(ArgumentError("unsupported angular-momentum label \"$label\" in legacy basis file"))
end

function _finalize_legacy_shell!(
    shells::Vector{Tuple{Int, Vector{Float64}, Vector{Float64}}},
    current_l::Union{Nothing, Int},
    current_exponents::Vector{Float64},
    current_coefficients::Vector{Float64},
)
    current_l === nothing && return
    push!(shells, (current_l, copy(current_exponents), copy(current_coefficients)))
    empty!(current_exponents)
    empty!(current_coefficients)
    return
end

function _legacy_basis_shells(
    atom::AbstractString,
    basis_name::AbstractString;
    basisfile::Union{Nothing, AbstractString} = nothing,
)
    path = _legacy_basisfile_path(; basisfile = basisfile)
    isfile(path) || throw(ArgumentError("legacy basis file not found at \"$path\""))

    atom_key = lowercase(String(atom))
    basis_key = lowercase(String(basis_name))
    shells = Tuple{Int, Vector{Float64}, Vector{Float64}}[]
    current_l = nothing
    current_exponents = Float64[]
    current_coefficients = Float64[]
    found_block = false
    in_target_block = false

    for line in eachline(path)
        basis_match = match(r"^\#BASIS\s+SET\:\s*([A-Za-z]+)\s+([0-9A-Za-z\-\+]+)\s*$", line)
        if basis_match !== nothing
            if in_target_block
                _finalize_legacy_shell!(shells, current_l, current_exponents, current_coefficients)
                break
            end
            line_atom = lowercase(basis_match.captures[1])
            line_basis = lowercase(basis_match.captures[2])
            if line_atom == atom_key && line_basis == basis_key
                found_block = true
                in_target_block = true
                current_l = nothing
            end
            continue
        end

        in_target_block || continue

        shell_match = match(r"^([A-Za-z]+)\s+([A-Za-z])\s*$", strip(line))
        if shell_match !== nothing
            _finalize_legacy_shell!(shells, current_l, current_exponents, current_coefficients)
            lowercase(shell_match.captures[1]) == atom_key || throw(
                ArgumentError("unexpected atom label inside legacy basis block for \"$atom\""),
            )
            current_l = _legacy_basis_angular_momentum(shell_match.captures[2])
            continue
        end

        if occursin(r"^\s*END\s*$"i, line)
            _finalize_legacy_shell!(shells, current_l, current_exponents, current_coefficients)
            break
        end

        data_match = match(r"^\s*([0-9\.\-\+eE]+)\s+([0-9\.\-\+eE]+)\s*$", line)
        if data_match !== nothing
            current_l === nothing && throw(
                ArgumentError("encountered primitive data before a shell header in legacy basis block for \"$atom\" / \"$basis_name\""),
            )
            push!(current_exponents, parse(Float64, data_match.captures[1]))
            push!(current_coefficients, parse(Float64, data_match.captures[2]))
        end
    end

    found_block || throw(ArgumentError("could not find legacy basis block for atom=\"$atom\", basis=\"$basis_name\" in \"$path\""))
    return shells, path
end

_zeta_to_width(zeta::Real) = inv(sqrt(2.0 * Float64(zeta)))

function _legacy_contracted_gaussian_representatives(
    primitive_gaussians::AbstractVector{<:Gaussian},
    contraction_matrix::AbstractMatrix{<:Real},
)
    primitive_set_1d = PrimitiveSet1D(
        AbstractPrimitiveFunction1D[gaussian for gaussian in primitive_gaussians];
        name = :legacy_supplement_primitives,
    )
    overlap = Matrix{Float64}(overlap_matrix(primitive_set_1d))
    position = Matrix{Float64}(position_matrix(primitive_set_1d))
    x2 = Matrix{Float64}(_x2_matrix(primitive_set_1d))
    coefficients = Matrix{Float64}(contraction_matrix)

    widths = Float64[]
    gaussians = Gaussian[]
    for column in axes(coefficients, 2)
        vector = view(coefficients, :, column)
        norm_value = Float64(dot(vector, overlap * vector))
        norm_value > 1.0e-12 || throw(
            ArgumentError("legacy contracted supplement requires nonzero contracted norm"),
        )
        center_value = Float64(dot(vector, position * vector) / norm_value)
        second_moment = Float64(dot(vector, x2 * vector) / norm_value)
        variance = second_moment - center_value^2
        variance > 1.0e-12 || throw(
            ArgumentError("legacy contracted supplement requires positive contracted variance"),
        )
        width_value = sqrt(2.0 * variance)
        push!(widths, width_value)
        push!(gaussians, Gaussian(center = center_value, width = width_value))
    end
    return widths, gaussians
end

function _legacy_filter_shell(
    shell::Tuple{Int, Vector{Float64}, Vector{Float64}},
    max_width::Union{Nothing, Real},
)
    shell_exponents = Float64[shell[2]...]
    shell_coefficients = Float64[shell[3]...]
    shell_widths = _zeta_to_width.(shell_exponents)
    if max_width !== nothing
        width_limit = Float64(max_width)
        keep = findall(width -> width <= width_limit, shell_widths)
        shell_exponents = shell_exponents[keep]
        shell_coefficients = shell_coefficients[keep]
        shell_widths = shell_widths[keep]
    end
    return shell_exponents, shell_widths, shell_coefficients
end

function _legacy_atomic_shell_data(
    atom::AbstractString,
    basis_name::AbstractString,
    shells_raw::Vector{Tuple{Int, Vector{Float64}, Vector{Float64}}},
    path::AbstractString;
    lmax::Int,
    center::Real,
    uncontracted::Bool,
    max_width::Union{Nothing, Real},
)
    filtered_shells = LegacyAtomicGaussianShell[]
    primitive_exponents = Float64[]
    primitive_widths = Float64[]

    for shell in shells_raw
        shell_l = shell[1]
        shell_l <= lmax || continue
        shell_exponents, shell_widths, shell_coefficients = _legacy_filter_shell(shell, max_width)
        isempty(shell_exponents) && continue
        push!(
            filtered_shells,
            LegacyAtomicGaussianShell(
                shell_l,
                shell_exponents,
                shell_widths,
                shell_coefficients,
            ),
        )
        if shell_l == 0
            append!(primitive_exponents, shell_exponents)
            append!(primitive_widths, shell_widths)
        end
    end

    isempty(filtered_shells) && throw(
        ArgumentError("legacy basis \"$basis_name\" for atom \"$atom\" contains no shells at or below lmax = $lmax after filtering"),
    )
    isempty(primitive_exponents) && throw(
        ArgumentError("current atomic legacy supplement needs at least one centered s shell in \"$basis_name\" for atom \"$atom\""),
    )

    s_shells = filter(shell -> shell.l == 0, filtered_shells)
    contraction_columns = Vector{Vector{Float64}}()
    for shell in s_shells
        if uncontracted
            for index in eachindex(shell.exponents)
                column = zeros(Float64, length(shell.exponents))
                column[index] = 1.0
                push!(contraction_columns, column)
            end
        else
            push!(contraction_columns, normalize(shell.coefficients))
        end
    end

    contraction_matrix = zeros(Float64, length(primitive_exponents), length(contraction_columns))
    primitive_offset = 0
    column_offset = 0
    for shell in s_shells
        nshell = length(shell.exponents)
        if uncontracted
            contraction_matrix[(primitive_offset + 1):(primitive_offset + nshell), (column_offset + 1):(column_offset + nshell)] .= Matrix{Float64}(I, nshell, nshell)
            column_offset += nshell
        else
            contraction_matrix[(primitive_offset + 1):(primitive_offset + nshell), column_offset + 1] .= normalize(shell.coefficients)
            column_offset += 1
        end
        primitive_offset += nshell
    end

    primitive_gaussians = Gaussian[Gaussian(center = center, width = width) for width in primitive_widths]
    widths, gaussians = _legacy_contracted_gaussian_representatives(
        primitive_gaussians,
        contraction_matrix,
    )

    return LegacyAtomicGaussianSupplement(
        String(atom),
        String(basis_name),
        String(path),
        lmax,
        filtered_shells,
        primitive_exponents,
        primitive_widths,
        primitive_gaussians,
        contraction_matrix,
        widths,
        gaussians,
        uncontracted,
        max_width === nothing ? nothing : Float64(max_width),
    )
end

"""
    legacy_atomic_gaussian_supplement(
        atom,
        basis_name;
        lmax = 0,
        basisfile = nothing,
        center = 0.0,
        uncontracted = false,
        max_width = nothing,
    )

Load one atom's legacy named basis in the old `getbasis(...; maxl = lmax)`
style and form the shared atomic supplement object used by the current
ordinary-QW and nested-QW atomic comparison path.

This pass is intentionally atomic-only. It does not place functions on multiple
atoms. It records all shell metadata up to `lmax` and also exposes the centered
`s` projection through the analytic 1D primitive route still used by the
one-dimensional hybrid builder.

The basis-file lookup order is:

1. explicit `basisfile = ...`
2. `GAUSSLETBASES_BASISSETS_PATH`
3. vendored repo copy at `data/legacy/BasisSets`
4. legacy fallback `~/BasisSets`

For atomic ordinary-QW and nested-QW, all shell content up to `lmax = 2`,
including pure `s` shells, is now consumed through an explicit atomic-centered
3D Cartesian shell supplement route. The one-dimensional hybrid builder remains
honestly `s`-only, and the separate two-center molecular shell route remains
narrower for now.
"""
function legacy_atomic_gaussian_supplement(
    atom::AbstractString,
    basis_name::AbstractString;
    lmax::Integer = 0,
    basisfile::Union{Nothing, AbstractString} = nothing,
    center::Real = 0.0,
    uncontracted::Bool = false,
    max_width::Union{Nothing, Real} = nothing,
)
    Int(lmax) >= 0 || throw(ArgumentError("legacy_atomic_gaussian_supplement requires lmax >= 0"))
    shells, path = _legacy_basis_shells(atom, basis_name; basisfile = basisfile)
    return _legacy_atomic_shell_data(
        atom,
        basis_name,
        shells,
        path;
        lmax = Int(lmax),
        center = center,
        uncontracted = uncontracted,
        max_width = max_width,
    )
end

"""
    legacy_bond_aligned_diatomic_gaussian_supplement(
        atom,
        basis_name,
        nuclei;
        lmax = 0,
        basisfile = nothing,
        uncontracted = false,
        max_width = nothing,
    )

Build the first narrow molecular supplement object for a bond-aligned
homonuclear diatomic.

The object reuses one atomic named-basis shell definition and places the
resulting explicit Cartesian shell content on the two supplied nuclear centers.
It is intended for the first honest molecular QW residual-Gaussian completion
and does not yet attempt arbitrary molecular placement.

When `max_width` is supplied, it is a molecular core/locality cutoff: primitives
wider than `max_width` are removed from each contracted shell before placement,
and a shell contributes no supplement orbital only if all of its primitives are
removed.
"""
function legacy_bond_aligned_diatomic_gaussian_supplement(
    atom::AbstractString,
    basis_name::AbstractString,
    nuclei::AbstractVector{<:NTuple{3,<:Real}};
    lmax::Integer = 0,
    basisfile::Union{Nothing,AbstractString} = nothing,
    uncontracted::Bool = false,
    max_width::Union{Nothing,Real} = nothing,
)
    length(nuclei) == 2 || throw(
        ArgumentError("legacy_bond_aligned_diatomic_gaussian_supplement currently expects exactly two nuclear centers"),
    )
    atomic_source = legacy_atomic_gaussian_supplement(
        atom,
        basis_name;
        lmax = lmax,
        basisfile = basisfile,
        center = 0.0,
        uncontracted = uncontracted,
        max_width = nothing,
    )
    return LegacyBondAlignedDiatomicGaussianSupplement(
        atomic_source,
        [(Float64(center[1]), Float64(center[2]), Float64(center[3])) for center in nuclei],
        max_width === nothing ? nothing : Float64(max_width),
    )
end

"""
    legacy_bond_aligned_heteronuclear_gaussian_supplement(
        atom_a,
        basis_name_a,
        atom_b,
        basis_name_b,
        nuclei;
        lmax = 0,
        basisfile = nothing,
        uncontracted = false,
        max_width = nothing,
    )

Build the first narrow mixed-species molecular supplement object for a
bond-aligned heteronuclear diatomic.

The object places one named atomic shell source on nucleus A and one distinct
named atomic shell source on nucleus B. It is intentionally limited to the
first honest two-center heteronuclear ordinary-QW line.

When `max_width` is supplied, it is a molecular core/locality cutoff: primitives
wider than `max_width` are removed from each contracted shell before placement,
and a shell contributes no supplement orbital only if all of its primitives are
removed.
"""
function legacy_bond_aligned_heteronuclear_gaussian_supplement(
    atom_a::AbstractString,
    basis_name_a::AbstractString,
    atom_b::AbstractString,
    basis_name_b::AbstractString,
    nuclei::AbstractVector{<:NTuple{3,<:Real}};
    lmax::Integer = 0,
    basisfile::Union{Nothing,AbstractString} = nothing,
    uncontracted::Bool = false,
    max_width::Union{Nothing,Real} = nothing,
)
    length(nuclei) == 2 || throw(
        ArgumentError("legacy_bond_aligned_heteronuclear_gaussian_supplement currently expects exactly two nuclear centers"),
    )
    source_a = legacy_atomic_gaussian_supplement(
        atom_a,
        basis_name_a;
        lmax = lmax,
        basisfile = basisfile,
        center = 0.0,
        uncontracted = uncontracted,
        max_width = nothing,
    )
    source_b = legacy_atomic_gaussian_supplement(
        atom_b,
        basis_name_b;
        lmax = lmax,
        basisfile = basisfile,
        center = 0.0,
        uncontracted = uncontracted,
        max_width = nothing,
    )
    return LegacyBondAlignedHeteronuclearGaussianSupplement(
        (source_a, source_b),
        [(Float64(center[1]), Float64(center[2]), Float64(center[3])) for center in nuclei],
        max_width === nothing ? nothing : Float64(max_width),
    )
end

"""
    legacy_s_gaussian_data(atom, basis_name; kwargs...)

Thin compatibility wrapper for the earlier He-`s` supplement helper.

This now forwards to `legacy_atomic_gaussian_supplement(...; lmax = 0)`.
"""
function legacy_s_gaussian_data(
    atom::AbstractString,
    basis_name::AbstractString;
    basisfile::Union{Nothing, AbstractString} = nothing,
    center::Real = 0.0,
    uncontracted::Bool = false,
    max_width::Union{Nothing, Real} = nothing,
)
    return legacy_atomic_gaussian_supplement(
        atom,
        basis_name;
        lmax = 0,
        basisfile = basisfile,
        center = center,
        uncontracted = uncontracted,
        max_width = max_width,
    )
end
