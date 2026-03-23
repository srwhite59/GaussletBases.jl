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
  supplement route

That mirrors the legacy `basisaddname` / `lmaxadd` intent while staying honest
about the current consumer model: the current separable Cartesian residual path
still consumes only the centered `s` channel analytically.

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

function _legacy_basisfile_path(; basisfile::Union{Nothing, AbstractString} = nothing)
    basisfile !== nothing && return String(basisfile)
    return get(ENV, "GAUSSLETBASES_BASISSETS_PATH", joinpath(homedir(), "BasisSets"))
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
atoms and it does not yet turn `l > 0` shells into a richer nonseparable
Cartesian supplement. Instead, it records all shell metadata up to `lmax` and
exposes the centered `s` projection through the same analytic primitive route
the current hybrid/QW code already consumes.
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
