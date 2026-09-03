"""
    CartesianBasisMetadata3D

Minimal public metadata bundle identifying a three-dimensional Cartesian basis.

This is the Cartesian analogue of `BasisMetadata1D`. It records:

- the high-level basis kind, such as `:direct_product`, `:nested_fixed_block`,
  or `:hybrid_residual`
- the axis-level one-dimensional metadata reused from the existing 1D contract
- the defining parent-space contract and final basis size
- light basis-level geometric data such as labels, centers, and working-box
  metadata when that construction metadata already exists

The object is basis-defining metadata only. It is not an operator cache.
"""
struct CartesianBasisMetadata3D{AT <: NamedTuple, RT <: NamedTuple}
    basis_kind::Symbol
    axis_sharing::Symbol
    axis_metadata::AT
    parent_kind::Symbol
    parent_axis_counts::NTuple{3,Int}
    parent_dimension::Int
    final_dimension::Int
    working_box::Union{Nothing,NTuple{3,UnitRange{Int}}}
    basis_labels::Vector{String}
    basis_centers::Matrix{Float64}
    route_metadata::RT
end

function Base.show(io::IO, metadata::CartesianBasisMetadata3D)
    print(
        io,
        "CartesianBasisMetadata3D(kind=:",
        metadata.basis_kind,
        ", nfinal=",
        metadata.final_dimension,
        ", parent_kind=:",
        metadata.parent_kind,
        ", axis_sharing=:",
        metadata.axis_sharing,
        ")",
    )
end

"""
    CartesianBasisRepresentation3D

Public in-memory three-dimensional Cartesian basis representation.

The representation stores:

- `metadata`
- explicit 1D axis representations reused from the existing 1D contract
- the contraction from the stored parent/raw basis to the final basis
- parent-space labels and centers
- support-space indexing data when the construction is support-local
- optional parent-data sidecars needed to make mixed raw-space contracts
  explicit, for example the Cartesian parent representation and supplement
  orbital metadata in the QW residual-Gaussian route

This first-pass contract is representation-only. It deliberately avoids any
operator cache beyond the basis-defining contraction data itself.
"""
struct CartesianBasisRepresentation3D{
    MT <: CartesianBasisMetadata3D,
    AT <: NamedTuple,
    PT <: NamedTuple,
}
    metadata::MT
    axis_representations::AT
    contraction_kind::Symbol
    coefficient_matrix::Union{Nothing,_CartesianCoefficientMap}
    parent_labels::Vector{String}
    parent_centers::Matrix{Float64}
    support_indices::Union{Nothing,Vector{Int}}
    support_states::Union{Nothing,Vector{NTuple{3,Int}}}
    parent_data::PT
end

function Base.show(io::IO, representation::CartesianBasisRepresentation3D)
    print(
        io,
        "CartesianBasisRepresentation3D(kind=:",
        representation.metadata.basis_kind,
        ", nfinal=",
        representation.metadata.final_dimension,
        ", parent_kind=:",
        representation.metadata.parent_kind,
        ", contraction=:",
        representation.contraction_kind,
        ")",
    )
end

basis_metadata(representation::CartesianBasisRepresentation3D) = representation.metadata
basis_representation(representation::CartesianBasisRepresentation3D) = representation

"""
    CartesianGaussianShellOrbitalRepresentation3D

Public basis-defining representation of one contracted 3D Cartesian Gaussian
shell orbital used on the supplement side of the hybrid residual-Gaussian
routes.

The stored data is exact enough for overlap:

- `label`
- `angular_powers`
- `center`
- `exponents`
- `coefficients`
- `primitive_normalization`

The current normalization contract is
`:axiswise_normalized_cartesian_gaussian`, meaning each primitive uses the same
axiswise normalized Cartesian-Gaussian prefactor already used in the live
QW/nested residual construction.
"""
struct CartesianGaussianShellOrbitalRepresentation3D
    label::String
    angular_powers::NTuple{3,Int}
    center::NTuple{3,Float64}
    exponents::Vector{Float64}
    coefficients::Vector{Float64}
    primitive_normalization::Symbol
end

function Base.show(io::IO, orbital::CartesianGaussianShellOrbitalRepresentation3D)
    print(
        io,
        "CartesianGaussianShellOrbitalRepresentation3D(label=",
        repr(orbital.label),
        ", l=(",
        orbital.angular_powers[1],
        ",",
        orbital.angular_powers[2],
        ",",
        orbital.angular_powers[3],
        "), nprimitive=",
        length(orbital.exponents),
        ")",
    )
end

"""
    CartesianGaussianShellSupplementRepresentation3D

Public representation of the explicit supplement-orbital sector used in the
hybrid residual-Gaussian Cartesian routes.

It stores:

- `supplement_kind`
- explicit orbital representations
- light source metadata such as basis name, `lmax`, and nuclei when that data
  already exists on the legacy supplement object
"""
struct CartesianGaussianShellSupplementRepresentation3D{MT <: NamedTuple}
    supplement_kind::Symbol
    orbitals::Vector{CartesianGaussianShellOrbitalRepresentation3D}
    metadata::MT
end

function Base.show(io::IO, supplement::CartesianGaussianShellSupplementRepresentation3D)
    print(
        io,
        "CartesianGaussianShellSupplementRepresentation3D(kind=:",
        supplement.supplement_kind,
        ", norbitals=",
        length(supplement.orbitals),
        ")",
    )
end

basis_representation(supplement::CartesianGaussianShellSupplementRepresentation3D) = supplement
basis_metadata(supplement::CartesianGaussianShellSupplementRepresentation3D) = supplement.metadata
