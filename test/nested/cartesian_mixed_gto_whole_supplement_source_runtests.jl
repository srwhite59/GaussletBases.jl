# Runtime role: whole-supplement source adapter fingerprint for CPB-local GTO work.
#
# This test validates existing legacy GTO supplement sources and their modern
# CartesianGaussianShellSupplementRepresentation3D conversion surface. It does
# not build multi-orbital CPB blocks, nuclear/Coulomb GTO terms, WL/PQS
# realization, route/global placement, Hamiltonian assembly, exports, artifacts,
# or IDA/MWG/PQS semantics.

using Test
using LinearAlgebra
using GaussletBases

function _check_whole_supplement_common(
    representation;
    expected_kind,
    expected_source_kind,
    expected_lmax,
    expected_nuclei,
    expected_uncontracted,
)
    @test representation isa CartesianGaussianShellSupplementRepresentation3D
    @test representation.supplement_kind == expected_kind
    @test representation.metadata.source_kind == expected_source_kind
    @test representation.metadata.lmax == expected_lmax
    @test representation.metadata.nuclei == expected_nuclei
    @test representation.metadata.uncontracted == expected_uncontracted
    @test representation.metadata.max_width == nothing
    @test !isempty(representation.orbitals)
    for orbital in representation.orbitals
        @test orbital.label isa String
        @test length(orbital.angular_powers) == 3
        @test all(power -> power >= 0, orbital.angular_powers)
        @test length(orbital.center) == 3
        @test !isempty(orbital.exponents)
        @test length(orbital.exponents) == length(orbital.coefficients)
        @test all(exponent -> isfinite(exponent) && exponent > 0.0, orbital.exponents)
        @test all(isfinite, orbital.coefficients)
        @test orbital.primitive_normalization ==
              :axiswise_normalized_cartesian_gaussian
    end
    return nothing
end

function _labels(representation)
    return [orbital.label for orbital in representation.orbitals]
end

@testset "Whole mixed GTO supplement source adapter" begin
    nuclei = [(0.0, 0.0, -0.7), (0.0, 0.0, 0.7)]

    atomic = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
    @test atomic isa LegacyAtomicGaussianSupplement
    @test atomic.atom == "He"
    @test atomic.basis_name == "cc-pVTZ"
    @test atomic.lmax == 1
    @test atomic.uncontracted == false
    @test endswith(atomic.basisfile, "data/legacy/BasisSets")

    atomic_representation = basis_representation(atomic)
    _check_whole_supplement_common(
        atomic_representation;
        expected_kind = :atomic_cartesian_shell,
        expected_source_kind = :legacy_atomic_gaussian_supplement,
        expected_lmax = 1,
        expected_nuclei = [(0.0, 0.0, 0.0)],
        expected_uncontracted = false,
    )
    @test atomic_representation.metadata.atom == "He"
    @test atomic_representation.metadata.basis_name == "cc-pVTZ"
    @test endswith(atomic_representation.metadata.basisfile, "data/legacy/BasisSets")
    @test length(atomic.shells) == 5
    @test length(atomic_representation.orbitals) == 9
    @test _labels(atomic_representation) == [
        "s1", "s2", "s3",
        "px1", "py1", "pz1",
        "px2", "py2", "pz2",
    ]
    @test atomic_representation.orbitals[1].angular_powers == (0, 0, 0)
    @test atomic_representation.orbitals[1].center == (0.0, 0.0, 0.0)
    @test isapprox(
        atomic_representation.orbitals[1].exponents,
        [234.0, 35.16, 7.989, 2.212];
        atol = 0.0,
        rtol = 0.0,
    )
    @test isapprox(
        atomic_representation.orbitals[1].coefficients,
        [0.002587, 0.019533, 0.090998, 0.27205] /
        norm([0.002587, 0.019533, 0.090998, 0.27205]);
        atol = 1.0e-14,
        rtol = 1.0e-14,
    )
    @test atomic_representation.orbitals[end].label == "pz2"
    @test atomic_representation.orbitals[end].angular_powers == (0, 0, 1)

    uncontracted_atomic =
        legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1, uncontracted = true)
    @test uncontracted_atomic isa LegacyAtomicGaussianSupplement
    uncontracted_representation = basis_representation(uncontracted_atomic)
    _check_whole_supplement_common(
        uncontracted_representation;
        expected_kind = :atomic_cartesian_shell,
        expected_source_kind = :legacy_atomic_gaussian_supplement,
        expected_lmax = 1,
        expected_nuclei = [(0.0, 0.0, 0.0)],
        expected_uncontracted = true,
    )
    @test length(uncontracted_representation.orbitals) == 12
    @test count(orbital -> orbital.angular_powers == (0, 0, 0), uncontracted_representation.orbitals) == 6
    @test count(sum(orbital.angular_powers) == 1 for orbital in uncontracted_representation.orbitals) == 6
    @test uncontracted_representation.orbitals[1].coefficients == [1.0, 0.0, 0.0, 0.0]

    diatomic = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H",
        "cc-pVTZ",
        nuclei;
        lmax = 1,
    )
    @test diatomic isa LegacyBondAlignedDiatomicGaussianSupplement
    @test diatomic.atomic_source.atom == "H"
    @test diatomic.atomic_source.basis_name == "cc-pVTZ"
    diatomic_representation = basis_representation(diatomic)
    _check_whole_supplement_common(
        diatomic_representation;
        expected_kind = :bond_aligned_diatomic_cartesian_shell,
        expected_source_kind = :legacy_bond_aligned_diatomic_gaussian_supplement,
        expected_lmax = 1,
        expected_nuclei = nuclei,
        expected_uncontracted = false,
    )
    @test diatomic_representation.metadata.atom == "H"
    @test diatomic_representation.metadata.basis_name == "cc-pVTZ"
    @test length(diatomic_representation.orbitals) == 18
    @test first(_labels(diatomic_representation)) == "a_s1"
    @test last(_labels(diatomic_representation)) == "b_pz2"
    @test diatomic_representation.orbitals[1].center == nuclei[1]
    @test diatomic_representation.orbitals[end].center == nuclei[2]

    heteronuclear = legacy_bond_aligned_heteronuclear_gaussian_supplement(
        "He",
        "cc-pVTZ",
        "H",
        "cc-pVTZ",
        nuclei;
        lmax = 0,
    )
    @test heteronuclear isa LegacyBondAlignedHeteronuclearGaussianSupplement
    @test heteronuclear.atomic_sources[1].atom == "He"
    @test heteronuclear.atomic_sources[2].atom == "H"
    heteronuclear_representation = basis_representation(heteronuclear)
    _check_whole_supplement_common(
        heteronuclear_representation;
        expected_kind = :bond_aligned_heteronuclear_cartesian_shell,
        expected_source_kind = :legacy_bond_aligned_heteronuclear_gaussian_supplement,
        expected_lmax = 0,
        expected_nuclei = nuclei,
        expected_uncontracted = false,
    )
    @test heteronuclear_representation.metadata.atom == nothing
    @test heteronuclear_representation.metadata.basis_name == nothing
    @test heteronuclear_representation.metadata.basisfile == nothing
    @test length(heteronuclear_representation.orbitals) == 6
    @test _labels(heteronuclear_representation) == [
        "a_s1", "a_s2", "a_s3",
        "b_s1", "b_s2", "b_s3",
    ]
    @test heteronuclear_representation.orbitals[1].center == nuclei[1]
    @test heteronuclear_representation.orbitals[end].center == nuclei[2]
end
